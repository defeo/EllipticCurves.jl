

######################################################################
# montgomerypoints.jl: xz-points for Montgomery curves
######################################################################

export XZPoint, xadd, xdouble, xinfinity, fixedtorsion, isfixedtorsion, XZzero

######################################################################
# Type definitions
######################################################################

"""
Concrete type for xz-points on Montgomery curves.
"""
mutable struct XZPoint{T<:Nemo.RingElem} <: ProjectivePoint{T}
    X::T
    Z::T
    curve::Montgomery{T}
end

base_ring(P::XZPoint) = Nemo.parent(P.X)

base_curve(P::XZPoint) = P.curve


"""
Get a description of an xz-point.
"""
show(io::IO, P::XZPoint) = print(io, "($(P.X):$(P.Z))")

"""
Create an xz-point from a projective point with 3 coordinates.
"""
function XZPoint(P::EllipticPoint{T}) where T<:Nemo.RingElem
    return P.Z == 0 ? xinfinity(P.curve) : XZPoint(P.X, P.Z, P.curve)
end


"""
Get a normalized xz-point from any xz-point.
"""
function normalized(P::XZPoint{T}) where T
    K = Nemo.parent(P.X)
    if isinfinity(P)
        return XZPoint(Nemo.one(K), 
                               Nemo.zero(K), 
                               P.curve)
    else
        return XZPoint(P.X // P.Z,
                               Nemo.one(K),
                               P.curve)
    end
end

"""
Normalizes a xz-point to a xz-point with Z-coordinate 1, in place.
"""
function normalize!(P::XZPoint{T}) where T<:Nemo.FieldElem
    K = Nemo.parent(P.X)
    if isinfinity(P)
        P.X = Nemo.one(K)
    else
        P.X = P.X // P.Z
        P.Z = Nemo.one(K)
    end
    return
end



"""
Decide whether two xz-points are given by the exact same coordinates.
"""
function samefields(P::XZPoint, Q::XZPoint)
	return (P.curve == Q.curve) & (P.X == Q.X) & (P.Z == Q.Z)
end

"""
Decide whether two xz-points are equal.
"""
==(P::XZPoint, Q::XZPoint) = P.X * Q.Z - P.Z * Q.X == 0


"""
Decide whether an xz-point is the point at infinity.
"""
isinfinity(P::XZPoint) = iszero(P.Z)

"""
Decide whether an xz-point is the fixed 2-torsion point of the Montgomery model (at the origin).
"""
isfixedtorsion(P::XZPoint) = iszero(P.X)



"""
Get the xz-point at infinity on a Montgomery curve.
"""
function xinfinity(E::Montgomery)
    K = Nemo.parent(E.A)
    return XZPoint(Nemo.one(K), Nemo.zero(K), E)
end

function XZzero(E::Montgomery)
	return xinfinity(E)
end


"""
Get the fixed xz- 2-torsion point on a Montgomery curve.
"""
function fixedtorsion(E::Montgomery)
    K = Nemo.parent(E.A)
    return XZPoint(Nemo.zero(K), Nemo.one(K), E)
end

function isvalid(P::XZPoint)
    return !(P.X == P.Z == 0) && issquare( P.curve.B * ( P.X^2 + P.curve.A*P.X*P.Z + P.Z^2 ) * P.X * P.Z )[1]
end




######################################################################
# Montgomery arithmetic
######################################################################


"""
Double any xz-point using the least possible field operations.
"""
function xdouble(P::XZPoint{T}) where T
    E = P.curve
    R = base_ring(E)

    v1 = P.X + P.Z
    v1 = v1^2
    v2 = P.X - P.Z
    v2 = v2^2
    
    X2 = v1 * v2
    
    v1 = v1 - v2
    v3 = ((E.A + 2) // 4) * v1
    v3 = v3 + v2
    
    Z2 = v1 * v3
    
    return XZPoint(X2, Z2, E)
end

"""
Differential addition on xz-points using the least possible field operations.

This function assumes the difference is not (0:0) or (0:1).
"""
function xadd(P::XZPoint{T}, Q::XZPoint{T}, Minus::XZPoint{T}) where T
    v0 = P.X + P.Z
    v1 = Q.X - Q.Z
    v1 = v1 * v0
    
    v0 = P.X - P.Z
    v2 = Q.X + Q.Z
    v2 = v2 * v0
    
    v3 = v1 + v2
    v3 = v3^2
    v4 = v1 - v2
    v4 = v4^2

    Xplus = Minus.Z * v3
    Zplus = Minus.X * v4
   
    return XZPoint(Xplus, Zplus, P.curve)
end

"""
Montgomery ladder to compute scalar multiplications of generic xz-points, using the least possible field operations.
"""
function xladder(k::Nemo.Integer, P::XZPoint)
    x0, x1 = P, xdouble(P)
    for b in Iterators.drop(Iterators.reverse(digits(Int8, k, base=2)), 1)
        if (b == 0)
        	x0, x1 = xdouble(x0), xadd(x0, x1, P)
        else
        	x0, x1 = xadd(x0, x1, P), xdouble(x1)
        end
    end
    return x0
end

"""
    Top-level function for scalar multiplications with xz-points on Montgomery curves
    """
function *(k::Nemo.Integer, P::XZPoint)
    E = P.curve
    if k == 0
	return xinfinity(E)
    elseif k<0
	return times(-k, P)
    else
	if isinfinity(P)
	    return xinfinity(E)
	elseif isfixedtorsion(P)
	    if k % 2 == 0
		return xinfinity(E)
	    else
		return fixedtorsion(E)
	    end
	else
	    return xladder(k, P)
	end
    end
end

*(k::Nemo.fmpz, P::XZPoint) = BigInt(k) * P

