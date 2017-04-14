module Montgomery

import Nemo
import Base.show
import ..EllipticCurves: EllipticCurve, ProjectivePoint, Weierstrass.WeierstrassCurve, Points.EllipticPoint, Weierstrass.AbstractWeierstrass


######################################################################
# Basic methods
######################################################################

"""
Concrete type for (twisted) Montgomery curves.
"""
immutable MontgomeryCurve{T<:Nemo.RingElem} <: AbstractWeierstrass{T}
    A::T
    B::T
end

function base_ring(E::MontgomeryCurve)
	return Nemo.parent(E.A)
end

"""
Get a description of a Montgomery curve.
"""
function show{T}(io::IO, E::MontgomeryCurve{T})
    print(io, "Elliptic Curve  $(E.B) y² = x³ + $(E.A) x² + x  over ")
    show(io, base_ring(E))
end

"""
Get the a-invariants of a Montgomery Curve.
"""

function a_invariants{T}(E::MontgomeryCurve{T})
	R = base_ring(E)
	zero = Nemo.zero(R)
	return (zero, E.A, zero, E.B, zero)
end


######################################################################
# Model changes
######################################################################

"""
Given an elliptic curve in Montgomery form, build an elliptic curve in long Weierstrass form isomorphic to it.

Returns an elliptic curve and two explicit maps givint the change of variables.
"""
function tolongWeierstrass{T<:Nemo.FieldElem}(E::MontgomeryCurve{T})
    zero = Nemo.zero(Nemo.parent(E.A))
    E1 = WeierstrassCurve(zero, E.A // E.B, zero, Nemo.inv(E.B ^ 2), zero)
	phi1 = ExplicitMap(E, E1,
		function(P::EllipticPoint)
			Xprime = P.X // E.B
			Yprime = P.Y // E.B
			Zprime = P.Z
			return EllipticPoint(Xprime, Yprime, Zprime, E1)
		end)
	phi2 = ExplicitMap(E1, E,
		function(P::EllipticPoint)
			Xprime = P.X * E.B
			Yprime = P.Y * E.B
			Zprime = P.Z
			return EllipticPoint(Xprime, Yprime, Zprime, E)
		end)
	return (E1, phi1, phi2)
end




######################################################################
# X-only points for Montgomery curves
######################################################################

"""
Concrete type for x-only points to speed up arithmetic.
"""
type XonlyPoint{T<:Nemo.RingElem, form<:MontgomeryCurve}
    X::T
    Z::T
    curve::form
end


"""
Get a description of an x-only point.
"""
show(io::IO, P::XonlyPoint) = print(io, "($(P.X):$(P.Z))")

"""
Create an x-only point from a projective point with 3 coordinates.
"""
function xonly{T<:Nemo.RingElem}(P::EllipticPoint{T, MontgomeryCurve{T}})
	return XonlyPoint(P.X, P.Z, P.curve)
end


"""
Get a normalized x-only point from any x-only point.
This requires the base ring to be a field.

Returns a new x-only point with Z-coordinate equal to 1, without changing the input.
"""
function normalized{T<:Nemo.FieldElem}(P::XonlyPoint{T, MontgomeryCurve{T}})
    K = Nemo.parent(P.X)
    if isinfinity(P)
        return XonlyPoint(Nemo.zero(K), 
                               Nemo.zero(K), 
                               P.curve)
    else
        return XonlyPoint(P.X // P.Z,
                               Nemo.one(K),
                               P.curve)
    end
end

"""
Normalizes a x-only point to a x-only point with Z-coordinate 1.

Does not create a new point, and changes the input.
"""
function normalize!{T<:Nemo.FieldElem}(P::XonlyPoint{T, MontgomeryCurve{T}})
    K = Nemo.parent(P.X)
    if isinfinity(P)
        P.X = Nemo.zero(K)
    else
        P.X = P.X // P.Z
        P.Z = Nemo.one(K)
    end
    return
end


"""
Decide whether an x-only point is the point at infinity.
"""
isinfinity(P::XonlyPoint) = Nemo.iszero(P.Z)

"""
Decide whether an x-only point is the fixed 2-torsion point of the Montgomery model (at the origin).
"""
isfixedtorsion(P::XonlyPoint) = Nemo.iszero(P.X)



"""
Get the x-only point at infinity on a Montgomery curve.
"""
function xinfinity(E::MontgomeryCurve)
    K = Nemo.parent(E.A)
    zero = Nemo.zero(K)
    return XonlyPoint(zero, zero, E)
end


"""
Get the fixed x-only 2-torsion point on a Montgomery curve.
"""
function fixedtorsion(E::MontgomeryCurve)
    K = Nemo.parent(E.A)
    return XonlyPoint(Nemo.zero(K), Nemo.one(K), E)
end






######################################################################
# Montgomery arithmetic
######################################################################

#taken from Ben Smith's review article

"""
Double any x-only point using the least possible field operations.
"""
function xdouble{T<:Nemo.RingElem}(P::XonlyPoint{T, MontgomeryCurve{T}})
	E = P.curve
	
    v1 = P.X + P.Z
    v1 = v1^2
    v2 = P.X - P.Z
    v2 = v2^2
    
    X2 = v1 * v2
    
    v1 = v1 - v2
    v3 = ((E.A + 2) // 4) * v1
    v3 = v3 + v2
    
    Z2 = v1 * v3
    
    if Z2 == 0
        X2 = 0
    end
    
    return XonlyPoint(X2, Z2, E)
end

"""
Differential addition on x-only points using the least possible field operations.

This function assumes the difference is not (0:0) or (0:1).
"""
function xadd{T<:Nemo.RingElem}(P::XonlyPoint{T, MontgomeryCurve{T}}, Q::XonlyPoint{T, MontgomeryCurve{T}}, Minus::XonlyPoint{T, MontgomeryCurve{T}})
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
   
    return XonlyPoint(Xplus, Zplus, P.curve)
end

"""
Montgomery ladder to compute scalar multiplications of generic x-only points, using the least possible field operations.
"""
function xladder{T<:Nemo.RingElem}(k::Nemo.Integer, P::XonlyPoint{T, MontgomeryCurve{T}})
	normalize!(P)
    x0, x1 = P, xdouble(P)
    for b in bin(k)[2:end]
        if (b == 0)
        	x0, x1 = xdouble(x0), xadd(x0, x1, P)
        else
        	x0, x1 = xadd(x0, x1, P), xdouble(x1)
        end
    end
    return x0
end

"""
Top-level function for scalar multiplications with x-only points on Montgomery curves
"""
function times{T<:Nemo.RingElem}(k::Nemo.Integer, P::XonlyPoint{T, MontgomeryCurve{T}})
	E = P.curve
	if k==0
		return infinity(E)
	elseif k<0
		return times(-k, P)
	else
		if isinfinity(P)
			return infinity(E)
		elseif isfixedtorsion(P)
			if k % 2 == 0
				return infinity(E)
			else
				return fixedtorsion(E)
			end
		else
			return xladder(k, P)
		end
	end
end



end # module





