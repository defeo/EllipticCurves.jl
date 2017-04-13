module Points

import Nemo
import Base.show

import ..EllipticCurves: EllipticCurve, ProjectivePoint, AbstractWeierstrass



######################################################################
# Basic methods for projective points
######################################################################


"""
Concrete types for projective points on elliptic curves.
"""

type EllipticPoint{T<:Nemo.RingElem, form<:EllipticCurve} <: ProjectivePoint{T}
	X::T
	Y::T
	Z::T
	curve::form
end

function basering(P::EllipticPoint)
	return Nemo.parent(P.X)
end

"""
Describes a projective point giving its X, Y and Z coordinates.
"""
function show(io::IO, P::EllipticPoint)
    print(io, "($(P.X):$(P.Y):$(P.Z))")
end


"""
Decides whether a projective point is at infinity.
"""
isinfinity(P::EllipticPoint) = Nemo.iszero(P.Z)


"""
Get a normalized projective point from any projective point.
This requires the base ring to be a field.

Returns a new projective point with Z-coordinate equal to 1, without changing the input.
"""
function normalized{T<:Nemo.FieldElem, form}(P::EllipticPoint{T, form})
    K = Nemo.parent(P.X)
    if isinfinity(P)
        return EllipticPoint(Nemo.zero(K),
                               Nemo.one(K),
                               Nemo.zero(K),
                               P.curve)
    else
        return EllipticPoint(P.X // P.Z,
                               P.Y // P.Z,
                               Nemo.one(K),
                               P.curve)
    end
end

"""
Normalizes a projective point to a projective point with Z-coordinate 1.

Does not create a new point, and changes the input.
"""
function normalize!{T<:Nemo.FieldElem, form}(P::EllipticPoint{T, form})
    K = Nemo.parent(P.X)
    if isinfinity(P)
        P.X = Nemo.zero(K)
        P.Y = Nemo.one(K)
        P.Z = Nemo.zero(K)
    else
        P.X = P.X // P.Z
        P.Y = P.Y // P.Z
        P.Z = Nemo.one(K)
    end
    return
end

"""
Decides whether two projective points are the same.
"""

function areequal(P::EllipticPoint, Q::EllipticPoint)
	return normalized(P) == normalized(Q) & P.curve == Q.curve
end


######################################################################
# Addition laws for projective points on Weierstrass curves
######################################################################


"""
Get the point at infinity on an elliptic curve in Weierstrass form.
"""
function infinity{T}(E::AbstractWeierstrass{T})
    R = basering(E)
    return EllipticPoint(Nemo.zero(R), Nemo.one(R), Nemo.zero(R), E)
end


"""
Get the opposite of a point on an elliptic curve in Weierstrass form.
"""
function minus{T<:Nemo.RingElem}(P::EllipticPoint{T, AbstractWeierstrass{T}})
	E = P.curve
	a1, _, a3, _, _ = a_invariants(E)
    return EllipticPoint(P.X, - P.Y - a1 * P.X - a3 * P.Z, P.Z, E)
end

"""
Get the sum of two normalized projective points on the same Weierstrass curve, assuming they are not equal and not inverse of each other.

This assumes the base ring is a field.
"""
function addgeneric{T<:Nemo.FieldElem}(P::EllipticPoint{T, AbstractWeierstrass{T}}, Q::EllipticPoint{T, AbstractWeierstrass{T}})
    E = P.curve
	a1, a2, a3, _, _ = a_invariants(E)
    denom = Q.X - P.X
    lambda = (Q.Y - P.Y) // denom
    nu = (P.Y * Q.X - P.X * Q.Y) // denom
    
    Xplus = lambda^2 + a1 * lambda - a2 - P.X - Q.X
    Yplus = -(lambda + a1) * Xplus - nu - a3
    Zplus = Nemo.one(basering(P))
    return EllipticPoint(Xplus, Yplus, Zplus, E)
end

"""
Get the sum of two normalized projective points on the same long Weierstrass curve, assuming they have equal x-coordinate.

This assumes the base ring is a field.
"""
function addequalx{T<:Nemo.FieldElem}(P::EllipticPoint{T, AbstractWeierstrass{T}}, Q::EllipticPoint{T, AbstractWeierstrass{T}})
    E = P.curve
	a1, a2, a3, a4, a6 = a_invariants(E)
    denom = P.Y + Q.Y + a1 * Q.X + a3
    if Nemo.iszero(denom)
	    return infinity(E)
    else
	    lambda = (3 * (P.X)^2 + 2 * a2 * P.X + a4 - a1 * P.Y) // denom
	    nu = (- (P.X)^3 + a4 * P.X + 2 * a6 - a3 * P.Y) // denom
	    Xplus = lambda^2 + a1 * lambda - a2 - P.X - Q.X
        Yplus = -(lambda + a1) * Xplus - nu - a3
        Zplus = Nemo.one(basering(P))
	    return EllipticPoint(Xplus, Yplus, Zplus, E)
    end
end

"""
Get the sum of two projective points on the same long Weierstrass curve.

This assumes the base ring is a field.
"""
function plus{T<:Nemo.FieldElem}(P::EllipticPoint{T, AbstractWeierstrass{T}}, Q::EllipticPoint{T, AbstractWeierstrass{T}})
    normalize!(P)
    normalize!(Q)
    if isinfinity(P)
        return Q
    elseif isinfinity(Q)
        return P
    elseif P.X == Q.X
		return addequalx(P,Q)
    else
		return addgeneric(P,Q)
    end
end





end #module
