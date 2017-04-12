module Points

import Nemo
import Base.show

import ..EllipticCurves: EllipticCurve, EllipticPoint

######################################################################
# Basic methods for projective points
######################################################################


"""
Concrete types for projective points on elliptic curves.
"""

type ProjectivePoint{T<:Nemo.RingElem, form<:EllipticCurve} <: EllipticPoint{T}
	X::T
	Y::T
	Z::T
	curve::form
end

function basering(P::ProjectivePoint)
	return Nemo.parent(P.X)
end

"""
Describes a projective point giving its X, Y and Z coordinates.
"""
function show(io::IO, P::ProjectivePoint)
    print(io, "($(P.X):$(P.Y):$(P.Z))")
end


"""
Decides whether a projective point is at infinity.
"""
isinfinity(P::ProjectivePoint) = Nemo.iszero(P.Z)

"""
Get a normalized projective point from any projective point.
This requires the base ring to be a field.

Returns a new projective point with Z-coordinate equal to 1, without changing the input.
"""
function normalized{T<:Nemo.FieldElem, form}(P::ProjectivePoint{T, form})
    K = Nemo.parent(P.X)
    if isinfinity(P)
        return ProjectivePoint(Nemo.zero(K),
                               Nemo.one(K),
                               Nemo.zero(K),
                               P.curve)
    else
        return ProjectivePoint(P.X // P.Z,
                               P.Y // P.Z,
                               Nemo.one(K),
                               P.curve)
    end
end

"""
Normalizes a projective point to a projective point with Z-coordinate 1.

Does not create a new point, and changes the input.
"""
function normalize!{T<:Nemo.FieldElem, form}(P::ProjectivePoint{T, form})
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

function areequal(P::ProjectivePoint, Q::ProjectivePoint)
	return normalized(P)==normalized(Q)
end


end #module
