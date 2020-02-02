
######################################################################
# points.jl: projective points on elliptic curves
######################################################################

export EllipticPoint, Point

export coordinates, base_curve, isinfinity, normalize!, normalized, samefields

"""
Concrete type for projective points on elliptic curves given by a projective planar equation.
"""

mutable struct EllipticPoint{T<:Nemo.RingElem} <: ProjectivePoint{T}
	X::T
	Y::T
	Z::T
	curve::EllipticCurve{T}
end

coordinates(P::EllipticPoint) = (P.X, P.Y, P.Z)

base_curve(P::EllipticPoint) = P.curve

######################################################################
# Calling with the 'Point' constructor
######################################################################

function Point(x::T, y::T, curve::EllipticCurve{T}) where T
	R = base_ring(x)
	return EllipticPoint(x, y, R(1), curve)
end

function Point(X::T, Y::T, Z::T, curve::EllipticCurve{T}) where T
	return EllipticPoint(X, Y, Z, curve)
end

######################################################################
# Basic methods for projective points
######################################################################



base_ring(P::EllipticPoint) = parent(P.X)

"""
Decide whether a given elliptic point is valid.
"""
function isvalid(P::EllipticPoint)
	Eq = projective_equation(P.curve)
	x, y, z = coordinates(P)
	return Eq(x, y, z) == 0 && (x, y, z) != (0, 0, 0)
end


"""
Decide whether two projective points are given by the exact same coordinates.
"""
samefields(P::EllipticPoint, Q::EllipticPoint) = (P.curve == Q.curve) & (P.X == Q.X) & (P.Y == Q.Y) & (P.Z == Q.Z)

"""
Describe a projective point giving its X, Y and Z coordinates, and the curve it lives on.
"""
show(io::IO, P::EllipticPoint) = print("Point (", P.X, ":", P.Y, ":", P.Z, ") on ", P.curve)


"""
Decide whether a projective point is at infinity.
"""
isinfinity(P::EllipticPoint) = iszero(P.Z)

"""
Get a new normalized projective point from any projective point.
"""
function normalized(P::EllipticPoint)
    K = base_ring(P)
    if isinfinity(P)
        return EllipticPoint(zero(K),
                               one(K),
                               zero(K),
                               P.curve)
    else
        return EllipticPoint(P.X // P.Z,
                               P.Y // P.Z,
                               one(K),
                               P.curve)
    end
end

"""
Normalize a projective point in place.
"""
function normalize!(P::EllipticPoint{T}) where T<:Nemo.FieldElem
    K = base_ring(P)
    if isinfinity(P)
        P.X = zero(K)
        P.Y = one(K)
        P.Z = zero(K)
    else
        P.X = P.X // P.Z
        P.Y = P.Y // P.Z
        P.Z = one(K)
    end
    return
end

"""
Decides whether two projective points are the same.
"""

==(P::EllipticPoint, Q::EllipticPoint) = samefields(normalized(P), normalized(Q))

