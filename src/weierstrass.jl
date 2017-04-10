module WeierstrassCurves

import Nemo
import Base.show
import ..EllipticCurves: EllipticCurve, AbstractWeierstrass, EllipticPoint, Map, Isogeny, BaseCurve, Domain, Image, Eval

######################################################################
# Basic methods
######################################################################

"""
Concrete types for elliptic curves in long Weierstrass form over a ring.

Immutable types with 5 fieds : a1, a2, a3, a4, a6.
"""
immutable WeierstrassCurve{T} <: AbstractWeierstrass{T}
    a1::T
    a2::T
    a3::T
    a4::T
    a6::T
end



"""
Concrete types for elliptic curves in short Weierstrass form over a ring.

Immutable type with 2 fields : a, b.
"""
immutable ShortWeierstrassCurve{T} <: AbstractWeierstrass{T}
    a::T
    b::T
end



"""
Get a description of an elliptic curve in long Weierstrass form.

Shows an explicit equation and the base ring.
"""
function show{T}(io::IO, E::WeierstrassCurve{T})
    print(io, "Elliptic Curve in long Weierstrass form y² + $(E.a1) xy + $(E.a3) y = x³ + $(E.a2) x² + $(E.a4) x + $(E.a6)  over ")
    show(io, Nemo.parent_type(T))
end



"""
Get a description of an elliptic curve in short Weierstrass form.

Shows an explicit equation and the base ring.
"""
function show{T}(io::IO, E::ShortWeierstrassCurve{T})
    print(io, "Elliptic Curve in short Weierstrass form y² = x³ + $(E.a) x + $(E.b)  over ")
    show(io, Nemo.parent_type(T))
end



"""
Get an elliptic curve in long Weierstrass form from an elliptic curve in short Weierstrass form.

Returns an elliptic curve in short Weierstrass curve with the same equation, and two maps which are the canonical isomorphisms between the curves.
"""

function ToLongWeierstrass{T}(E::ShortWeierstrassCurve{T})

end


######################################################################
# Invariants
######################################################################



"""
Get the a-invariants (a1, a2, a3, a4, a6) of an elliptic curve in long Weierstrass form.

Returns a tuple of five elements in the base ring.
"""
function a_invariants(E::WeierstrassCurve)
    return (E.a1, E.a2, E.a3, E.a4, E.a6)
end



"""
Get the a-invariants (a1, a2, a3, a4, a6) of an elliptic curve in short Weierstrass form.

Returns a tuple of five elements in the base ring, the first three being zeroes.
"""
function a_invariants(E::ShortWeierstrassCurve)
    zero = Nemo.zero(Nemo.parent(E.a))
    return (zero, zero, zero, E.a, E.b)
end



"""
Get the b-invariants (b2, b4, b6, b8) of an elliptic curve in Weierstrass form.

Returns an array of four elements in the base ring.
"""
function b_invariants{T}(E::AbstractWeierstrass{T})
    a1, a2, a3, a4, a6 = a_invariants(E)
    return T[a1*a1 + 4*a2,
             a1*a3 + 2*a4,
             a3^2 + 4*a6,
             a1^2 * a6 + 4*a2*a6 - a1*a3*a4 + a2*a3^2 - a4^2]
end



"""
Get the c-invariants (c4, c6) of an elliptic curve in Weierstrass form.

Returns a tuple of two elements in the base ring.
"""
function c_invariants{T}(E::AbstractWeierstrass{T})
    b2, b4, b6, b8 = b_invariants(E)
    return (b2^2 - 24*b4, -b2^3 + 36*b2*b4 - 216*b6)
end




"""
Get the discriminant of an elliptic curve in long Weierstrass form.

Returns an element in the base ring.
"""
function discriminant{T}(E::WeierstrassCurve{T})
    b2, b4, b6, b8 = b_invariants(E)
    return -b2^2*b8 - 8*b4^3 - 27*b6^2 + 9*b2*b4*b6
end



"""
Get the discriminant of an elliptic curve in short Weierstrass form.

Returns an element in the base ring.
"""
function discriminant{T}(E::ShortWeierstrassCurve{T})
    return -16*(4*(E.a)^3 + 27*(E.b)^3)
end

"""
Get the j-invariant of an elliptic curve in long Weierstrass form.
This requires the base ring to be a field.

Returns an element in the base field.
"""
function j_invariant{T<:Nemo.FieldElem}(E::WeierstrassCurve{T})
    c4, _ = c_invariants(E)
    disc = discriminant(E)
    return c4^3 // disc
end



"""
Get the j-invariant of an elliptic curve in short Weierstrass form.
This requires the base ring to be a field.

Returns an element in the base field.
"""
function j_invariant{T<:Nemo.FieldElem}(E::ShortWeierstrassCurve{T})
    disc = discriminant(E)
    return -1728 * (4 * (E.a))^3 // disc
end


######################################################################
# Basic methods for projective points
######################################################################


"""
Concrete types for projective points on elliptic curves in Weierstrass form.

Mutable type with 3 fields : X, Y, Z.
"""
type ProjectivePoint{T, E <: AbstractWeierstrass{T}} <: EllipticPoint{T}
    X::T
    Y::T
    Z::T
    curve::E
end

"""
Describes a projective point giving its X, Y and Z coordinates.
"""
function show(io::IO, P::ProjectivePoint)
    print(io, "($(P.X):$(P.Y):$(P.Z))")
end
"""
Decides whether a projective point represents the point at infinity on a Weierstrass curve.
"""
is_identity(P::ProjectivePoint) = Nemo.iszero(P.Z)

"""
Get a normalized projective point from any projective point.
This requires the base ring to be a field.

Returns a new projective point with Z-coordinate equal to 1, without changing the input.
"""
function normalized{T<:Nemo.FieldElem, E}(P::ProjectivePoint{T,E})
    K = Nemo.parent(P.X)
    if is_identity(P)
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
function normalize!{T<:Nemo.FieldElem, E}(P::ProjectivePoint{T,E})
    K = Nemo.parent(P.X)
    if is_identity(P)
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



######################################################################
# Addition laws for projective points
######################################################################

"""
Get the opposite of a point on an elliptic curve in short Weierstrass form.
"""
function -{T, E::ShortWeierstrassCurve{T}}(P::ProjectivePoint{T,E})
    return ProjectivePoint(P.X, -P.Y, P.Z, E)
end

"""
Get the opposite of a point on an elliptic curve in general Weierstrass form.
"""
function -{T, E::WeierstrassCurve{T}}(P::ProjectivePoint{T,E})
    return ProjectivePoint(
		P.X, 
		- P.Y - (E.a1) * P.X - (E.a3) * P.Z, 
		P.Z, 
		E)
end

"""
Get the sum of two projective points on the same long Weierstrass curve, assuming they are not equal and not inverse of each other.

This does not assume the base ring is a field.
"""

function addgeneric{T, E::WeierstrassCurve{T}}(P::ProjectivePoint{T,E}, Q::ProjectivePoint{T,E})
    lambda = P.Z * Q.Y - P.Y * Q.Z
    nu = P.Y * Q.X - P.X * Q.Y
    denom = P.Z * Q.X - P.X * Q.Z
    Xplus = lambda^2 + (E.a1) * lambda * denom - (E.a2) * denom^2 - 
    Yplus =
    Zplus =
    return ProjectivePoint(Xplus, Yplus, Zplus, E)
end

######################################################################
# Affine points, just for laughs
######################################################################

type AffinePoint{T<:Nemo.FieldElem, E<:EllipticCurve} <: EllipticPoint{T}
    proj::ProjectivePoint{T,E}
end

function show(io::IO, P::AffinePoint)
    if is_identity(P.proj)
        print(io, "𝒪")
    else
        Q = normalized(P.proj)
        print(io, "($(Q.X), $(Q.Y))")
    end
end

end
