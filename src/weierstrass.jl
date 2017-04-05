module WeierstrassCurves

import Nemo
import Base.show
import ..EllipticCurves: EllipticCurve, a_invariants, EllipticPoint

immutable WeierstrassCurve{T} <: EllipticCurve{T}
    a1::T
    a2::T
    a3::T
    a4::T
    a6::T
end

function a_invariants(E::WeierstrassCurve)
    return (E.a1, E.a2, E.a3, E.a4, E.a6)
end

immutable SimplifiedWeierstrassCurve{T} <: EllipticCurve{T}
    a::T
    b::T
end

function a_invariants(E::SimplifiedWeierstrassCurve)
    zero = Nemo.zero(Nemo.parent(E.a))
    return (zero, zero, zero, E.a, E.b)
end

function show{T}(io::IO, E::SimplifiedWeierstrassCurve{T})
    print(io, "Elliptic Curve  y² = x³ + $(E.a) x + $(E.b)  over ")
    show(io, Nemo.parent_type(T))
end

######################################################################
# Elliptic curve points
######################################################################

type ProjectivePoint{T, E<:EllipticCurve} <: EllipticPoint{T}
    X::T
    Y::T
    Z::T
    curve::E
end

show(io::IO, P::ProjectivePoint) = print(io, "($(P.X):$(P.Y):$(P.Z))")

is_identity(P::ProjectivePoint) = Nemo.iszero(P.Z)

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
