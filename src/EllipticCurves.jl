module EllipticCurves

import Nemo
import Base.show

export a_invariants, b_invariants, c_invariants, discriminant, j_invariant

######################################################################
# Abstract types
######################################################################

"""
Abstract classes for elliptic curves.

Every elliptic curve inherits from this.
"""
abstract EllipticCurve{T<:Nemo.RingElem}

"""
Get the a-invariants of an elliptic curve.

Every instance of EllipticCurve must implement this method.
"""
function a_invariants
end

"""
Abstract classes for elliptic curve points.

Every elliptic curve point inherits from this.
"""
abstract EllipticPoint{T<:Nemo.RingElem}

######################################################################
# Generic methods
######################################################################

function show{T}(io::IO, E::EllipticCurve{T})
    a1, a2, a3, a4, a6 = a_invariants(E)
    print(io, "Elliptic Curve  y² + $a1 xy + $a3 y = x³ + $a2 x² + $a4 x + $a6  over ")
    show(io, Nemo.parent_type(T))
end

function b_invariants{T}(E::EllipticCurve{T})
    a1, a2, a3, a4, a6 = a_invariants(E)
    return T[a1*a1 + 4*a2,
             a1*a3 + 2*a4,
             a3^2 + 4*a6,
             a1^2 * a6 + 4*a2*a6 - a1*a3*a4 + a2*a3^2 - a4^2]
end

function c_invariants{T}(E::EllipticCurve{T})
    b2, b4, b6, b8 = b_invariants(E)
    return (b2^2 - 24*b4, -b2^3 + 36*b2*b4 - 216*b6)
end

function discriminant{T}(E::EllipticCurve{T})
    b2, b4, b6, b8 = b_invariants(E)
    return -b2^2*b8 - 8*b4^3 - 27*b6^2 + 9*b2*b4*b6
end

function j_invariant{T<:Nemo.FieldElem}(E::EllipticCurve{T})
    c4, _ = c_invariants(E)
    disc = discriminant(E)
    return c4^3 // disc
end

######################################################################
# Weierstrass curves
######################################################################
include("weierstrass.jl")

end # module
