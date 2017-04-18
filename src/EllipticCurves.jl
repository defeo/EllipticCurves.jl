module EllipticCurves

import Nemo
import Nemo: base_ring, discriminant, degree

import Base: show, normalize!, isvalid, ==

export EllipticCurve, AbstractWeierstrass, ProjectivePoint, Map, base_curve, Eval, base_ring, a_invariants, areequal, EllipticPoint, base_ring, isinfinity, normalized, infinity, minus, plus, isvalid, ExplicitMap, Isogeny, WeierstrassCurve, ShortWeierstrassCurve, b_invariants, c_invariants, j_invariant, tolongWeierstrass, toshortWeierstrass, MontgomeryCurve, XonlyPoint, xonly, isfixedtorsion, xinfinity, fixedtorsion, times, codomain, domain, kernel, divisionpolynomial

######################################################################
# Abstract types
######################################################################

"""
Abstract classes for elliptic curves.

Every elliptic curve inherits from this.
"""
abstract EllipticCurve{T<:Nemo.RingElem}

"""
Abstract classes for elliptic curves in Weierstrass form.

Every elliptic curve in Weierstrass form inherits from this.
"""
abstract AbstractWeierstrass{T<:Nemo.RingElem} <: EllipticCurve{T}


"""
Abstract classes for elliptic curve points.

Every elliptic curve point inherits from this.
"""
abstract ProjectivePoint{T<:Nemo.RingElem}

"""
Abstract class for maps.

Every map between elliptic curves inherits from this.
"""
abstract Map{T<:Nemo.RingElem}


######################################################################
# Abstract functions
######################################################################


function a_invariants(E::AbstractWeierstrass)
end


function ==(E1::AbstractWeierstrass, E2::AbstractWeierstrass)
	return (a_invariants(E1) == a_invariants(E2))
end

function areequal(P::ProjectivePoint, Q::ProjectivePoint)
end

function ispoint{T}(a::T, b::T, c::T, E::EllipticCurve{T})
end

######################################################################
# Points on elliptic curves
######################################################################

include("points.jl")

using .Points


######################################################################
# Maps
######################################################################

include("maps.jl")

using .Maps



######################################################################
# Weierstrass curves
######################################################################

include("weierstrass.jl")

using .Weierstrass

######################################################################
# Montgomery curves
######################################################################

include("montgomery.jl")

using .Montgomery



end # module
