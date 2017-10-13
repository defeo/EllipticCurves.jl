module EllipticCurves

import Nemo
import Nemo: base_ring, discriminant, degree, PolynomialRing, order, FinField, FinFieldElem, PolyElem, FiniteField, roots, issquare, parent, convert, trace, degree, coeff, characteristic, jacobi, Integer, AbsSeriesElem, SeriesRing, PowerSeriesRing, PolyRing, PolynomialRing, PolyElem, FieldElem, gen, shift_left, shift_right, truncate, inv, sqrt, set_prec!, divrem, compose, setcoeff!, sqrt, isprime, ZZ

import Base: show, normalize!, isvalid, ==, rand

#Exporting types

export EllipticCurve, AbstractWeierstrass, ProjectivePoint, Map, ExplicitMap, Isogeny, WeierstrassCurve, ShortWeierstrassCurve, EllipticPoint, MontgomeryCurve, XonlyPoint

#Exporting functions

export Eval, base_curve, base_ring, a_invariants, b_invariants, c_invariants, j_invariant, areequal, isinfinity, normalized, infinity, minus, plus, isvalid, tolongWeierstrass, toshortWeierstrass, xonly, isfixedtorsion, xinfinity, fixedtorsion, times, codomain, domain, kernel, divisionpolynomial, random, cardinality, frobeniuspolynomial, any_root, randomxonly, torsionpoint, torsionxonly, base_extend, isgoodprime, goodprimes, first_isogeny, first_isogeny_x

#Functions we may want to turn internal

export xdouble, xadd, addequalx, addgeneric, xladder, subgrouppolynomial

######################################################################
# Abstract types
######################################################################

"""
Abstract classes for elliptic curves.

Every elliptic curve inherits from this.
"""
abstract type EllipticCurve{T<:Nemo.RingElem} end

"""
Abstract classes for elliptic curves in Weierstrass form.

Every elliptic curve in Weierstrass form inherits from this.
"""
abstract type AbstractWeierstrass{T<:Nemo.RingElem} <: EllipticCurve{T} end


"""
Abstract classes for elliptic curve points.

Every elliptic curve point inherits from this.
"""
abstract type ProjectivePoint{T<:Nemo.RingElem} end

"""
Abstract class for maps.

Every map between elliptic curves inherits from this.
"""
abstract type Map{T<:Nemo.RingElem} end


######################################################################
# Points on elliptic curves
######################################################################

include("points.jl")


######################################################################
# Maps
######################################################################

include("maps.jl")


######################################################################
# Weierstrass curves
######################################################################

include("weierstrass.jl")


######################################################################
# Montgomery curves
######################################################################

include("montgomery.jl")


######################################################################
# Over finite fields
######################################################################

include("finfields.jl")


end # module
