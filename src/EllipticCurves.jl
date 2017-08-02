
######################################################################
# The EllipticCurves module for Julia/Nemo
######################################################################

module EllipticCurves

######################################################################
# Julia/Nemo imports
######################################################################

import Nemo
import Nemo: base_ring, discriminant, degree, PolynomialRing, order, FinField, FinFieldElem, PolyElem, FiniteField, roots, issquare, parent, convert, trace, degree, coeff, characteristic, jacobi, Integer, AbsSeriesElem, SeriesRing, PowerSeriesRing, PolyRing, PolynomialRing, PolyElem, FieldElem, gen, shift_left, shift_right, truncate, inv, sqrt, set_prec!, divrem, compose, setcoeff!, sqrt, isprime, ZZ, GenFrac, GenMPoly, evaluate, divexact

import Base: show, normalize!, isvalid, ==, rand, *, +, -

######################################################################
# Julia/Nemo exports
######################################################################

export EllipticCurve, AbstractWeierstrass, ProjectivePoint, Map, ExplicitMap, Isogeny, Weierstrass, ShortWeierstrass, EllipticPoint, MontgomeryCurve, XonlyPoint

#Exporting functions

export Eval, base_curve, base_ring, a_invariants, b_invariants, c_invariants, j_invariant, areequal, isinfinity, normalized, infinity, minus, plus, isvalid, tolongWeierstrass, toshortWeierstrass, xonly, isfixedtorsion, xinfinity, fixedtorsion, times, codomain, domain, kernel, divisionpolynomial, random, cardinality, frobeniuspolynomial, any_root, randomxonly, torsionpoint, torsionxonly, base_extend, isgoodprime, goodprimes, first_isogeny, first_isogeny_x, image, squarefree_kernel

#Functions we may want to turn internal

export xdouble, xadd, addequalx, addgeneric, xladder, subgrouppolynomial

######################################################################
# Abstract types
######################################################################

"""
Abstract class for elliptic curves.

Every elliptic curve concrete type inherits from this.
"""
abstract type EllipticCurve{T<:Nemo.RingElem} end

"""
Abstract class for elliptic curves in Weierstrass form.

Every concrete type for elliptic curves in Weierstrass form inherits from this.
"""
abstract type AbstractWeierstrass{T<:Nemo.RingElem} <: EllipticCurve{T} end


"""
Abstract class for elliptic curve points.

Every concrete type for points inherits from this.
"""
abstract type ProjectivePoint{T<:Nemo.RingElem} end


######################################################################
# Points on elliptic curves
######################################################################

include("tools.jl")


######################################################################
# Points on elliptic curves
######################################################################

include("points.jl")


######################################################################
# Weierstrass curves
######################################################################

include("weierstrass.jl")


######################################################################
# Points on Weierstrass curves
######################################################################

include("weierstrasspoints.jl")


######################################################################
# Montgomery curves
######################################################################

include("montgomery.jl")


######################################################################
# Points on Montgomery curves
######################################################################

include("montgomerypoints.jl")


######################################################################
# Isogenies
######################################################################

include("isogenies.jl")



######################################################################
# Over finite fields
######################################################################

include("finfields.jl")


end # module
