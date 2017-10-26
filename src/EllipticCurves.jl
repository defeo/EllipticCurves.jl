
######################################################################
# The EllipticCurves module for Julia/Nemo
######################################################################

module EllipticCurves


######################################################################
# Conventions
######################################################################

#Polynomial rings have variables X, Y, Z, U, V...


######################################################################
# Julia/Nemo imports
######################################################################

import Nemo
import Nemo: base_ring, discriminant, degree, PolynomialRing, order, FinField, FinFieldElem, PolyElem, FiniteField, roots, issquare, parent, convert, trace, degree, coeff, characteristic, jacobi, Integer, AbsSeriesElem, SeriesRing, PowerSeriesRing, PolyRing, PolynomialRing, PolyElem, FieldElem, shift_left, shift_right, truncate, inv, sqrt, set_prec!, divrem, compose, setcoeff!, sqrt, isprime, ZZ, GenFrac, GenMPoly, evaluate, divexact, gen, derivative, gens, zero, one, ResidueRing, ResRing, div, FieldElem, rem, FractionField, numerator, denominator, sqrt, gcd

import Base: show, normalize!, isvalid, ==, rand, *, +, -

######################################################################
# Julia/Nemo exports
######################################################################

export EllipticCurve, AbstractWeierstrass, ProjectivePoint

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
# Tools
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
# Edwards curves
######################################################################

include("edwards.jl")


######################################################################
# Points on Edwards curves
######################################################################

include("edwardspoints.jl")


######################################################################
# Isogenies
######################################################################

include("isogenies.jl")

#=
######################################################################
# Over finite fields
######################################################################

include("finfields.jl")

=#
end # module
