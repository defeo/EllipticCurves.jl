module Weierstrass


import Nemo

import ..EllipticCurves: EllipticCurve, AbstractWeierstrass, ProjectivePoint, Map, ExplicitMap, Eval, Isogeny, EllipticPoint, base_ring, normalize!, isinfinity, isvalid, a_invariants, show, ==, discriminant, ispoint

export WeierstrassCurve, ShortWeierstrassCurve, b_invariants, c_invariants, j_invariant, tolongWeierstrass, toshortWeierstrass

######################################################################
# Basic methods
######################################################################

"""
Concrete types for elliptic curves in long Weierstrass form over a ring.
"""
immutable WeierstrassCurve{T} <: AbstractWeierstrass{T}
    a1::T
    a2::T
    a3::T
    a4::T
    a6::T
end

function base_ring(E::WeierstrassCurve)
    return Nemo.parent(E.a1)
end




"""
Concrete types for elliptic curves in short Weierstrass form over a ring.
"""
immutable ShortWeierstrassCurve{T} <: AbstractWeierstrass{T}
    a::T
    b::T
end

function base_ring(E::ShortWeierstrassCurve)
    return Nemo.parent(E.a)
end

"""
Get a description of an elliptic curve in long Weierstrass form.

Shows an explicit equation and the base ring.
"""
function show(io::IO, E::WeierstrassCurve)
    print(io, "Elliptic Curve in long Weierstrass form y² + $(E.a1) xy + $(E.a3) y = x³ + $(E.a2) x² + $(E.a4) x + $(E.a6)  over ")
    show(io, base_ring(E))
end



"""
Get a description of an elliptic curve in short Weierstrass form.

Shows an explicit equation and the base ring.
"""
function show(io::IO, E::ShortWeierstrassCurve)
    print(io, "Elliptic Curve in short Weierstrass form y² = x³ + $(E.a) x + $(E.b)  over ")
    show(io, base_ring(E))
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
    zero = Nemo.zero(base_ring(E))
    return (zero, zero, zero, E.a, E.b)
end



"""
Get the b-invariants (b2, b4, b6, b8) of an elliptic curve in Weierstrass form.

Returns a tuple of four elements in the base ring.
"""
function b_invariants(E::AbstractWeierstrass)
    a1, a2, a3, a4, a6 = a_invariants(E)
    return (a1*a1 + 4*a2,
             a1*a3 + 2*a4,
             a3^2 + 4*a6,
             a1^2 * a6 + 4*a2*a6 - a1*a3*a4 + a2*a3^2 - a4^2)
end



"""
Get the c-invariants (c4, c6) of an elliptic curve in Weierstrass form.

Returns a tuple of two elements in the base ring.
"""
function c_invariants(E::AbstractWeierstrass)
    b2, b4, b6, b8 = b_invariants(E)
    return (b2^2 - 24*b4, -b2^3 + 36*b2*b4 - 216*b6)
end




"""
Get the discriminant of an elliptic curve in long Weierstrass form.

Returns an element in the base ring.
"""
function discriminant(E::AbstractWeierstrass)
    b2, b4, b6, b8 = b_invariants(E)
    return -b2^2*b8 - 8*b4^3 - 27*b6^2 + 9*b2*b4*b6
end

"""
Check if a given instance of AbstractWeierstrass is indeed a nonsingular curve.
"""
function isvalid(E::AbstractWeierstrass)
	return !Nemo.iszero(discriminant(E)) #isunit ?
end


"""
Get the discriminant of an elliptic curve in short Weierstrass form.

Returns an element in the base ring.
"""
function discriminant(E::ShortWeierstrassCurve)
    return -16*(4*(E.a)^3 + 27*(E.b)^2)
end

"""
Get the j-invariant of an elliptic curve in long Weierstrass form.
This requires the base ring to be a field.

Returns an element in the base field.
"""
function j_invariant{T<:Nemo.FieldElem}(E::AbstractWeierstrass{T})
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

"""
Decide if a triple satisfies the equation of the curve.
"""

function ispoint{T}(x::T, y::T, z::T, E::AbstractWeierstrass{T})
	a1, a2, a3, a4, a6 = a_invariants(E)
	return z * y^2 + a1 * x * y * z + a3 * y * z^2 == x^3 + a2 * x^2 * z + a4 * x * z^2 + a6 * z^3
end


######################################################################
# Model changes
######################################################################


"""
Get an elliptic curve in long Weierstrass form from an elliptic curve in short Weierstrass form.

Returns an elliptic curve in long Weierstrass form with the same equation, and two maps which are the canonical isomorphisms between the curves.
"""
function tolongWeierstrass(E::ShortWeierstrassCurve)
	zero = Nemo.zero(base_ring(E))
	E2 = WeierstrassCurve(zero, zero, zero, E.a, E.b)
	phi1 = ExplicitMap(E, E2,
		function(P::EllipticPoint)
			return EllipticPoint(P.X, P.Y, P.Z, E2)
		end)
	phi2 = ExplicitMap(E2, E,
		function(P::EllipticPoint)
			return EllipticPoint(P.X, P.Y, P.Z, E)
		end)
	return E2, phi1, phi2
end


"""
Get an elliptic curve in short Weierstrass form from an elliptic curve in long Weierstrass form. This reduction is not canonical.
This assumes 2 and 3 are invertible in the base ring.

Returns an elliptic curve in short Weierstrass form, and two explicit maps giving the change of variables.
"""
function toshortWeierstrass(E::WeierstrassCurve)
	b2, _, _, _ = b_invariants(E)
	c4, c6 = c_invariants(E)
	E2 = ShortWeierstrassCurve(-27 * c4, -54 * c6)
	phi1 = ExplicitMap(E, E2,
		function(P::EllipticPoint)
			Yprime = (P.Y - E.a1 * P.X - E.a3 * P.Z) // 216
			Xprime = (P.X - 3 * b2 * P.Z) // 36
			Zprime = P.Z
			return EllipticPoint(Xprime, Yprime, Zprime, E2)
		end)
	phi2 = ExplicitMap(E2, E,
		function(P::EllipticPoint)
			Xprime = 36 * P.X + 3 * b2 * P.Z
			Yprime = 216 * P.Y + E.a3 * P.Z + E.a1 * Xprime
			Zprime = P.Z
			return EllipticPoint(Xprime, Yprime, Zprime, E)
		end)
	return E2, phi1, phi2
end

######################################################################
# Modular polynomials
######################################################################

include("modularpoly.jl")


######################################################################
# BMSS algorithm
######################################################################

include("bmss.jl")


######################################################################
# Isogenies between Weierstrass curves
######################################################################


"""
Get the polynomial associated to an l-torsion rational point, l being an odd prime.

The input is not checked.
"""
function subgrouppoly{T<:Nemo.FieldElem}(Q::EllipticPoint{T}, l::Nemo.Integer)
	normalize!(Q)
	K = base_ring(Q)
	R, x = PolynomialRing(K, "x")
	poly = Nemo.one(R)
	point = Q
	for k = 1 : ((l-1) // 2)
		poly *= (x - point.X)
		point += Q
	end
	return poly
end

"""
Build an isogeny given its domain and the polynomial defining its kernel.

This isogeny is separable and normalized.
"""
function Isogeny{T}(E::WeierstrassCurve{T}, poly::Nemo.PolyElem{T})
	a1, a2, a3, a4, a6 = a_invariants(E)
    b2, b4, b6, b8 = b_invariants(E)
    n = degree(poly)
    
    s1 = - poly[n - 1]
    s2 = poly[n - 2]
    s3 = - poly[n - 3]
    t = 6 * (s1^2 - 2 * s2) + b2 * s1 + n * b4
    w = 10 * (s1^3 - 3 * s1 * s2 + 3 * s3) + 2 * b2 * (s1^2 - 2 * s2) + 3 * b4 * s1 + n * b6
    
    E1 = WeierstrassCurve(a1, a2, a3, a4 - 5 * t, a6 - b2 * t - 7 * w)
    
    return Isogeny(E, poly, E1)
end


function Isogeny{T}(E::ShortWeierstrassCurve{T}, poly::Nemo.PolyElem{T})
	a1, a2, a3, a4, a6 = a_invariants(E)
	b2, b4, b6, b8 = b_invariants(E)
	n = degree(poly)
	
	s1 = - poly[n - 1]
	s2 = poly[n - 2]
	s3 = - poly[n - 3]
	t = 6 * (s1^2 - 2 * s2) + b2 * s1 + n * b4
	w = 10 * (s1^3 - 3 * s1 * s2 + 3 * s3) + 2 * b2 * (s1^2 - 2 * s2) + 3 * b4 * s1 + n * b6
	
	E1 = ShortWeierstrassCurve(a4 - 5 * t, a6 - b2 * t - 7 * w)
	
	return Isogeny(E, poly, E1)
end

"""
Build an isogeny given an elliptic curve in short Weierstrass form, the j-invariant of the targetted elliptic curve, and the degree.

This isogeny is separable and normalized.
"""

function Isogeny{T}(E::ShortWeierstrassCurve{T}, degree::Nemo.Integer, jprime::T)
	K = base_ring(E)
	Phi_l = modularpoly(degree)
	R, x = PolynomialRing(K, "x")
	Phi_l = R(Phi_l)
	j = j_invariant(E)
	l = K(degree)
	
	derx = derivative(Phi_l)(j)(jprime)
	dery = derivative(Phi_l)(jprime)(j) #works because Phi_l is symmetric
	
	J = - K(18) // l * (E.b // E.a) * (derx // dery) * j
	jj = jprime * (jprime - 1728)
	
	aprime = - J^2 // (48 * l^4 * jj)
	bprime = - J^3 // (864 * l^6 * jprime * jj)
	
	Eprime = ShortWeierstrassCurve(aprime, bprime)
	
	poly = kernelpoly(E, Eprime, degree)
	
	return Isogeny(E, poly, Eprime)
end


	
end # module


