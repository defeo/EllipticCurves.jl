
######################################################################
# weierstrass.jl : curves given by weierstrass equations
######################################################################


export Weierstrass, ShortWeierstrass, SeparatedWeierstrass

export curvetype, a_invariants, b_invariants, c_invariants, discriminant, j_invariant, equation, projective_equation, divisionpolynomial


######################################################################
# Type definitions
######################################################################

"""
Concrete type for elliptic curves in long Weierstrass form over a ring.
"""
immutable Weierstrass{T} <: AbstractWeierstrass{T}
    a1::T
    a2::T
    a3::T
    a4::T
    a6::T
end


"""
Concrete types for elliptic curves in short Weierstrass form over a ring.
"""
immutable ShortWeierstrass{T} <: AbstractWeierstrass{T}
    a::T
    b::T
end


immutable SeparatedWeierstrass{T} <: AbstractWeierstrass{T}
	a2::T
	a4::T
	a6::T
end

curvetype(E::SeparatedWeierstrass) = SeparatedWeierstrass

curvetype(E::ShortWeierstrass) = ShortWeierstrass

curvetype(E::Weierstrass) = Weierstrass

######################################################################
# Calling with six arguments
######################################################################

function ShortWeierstrass{T}(a1::T, a2::T, a3::T, a4::T, a6::T)
	a1 == 0 && a2 == 0 && a3 == 0 || throw(ArgumentError("Tried to call ShortWeierstrass with nonzero coefficients"))
	return ShortWeierstrass(a4, a6)
end

function SeparatedWeierstrass{T}(a1::T, a2::T, a3::T, a4::T, a6::T)
	a1 == 0 && a3 == 0 || throw(ArgumentError("Tried to call SeparatedWeierstrass with nonzero coefficients"))
	return SeparatedWeierstrass(a2, a4, a6)
end

######################################################################
# Basic methods
######################################################################



base_ring(E::Weierstrass) = Nemo.parent(E.a1)

base_ring(E::ShortWeierstrass) = Nemo.parent(E.a)

base_ring(E::SeparatedWeierstrass) = Nemo.parent(E.a2)

"""
Get a description of an elliptic curve in long Weierstrass form.
"""
function show(io::IO, E::Weierstrass)
    print(io, "Elliptic Curve in long Weierstrass form y² + $(E.a1) xy + $(E.a3) y = x³ + $(E.a2) x² + $(E.a4) x + $(E.a6)  over ")
    show(io, base_ring(E))
end

"""
Get a description of an elliptic curve in short Weierstrass form.
"""
function show(io::IO, E::ShortWeierstrass)
    print(io, "Elliptic Curve in short Weierstrass form y² = x³ + $(E.a) x + $(E.b)  over ")
    show(io, base_ring(E))
end

"""
Get a description of an elliptic curve in short Weierstrass form.
"""
function show(io::IO, E::SeparatedWeierstrass)
    print(io, "Elliptic Curve in separated Weierstrass form y² = x³ $(E.a2) + $(E.a4) x + $(E.a6)  over ")
    show(io, base_ring(E))
end

######################################################################
# Invariants and equations
######################################################################

"""
Get the a-invariants (a1, a2, a3, a4, a6) of an elliptic curve in long Weierstrass form.
"""
a_invariants(E::Weierstrass) = (E.a1, E.a2, E.a3, E.a4, E.a6)

"""
Get the a-invariants (a1, a2, a3, a4, a6) of an elliptic curve in short Weierstrass form.
"""
function a_invariants(E::ShortWeierstrass)
    zero = Nemo.zero(base_ring(E))
    return (zero, zero, zero, E.a, E.b)
end

"""
Get the a-invariants (a1, a2, a3, a4, a6) of an elliptic curve in separated Weierstrass form.
"""
function a_invariants(E::SeparatedWeierstrass)
    zero = Nemo.zero(base_ring(E))
    return (zero, E.a2, zero, E.a4, E.a6)
end

==(E1::AbstractWeierstrass, E2::AbstractWeierstrass) = a_invariants(E1) == a_invariants(E2)



"""
Get the b-invariants (b2, b4, b6, b8) of an elliptic curve in Weierstrass form.
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
"""
function c_invariants(E::AbstractWeierstrass)
    b2, b4, b6, b8 = b_invariants(E)
    return (b2^2 - 24*b4, -b2^3 + 36*b2*b4 - 216*b6)
end


"""
Get the discriminant of an elliptic curve in Weierstrass form.
"""
function discriminant(E::AbstractWeierstrass)
    b2, b4, b6, b8 = b_invariants(E)
    return -b2^2*b8 - 8*b4^3 - 27*b6^2 + 9*b2*b4*b6
end

"""
Check if a given instance of AbstractWeierstrass is indeed a nonsingular curve.
"""
isvalid(E::AbstractWeierstrass) = !iszero(discriminant(E))



"""
Get the j-invariant of an elliptic curve in Weierstrass form.
This may fail if the base ring is not a field.
"""
function j_invariant(E::AbstractWeierstrass)
    c4, _ = c_invariants(E)
    disc = discriminant(E)
	return c4^3 // disc
end

function equation(E::AbstractWeierstrass)
	K = base_ring(E)
	R, (X, Y) = PolynomialRing(K, ["X", "Y"])
	a1, a2, a3, a4, a6 = a_invariants(E)
	return Y^2 + a1 * X * Y + a3 * Y - X^3 - a2 * X^2 - a4 * X - a6
end

function projective_equation(E::AbstractWeierstrass)
	K = base_ring(E)
	R, (X, Y, Z) = PolynomialRing(K, ["X", "Y", "Z"])
	a1, a2, a3, a4, a6 = a_invariants(E)
	return Z * Y^2 + a1 * X * Y * Z + a3 * Y * Z^2 - X^3 - a2 * X^2 * Z - a4 * X * Z^2 - a6 * Z^3
end


######################################################################
# Division polynomials for Weierstrass curves
######################################################################

"""
Internal function to compute the 2-division polynomial of an elliptic curve in Weierstrass form.
"""
_2divpoly{T}(E::AbstractWeierstrass{T}, x::PolyElem{T}, D::Dict{Int, PolyElem{T}}) = _divpoly(E, -1, x, D)


"""
Internal function to compute the mth division polynomial of an elliptic curve in Weierstrass form, using the well-known induction formulas. If the index is even, *the usual factor 2*y is dropped*.
"""
function _divpoly{T}(E::AbstractWeierstrass{T}, m::Integer, x::PolyElem{T}, D::Dict{Int, PolyElem{T}})
	b2, b4, b6, b8 = b_invariants(E)
	R = parent(x)
	if m in keys(D)
		return D[m]
	elseif m == -2
		D[-2] = _divpoly(E, - 1, x, D)^2
		return D[-2]
	elseif m == -1
		D[-1] = 4 * x^3 + b2 * x^2 + 2 * b4  *x + b6
		return D[-1]
	elseif m <= 0
		throw(ArgumentError("Non-positive integers in _divpoly must be either -1 or -2, got $m"))
	elseif m == 1
		D[1] = R(1)
		return D[1]
	elseif m == 2
		D[2] = R(1)
		return D[2]
	elseif m == 3
		D[3] = 3 * x^4 + b2 * x^3 + 3 * b4 * x^2 + 3 * b6 * x + b8
		return D[3]
	elseif m == 4
		pm2 = _divpoly(E, -2, x, D)
		p3 = _divpoly(E, 3, x, D)
		D[4] = - pm2 + (6*x^2 + b2*x + b4) * p3
		return D[4]
	elseif isodd(m)
		k = div(m - 1, 2)
		m1 = _divpoly(E, k - 1, x, D)
		eq = _divpoly(E, k, x, D)
		p1 = _divpoly(E, k + 1, x, D)
		p2 = _divpoly(E, k + 2, x, D)
		if iseven(k)
			D[m] = _divpoly(E, -2, x, D) * p2 * eq^3 - m1 * p1^3
			return D[m]
		else #k is odd
			D[m] = p2 * eq^3 - _divpoly(E, -2, x, D) * m1 * p1^3
			return D[m]
		end
	else  #m is even
		k = div(m, 2)
		m2 = _divpoly(E, k - 2, x, D)
		m1 = _divpoly(E, k - 1, x, D)
		eq = _divpoly(E, k, x, D)
		p1 = _divpoly(E, k + 1, x, D)
		p2 = _divpoly(E, k + 2, x, D)
		D[m] = eq * (p2 * m1^2 - m2 * p1^2)
		return D[m]
	end
end

"""
Compute the division polynomial of an elliptic curve in Weierstrass form.

If m is even, you may specify an additional parameter two_torsion:
* 0 simply ignores the points of order two
* 1 returns the exact bivariate division polynomial
* 2 (default) adds extra roots at points of two-torsion in order to get a univariate polynomial.
"""

function divisionpolynomial{T}(E::AbstractWeierstrass{T}, m::Nemo.Integer, two_torsion::Integer = 2)
	m <= -3 && throw(ArgumentError("m must be positive or -1 or -2"))
	
	if m<=0
		K = base_ring(E)
		R, x = Nemo.PolynomialRing(K, "x")
		D = Dict{Int, Nemo.PolyElem{T}}()
		poly = _divpoly(E, m, x, D)
		return poly
	elseif m%2 == 0
		if two_torsion == 1
			K = base_ring(E)
			a1, _, a3, _, _ = a_invariants(E)
			R, (x, y) = Nemo.PolynomialRing(K, ["x", "y"])
			D = Dict{Int, Nemo.PolyElem{T}}()
			poly = _divpoly(E, m, x, D)
			return (2 * y + a1 * x + a3) * poly
		elseif two_torsion == 2
			K = base_ring(E)
			R, x = Nemo.PolynomialRing(K, "x")
			D = Dict{Int, Nemo.PolyElem{T}}()
			poly = _divpoly(E, m, x, D)
			pm1 = _divpoly(E, -1, x, D)
			return pm1 * poly
		else #two_torsion == 0
			K = base_ring(E)
			R, x = Nemo.PolynomialRing(K, "x")
			D = Dict{Int, Nemo.PolyElem{T}}()
			poly = _divpoly(E, m, x, D)
			return poly
		end
	else # m is odd
		K = base_ring(E)
		R, x = Nemo.PolynomialRing(K, "x")
		D = Dict{Int, Nemo.PolyElem{T}}()
		poly = _divpoly(E, m, x, D)
		return poly
	end
end




