
######################################################################
# isogenies.jl: Isogenies between elliptic curves
######################################################################

#Only odd-degree isogenies are implemented; create a special class for degree 2 isogenies

export Isogeny

export domain, image, kernel, degree, rational_fractions, areisomorphic, compose, squarefree_kernel, isomorphism_coefficients, subgroup, dual, isomorphism, find_isomorphism, multiplication_isogeny

######################################################################
# Type definitions
######################################################################

"""
Concrete type for isogenies between Weierstrass or Montgomery elliptic curves.
Only odd degree isogenies are currently implemented.
The kernel polynomial of an l-isogeny is a monic, squarefree polynomial of degree (l-1)/2.
In the rational fractions, the usual factor y in the ordinate is dropped in order to use only univariate polynomials.
"""
type Isogeny{T}
	domain::EllipticCurve{T}
	degree::Integer
	kernel::PolyElem{T}
	u::T
	r::T
	s::T
	t::T
	image::EllipticCurve{T}
	x::Nullable{GenFrac{U}} where U <: PolyElem{T}
	y::Nullable{GenFrac{U}} where U <: PolyElem{T}
end

base_ring(phi::Isogeny) = base_ring(domain(phi))

domain(phi::Isogeny) = phi.domain

image(phi::Isogeny) = phi.image

kernel(phi::Isogeny) = phi.kernel

degree(phi::Isogeny) = phi.degree

"""
Get the u, r, s, t coefficients of the isomophism coming after the normalized part of the isogeny
"""
isomorphism_coefficients(phi::Isogeny) = (phi.u, phi.r, phi.s, phi.t)

"""
Get the rational fractions defining the isogeny. They are computed at the first call.
"""
function rational_fractions(phi::Isogeny)
	try
		phix = get(phi.rationalx)
		phiy = get(phi.rationaly)
		return (phix, phiy)
	catch
		fx, fy = _compute_rational_fractions(phi)
		phi.x = Nullable(fx)
		phi.y = Nullable(fy)
		return fx, fy
	end
end

"""
Decide if two isogenies differ by postcomposition by an isomorphism.
"""
areisomorphic(phi::Isogeny, psi::Isogeny) = domain(phi) == domain(psi) && degree(phi) == degree(psi) && kernel(phi) == kernel(psi)


"""
Decide if two isogenies are the same.
"""
==(phi::Isogeny, psi::Isogeny) =  areisomorphic(phi, psi) && isomorphism_coefficients(phi) == isomorphism_coefficients(psi)

"""
Shows a description of an isogeny between elliptic curves.
"""
function show(io::IO, phi::Isogeny)
	print("Isogeny of degree ")
	show(phi.degree)
	print(" between ")
	show(io, phi.domain)
	print(" and ")
	show(io, phi.image)
	print(" with kernel ")
	show(io, phi.kernel)
end

"""
Decide if an isogeny has valid formulas.
"""
function isvalid(phi::Isogeny)
	E1 = domain(phi)
	E2 = image(phi)
	a11, a21, a31, a41, a61 = a_invariants(E1) #we have a1 = a3 = 0
	a12, a22, a32, a42, a62 = a_invariants(E2) #we have a1 = a3 = 0
	phix, phiy = rational_fractions(phi)
	X = gen(base_ring(parent(phix)))
	lhs = phiy^2 * (X^3 + a21 * X^2 + a41 * X + a61)
	rhs = (phix^3 + a22 * phix^2 + a42 * phix + a62)
	test = lhs - rhs
	return test==0
end
	

######################################################################
# Rational fractions and evaluation
######################################################################


"""
Get the evaluation of an isogeny on a point.
"""
function (phi::Isogeny)(P::EllipticPoint)
	Q = normalized(P)
	x1 = Q.X
	y1 = Q.Y
	phix, phiy = rational_fractions(phi)
	try
		x2 = phix(x1)
		y2 = y1 * phiy(x1)
		return EllipticPoint(x2, y2, one(base_ring(Q)), image(phi))
	catch
		return infinity(image(phi))
	end
end


#DOES NOT WORK if u or t != 0...
"""
Compute the rational functions defining an isogeny
"""
function _compute_rational_fractions(phi::Isogeny)
	phixn, phiyn = _normalized_rational_fractions(phi)
	u, r, s, t = isomorphism_coefficients(phi)
	return (u^2 * phixn + r, u^3 * phiyn + s * u^2 * phixn + t)
end

function _normalized_rational_fractions(phi::Isogeny)
	l = degree(phi)
	E = domain(phi)
	K = base_ring(E)
	b2, b4, b6, b8 = b_invariants(E)
	a1, a2, a3, a4, a6 = a_invariants(E)
	
	if l==1 #phi is an isomorphism
		A, x = PolynomialRing(K, "x")
		Ffield = FractionField(A)
		return (Ffield(x), Ffield(1))
	
	else #remember l is odd
		
		n = div(l - 1, 2)
		kpoly = kernel(phi)
		s1 = - coeff(kpoly, n-1)
		kpolysquare = kpoly^2
		dkpoly = derivative(kpoly)
		d2kpoly = derivative(dkpoly)
		
		FX = parent(kpoly)
		psi2 = FX(divisionpolynomial(E, -1))
		X = gen(FX)
		N = psi2 * (dkpoly^2 - d2kpoly * kpoly) - (6 * X^2 + b2 * X + b4) * dkpoly * kpoly + (l * X - 2 * s1) * kpolysquare
		
		if K(2) != 0
			dN = derivative(N)
			omega = dN * kpoly - 2 * N * dkpoly    #+ (1 // K(2)) * ((a1 * N + a3 * kpolysquare) * kpoly)(x): we force to have a1 = a3 = 0
			return (N // kpolysquare , omega // (kpoly * kpolysquare))
			
		else #K has characteristic 2
			throw(NotImplementedError("Rational fractions not yet implemented in characteristic 2"))
		end
	end	
end

######################################################################
# Composition
######################################################################

function compose(phi1::Isogeny, phi2::Isogeny)
	E1 = domain(phi1)
	E2 = image(phi1)
	E3 = image(phi2)
	@assert E2 == domain(phi2)
	
	phi1x, phi1y = rational_fractions(phi1)
	phi2x, phi2y = rational_fractions(phi2)
	
	resx = phi2x(phi1x)
	resy = phi1y * phi2y(phi1x)
	
	return Isogeny(E1, E3, resx, resy)
end


######################################################################
# Defining isogenies with various inputs
######################################################################

"""
Build an isogeny given its domain, kernel, degree, and image
"""
function Isogeny{T}(E::AbstractWeierstrass{T}, poly::PolyElem{T}, Eprime::AbstractWeierstrass{T})
	phi1 = Isogeny(E, poly)
	postisomorphism = find_isomorphism(image(phi1), Eprime)
	phi = compose(phi1, postisomorphism)
	return phi
end

function Isogeny{T, U<:PolyElem{T}}(E1::AbstractWeierstrass{T}, E2::AbstractWeierstrass{T}, phix::GenFrac{U}, phiy::GenFrac{U})
	K = sqrt(denominator(phix))
	F = base_ring(K)
	normalizedphix, normalizedphiy = rational_fractions(Isogeny(E1, K))
	usquare = phix//normalizedphix
	ucube = phiy//normalizedphiy
	
	#convert to base field
	num = numerator(usquare)
	den = denominator(usquare)
	@assert den == 1
	@assert degree(num) == 0
	usquare = coeff(num, 0)
	num = numerator(ucube)
	den = denominator(ucube)
	@assert den == 1
	@assert degree(num) == 0
	ucube = coeff(num, 0)
	
	u = ucube // usquare
	testphiy = u^3 * normalizedphiy
	if testphiy == - phiy
		u = -u
	end
	@assert degree(numerator(phix)) == 2 * degree(K) + 1
	res = Isogeny(E1, 2 * degree(K) + 1, K, u, F(0), F(0), F(0), E2, Nullable(phix), Nullable(phiy))
	@assert isvalid(res)
	return res
end

"""
Build an isogeny given its domain and its kernel polynomial.
"""
function Isogeny{T}(E::AbstractWeierstrass{T}, poly::PolyElem{T})
	
	CurveType = curvetype(E)
	
	K = base_ring(E)
	a1, a2, a3, a4, a6 = a_invariants(E)
	b2, b4, b6, b8 = b_invariants(E)
	n = Nemo.degree(poly)
	
	coeff(poly, n) == 1 || throw(ArgumentError("kernel polynomial must be monic"))
	
	s1, s2, s3 = zero(K), zero(K), zero(K)

	n > 0 && (s1 = - coeff(poly, n - 1))
	n > 1 && (s2 = coeff(poly, n - 2))
	n > 2 && (s3 = - coeff(poly, n - 3))
	t = 6 * (s1^2 - 2 * s2) + b2 * s1 + n * b4
	w = 10 * (s1^3 - 3 * s1 * s2 + 3 * s3) + 2 * b2 * (s1^2 - 2 * s2) + 3 * b4 * s1 + n * b6
	E1 = CurveType(a1, a2, a3, a4 - 5 * t, a6 - b2 * t - 7 * w)

	return Isogeny(E, 2 * n + 1, poly, K(1), K(0), K(0), K(0), E1, Nullable{GenFrac{PolyElem{T}}}(), Nullable{GenFrac{PolyElem{T}}}())
end


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
function subgroup{T<:Nemo.FieldElem}(Q::EllipticPoint{T}, l::Nemo.Integer)
	normalize!(Q)
	K = base_ring(Q)
	R, x = PolynomialRing(K, "x")
	poly = Nemo.one(R)
	point = Q
	for k = 1 : div(l-1, 2)
		poly *= (x - point.X)
		point += Q
	end
	return poly
end



"""
Build an isogeny given an elliptic curve in short Weierstrass form, the j-invariant of the targetted elliptic curve, and the degree, assuming it exists.

Returns a normalized isogeny
"""

function Isogeny{T}(E::ShortWeierstrass{T}, degree::Nemo.Integer, jprime::T)
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
	
	Eprime = ShortWeierstrass(aprime, bprime)
	
	poly = kernelpoly(E, Eprime, degree)
	
	return Isogeny(E, degree, poly, Eprime)
end


"""
Build an isogeny given an odd integer l and a rational torsion point of this order. The input is not checked.
"""

function Isogeny{T<:Nemo.FieldElem}(E::AbstractWeierstrass{T}, Q::EllipticPoint{T}, l::Nemo.Integer)
	poly = subgroup(Q, l)
	return Isogeny(E, poly)
end

"""
Build an isogeny of given degree between two curves, assuming a normalized separable isogeny of this degree exists.
"""

function Isogeny{T<:Nemo.FieldElem}(E1::AbstractWeierstrass{T}, degree::Nemo.Integer, E2::AbstractWeierstrass{T})
	poly = kernelpoly(E1, E2, degree)
	return Isogeny(E1, poly)
end

"""
Compute the scalar multiplication by m on a weierstrass curve as a normalized isogeny.
"""
function multiplication_isogeny(E::AbstractWeierstrass, m::Nemo.Integer)
	@assert isodd(m)
	K = base_ring(E)
	psim = divisionpolynomial(E, m)
	n = degree(psim)
	phi1 = Isogeny(E, 1 // coeff(psim, n) * psim)
	E2 = image(phi1)
	postisomorphism = isomorphism(E2, 1 // K(m))
	@assert image(postisomorphism) == E
	return compose(phi1, postisomorphism)
end



######################################################################
# Dual isogeny
######################################################################

function dual(phi::Isogeny)
	K = base_ring(phi)
	E1 = domain(phi)
	E2 = image(phi)
	d = degree(phi)
	E1_0 = image(isomorphism(E1, K(d)))
	phi_hat = Isogeny(E2, d, E1_0) #use BMSS
	#de-normalize !
	postisomorphism = isomorphism(image(phi_hat), 1//K(d))
	phi_hat = compose(phi_hat, postisomorphism)
	@assert isvalid(phi_hat)
	@assert image(phi_hat) == E1
	return phi_hat
end

function isomorphism{T}(E::AbstractWeierstrass{T}, u::T)
	a1, a2, a3, a4, a6 = a_invariants(E) #we have a1 = a3 = 0
	a2prime = a2 * u^2
	a4prime = a4 * u^4
	a6prime = a6 * u^6
	K = base_ring(E)
	R, X = PolynomialRing(base_ring(E), "X")
	phi = Isogeny(E, 1, R(1), u, K(0), K(0), K(0), EllipticCurve(a2prime, a4prime, a6prime), Nullable{GenFrac{PolyElem{T}}}(), Nullable{GenFrac{PolyElem{T}}}())
	return phi
end

function find_isomorphism{T}(E1::AbstractWeierstrass{T}, E2::AbstractWeierstrass{T})
	a11, a21, a31, a41, a61 = a_invariants(E1)
	a12, a22, a32, a42, a62 = a_invariants(E2)
	u = sqrt(a22 // a21)
	@assert a42 = a41 * u^4
	@assert a62 = a61 * u^6
	return isomophism(E1, u)
end

######################################################################
# Operations
######################################################################

function +(phi::Isogeny, psi::Isogeny)
	E = domain(phi)
	@assert E == domain(psi)
	Eprime = image(phi)
	@assert Eprime == image(psi)
	phix, phiy = rational_fractions(phi)
	psix, psiy = rational_fractions(phi)
	return Isogeny(E, Eprime, phix + psix, phiy + psiy)
end

function -(phi::Isogeny)
	pix, phiy = rational_fractions(phi)
	return Isogeny(domain(phi), image(phi), phix, -phiy)
end

function -(phi::Isogeny, psi::Isogeny)
	return phi + (-psi)
end







