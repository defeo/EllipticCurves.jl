
######################################################################
# isogenies.jl: Isogenies between elliptic curves
######################################################################

export Isogeny

export domain, image, kernel, degree, rational_fractions, areisomorphic, compose, squarefree_kernel

######################################################################
# Type definitions
######################################################################

immutable IsogenyWeierstrass{T}
	domain::AbstractWeierstrass{T}
	degree::Integer
	two_tors_kernel::PolyElem{T}
	sqfkernel::PolyElem{T}
	image::AbstractWeierstrass{T}
end

immutable IsogenyMontgomery{T}
	domain::Montgomery{T}
	degree::Integer
	two_tors_kernel::PolyElem{T}
	sqfkernel::PolyElem{T}
	image::Montgomery{T}
end

"""
Concrete type for isogenies between Weierstrass or Montgomery elliptic curves.
In the case of odd degree, the kernel is assumed to be squarefree.
"""
type Isogeny{T}
	domain::EllipticCurve{T}
	degree::Integer
	two_torsion_kernel::PolyElem{T}
	sqfoddkernel::PolyElem{T}
	u::T
	r::T
	s::T
	t::T
	image::EllipticCurve{T}
	rationalx::Nullable{Frac{PolyElem{T}}}
	rationaly::Nullable{Frac{PolyElem{T}}}
end

domain(phi::Isogeny) = phi.domain

image(phi::Isogeny) = phi.image

two_torsion_kernel(phi::Isogeny) = phi.two_torsion_kernel

kernel(phi::Isogeny) = phi.two_torsion_kernel * (phi.sqfoddkernel)^2

squarefree_kernel(phi::Isogeny) = phi.two_torsion_kernel * phi.sqfoddkernel

squarefree_odd_kernel(phi::Isogeny) = phi.sqfoddkernel

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
		phi.rational_x = Nullable(fx)
		phi.rational_y = Nullable(fy)
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
	R1 = coordinate_ring(E1)
	P2 = equation(E2)
	X1, Y1 = gens(R1)
	phix, phiy = rational_fractions(phi)
	return P2(phix(X1, Y1), phiy(X1, Y1)) == 0
end
	

######################################################################
# Rational fractions and evaluation
######################################################################


"""
Get the evaluation of an isogeny on a point.
"""
function (phi::Isogeny)(P::EllipticPoint)
	Q = normalize(P)
	x1 = Q.X
	y1 = Q.Y
	phix, phiy = rational_fractions(phi)
	try
		x2 = phix(x1, y1)
		y2 = phiy(x1, y1)
		return EllipticPoint(x2, y2, one(base_ring(Q)), image(phi))
	catch
		return infinity(image(phi))
	end
end

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
	
	if l==1
		A, (x, y) = PolynomialRing(K, ["x", "y"])
		Ffield = FractionField(A)
		return (Ffield(x), Ffield(y))
	
	elseif isodd(l)
		s1 = - coeff(kpoly, n-1)
		n = div(l - 1, 2)
		kpoly = squarefree_kernel(phi)
		kpolysquare = kernel(phi)
		dkpoly = derivative(kpoly)
		d2kpoly = derivative(dkpoly)
		
		psi2 = parent(kpoly)(divisionpolynomial(E, -1))
		X = gen(parent(kpoly))
		phi = psi2 * (dkpoly^2 - d2kpoly * kpoly) - (6 * X^2 + b2 * X + b4) * dkpoly * kpoly + (l * X - 2 * s1) * kpolysquare
		
		if K(2) != 0
			dphi = derivative(phi)
			psi2 = divisionpolynomial(E, 2, 1) #bivariate polynomial
			A = parent(phi)
			x, y = gens(A)
			
			omega = (1 // K(2)) * (dphi * kpoly)(x) * psi2 - (phi * dkpoly)(x) * psi2 + (1 // K(2)) * ((a1 * phi + a3 * kpolysquare) * kpoly)(x)
			
			return (phi(x) // kpolysquare(x) , omega // (kpoly^3)(x))
			
		else #K has characteristic 2
			throw(NotImplementedError("Rational fractions not yet implemented in characteristic 2"))
		end
	
	else #l is even : we decompose "odd" and 2-torsion parts
		two_kernel = two_torsion_kernel(phi)
		oddkernel = squarefree_odd_kernel(phi)
		
		if degree(oddkernel) > 0
			(phi2, phiodd) = decompose_two_torsion(phi)
			phi2x, phi2y = rational_fractions(phi2)
			phioddx, phioddy = rational_fractions(phiodd)
			
			return (phioddx(phi2x, phi2y), phioddy(phi2x, phi2y))
			
		else #kernel of phi is contained in the two-torsion subgroup
			if degree(two_kernel) == 1
				x0 = - coeff(two_kernel, 0)
				if K(2) != 0
					y0 = (- a1 * x0 + a3) // 2
				else #K has characteristic 2
					throw(NotImplementedError("Rational fractions not yet implemented in characteristic 2"))
				end
				
				c = 3 * x0^2 + 2 * a2 * x0 + a4 - a1 * y0
				A, (x, y) = PolynomialRing(K, ["x", "y"])
				Ffield = FractionField(A)
				return (x + Ffield(c)//(x - x0) , y - c * (a1 * (x - x0) + (y - y0)) // (x - x0)^2)
			
			else #two_kernel is the full 2-torsion, and car(K) != 2 since isogenies are assumed to be separable here
				s1 = -coeff(two_kernel, 2)
				dtwo_kernel = derivative(two_kernel)
				d2two_kernel = derivative(dtwo_kernel)
				phi = dtwo_kernel^2 - 2 * d2two_kernel * two_kernel + (4 * gen(parent(two_kernel)) - s1) * two_kernel^2
				
				psi2 = divisionpolynomial(E, 2, 1)
				A = parent(psi2)
				x, y = gens(A)
				
				omega = (1 // K(2)) * psi2 * (derivative(phi) * two_kernel - phi * dtwo_kernel)(x) - (1 // K(2)) * (a1 * phi + a3 * two_kernel)(x)
				return (phi(x) // (two_kernel^2)(x) , omega // (two_kernel^3)(x))
			
			end
		end
	end	
end

function decompose_two_torsion(phi::Isogeny)
	throw(NotImplementedError)
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
	
	resx = phi2x(phi1x, phi1y)
	resy = phi2y(phi2x, phi2y)
	
	return Isogeny(E1, E3, resx, resy)
end


######################################################################
# Defining isogenies with various inputs
######################################################################

"""
Build an isogeny given its domain, kernel, degree, and image
"""
function Isogeny{T}(E::AbstractWeierstrass{T}, degree::Integer, poly::PolyElem{T}, image)
	
end

function Isogeny{T}(E1::AbstractWeierstrass{T}, E2::AbstractWeierstrass{T}, phix::Frac{PolyElem{T}}, phiy::Frac{PolyElem{T}})
	
end

"""
Build an isogeny given its domain and its kernel polynomial.
"""
function Isogeny{T}(E::AbstractWeierstrass{T}, poly::PolyElem{T})
	
	CurveType = curvetype(E)
	psi2 = divisionpolynomial(E, 2, 2)(gen(parent(poly)))
	two_torsion_part = gcd(psi2, poly)
	odd_part = divexact(poly, two_torsion_part)
	
	K = base_ring(E)
	a1, a2, a3, a4, a6 = a_invariants(E)
	b2, b4, b6, b8 = b_invariants(E)
	n = Nemo.degree(poly)
	
	coeff(poly, n) == 1 || throw(ArgumentError("kernel polynomial must be monic"))
	
	if two_torsion_part == 1
	
		s1, s2, s3 = zero(K), zero(K), zero(K)
	
		n > 0 && (s1 = - coeff(poly, n - 1))
		n > 1 && (s2 = coeff(poly, n - 2))
		n > 2 && (s3 = - coeff(poly, n - 3))
		t = 6 * (s1^2 - 2 * s2) + b2 * s1 + n * b4
		w = 10 * (s1^3 - 3 * s1 * s2 + 3 * s3) + 2 * b2 * (s1^2 - 2 * s2) + 3 * b4 * s1 + n * b6
		E1 = CurveType(a1, a2, a3, a4 - 5 * t, a6 - b2 * t - 7 * w)
	
		return Isogeny(E, 2 * n + 1, one(parent(poly)), poly, K(1), K(0), K(0), K(0), E1, Nullable{Frac{PolyElem{T}}}(), Nullable{Frac{PolyElem{T}}}())
	
	else #there is a two-torsion part
		throw(NotImplementedError)
	end
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
function _subgrouppoly{T<:Nemo.FieldElem}(Q::EllipticPoint{T}, l::Nemo.Integer)
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

function Isogeny{T<:Nemo.FieldElem}(E::Weierstrass{T}, Q::EllipticPoint{T}, l::Nemo.Integer)
	poly = _subgrouppoly(Q, l)
	return Isogeny(E, poly)
end

"""
Build an isogeny of given degree between two curves, assuming a normalized separable isogeny of this degree exists.
"""

function Isogeny{T<:Nemo.FieldElem}(E1::AbstractWeierstrass{T}, E2::AbstractWeierstrass{T},
	degree::Nemo.Integer)
	poly = kernelpoly(E1, E2, degree)
	return Isogeny(E1, poly)
end

"""
Compute the scalar multiplication by m on a weierstrass curve as a normalized isogeny.
"""
function multiplication_isogeny(E::AbstractWeierstrass, m::Nemo.Integer)
	return Isogeny(E, divisionpolynomial(E, m))
end






