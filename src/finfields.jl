

export has_montgomery, isgoodprime, first_isogeny_x, first_isogeny, rand

######################################################################
# Elliptic curves over finite fields
######################################################################


#Trace, cardinality and Frobenius polynomial

function trace{T<:FinFieldElem}(E::EllipticCurve{T})
	#Call PARI
end

function cardinality{T<:FinFieldElem}(E::EllipticCurve{T})
	t = trace(E)
	q = order(base_ring(E))
	return q + 1 - t
end

function frobeniuspolynomial{T<:FinFieldElem}(E::EllipticCurve{T}, Card::Nemo.fmpz)
	q = order(base_ring(E))
	t = q + 1 - Card
	A, X = PolynomialRing(Nemo.ZZ, "X")
	return X^2 - t * X + q
end

#rand points on elliptic curves

function rand(E::AbstractWeierstrass)
	a1, a2, a3, a4, a6 = a_invariants(E)
	K = base_ring(E)
	A, Y = PolynomialRing(K, "Y")
	x = rand(K)
	poly = Y^2 + a1 * x * Y + a3 * Y - x^3 - a2 * x^2 - a4 * x - a6
	bool, y = any_root(poly)
	while bool == false
		x = rand(K)
		poly = Y^2 + a1 * x * Y + a3 * Y - x^3 - a2 * x^2 - a4 * x - a6
		(bool, y) = any_root(poly)
	end
	return EllipticPoint(x, y, Nemo.one(K), E)
end

function rand{T<:FinFieldElem}(E::Montgomery{T})
	K = base_ring(E)
	A, Y = PolynomialRing(K, "Y")
	x = rand(K)
	poly = E.B * Y^2 - x^3 - E.A * x^2 - x
	bool, y = any_root(poly)
	while bool == false
		x = rand(K)
		poly = E.B * Y^2 - x^3 - E.A * x^2 - x
		bool, y = any_root(poly)
	end
	return EllipticPoint(x, y, Nemo.one(K), E)
end

function randxonly{T<:FinFieldElem}(E::Montgomery{T})
	K = base_ring(E)
	A, Y = PolynomialRing(K, "Y")
	x = rand(K)
	poly = E.B * Y^2 - x^3 - E.A * x^2 - x
	bool, y = any_root(poly)
	while bool == false
		x = rand(K)
		poly = E.B * Y^2 - x^3 - E.A * x^2 - x
		bool, y = any_root(poly)
	end
	return XonlyPoint(x, Nemo.one(K), E)
end

#Cardinality over extensions

function card_ext(Card::Nemo.fmpz, p::Nemo.fmpz, r::Int)
	t = p + 1 - Card
	A, X = PolynomialRing(Nemo.QQ, "X")
	poly = X^2 - t * X + p
	_, z = Nemo.NumberField(poly, "z")
	poly(z) == 0 || throw(ArgumentError("NumberField does not work as planned"))
	alpha = z
	beta = t - z
	poly(beta) == 0 || throw(ArgumentError("NumberField does not work as planned"))
	# print(alpha * beta, "\n", z * t - z^2, "\n")
	# return num(coeff(p^r + 1 - alpha^r - beta^r, 0))
	# print(p^r + 1 - alpha^r - beta^r,"\n")
	return num(trace(p^r + 1 - alpha^r - beta^r) // 2)
end

#Montgomery models

"""
Decide whether a given curve in short Weierstrass form has a Mongomery rational model, and return it if so.
"""
function has_montgomery{T<:FinFieldElem}(E::ShortWeierstrass{T})
	K = base_ring(E)
	A, X = PolynomialRing(K, "X")
	poly = X^3 + E.a * X + E.b
	r = roots(poly)
	n = length(r)
	for i = 1:n
		alpha, _ = r[i]
		(test, beta) = issquare(3*alpha^2 + E.a)
		test && return (true, Montgomery(3 * alpha // beta, 1 // beta))
	end
	return (false, Montgomery(K(1), K(1)))
end

function has_montgomery{T<:FinFieldElem}(E::Weierstrass{T})
	E2, _, _ = toshortWeierstrass(E)
	return has_montgomery(E2)
end

#Vélu's formulae for Montgomery curves

function Isogeny{T<:FinFieldElem}(E::Montgomery{T}, poly::PolyElem{T})
	K = base_ring(E)
	E2 = Weierstrass(zero(K), one(K), zero(K), E.A, zero(K)) #tolongWeierstrass(E)
	X = Nemo.gen(parent(poly))
	poly2 = poly(E.B * X)
	E3 = image(Isogeny(E2, poly2))
	bool, E4 = has_montgomery(E3)
	bool || throw(ArgumentError("Target curve has no rational Montgomery model"))
	return Isogeny(E, 2 * degree(poly) + 1, poly, E4)
end

#Base extensions

"""
Given an elliptic curve over a prime finite field, compute the curve obtained after base extension.
"""
function base_extend{T<:FinFieldElem}(E::EllipticCurve{T}, K::FinField)
end

function base_extend{T<:FinFieldElem}(E::ShortWeierstrass{T}, K::FinField)
	K1 = base_ring(E)
	p1 = characteristic(K1)
	p = characteristic(K)
	((degree(K1) == 1) & (p1 == p)) || throw(ArgumentError("Invalid field extension"))
	return ShortWeierstrass(convert(E.a, K), convert(E.b, K))
end

function base_extend{T<:FinFieldElem}(E::Weierstrass{T}, K::FinField)
	K1 = base_ring(E)
	p1 = characteristic(K1)
	p = characteristic(K)
	((degree(K1) == 1) & (p1 == p)) || throw(ArgumentError("Invalid field extension"))
	return Weierstrass(convert(E.a1, K), convert(E.a2, K), convert(E.a3, K), convert(E.a4, K), convert(E.a6, K))
end

function base_extend{T<:FinFieldElem}(E::Montgomery{T}, K::FinField)
	K1 = base_ring(E)
	p1 = characteristic(K1)
	p = characteristic(K)
	((degree(K1) == 1) & (p1 == p)) || throw(ArgumentError("Invalid field extension"))
	return Montgomery(convert(E.A, K), convert(E.B, K))
end

######################################################################
# More specific to the isogeny problem
######################################################################

#Computing torsion points


function torsionpoint{T<:FinFieldElem}(E::EllipticCurve{T}, l::Int, Card::Nemo.fmpz)
	cofactor = Nemo.divexact(Card, l)
	P = rand(E)
	Q = cofactor * P
	while isinfinity(Q)
		P = rand(E)
		Q = cofactor * P
	end
	isinfinity(l * Q) || throw(ArgumentError("Given curve has no such torsion rational points"))
	return Q
end

function torsionxonly{T<:FinFieldElem}(E::Montgomery{T}, l::Int, Card::Nemo.fmpz)
	cofactor = Nemo.divexact(Card, l)
	P = randxonly(E)
	Q = times(cofactor, P)
	while isinfinity(Q)
		P = randxonly(E)
		Q = cofactor * P
	end
	isinfinity(times(l, Q)) || throw(ArgumentError("Given curve has no such torsion rational points"))
	return Q
end

#Computing good primes for a given curve : a prime ell is good if the frobenius eigenvalues over E[ell](kbar) have different multiplicative orders

function jacobi(D::Nemo.fmpz, l::Int)
	r = rem(D, l)
	if r<0
		return jacobi(Nemo.ZZ(r+l), Nemo.ZZ(l))
	else
		return jacobi(Nemo.ZZ(r), Nemo.ZZ(l))
	end
end

function isgoodprime{T<:FinFieldElem}(E::EllipticCurve{T}, l::Int, Card::Nemo.fmpz)
	poly = frobeniuspolynomial(E, Card)
	discr = Nemo.discriminant(poly)
	if !isprime(ZZ(l)) | jacobi(discr, l) != 1
		return (false, 0)
	else
		Fl, _ = Nemo.FiniteField(l, 1, "alpha")
		A, _ = Nemo.PolynomialRing(Fl, "X")
		polymod = A(poly)
		r = roots(polymod)
		r1 = r[1][1]
		r2 = r[2][1]
		o1 = order(r1)
		o2 = order(r2)
		return (o1 != o2, min(o1, o2))
	end
end

function goodprimes{T<:FinFieldElem}(E::EllipticCurve{T}, Card::Nemo.fmpz, Bound::Int)
	res = Dict{Int, Pair{Bool, Int}}()
	for l in 3:Bound
		bool, ord = isgoodprime(l)
		bool && (res[l] = ord)
	end
	return res
end

#Computing the rational isogeny corresponding to a good prime

function subgrouppolynomial{T<:FinFieldElem}(E::EllipticCurve{T}, l::Int, Q::EllipticPoint{T})
	K = base_ring(E)
	A, Y = PolynomialRing(K, "Y")
	poly = Nemo.one(A)
	P = Q
	for i = 1 : div(l-1, 2)
		poly *= (Y - P.X)
		P = P + Q
	end
	return poly
end

function first_isogeny{T<:FinFieldElem}(E::AbstractWeierstrass{T}, l::Int, Cards::Dict{Int, Nemo.fmpz})
	Card = Cards[1]
	(bool, order) = isgoodprime(E, l, Card)
	bool || throw(ArgumentError("Given degree is not a good prime for the chosen curve"))
	K = base_ring(E)
	p = characteristic(K)
	(degree(K) == 1) || throw(ArgumentError("Base field must be prime"))
	
	if order == 1
		Q = torsionpoint(E, l, Card)
		poly = subgrouppolynomial(E, l, Q)
	else #Use an extended curve
		haskey(Cards, order) || throw(ArgumentError("Cardinal for this extension degree was not computed"))
		Cprime = Cards[order]
		Kprime, _ = FiniteField(p, order, "alpha")
		Eprime = base_extend(E, Kprime)
		Q = torsionpoint(Eprime, l, Cprime)
		poly = subgrouppolynomial(Eprime, l, Q)
		poly = convert(poly, K)
	end
	return Isogeny(E, poly)
end

#With Montgomery curves

function subgrouppolynomial{T<:FinFieldElem}(E::Montgomery{T}, l::Int, Q::XonlyPoint{T})
	K = base_ring(E)
	A, Y = PolynomialRing(K, "Y")
	poly = Nemo.one(A)
	P = xinfinity(E)
	for i = 1 : div(l-1, 2)
		P = times(i, Q)
		poly *= (Y - P.X // P.Z)
	end
	return poly
end

function first_isogeny_x{T<:FinFieldElem}(E::Montgomery{T}, l::Int, Cards::Dict{Int, Nemo.fmpz})
	Card = Cards[1]
	(bool, order) = isgoodprime(E, l, Card)
	bool || throw(ArgumentError("Given degree is not a good prime for the chosen curve"))
	K = base_ring(E)
	p = characteristic(K)
	(degree(K) == 1) || throw(ArgumentError("Base field must be prime"))
	
	if order == 1
		Q = torsionxonly(E, l, Card)
		poly = subgrouppolynomial(E, l, Q)
	else #Use an extended curve
		haskey(Cards, order) || throw(ArgumentError("Cardinal for this extension degree was not computed"))
		Cprime = Cards[order]
		Kprime, _ = FiniteField(p, order, "alpha")
		Eprime = base_extend(E, Kprime)
		Q = torsionxonly(Eprime, l, Cprime)
		poly = subgrouppolynomial(Eprime, l, Q)
		poly = convert(poly, K)
	end
	return Isogeny(E, poly)
end



	
