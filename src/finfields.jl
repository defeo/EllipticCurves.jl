
######################################################################
# finfields.jl: Elliptic curves over finite fields
######################################################################

export frobeniustrace, cardinality, frobeniuspolynomial, randXZ, has_montgomery, torsionpoint, torsionXZ, card_over_extension, base_extend


######################################################################
# Trace, cardinality and Frobenius polynomial
######################################################################


function frobeniustrace{T<:FinFieldElem}(E::EllipticCurve{T})
	#Call PARI
	throw(NotImplementedError("Link to PARI to compute cardinalities of elliptic curves over finite fields"))
end

function cardinality{T<:FinFieldElem}(E::EllipticCurve{T})
	t = trace(E)
	q = order(base_ring(E))
	return q + 1 - t
end

function frobeniuspolynomial{T<:FinFieldElem}(E::EllipticCurve{T}, Card)
	q = order(base_ring(E))
	t = q + 1 - Card
	A, X = PolynomialRing(Nemo.ZZ, "X")
	return X^2 - t * X + q
end


######################################################################
# Random points on elliptic curves
######################################################################

function rand{T<:FinFieldElem}(E::EllipticCurve{T})
	K = base_ring(E)
	A, Y = PolynomialRing(K, "Y")
	x = rand(K)
	poly = equation(E, x, Y)
	bool, y = any_root(poly)
	while bool == false
		x = rand(K)
		poly = equation(E, x, Y)
		(bool, y) = any_root(poly)
	end
	return Point(x, y, Nemo.one(K), E)
end

function randXZ{T<:FinFieldElem}(E::Montgomery{T})
	P = rand(E)
	return XZPoint(P)
end


######################################################################
# Torsion points
######################################################################

function torsionpoint{T<:FinFieldElem}(E::EllipticCurve{T}, l, Card)
	cofactor = div(Card, l)
	@assert cofactor * l == Card
	P = rand(E)
	Q = cofactor * P
	while isinfinity(Q)
		P = rand(E)
		Q = cofactor * P
	end
	isinfinity(l * Q) || throw(ArgumentError("Given curve has no such rational torsion points"))
	return Q
end

function torsionXZ{T<:FinFieldElem}(E::Montgomery{T}, l, Card)
	cofactor = div(Card, l)
	@assert cofactor * l == Card
	P = randXZ(E)
	Q = cofactor * P
	while isinfinity(Q)
		P = randXZ(E)
		Q = cofactor * P
	end
	isinfinity(l * Q) || throw(ArgumentError("Given curve has no such torsion rational points"))
	return Q
end


######################################################################
# Base extensions
######################################################################


function base_extend{T<:FinFieldElem}(E::AbstractWeierstrass{T}, K::FinField)
	CurveType = curvetype(E)
	K1 = base_ring(E)
	p1 = characteristic(K1)
	p = characteristic(K)
	a1, a2, a3, a4, a6 = a_invariants(E)
	((degree(K1) == 1) & (p1 == p)) || throw(ArgumentError("Invalid field extension"))
	return CurveType(convert(a1, K), convert(a2, K), convert(a3, K), convert(a4, K), convert(a6, K))
end

function base_extend{T<:FinFieldElem}(E::Montgomery{T}, K::FinField)
	K1 = base_ring(E)
	p1 = characteristic(K1)
	p = characteristic(K)
	((degree(K1) == 1) & (p1 == p)) || throw(ArgumentError("Invalid field extension"))
	return Montgomery(convert(E.A, K), convert(E.B, K))
end

function card_over_extension(Card, p, r)
	t = p + 1 - Card
	A, X = PolynomialRing(Nemo.QQ, "X")
	poly = X^2 - t * X + p
	_, z = Nemo.NumberField(poly, "z")
	poly(z) == 0 || throw(ArgumentError("NumberField does not work as planned"))
	alpha = z
	beta = t - z
	poly(beta) == 0 || throw(ArgumentError("NumberField does not work as planned"))
	return num(trace(p^r + 1 - alpha^r - beta^r) // 2)
end


######################################################################
# Existence of a Montgomery model
######################################################################

"""
Decide whether a given curve in short Weierstrass form has a Mongomery rational model.
"""
function has_montgomery{T<:FinFieldElem}(E::ShortWeierstrass{T})
	K = base_ring(E)
	A, X = PolynomialRing(K, "X")
	_, _, _, a, b = a_invariants(E)
	poly = X^3 + a * X + b
	r = roots(poly)
	n = length(r)
	for i = 1:n
		alpha, _ = r[i]
		(test, beta) = issquare(3*alpha^2 + E.a)
		test && return (true, Montgomery(3 * alpha // beta, 1 // beta))
	end
	return (false, Montgomery(K(1), K(1)))
end

	
