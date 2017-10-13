
export order, roots, any_root, issquare

######################################################################
# tools.jl: useful tools for EllipticCurves
######################################################################

function (P::MPoly)(x, y, z)
	return evaluate(P, [x, y, z])
end

function (P::Frac)(args...)
	n = P.num
	d = P.den
	return n(args...) // d(args...)
end

######################################################################
# Useful functions for finite fields
######################################################################

#Multiplicative order in (small) finite fields

function order(x::FinFieldElem)
	K = parent(x)
	if x == Nemo.zero(K)
		throw(ArgumentError("Argument must be non-zero"))
	else
		k = 1
		z = x
		while z != Nemo.one(K)
			z *= x
			k += 1
		end
		return k
	end
end


#Computing roots of polynomials over finite fields

function roots{T<:FinFieldElem}(f::PolyElem{T})
	fac = Nemo.factor(f)
	res = Dict{T, Int}()
	for index in fac
		P, exp = index
		if degree(P) == 1
			res[- coeff(P, 0)] = exp
		end
	end
	res = collect(res)
	return res
end
	
function any_root{T<:FinFieldElem}(f::PolyElem{T})
	r = roots(f)
	if any(_->true, r) #tests whether r contains any element
		y0, _ = r[1]
		return (true, y0)
	else
		return (false, Nemo.zero(base_ring(f)))
	end
end

function issquare(x::FinFieldElem)
	K = parent(x)
	A, y = PolynomialRing(K, "y")
	return any_root(y^2 - x)
end

#Conversions between finite fields

"""
Convert an element of a finite field to another one.

No control is made on the input.
"""
function convert(x::FinFieldElem, K::FinField)
	K1 = parent(x)
	p1 = characteristic(K1)
	p = characteristic(K)
	(p1 == p) || throw(ArgumentError("Fields must have the same characteristic"))
	y = deepcopy(x)
	y.parent = K
	return y
end

"""
When used with polynomials, converts each coefficient.
"""
function convert{T<:FinFieldElem}(P::PolyElem{T}, K::FinField)
	A, Y = PolynomialRing(K, "Y")
	poly = Nemo.zero(A)
	for i = 0:degree(P)
		Nemo.setcoeff!(poly, i, convert(coeff(P, i), K))
	end
	return poly
end





