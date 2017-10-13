
######################################################################
# Computing with polynomials
######################################################################


"""
Decide whether a given monic polynomial is a square, and if so compute a square root.

This test can fail in positive characteristic.
"""

function issquare(P::Nemo.PolyElem)
	d = degree(P)
	if isodd(d)
		return (false, P)
	else
		dP = Nemo.derivative(P)
		sqrt = Nemo.gcd(dP, P)
		if sqrt^2 == P
			return (true, sqrt)
		else
			return (false, sqrt)
		end
	end
end


######################################################################
# Computing with power series
######################################################################

"""
Get the polynomial associated with a power series.
"""
function convert{T}(R::PolyRing{T}, F::AbsSeriesElem{T})
	d = precision(F)
	pol = zero(R)
	for i = 0: d-1
		setcoeff!(pol, i, coeff(F, i))
	end
	return pol
end


"""
Get the power series associated with a polynomial.
"""
function convert{T}(S::SeriesRing{T}, P::PolyElem{T})
	d = min(degree(P), Nemo.max_precision(S)-1)
	series = zero(S)
	for i = 0:d
		setcoeff!(series, i, coeff(P, i))
	end
	return series
end


"""
Transform the coefficients of an absolute power series.
"""
function map(F::AbsSeriesElem, f::Function)
	R = parent(F)
	G = zero(R)
	d = precision(F)
	set_prec!(G, d)
	for i = 0: d-1
		setcoeff!(G, i, f(i, coeff(F, i)))
	end
	return G	
end


"""
Compute the derivative of an absolute capped power series.
"""
function derivative(F::AbsSeriesElem)
	return shift_right(map(F, (i, c)-> i * c), 1)
end


"""
Compute the integral of an absolute capped power series, with constant term zero.

This does not change the maximum precision of the parent power series ring, however.
"""
function integral{T<:FieldElem}(F::AbsSeriesElem{T})
	return shift_left(map(F, (i, c)-> c // (i+1)), 1)
end


"""
Compose two absolute capped power series.
"""
function compose{T}(F::AbsSeriesElem{T}, G::AbsSeriesElem{T})
	R = parent(F)
	K = base_ring(R)
	A, _ = PolynomialRing(K, "x")
	Fpol = convert(A, F)
	Gpol = convert(A, G)
	comp = compose(Fpol, Gpol)
	return convert(R, comp)
end


######################################################################
# Berlekamp-Massey
######################################################################


"""
Berlekamp-Massey algorithm to compute the rational reconstruction of a given polynomial of degree 2M - 1 (or 2M - 2) mod X^2M.
"""

function berlekamp_massey{T<:FieldElem}(a::PolyElem{T})
    M = 1 + div(degree(a), 2)
    A = parent(a)
    x = gen(A)
	
	#Extended euclidean algorithm : we stop when the remainder has degree < M
    fjm2 = a
    fjm1 = x^(2*M)
	sjm1 = 0
    sjm2 = 1
    j = 0
    while degree(fjm1) >= M
        j += 1
        q, f = divrem(fjm2, fjm1)
        s = sjm2 - q * sjm1
		fjm1, fjm2 = f, fjm1
		sjm1, sjm2 = s, sjm1
	end

	t = sjm1
	u = fjm1
    den = inv(Nemo.lead(t)) * t  # make monic
	num = inv(Nemo.lead(u)) * u
    return num, den
end

######################################################################
# BMSS algorithm
######################################################################



"""
Compute the rational fraction of the normalized rational isogeny between elliptic curves of the form y^2 = f(x).
``E1`` and ``E2`` of *odd* degree ``deg``.

Assuming a rational normalized separable isogeny of degree ``deg`` exists between
``E1`` and
``E2``, kernelpoly returns the squarefree polynomial vanishing on the
abscissae of its kernel. unsafe_kernelpoly returns its square if the isogeny has odd degree, with some correction made on two-torsion points otherwise.

This algorithm works only if the characteristic is zero or greater than 4*deg - 1.
"""

function unsafe_kernelpoly{T<:FieldElem}(E1::AbstractWeierstrass{T}, E2::AbstractWeierstrass{T},
	deg::Integer)

	K = base_ring(E1)
	#We hope no one will use this with function fields...
	if isa(K, Nemo.FinField)
		p = Nemo.characteristic(K)
	else
		p = 0
	end

	((p > 0) & (p <= 4*deg-1)) &&
		throw(ArgumentError("BMSS algorithm only works for characteristic 0 or greater than 4*deg - 1."))

    (a1, a2, a3, a4, a6) = a_invariants(E1)
    (b1, b2, b3, b4, b6) = a_invariants(E2)

    ((a1 != 0) | (a3 != 0) | (b1 != 0) | (b3 != 0)) &&
        throw(ArgumentError("Curves must have a model of the form y^2 = f(x)."))

	#We need precision 2*deg + 1
	R, x = PowerSeriesRing(K, 2*deg +1, "x", model=:capped_absolute)
    G = a6 * x^3 + a4 * x^2 + a2 * x + one(R)
    H = b6 * x^3 + b4 * x^2 + b2 * x + one(R)

    # solve the differential equation G(x) T'^2 = (T/x) H(T)
    sol = _BMSS_diffeq(G, H, 2*deg + 1)

    # We recover the rational fraction using the relation T = D(1/x) / N(1/x)
    U = shift_right(sol, 1)   #why ?
	A, _ = PolynomialRing(K, "x")
	Upol = convert(A, U)
    num, den = berlekamp_massey(Upol)
	N = reverse(den) #maybe not the right numerator
	D = reverse(num)

    # If the points of abscissa 0 are in the kernel, correct the degree of D
    gap = degree(N) - degree(D) - 1
    (gap > 0) && (D = shift_left(D, gap))

	#Making denominator monic
	lambda = inv(Nemo.lead(D))
	D = lambda * D
	N = lambda * N
    return (N, D)
end

"""
    Compute a power-series solution to the differential equation used by
    the BMSS algorithm. The differential equation is

    .. math::  P(T,x) = (T/x) H(T) - G(x) {T'}^2 = 0

    with initial conditions `T = O(x)`. The output is truncated to
    precision ``prec``.

    :param G: a power series with non-zero constant coefficient.
    :param H: a power series with non-zero constant coefficient.
    :param prec: an integer denoting the truncation order of the solution.

    :returns: the solution to the differential equation.
    :rtype: Power series

    :raises ZeroDivisionError: when the characteristic is smaller than
        ``2*prec-2`` and not 0.
    :raises ValueError: when ``G`` or ``H`` have constant term \neq 1.

    ALGORITHM:

    It is a Newton iteration. After some substitution we get

    .. math::

       T_0 = ax + O(x^2)\\
       T_{i+1} = T_i + T_i' \sqrt{G} \sqrt{x} \int \frac{k_i(x)}{2\sqrt{x}}
       + O(x^{2^i+1}),

    where

    .. math::

       k_i(x) = \frac{P(T_i, x)} {{T_i'}^2 G \sqrt{G}}.
"""

function _BMSS_diffeq{T<:FieldElem}(G::AbsSeriesElem{T}, H::AbsSeriesElem{T}, prec::Integer)
    
    # The power series ring
    R = parent(G)
	one = Nemo.one(R)
	x = gen(R)
	K = base_ring(R)
	#We hope no one will use this with function fields...
	if isa(K, Nemo.FinField)
    	p = Nemo.characteristic(base_ring(R))
	else
		p = 0
	end

	#Checking the arguments are sane
    ((coeff(G,0) != 1) | (coeff(H,0) != 1)) &&
        throw(ArgumentError("Arguments must have constant coefficient one."))

    (0 < p & p <= 2*prec-3) &&
        throw(DivideError("Characteristic must be greater than 2*prec - 3 in order to compute a solution to precision 'prec'."))

    # the precision to which the solution T is known
    d = 1
    # 1/(T'^2 G) to precision d ; at the beginning we set T = x * H/G
    diffT2G = R(1) + Nemo.O(x)
    # Sqrt(G) and Sqrt(1/G) to precision d
    sqG = R(1) + Nemo.O(x)
	invsqG = R(1) + Nemo.O(x)

    sol = x + Nemo.O(x^2)

    while d < prec-1
        # update diffT, sqG and invsqG to precision d
        # (nothing changes in the first iteration)
        diffT2G = diffT2G * (2 - G * derivative(sol)^2 * diffT2G)
        sqG = inv(K(2)) * (sqG + G * invsqG * (2 - sqG * invsqG))
        invsqG = invsqG * (2 - sqG * invsqG)

        # double the current precision
        d = min(2*d, prec-1)
        set_prec!(diffT2G, d)
        set_prec!(sqG, d)
        set_prec!(invsqG, d)
        set_prec!(sol, d+1)

        f = (compose(H, sol) * shift_right(sol, 1) - G * derivative(sol)^2) * diffT2G * invsqG
        # F = 1/2 sqrt(x) integral(k/sqrt(x))
        F = shift_left( map(f, (i, c) -> c // (2*i+1)), 1)

        # update the solution
        sol += derivative(sol) * sqG * F
	end
    return sol
end


"""
Sanity checks
"""

function kernelpoly{T<:FieldElem}(E1::AbstractWeierstrass{T}, E2::AbstractWeierstrass{T},
	deg::Integer)
	N, D = unsafe_kernelpoly(E1, E2, deg)

	evenpart = Nemo.gcd(divisionpolynomial(E1, 2), D)
    oddpart = Nemo.divexact(D, evenpart)
    check, sqrtoddpart = issquare(oddpart)
    if !check
        throw(ArgumentError("The two curves are not linked by a rational normalized isogeny of this degree"))
	end
    return evenpart * sqrtoddpart
end

