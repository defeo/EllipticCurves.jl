
######################################################################
# montgomery.jl: Montgomery models
######################################################################



######################################################################
# Basic methods
######################################################################

"""
Concrete type for (twisted) Montgomery curves.
"""
immutable MontgomeryCurve{T<:Nemo.RingElem} <: EllipticCurve{T}
    A::T
    B::T
end

function base_ring(E::MontgomeryCurve)
	return Nemo.parent(E.A)
end

"""
Get a description of a Montgomery curve.
"""
function show{T}(io::IO, E::MontgomeryCurve{T})
    print(io, "Elliptic Curve  $(E.B) y² = x³ + $(E.A) x² + x  over ")
    show(io, base_ring(E))
end

"""
Get the j-invariant of a Montgomery Curve.
"""

function j_invariant{T<:Nemo.FieldElem}(E::MontgomeryCurve{T})
	K = base_ring(E)
	zero = Nemo.zero(K)
	return j_invariant(Weierstrass(zero, E.A, zero, Nemo.one(K), zero))
end

function projective_equation(E::MontgomeryCurve)
	K = base_ring(E)
	A, (X, Y, Z) = PolynomialRing(K, ["X", "Y", "Z"])
	return E.B * Z * Y^2 - X^3 - E.A * X^2 * Z - X * Z^2 
end

function isvalid(E::MontgomeryCurve)
	return (E.B != 0) & (E.A != 0) & (E.A != 2) & (E.A != -2)
end








