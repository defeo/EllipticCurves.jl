
######################################################################
# montgomery.jl: Montgomery models
######################################################################

export Montgomery

######################################################################
# Basic methods
######################################################################

"""
Concrete type for (twisted) Montgomery curves.
"""
struct Montgomery{T<:Nemo.RingElem} <: EllipticCurve{T}
    A::T
    B::T
end

function base_ring(E::Montgomery)
	return Nemo.parent(E.A)
end

"""
Get a description of a Montgomery curve.
"""
function show(io::IO, E::Montgomery{T}) where T
    print(io, "Elliptic Curve  $(E.B) y² = x³ + $(E.A) x² + x  over ")
    show(io, base_ring(E))
end

"""
Get the j-invariant of a Montgomery Curve.
"""

function j_invariant(E::Montgomery{T}) where T<:Nemo.FieldElem
	K = base_ring(E)
	zero = Nemo.zero(K)
	return j_invariant(Weierstrass(zero, E.A, zero, Nemo.one(K), zero))
end

function projective_equation(E::Montgomery)
	K = base_ring(E)
	A, (X, Y, Z) = PolynomialRing(K, ["X", "Y", "Z"])
	return E.B * Z * Y^2 - X^3 - E.A * X^2 * Z - X * Z^2 
end

function equation(E::Montgomery)
	K = base_ring(E)
	A, (X, Y) = PolynomialRing(K, ["X", "Y"])
	return E.B * Y^2 - X^3 - E.A * X^2 - X
end

function equation(E::Montgomery, X, Y)
	K = base_ring(E)
	return E.B * Y^2 - X^3 - E.A * X^2 - X
end

function isvalid(E::Montgomery)
	return (E.B != 0) & (E.A != 0) & (E.A != 2) & (E.A != -2)
end



