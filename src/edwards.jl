
######################################################################
# edwards.jl: Edwards models
######################################################################

#export Edwards, isvalid

######################################################################
# Basic methods
######################################################################

"""
Concrete type for (twisted) Edwards curves.
"""
struct Edwards{T<:Nemo.RingElem} <: EllipticCurve{T}
    a::T
    d::T
end

Edwards(d) = Edwards(one(parent(d)), d)

function base_ring(E::Edwards)
	return Nemo.parent(E.a)
end

"""
Get a description of a Edwards curve.
"""
function show(io::IO, E::Edwards{T}) where T
	coef = (E.c)^2
    print(io, "Elliptic Curve in Edwards form $E.a x² + y² = (1 + $E.d x²y²)  over ")
    show(io, base_ring(E))
end

function projective_equation(E::Edwards)
	K = base_ring(E)
	A, (X, Y, Z) = PolynomialRing(K, ["X", "Y", "Z"])
	return E.a * X^2 * Z^2 + Y^2 * Z^2 - (1 + E.d * X^2 * Y^2)
end

function equation(E::Edwards)
	K = base_ring(E)
	A, (X, Y) = PolynomialRing(K, ["X", "Y"])
	return E.a * X^2 + Y^2 - (1 + E.D * X^2 * Y^2)
end

function isvalid(E::Edwards)
	return (E.a =! 0) & (E.d != 0) & (E.a != E.d)
end
