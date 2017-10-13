######################################################################
# Concrete types for maps
######################################################################

"""
Concrete type for explicit (simple) maps between ellitpic curves, e.g. isomorphisms or changes of variables.

The 'map' field must contain a function sending a point on the domain curve to a point on the image curve.
"""
immutable ExplicitMap{T} <: Map{T}
	domain::EllipticCurve{T}
	image::EllipticCurve{T}
	map::Function
end

function domain(phi::ExplicitMap)
	return phi.domain
end

function codomain(phi::ExplicitMap)
	return phi.image
end

"""
Concrete type for isogenies between Weierstrass or Montgomery elliptic curves. It is assumed that the kernel may be descibed by an univariate polynomial.
"""
immutable Isogeny{T} <: Map{T}
	domain::EllipticCurve{T}
	degree::Nemo.Integer
	kernel::Nemo.PolyElem{T}
	image::EllipticCurve{T}
end

function domain(phi::Isogeny)
	return phi.domain
end

function codomain(phi::Isogeny)
	return phi.image
end

function kernel(phi::Isogeny)
	return phi.kernel
end

function degree(phi::Isogeny)
	return phi.degree
end

######################################################################
# Basic methods
######################################################################


"""
Get the evaluation of an explicit map on a point.
"""
function Eval{T}(phi::ExplicitMap{T}, P::EllipticPoint{T})
	return phi.map(P)
end

"""
Shows a description of an explicit map. Since the map itself can be anything, only the domain and image curves are given.
"""
function show(io::IO, phi::ExplicitMap)
	print(io, "Explicit map between ")
	show(io, phi.domain)
	print(io, " and ")
	show(io, phi.image)
end

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
end

