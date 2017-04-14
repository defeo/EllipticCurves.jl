module Maps

import Nemo
import Base.show
import ..EllipticCurves: EllipticCurve, AbstractWeierstrass, ProjectivePoint, Map, basecurve, domain, image, Eval, Points.EllipticPoint

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

"""
Concrete type for isogenies between Weierstrass or Montgomery elliptic curves. It is assumed that the kernel may be descibed by an univariate polynomial.
"""
immutable Isogeny{T} <: Map{T}
	domain::EllipticCurve{T}
	kernel::Nemo.PolyElem{T}
	image::EllipticCurve{T}
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



end # module
