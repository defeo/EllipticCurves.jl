module EllipticCurves

import Nemo

export EllipticCurve

######################################################################
# Abstract types
######################################################################

"""
Abstract classes for elliptic curves.

Every elliptic curve inherits from this.
"""
abstract EllipticCurve{T<:Nemo.RingElem}

"""
Abstract classes for elliptic curves in Weierstrass form.

Every elliptic curve in Weierstrass form inherits from this.
"""
abstract AbstractWeierstrass{T<:Nemo.RingElem} <: EllipticCurve{T}


"""
Abstract classes for elliptic curve points.

Every elliptic curve point inherits from this.
"""
abstract EllipticPoint{T<:Nemo.RingElem}

"""
Abstract class for maps.

Every map between elliptic curves inherits from this.
"""
abstract Map{T<:Nemo.RingElem}



######################################################################
# Abstract functions
######################################################################

"""
Abstract function to get the base curve of an EllipticPoint.

Every EllipticPoint must implement this method.

Returns an elliptic curve.
"""
function BaseCurve{T}(P::EllipticPoint{T})
end

"""
Abstract function to get the domain curve of a Map.

Every map must implement this method.

Returns an elliptic curve.
"""
function Domain{T}(phi::Map{T})
end

"""
Abstract function to get the image curve of a Map.

Every map must implement this method.

Returns an elliptic curve.
"""
function Image{T}(phi::Map{T})
end

"""
Abstract function to evaluate maps at points.

Every map between elliptic curves must implement this method with the points of the domain curve.

Returns a point on the image curve, and throws an exception if the curves do not match.
"""
function Eval{T}(phi::Map{T}, P::EllipticPoint{T})
end

"""
Get the base ring of an elliptic curve.
"""
function basering{T}(E::EllipticCurve{T})
    return Nemo.parent_type(T)
end

function A{T}(a::Function{T,T}, b::T) a(b) end
######################################################################
# Maps
######################################################################

include("maps.jl")


######################################################################
# Modular polynomials
######################################################################

include("modularpoly.jl")


######################################################################
# Weierstrass curves
######################################################################

include("weierstrass.jl")

#=
######################################################################
# Montgomery curves
######################################################################

include("montgomery.jl")

=#






end # module
