module EllipticCurves

import Nemo

export EllipticCurve, AbstractWeierstrass, ProjectivePoint, Map, ProjectivePoint, WeierstrassCurve, ShortWeierstrassCurve, MontgomeryCurve, XonlyPoint, ExplicitMap, Isogeny, basecurve, domain, image, Eval, base_ring, normalize!, normalized, a_invariants, b_invariants, c_invariants, j_invariant, discriminant, infinity, Eval, xonly, isinfinity, isfixedtorsion, fixedtorsion

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
abstract ProjectivePoint{T<:Nemo.RingElem}

"""
Abstract class for maps.

Every map between elliptic curves inherits from this.
"""
abstract Map{T<:Nemo.RingElem}



######################################################################
# Abstract functions
######################################################################

"""
Abstract function to get the base curve of an ProjectivePoint.

Every ProjectivePoint must implement this method.

Returns an elliptic curve.
"""
function basecurve{T}(P::ProjectivePoint{T})
end

"""
Abstract function to get the domain curve of a Map.

Every map must implement this method.

Returns an elliptic curve.
"""
function domain{T}(phi::Map{T})
end

"""
Abstract function to get the image curve of a Map.

Every map must implement this method.

Returns an elliptic curve.
"""
function image{T}(phi::Map{T})
end

"""
Abstract function to evaluate maps at points.

Every map between elliptic curves must implement this method with the points of the domain curve.

Returns a point on the image curve, and throws an exception if the curves do not match.
"""
function Eval{T}(phi::Map{T}, P::ProjectivePoint{T})
end

"""
Abstract function to get the base ring of an elliptic curve.

Every elliptic curve must implement this method.
"""
function base_ring{T}(E::EllipticCurve{T})
end



######################################################################
# Points on elliptic curves
######################################################################

include("points.jl")


######################################################################
# Maps
######################################################################

include("maps.jl")





######################################################################
# Weierstrass curves
######################################################################

include("weierstrass.jl")


######################################################################
# Montgomery curves
######################################################################

include("montgomery.jl")




end # module
