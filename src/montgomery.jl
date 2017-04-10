module MontgomeryCurves

import Nemo
import Base.show
import ..EllipticCurves: EllipticCurve, EllipticPoint
import ..WeierstrassCurves: WeierstrassCurve


immutable MontgomeryCurve{T<:Nemo.FieldElem} <: EllipticCurve{T}
    A::T
    B::T
end

function show{T}(io::IO, E::MontgomeryCurve{T})
    print(io, "Elliptic Curve  $(E.B) y² = x³ + $(E.A) x² + x  over ")
    show(io, Nemo.parent_type(T))
end

function weierstrass_form(E::MontgomeryCurve)
    zero = Nemo.zero(Nemo.parent(E.A))
    return WeierstrassCurve(zero, E.A // E.B, zero, Nemo.inv(E.B ^ 2), zero)
end


######################################################################
# Montgomery curve points
######################################################################


type MontgomeryPoint{T, E<:MontgomeryCurve} <: EllipticPoint{T}
    X::T
    Z::T
    curve::E
end

show(io::IO, P::MontgomeryPoint) = print(io, "($(P.X):$(P.Z))")

isidentity(P::MontgomeryPoint) = Nemo.iszero(P.Z)

isfixedtorsion(P::MontgomeryPoint) = Nemo.iszero(P.X)

function normalized{T<:Nemo.FieldElem, E}(P::MontgomeryPoint{T,E})
    K = Nemo.parent(P.X)
    if isidentity(P)
        return MontgomeryPoint(Nemo.zero(K), 
                               Nemo.zero(K), 
                               P.curve)
    else
        return MontgomeryPoint(P.X // P.Z,
                               Nemo.one(K),
                               P.curve)
    end
end

function normalize!{T<:Nemo.FieldElem, E}(P::MontgomeryPoint{T,E})
    K = Nemo.parent(P.X)
    if isidentity(P)
        P.X = Nemo.zero(K)
    else
        P.X = P.X // P.Z
        P.Z = Nemo.one(K)
    end
    return
end

function identity(E::MontgomeryCurve)
    K = Nemo.parent(E.A)
    zero = Nemo.zero(K)
    return MontgomeryPoint(zero, zero, E)
end

function fixedtorsion(E::MontgomeryCurve)
    K = Nemo.parent(E.A)
    return MontgomeryPoint(Nemo.zero(K), Nemo.one(K), E)
end


######################################################################
# Montgomery arithmetic
######################################################################

#taken from Ben Smith's review article

function xdouble{T,E}(P::MontgomeryPoint{T,E})
    v1 = P.X + P.Z
    v1 = v1^2
    v2 = P.X - P.Z
    v2 = v2^2
    
    X2 = v1 * v2
    
    v1 = v1 - v2
    v3 = ((E.A + 2) // 4) * v1
    v3 = v3 + v2
    
    Z2 = v1 * v3
    
    if Z2 == 0
        X2 = 0
    end
    
    return MontgomeryPoint(X2, Z2, E)
end

function xadd{T,E}(P::MontgomeryPoint{T,E}, Q::MontgomeryPoint{T,E}, Minus::MontgomeryPoint{T,E})
    v0 = P.X + P.Z
    v1 = Q.X - Q.Z
    v1 = v1 * v0
    
    v0 = P.X - P.Z
    v2 = Q.X + Q.Z
    v2 = v2 * v0
    
    v3 = v1 + v2
    v3 = v3^2
    v4 = v1 - v2
    v4 = v4^2
    
    Xplus = Minus.Z * v3
    Zplus = Minus.X * v4
   
    return MontgomeryPoint(Xplus, Zplus, P.curve) #valid if Minus is not (0:0) or (1:0)
end

"""
function xladder{T}(k::Nemo.Integer, P::MontgomeryPoint{T})
   normalize!(P)
    x0, x1 = P, xdouble(P)
    for #b running down k's bits, except the highest-weight one
        if b == 0
           x0, x1 = xdouble(x0), xadd(x0, x1, P)
        else
            x0, x1 = xadd(x0, x1, P), xdouble(x1)
        end
    end
    return x0
end
"""

function *{T,E}(k::Nemo.Integer, P::MontgomeryPoint{T,E})
    if isidentity(P)
        return P
    elseif isfixedtorsion(P)
        if k % 2 == 0
            return identity(E)
        else
            return fixedtorsion(E)
	end
    else
        return xladder(k, P)
    end
end



end


































