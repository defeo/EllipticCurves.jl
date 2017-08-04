
######################################################################
# weierstrasspoints.jl: Addition laws for projective points on Weierstrass curves
######################################################################

export infinity


######################################################################
# Basic methods
######################################################################

"""
Get the point at infinity on an elliptic curve in Weierstrass form.
"""
function infinity(E::AbstractWeierstrass)
    R = base_ring(E)
    return EllipticPoint(Nemo.zero(R), Nemo.one(R), Nemo.zero(R), E)
end


"""
Get the opposite of a point on an elliptic curve in Weierstrass form.
"""
function -(P::EllipticPoint)
	E = P.curve
	x, y, z = coordinates(P)
	a1, _, a3, _, _ = a_invariants(E)
    return EllipticPoint(x, - y - a1 * x - a3 * z, z, E)
end

"""
Get the sum of two normalized projective points on the same Weierstrass curve, assuming they are not equal and not inverse of each other.
"""
function _addgeneric{T<:Nemo.FieldElem}(P::EllipticPoint{T}, Q::EllipticPoint{T})
    E = P.curve
	a1, a2, a3, _, _ = a_invariants(E)
	
	xq, yq, zq = coordinates(Q)
	@assert zq == 1
	xp, yp, zp = coordinates(P)
	@assert zp == 1
	
	#sanity check
	@assert xp != xq
	
    denom = xq - xp
	inverse = 1//denom
    lambda = (yq - yp) * inverse
    nu = (yp * xq - xp * yq) * inverse
    
    Xplus = lambda^2 + a1 * lambda - a2 - xp - xq
    Yplus = -(lambda + a1) * Xplus - nu - a3
    Zplus = Nemo.one(base_ring(P))
    return EllipticPoint(Xplus, Yplus, Zplus, E)
end

"""
Get the sum of two normalized projective points on the same Weierstrass curve, assuming they have equal x-coordinate.
"""
function _addequalx{T<:Nemo.FieldElem}(P::EllipticPoint{T}, Q::EllipticPoint{T})
    E = P.curve
	a1, a2, a3, a4, a6 = a_invariants(E)
	
	xp, yp, zp = coordinates(P)
	@assert zp == 1
	xq, yq, zq = coordinates(Q)
	@assert zq == 1
	
	#sanity check
	@assert xp == xq
	
    denom = yp + yq + a1 * xq + a3
    if iszero(denom)
	    return infinity(E)
    else
		inverse = 1//denom
	    lambda = (3 * xp^2 + 2 * a2 * xp + a4 - a1 * yp) * inverse
	    nu = (- xp^3 + a4 * xp + 2 * a6 - a3 * yp) * inverse
	    Xplus = lambda^2 + a1 * lambda - a2 - xp - xq
        Yplus = -(lambda + a1) * Xplus - nu - a3
        Zplus = one(base_ring(P))
	    return EllipticPoint(Xplus, Yplus, Zplus, E)
    end
end

"""
Get the double of a normalized point.
"""
function _double(P::EllipticPoint)
	if isinfinity(P)
		return infinity(P.curve)
	else
		return _addequalx(P, P)
	end
end


"""
Get the sum of two projective points on the same Weierstrass curve.

"""
function +{T<:Nemo.FieldElem}(P::EllipticPoint{T}, Q::EllipticPoint{T})
    P = normalized(P)
    Q = normalized(Q)
	xp, _, _ = coordinates(P)
	xq, _, _ = coordinates(Q)
    if isinfinity(P)
        return Q
    elseif isinfinity(Q)
        return P
    elseif xp == xq
		return _addequalx(P,Q)
    else
		return _addgeneric(P,Q)
    end
end

function -{T<:FieldElem}(P::EllipticPoint{T}, Q::EllipticPoint{T})
	return P + (-Q)
end

"""
Get a scalar multiple of a point on a Weierstrass curve.
"""

function *(k::Nemo.Integer, P::EllipticPoint)
	E = P.curve
	P = normalized(P)
	if k == 0
		return infinity(E)
	elseif k<0
		return (-k) * (-P)
	else
		if isinfinity(P)
			return infinity(E)
		else
			return _ladder(k, P)
		end
	end
end

function *(k::Nemo.fmpz, P::EllipticPoint)
	return BigInt(k) * P
end

function _ladder(m::Nemo.Integer, P::EllipticPoint)
	p0 = P
	for b in bin(m)[2:end]
        if (b == '0')
        	p0 = _double(p0)
        else
        	p0 = P + _double(p0)
        end
    end
	return p0
end

