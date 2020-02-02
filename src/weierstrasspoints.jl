
######################################################################
# weierstrasspoints.jl: Addition laws for projective points on Weierstrass curves
######################################################################

export infinity, projective_add, projective_scalar_mul


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

function zero(E::AbstractWeierstrass)
	return infinity(E)
end


######################################################################
# Addition law
######################################################################

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
function _addgeneric(P::EllipticPoint{T}, Q::EllipticPoint{T}) where T<:Nemo.FieldElem
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
function _addequalx(P::EllipticPoint{T}, Q::EllipticPoint{T}) where T<:Nemo.FieldElem
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
function +(P::EllipticPoint{T}, Q::EllipticPoint{T}) where T<:Nemo.FieldElem
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

function -(P::EllipticPoint{T}, Q::EllipticPoint{T}) where T<:FieldElem
	return P + (-Q)
end

"""
Get a scalar multiple of a point on a Weierstrass curve.
"""

function *(k::Integer, P::EllipticPoint)
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

# function *(k::Integer, P::EllipticPoint)
# 	return k * P
# end

function _ladder(m::Integer, P::EllipticPoint)
    p0 = P
    for b in Iterators.drop(Iterators.reverse(digits(Int8, m, base=2)), 1)
        if (b == 0)
            p0 = _double(p0)
        else
            p0 = P + _double(p0)
        end
    end
    return p0
end

######################################################################
# Projective addition law for Short Weierstrass curves
######################################################################

# P = (x1, y1, z1)
# Q = (x2, y2, z2)
# P + Q = (x3, y3, z3)
# if P != Q then
# x3 = (x2 z1 - x1 z2) [ (y2 z1 - y1 z2)^2 z1 z2 - (x2 z1 - x1 z2)^2 (x2 z1 + x1 z2) ]
# y3 = (y2 z1 - y1 z2) [ (x2 z1 - x1 z2)^2 (x2 z1 + 2 x1 z2) - (y2 z1 - y1 z2)^2 z1 z2 ] - (x2 z1 - x1 z2)^3 y1 z2
# z3 = (x2 z1 - x1 z2)^3 z1 z2
# if P = Q then
# x3 = 2 y1 z1 [ (a z1^2 + 3 x1^2)^2 - 8 x1 y1^2 z1 ]
# y3 = (a z1^2 + 3 x1^2 ) [ 12 x1 y1^2 z1 - a^2 (z1^2 + 3 x1^2)^2 ] - 8 y1^4 z1^2
# z3 = (2 y1 z1)^3

#This function is only to be used with distinct points on the same short Weierstrass curve
#xdet is x2 z1 - x1 z2, and ydet is y2 z1 - y1 z2
function _projective_add_neq(P::EllipticPoint{T}, Q::EllipticPoint{T}, xdet::T, ydet::T) where T
	x1, y1, z1 = coordinates(P)
	x2, y2, z2 = coordinates(Q)
	xdet2 = xdet^2
	xdet3 = xdet2 * xdet
	ydet2 = ydet^2
	z1z2 = z1 * z2
	#we could save a few multiplications in what follows
	x3 = xdet * ( ydet2 * z1z2 - xdet2 * (x2 * z1 + x1 * z2) )
	y3 = ydet * ( xdet2 * (x2 * z1 + 2 * x1 * z2) - ydet2 * z1z2 ) - xdet3 * y1 * z2
	z3 = xdet3 * z1z2
	res = Point(x3, y3, z3, base_curve(P))
	#@assert isvalid(res) #sanity check
	return res
end

#This function is only to be used with a point on a short Weierstrass curve
function _projective_dbl(P::EllipticPoint{T}) where T
	x, y, z = coordinates(P)
	_, _, _, a, b = a_invariants(base_curve(P))
	factor = a * z^2 + 3 * x^2
	factor2 = factor^2
	yz = y * z
	xy2z = x * y * yz
	#we could save again a few multiplications in what follows
	xprime = 2 * yz * ( factor2 - 8 * xy2z )
	yprime = factor * ( 12 * xy2z - factor2 ) - 8 * (yz)^2 * y^2
	zprime = 8 * (yz)^3
	res = Point(xprime, yprime, zprime, base_curve(P))
	#@assert isvalid(res) #sanity check
	return res
end

#This function is only to be used with two points on the same short Weierstrass curve
#Compute xdet and ydet, if they are both zero go to _projective_dbl
function projective_add(P::EllipticPoint{T}, Q::EllipticPoint{T}) where T
	x1, y1, z1 = coordinates(P)
	x2, y2, z2 = coordinates(Q)
	xdet = x2 * z1 - x1 * z2
	ydet = y2 * z1 - y1 * z2
	if ((xdet == 0) & (ydet == 0))
		return _projective_dbl(P)
	else
		return _projective_add_neq(P, Q, xdet, ydet)
	end
end

#Here we assume v is a positive integer and P lives on a short Weierstrass curve
function _projective_scalar_mul(P::EllipticPoint, v::Integer)
    P0 = P
    for b in Iterators.drop(Iterators.reverse(digits(Int8, v, base=2)), 1)
        P0 = _projective_dbl(P0)
        if (b == 1)
            P0 = projective_add(P0, P)
        end
    end
    return P0
end

function projective_scalar_mul(P::EllipticPoint, v::Integer)
	if v == 0
		return infinity(base_curve(E))
	elseif v < 0
		return - _projective_scalar_mul(P, -v)
	else
		return _projective_scalar_mul(P, v)
	end
end



