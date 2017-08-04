
######################################################################
# edwardspoints.jl: Projective and extended points for Edwards curves
######################################################################

export EdwardsPoint, ExtendedEdwardsPoint, neutral, isneutral, extended_neutral

######################################################################
# Type definitions
######################################################################

"""
Concrete type for classical projective points on Edwards curves.
"""
type EdwardsPoint{T<:Nemo.RingElem} <: ProjectivePoint{T}
    X::T
	Y::T
    Z::T
    curve::Edwards{T}
end

"""
Concrete type for extended projective points on Edwards curves.
"""
type ExtendedEdwardsPoint{T<:Nemo.RingElem} <: ProjectivePoint{T}
    X::T
	Y::T
    Z::T
	U::T
    curve::Edwards{T}
end

base_ring(P::EdwardsPoint) = parent(P.X)

base_curve(P::EdwardsPoint) = P.curve

base_ring(P::ExtendedEdwardsPoint) = parent(P.X)

base_curve(P::ExtendedEdwardsPoint) = P.curve


"""
Get a description of a classical projective point on Edwards curves.
"""
function show(io::IO, P::EdwardsPoint)
	print(io, "Point ($(P.X):&(P.Y):$(P.Z)) on ")
	show(io, P.curve)
end

"""
Get a description of an extended projective point on Edwards curves.
"""
function show(io::IO, P::ExtendedEdwardsPoint)
	print(io, "Extended point ($(P.X):&(P.Y):$(P.Z):$P.U) on ")
	show(io, P.curve)
end

######################################################################
# Conversions
######################################################################

EdwardsPoint(x, y, E::Edwards) = EdwardsPoint(x, y, one(parent(X)), E)

ExtendedEdwardsPoint(P::EdwardsPoint) = ExtendedEdwardsPoint(P.X * P.Z, P.Y * P.Z, P.Z^2, P.X * P.Y, P.curve)

EdwardsPoint(P::ExtendedEdwardsPoint) = EdwardsPoint(P.X, P.Y, P.Z, P.curve)


######################################################################
# Normalization, equality and validity
######################################################################

normalized(P::EdwardsPoint) = EdwardsPoint(P.X // P.Z, P.Y // P.Z, one(base_ring(P)), P.curve)

function isvalid(P::EdwardsPoint)
	Eq = projective_equation(P.curve)
	return (P.Z != 0) && Eq(P.X, P.Y, P.Z) == 0
end

function ==(P::EdwardsPoint, Q::EdwardsPoint)
	return (P.X == Q.X) && (P.Y == Q.Y) && (P.Z == Q.Z) && (P.curve == Q.curve)
end

function areequal(P::EdwardsPoint, Q::EdwardsPoint)
	c = P.Z // Q.Z
	return (P.X = c * Q.X) && (P.Y = c * Q.Y) && (P.curve == Q.curve)
end

function isvalid(P::ExtendedEdwardsPoint)
	Q = EdwardsPoint(P)
	return isvalid(Q) && (P.U * P.Z == P.X * P.Y)
end
	
function ==(P::ExtendedEdwardsPoint, Q::ExtendedEdwardsPoint)
	return (P.X == Q.X) && (P.Y == Q.Y) && (P.Z == Q.Z) && (P.U == Q.U) && (P.curve == Q.curve)
end

function areequal(P::ExtendedEdwardsPoint, Q::ExtendedEdwardsPoint)
	c = P.Z // Q.Z
	return (P.X = c * Q.X) && (P.Y = c * Q.Y) && (P.U == Q.U) && (P.curve == Q.curve)
end


######################################################################
# Special points
######################################################################

function neutral(E::Edwards)
	K = base_ring(E)
	return EdwardsPoint(zero(K), one(K), one(K))
end

isneutral(P::EdwardsPoint) = areequal(P, neutral(P.curve))

extended_neutral(E::Edwards) = ExtendedEdwardsPoint(neutral(E))

isneutral(E::ExtendedEdwardsPoint) = isneutral(EdwardsPoint(P))


######################################################################
# Arithmetic of ordinary projective points
######################################################################

complete_formulas(E::Edwards) = issquare(E.a) && !issquare(E.d)

#For now addition formulae *assume* that the formulas are complete.

function +(P::EdwardsPoint, Q::EdwardsPoint)
	Ed = P.curve
	A = P.Z * Q.Z
	B = A^2
	C = P.X * Q.X
	D = P.Y * Q.Y
	E = Ed.d * C * D
	F = B - E
	G = B + E
	resx = A * F * ((P.X + Q.X) * (P.Y + Q.Y) - C - D)
	resy = A * G * (D - Ed.a * C)
	resz = F * G
	return EdwardsPoint(resx, resy, resz, Ed)
end

function double(P::EdwardsPoint)
	Ed = P.curve
	B = (P.X + P.Y)^2
	C = P.X^2
	D = P.Y^2
	E = Ed.a * C
	F = E + D
	H = P.Z^2
	J = F - 2*H
	resx = (B - C - D) * J
	resy = F * (E - D)
	resz = F * J
	return EdwardsPoint(resx, resy, resz, Ed)
end

function *(k::Nemo.Integer, P::EdwardsPoint)
	E = P.curve
	P = normalized(P)
	if k == 0
		return neutral(E)
	elseif k<0
		return (-k) * (-P)
	else
		if isneutral(P)
			return neutral(E)
		else
			return ladder(k, P)
		end
	end
end

*(k::Nemo.fmpz, P::EdwardsPoint) = BigInt(k) * P

function ladder(m::Integer, P::EdwardsPoint)
	p0 = P
	for b in bin(m)[2:end]
        if (b == '0')
        	p0 = double(p0)
        else
        	p0 = P + double(p0)
        end
    end
	return p0
end

 -(P::EdwardsPoint) = EdwardsPoint( - P.X, P.Y, P.Z, P.curve)
 
 -(P::EdwardsPoint, Q::EdwardsPoint) = P + (-Q)

######################################################################
# Arithmetic of extended projective points
######################################################################

function +(P::ExtendedEdwardsPoint, Q::ExtendedEdwardsPoint)
	Ed = P.curve
	A = P.X * Q.X
	B = P.Y * Q.Y
	C = Ed.d * P.U * Q.U
	D = P.Z * Q.Z
	E = (P.X + Q.X) * (P.Y + Q.Y) - A - B
	F = D - C
	G = D + C
	H = B - Ed.a * A
	resx = E * F
	resy = G * H
	resz = E * H
	resu = F * G
	return ExtendedEdwardsPoint(resx, resy, resz, resu, Ed)
end

function double(P::ExtendedEdwardsPoint)
	Ed = P.curve
	A = P.X^2
	B = P.Y^2
	C = 2 * P.Z^2
	D = Ed.a * A
	E = (P.X + P.Y)^2 - A - B
	G = D + B
	F = G - C
	H = D - B
	resx = E * F
	resy = G * H
	resz = F * G
	resu = E * H
	return ExtendedEdwardsPoint(resx, resy, resz, resu, Ed)
end

function *(k::Integer, P::ExtendedEdwardsPoint)
	E = P.curve
	P = normalized(P)
	if k == 0
		return neutral(E)
	elseif k<0
		return (-k) * (-P)
	else
		if isneutral(P)
			return neutral(E)
		else
			return ladder(k, P)
		end
	end
end

*(k::Nemo.fmpz, P::ExtendedEdwardsPoint) = BigInt(k) * P

function ladder(m::Nemo.Integer, P::ExtendedEdwardsPoint)
	p0 = P
	for b in bin(m)[2:end]
        if (b == '0')
        	p0 = double(p0)
        else
        	p0 = P + double(p0)
        end
    end
	return p0
end

 -(P::ExtendedEdwardsPoint) = ExtendedEdwardsPoint( - P.X, P.Y, P.Z, -P.U, P.curve)
 
 -(P::ExtendedEdwardsPoint, Q::ExtendedEdwardsPoint) = P + (-Q)


