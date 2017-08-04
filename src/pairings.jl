
######################################################################
# pairings.jl: Pairings on elliptic curves
######################################################################

export




######################################################################
# Miller's algorithm
######################################################################

function tangent(P::EllipticPoint)
end

function line(P::EllipticPoint, Q::EllipticPoint)
end

function miller(P::EllipticPoint, l::Integer)
	isneutral(l * P) || throw(ArgumentError("Given point in 'miller' is not a torsion point"))
	K = base_ring(P.curve)
	A, (x, y) = PolynomialRing(K, ["X", "Y"])
	E = P.curve
	O = neutral(E)
	
	f = 1
	T = P
	
	for b in bin(l)[2:end]
		D = tangent(P)
		V = line(2 * T, 0)
		f = f^2 * D // V
		if b == '1'
			...
		end
	end
end


######################################################################
# The Weil pairing
######################################################################







######################################################################
# The Tate pairing
######################################################################


