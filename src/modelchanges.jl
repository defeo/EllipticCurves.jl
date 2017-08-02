
######################################################################
# Model changes
######################################################################


"""
Get an elliptic curve in long Weierstrass form from an elliptic curve in short Weierstrass form.

Returns an elliptic curve in long Weierstrass form with the same equation, and two maps which are the canonical isomorphisms between the curves.
"""
function tolongWeierstrass(E::ShortWeierstrassCurve)
	zero = Nemo.zero(base_ring(E))
	E2 = WeierstrassCurve(zero, zero, zero, E.a, E.b)
	phi1 = ExplicitMap(E, E2, P::EllipticPoint -> EllipticPoint(P.X, P.Y, P.Z, E2))
	phi2 = ExplicitMap(E2, E, P::EllipticPoint -> EllipticPoint(P.X, P.Y, P.Z, E))
	return E2, phi1, phi2
end


"""
Get an elliptic curve in short Weierstrass form from an elliptic curve in long Weierstrass form. This reduction is not canonical.
This assumes 2 and 3 are invertible in the base ring.

Returns an elliptic curve in short Weierstrass form, and two explicit maps giving the change of variables.
"""
function toshortWeierstrass(E::WeierstrassCurve)
	b2, _, _, _ = b_invariants(E)
	c4, c6 = c_invariants(E)
	E2 = ShortWeierstrassCurve(-27 * c4, -54 * c6)
	phi1 = ExplicitMap(E, E2,
		function(P::EllipticPoint)
			Yprime = (P.Y - E.a1 * P.X - E.a3 * P.Z) // 216
			Xprime = (P.X - 3 * b2 * P.Z) // 36
			Zprime = P.Z
			return EllipticPoint(Xprime, Yprime, Zprime, E2)
		end)
	phi2 = ExplicitMap(E2, E,
		function(P::EllipticPoint)
			Xprime = 36 * P.X + 3 * b2 * P.Z
			Yprime = 216 * P.Y + E.a3 * P.Z + E.a1 * Xprime
			Zprime = P.Z
			return EllipticPoint(Xprime, Yprime, Zprime, E)
		end)
	return E2, phi1, phi2
end



######################################################################
# Model changes
######################################################################

"""
Given an elliptic curve in Montgomery form, build an elliptic curve in long Weierstrass form isomorphic to it.

Returns an elliptic curve and two explicit maps givint the change of variables.
"""
function tolongWeierstrass{T<:Nemo.FieldElem}(E::MontgomeryCurve{T})
    zero = Nemo.zero(Nemo.parent(E.A))
    E1 = WeierstrassCurve(zero, E.A // E.B, zero, Nemo.inv(E.B ^ 2), zero)
	phi1 = ExplicitMap(E, E1,
		function(P::EllipticPoint)
			Xprime = P.X // E.B
			Yprime = P.Y // E.B
			Zprime = P.Z
			return EllipticPoint(Xprime, Yprime, Zprime, E1)
		end)
	phi2 = ExplicitMap(E1, E,
		function(P::EllipticPoint)
			Xprime = P.X * E.B
			Yprime = P.Y * E.B
			Zprime = P.Z
			return EllipticPoint(Xprime, Yprime, Zprime, E)
		end)
	return (E1, phi1, phi2)
end