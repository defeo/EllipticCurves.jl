using EllipticCurves
using Base.Test

import Nemo

import ..EllipticCurves: Weierstrass.ShortWeierstrassCurve, Weierstrass.WeierstrassCurve, Weierstrass.a_invariants, Weierstrass.b_invariants, Weierstrass.c_invariants, Weierstrass.j_invariant, Weierstrass.discriminant, Weierstrass.infinity, Points.isinfinity, Weierstrass.+, Weierstrass.-

# write your own tests here
@test 1 != 2


"""
Testing Weierstrass equations and invariants
"""

function TestShortW(a,b)
	E = ShortWeierstrassCurve(a,b)
	show(E)
	print("\n")
	print("a's : ",a_invariants(E),"\n")
	print("b's : ",b_invariants(E),"\n")
	print("c's : ",c_invariants(E),"\n")
	print("Delta : ",discriminant(E),"\n")
	print("j : ",j_invariant(E))
end

function TestShortW(a1,a2,a3,a4,a6)
	E = WeierstrassCurve(a1,a2,a3,a4,a6)
	show(E)
	print("\n")
	print("a's : ",a_invariants(E),"\n")
	print("b's : ",b_invariants(E),"\n")
	print("c's : ",c_invariants(E),"\n")
	print("Delta : ",discriminant(E),"\n")
	print("j : ",j_invariant(E))
end

"""
Testing projective points on Weierstrass curves
"""

E = ShortWeierstrassCurve(Nemo.QQ(1), Nemo.QQ(2))
P = infinity(E)

@test isinfinity(P)






