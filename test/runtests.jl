using EllipticCurves
using Base.Test

# write your own tests here
@test 1 == 2


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
	E = ShortWeierstrassCurve(a1,a2,a3,a4,a6)
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

E = ShortWeierstrassCurve(QQ(1), QQ(2))
P = infinity(E)

@test is_identity(P)

@test -P == P






