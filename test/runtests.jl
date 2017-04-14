using EllipticCurves
using Base.Test

# write your own tests here
@test 1 == 2


######################################################################
# Testing weierstrass.jl
######################################################################

"""
Example is taken from Silverman 1 p.59. Testing basic methods on projective points
"""

Q = Nemo.QQ
E = ShortWeierstrassCurve(Q(0), Q(17))
Eprime = ShortWeierstrassCurve(Q(0), Q(16))
P1 = EllipticPoint(Q(-2), Q(3), Q(1), E)
P3 = EllipticPoint(Q(2), Q(5), Q(1), E)

@test isvalid(P1)
@test isvalid(P3)

P2 = EllipticPoint(QQ(-4), QQ(6), QQ(2), E)

@test areequal(P1, P2)
@test !areequal(P1, P3)

P2 = EllipticPoint(Q(-2), Q(3), Q(1), Eprime)

@test !isvalid(P2)
@test !areequal(P1, P2)

P2 = infinity(E)

@test isvalid(P2)
@test isinfinity(P2)
@test !areequal(P1, P2)

@test base_ring(P1) = Q

"""
Testing invariants; values are given by Sage
"""

@test a_invariants(E) == (0, 0, 0, 0, 17)
@test b_invariants(E) == (0, 0, 68, 0)
@test c_invariants(E) == (0, -14688)
@test discriminant(E) == -124848
@test j_invariant(E) == 0
@test isvalid(E)

E2 = WeierstrassCurve(Q(1), Q(2), Q(3), Q(4), Q(6))

@test a_invariants(E2) == (1, 2, 3, 4, 6)
@test b_invariants(E2) == (9, 11, 33, 44)
@test c_invariants(E2) == (-183, -4293)
@test discriminant(E2) == -14212
@test j_invariant(E2) == 6128487//14212
@test isvalid(E2)



"""
Testing changes of variables
"""

E2 = WeierstrassCurve(Q(0), Q(0), Q(0), Q(0), Q(17))
P2 = EllipticPoint(Q(-2), Q(3), Q(1), E2)
E3, phi, psi = tolongWeierstrass(E)

@test E2 == E3
@test Eval(phi, P1) == P2
@test Eval(psi, P2) == P1

E3, phi, psi = toshortWeierstrass(E2)

@test E1 == E3
@test Eval(phi, P2) == P1
@test Eval(psi, P1) == P2

E2 = WeierstrassCurve(Q(1), Q(2), Q(3), Q(4), Q(6))
E3, _, _ = toshortWeierstrass(E2)

@test j_invariant(E2) == j_invariant(E3)

#=
"""
Testing isogenies between Weierstrass curves
"""
=#

######################################################################
# Testing points.jl
######################################################################

"""
Addition laws for points on Weierstrass curves. Example is from Silverman 1 p.59
"""

E = ShortWeierstrassCurve(Q(0), Q(17))
P1 = EllipticPoint(Q(-2), Q(3), Q(1), E)
P2 = EllipticPoint(Q(-1), Q(4), Q(1), E)
P3 = EllipticPoint(Q(2), Q(5), Q(1), E)
Pinf = infinity(E)

@test plus(P1, P1) == addequalx(P1, P1) == EllipticPoint(Q(8), Q(-23), Q(1), E)
@test plus(Pinf, Pinf) == Pinf
@test plus(P1, Pinf) == plus(Pinf, P1) == P1
@test plus(P2, P2) == addequalx(P2, P2) == EllipticPoint(Q(137//64), Q(-2651//512), Q(1), E)
@test plus(P2, P3) == addgeneric(P2, P3) == EllipticPoint(Q(-8//9), Q(-109//27), Q(1), E)

mP3 = minus(P3)

@test mP3 == EllipticPoint(Q(2), Q(-5), Q(1), E)
@test plus(P3, mP3) == Pinf
@test plus(P1, mP3) == EllipticPoint(Q(4), Q(9), Q(1), E)


######################################################################
# Testing bmss.jl
######################################################################




######################################################################
# Testing modular.jl
######################################################################




######################################################################
# Testing montgomery.jl
######################################################################


"""
Testing basic functions
"""

E = MontgomeryCurve(Q(7), Q(1))
P1 = EllipticPoint(Q(1), Q(3), Q(1), E)
P0 = EllipticPoint(Q(0), Q(0), Q(1), E)

@test base_ring(E) == Q
@test a_invariants(E) == (0, 7, 0, 1, 0)
@test j_invariant(E) == 24918016//45

@test isvalid(P1) & isvalid(P2)
@test isvalid(E)

"""
Testing model changes
"""

Z = Nemo.ZZ
E2, _, _ = tolongWeierstrass(E)

@test E2 == WeierstrassCurve(a_invariants(E))

"""
Testing x-only arithmetic
"""

xP1 = xonly(P1)
xP0 = xonly(P0)

@test isfixedtorsion(xP0)
@test !isfixedtorsion(xP1)

@test plus(P0, P0) == infinity(E)
@test plus(P1, P1) == P0

@test xdouble(xP0) == xinfinity(E)
@test xdouble(xP1) == xP0
@test xdouble(xinfinity(E)) == xinfinity(E)

@test xadd(xP1, xP0) == xP1

@test times(Z(0), xP1) == xinfinity(E)
@test times(Z(1), xP1) == xP1
@test times(Z(2), xP1) == xP0
@test times(Z(1), xP0) == xP0
@test times(Z(2), xP0) == xinfinity(E)
@test times(Z(45), xinfinity(E)) == xinfinity(E)
@test times(Z(4), xP1) == xinfinity(E)

@test xladder(Z(8), xP1) == xinfinity(E)

