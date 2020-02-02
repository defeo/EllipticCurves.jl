using Nemo
using EllipticCurves
using ClassPolynomials

using Test

######################################################################
# Testing tools.jl
######################################################################

#This file should be ok if the tests below pass


######################################################################
# Testing weierstrass.jl
######################################################################

print("Testing weierstrass.jl... ")

E1 = Weierstrass(QQ(0), QQ(0), QQ(0), QQ(4), QQ(5))
E2 = ShortWeierstrass(QQ(0), QQ(0), QQ(0), QQ(4), QQ(5))
E3 = SeparatedWeierstrass(QQ(0), QQ(0), QQ(0), QQ(4), QQ(5))
E4 = EllipticCurve(QQ(4), QQ(5))

@test base_ring(E1) == QQ
@test base_ring(E2) == QQ
@test base_ring(E3) == QQ
@test E1 == E2 == E3 == E4

E = EllipticCurve(QQ(0), QQ(17))
@test a_invariants(E) == (QQ(0), QQ(0), QQ(0), QQ(0), QQ(17))
@test b_invariants(E) == (QQ(0), QQ(0), QQ(68), QQ(0))
@test c_invariants(E) == (QQ(0), QQ(-14688))
@test discriminant(E) == QQ(-124848)
@test j_invariant(E) == QQ(0)
@test isvalid(E)

eq = equation(E)
Y = gen(parent(eq))
X = gen(base_ring(eq))
@test eq == Y^2 - X^3 - 17

preq = projective_equation(E)
X, Y, Z = gens(parent(preq))
@test preq == Z * Y^2 - X^3 - 17 * Z^3

f3 = divisionpolynomial(E, 3)
f9 = divisionpolynomial(E, 9)
@test f9 % f3 == 0

print("done\n")


######################################################################
# Testing points.jl
######################################################################

print("Testing points.jl... ")

E = EllipticCurve(QQ(0), QQ(17))
P1 = Point(QQ(-2), QQ(3), QQ(1), E)
P2 = Point(QQ(-4), QQ(6), QQ(2), E)
P3 = Point(QQ(2), QQ(5), QQ(1), E)

@test coordinates(P1) == (-2, 3, 1)
@test base_curve(P1) == E
@test base_ring(P1) == QQ
@test isvalid(P1)
@test isvalid(P3)

@test samefields(normalized(P2), P1)
@test P1 == P2

print("done\n")


######################################################################
# Testing weierstrasspoints.jl
######################################################################

print("Testing weierstrasspoints.jl... ")

#Example is taken from Silverman 1 p.59.

E = EllipticCurve(QQ(0), QQ(17))
P1 = Point(QQ(-2), QQ(3), QQ(1), E)
P2 = Point(QQ(-1), QQ(4), QQ(1), E)
P3 = Point(QQ(2), QQ(5), QQ(1), E)
Pinf = zero(E)

@test P1 + P1 == 2 * P1 == Point(QQ(8), QQ(-23), QQ(1), E)
@test Pinf + Pinf == 2 * Pinf == Pinf
@test P1 + Pinf == Pinf + P1 == P1
@test P2 + P2 == 2 * P2 == Point(QQ(137//64), QQ(-2651//512), QQ(1), E)
@test P2 + P3 == Point(QQ(-8//9), QQ(-109//27), QQ(1), E)
@test -P3 == Point(QQ(2), QQ(-5), QQ(1), E)
@test P3 - P3 == 0 * P3 == Pinf
@test P1 - P3 == 1 * P1 + (-1) * P3 == Point(QQ(4), QQ(9), QQ(1), E)

@test projective_add(P2, P3) == P2 + P3
@test projective_scalar_mul(P2, 7) == 7 * P2

E = EllipticCurve(QQ(2), QQ(1))
P = Point(QQ(1), QQ(2), QQ(1), E)
@assert projective_add(P, P) == 2 * P
@assert projective_scalar_mul(P, 17) == 17 * P
print("done\n")


######################################################################
# Testing montgomery.jl
######################################################################

#Really trivial code, no need for testing

######################################################################
# Testing montgomerypoints.jl and curves
######################################################################

print("Testing montgomerypoints.jl... ")

E = Montgomery(QQ(7), QQ(1))
@test base_ring(E) == QQ
@test j_invariant(E) == 24918016//45
@test isvalid(E)

P1 = Point(QQ(1), QQ(3), QQ(1), E)
P0 = Point(QQ(0), QQ(0), QQ(1), E)
@test isvalid(P1) & isvalid(P2)

xP1 = XZPoint(P1)
xP0 = XZPoint(P0)
xinf = XZzero(E)

@test isfixedtorsion(xP0)
@test !isfixedtorsion(xP1)
@test fixedtorsion(E) == xP0

@test xdouble(xP0) == xinf
@test xdouble(xP1) == xP0
@test xdouble(xinf) == xinf

@test xadd(xP1, xP0, xP1) == xP1

#Scalar multiplications

@test 0 * xP1 == xinf
@test 1 * xP1 == xP1
@test 2 * xP1 == xP0
@test 1 * xP0 == xP0
@test 2 * xP0 == xinf
@test 45 * xinf == xinf
@test 4 * xP1 == xinf

@test 8 * xP1 == xinf

print("done\n")

######################################################################
# Testing edwards.jl
######################################################################

#Do we really want to keep this code ?

######################################################################
# Testing edwardspoints.jl
######################################################################

#Do we really want to keep this code ?

######################################################################
# Testing pairings.jl
######################################################################

#Do we really want to keep this code ?

######################################################################
# Testing isogenies.jl
######################################################################

print("Testing isogenies.jl... ")

F, x = FiniteField(101, 1, "x")
E = EllipticCurve(F(1), F(3))
P = Point(F(3), F(29), F(1), E)
Q = Point(F(4), F(24), F(1), E) #another point

#P is a point of three-torsion
@test isvalid(P)
@test 3 * P == zero(E)
@test 3 * Q != zero(E)

f3 = divisionpolynomial(E, 3)
@test f3(P.X) == 0

phi = Isogeny(E, P, 3)
Eprime = image(phi)
@test domain(phi) == E
@test Eprime == EllipticCurve(F(24), F(24)) #given by Sage
@test phi(P) == zero(Eprime)
@test isvalid(phi(Q))
@test phi(Q) == Point(F(91), F(20), F(1), Eprime)

K = subgroup(P, 3)

@test kernel(phi) == K
@test Isogeny(E, K) == phi
@test Isogeny(E, 3, Eprime) == phi

@test isvalid(phi)
psi = dual(phi)


@test compose(phi, psi) == multiplication_isogeny(E, 3)
@test multiplication_isogeny(E, 3)(P) == 3 * P

id = Isogeny(E, 1, E)
@test kernel(id) == 1
@test id(P) == P

@test isvalid(psi)
@test isvalid(id)

phix, phiy = rational_fractions(phi)
@test Isogeny(E, Eprime, phix, phiy) == phi

jprime = j_invariant(Eprime)
@test Isogeny(E, 3, jprime) == phi

#With modular equations

ell = 5
FX, X = PolynomialRing(F, "X")
bb, j1 = any_root(ClassicalModularPolynomial(ell, j_invariant(E), X))
@test bb
@test isvalid(Isogeny(E, ell, j1))

@test isvalid(Isogeny(E, ell, j1; USE_ATKIN = true))
print("done\n")

######################################################################
# Testing finfields.jl
######################################################################

print("Testing finfields.jl... ")

#Conversions

F, u = FiniteField(233, 1, "u")
G, z = FiniteField(233, 3, "z")
@test convert(G(1), F) == F(1)
@test convert(F(1), G) == G(1)

#Roots of polynomials

R, X = PolynomialRing(F, "X")
(bool, root) = any_root(X^2 - 1)
@test bool
@test root^2 == 1

#Cardinality

E = EllipticCurve(F(3), F(4))
#@test cardinality(E) == 260: no point counting algorithm yet
#@test frobeniustrace(E) == -26
E2 = Montgomery(F(3), F(1))
#@test cardinality(E2) == 240

#Random points

P = rand(E)
@test isvalid(P)
@test 260 * P == zero(E)
P2 = randXZ(E2)
@test isvalid(P2)
@test 240 * P2 == XZzero(E2)

#Torsion points

Q = torsionpoint(E, 13, 260)
@test Q !== zero(E)
@test 13 * Q == zero(E)
Q2 = torsionXZ(E2, 5, 240)
@test Q2 !== XZzero(E2)
@test 5 * Q2 == XZzero(E2)

#Base extensions

P = rand(base_extend(E, G))
@test isvalid(P)

@test card_over_extension(210, 233, 3) == 12652290

#Montgomery models
E = EllipticCurve(F(204), F(56))
(bool, Eprime) = has_montgomery(E)
@test bool
@test j_invariant(E) == j_invariant(Eprime)

print("done\n")

print("\n")



