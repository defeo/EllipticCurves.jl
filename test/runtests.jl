using Nemo
using Primes
using EllipticCurves

using Base.Test

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

print("done\n")


######################################################################
# Testing montgomery.jl
######################################################################

#Really trivial code, no need for testing

######################################################################
# Testing montgomerypoints.jl
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

@test domain(phi) == E
@test image(phi) == EllipticCurve(F(24), F(24)) #given by Sage
@test phi(P) == zero(image(phi))
@test isvalid(phi(Q))
@test phi(Q) == Point(F(91), F(20), F(1), image(phi))

K = subgroup(P, 3)

@test kernel(phi) == K
@test Isogeny(E, K) == phi
@test Isogeny(E, 3, image(phi)) == phi

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
@test Isogeny(E, image(phi), phix, phiy) == phi

print("done\n")

#Isogeny(E, d, jprime) is not tested, since there are no modular polynomials yet

#=

######################################################################
# Testing finite.jl
######################################################################

print("Testing isogeny computation over finite fields...\n")

F, _ = FiniteField(233, 1, "u")
Cards = Dict{Int,Nemo.fmpz}()

#Roots of polynomials

R, X = PolynomialRing(F, "X")
(bool, root) = any_root(X^2 - 1)
@test bool
@test root^2 == 1

#Conversions

G, z = FiniteField(233, 3, "z")
@test convert(G(1), F) == F(1)
@test convert(F(1), G) == G(1)

#Testing isogeny computations : cardinality is given by Sage

E = Weierstrass(F(1), F(2), F(3), F(4), F(5))
Cards[1] = ZZ(224)
Card = Cards[1]

(bool, order) = isgoodprime(E, 19, Card)
@test !bool
(bool, order) = isgoodprime(E, 7, Card)
@test bool
@test order == 1

phi = first_isogeny(E, 7, Cards)
@test j_invariant(image(phi)) == F(52) #Given by Sage


#With Montgomery curves

E = Montgomery(F(3), F(1))
Cards[1] = ZZ(240)
Card = Cards[1]

(bool, order) = isgoodprime(E, 17, Card)
@test !bool
(bool, order) = isgoodprime(E, 5, Card)
@test bool
@test order == 1


phi = first_isogeny_x(E, 5, Cards)
@test j_invariant(image(phi)) == F(76) #Given by Sage

#Testing over extensions : p must be big

F, _ = FiniteField(1267650600228229401496703205653, 1, "u")
E = Montgomery(F(3), F(1))
Cards[1] = ZZ(1267650600228228574627063707104)
Card = Cards[1]
Cards[4] = ZZ(2582249878086908589655919174260047736108316376414661336328042639672071392914646520434415276111057972712143504852274942464)

(bool, order) = isgoodprime(E, 109, Card)
@test bool
@test order == 4

phi = first_isogeny_x(E, 109, Cards)
=#


