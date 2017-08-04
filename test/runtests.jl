using Nemo
using Primes
using EllipticCurves

using Base.Test


######################################################################
# Testing weierstrass.jl
######################################################################


#Example is taken from Silverman 1 p.59.

print("Testing basic methods of elliptic curves and projective points...\n")

E = ShortWeierstrass(QQ(0), QQ(17))
Eprime = ShortWeierstrass(QQ(0), QQ(16))
P1 = EllipticPoint(QQ(-2), QQ(3), QQ(1), E)
P3 = EllipticPoint(QQ(2), QQ(5), QQ(1), E)

@test isvalid(P1)
@test isvalid(P3)

P2 = EllipticPoint(QQ(-4), QQ(6), QQ(2), E)

@test areequal(P1, P2)
@test !areequal(P1, P3)

P2 = EllipticPoint(QQ(-2), QQ(3), QQ(1), Eprime)

@test !isvalid(P2)
@test !areequal(P1, P2)

P2 = infinity(E)

@test isvalid(P2)
@test isinfinity(P2)
@test !areequal(P1, P2)



#Testing invariants; values are given by Sage

print("Testing invariants and model changes...\n")


@test a_invariants(E) == (QQ(0), QQ(0), QQ(0), QQ(0), QQ(17))
@test b_invariants(E) == (QQ(0), QQ(0), QQ(68), QQ(0))
@test c_invariants(E) == (QQ(0), QQ(-14688))
@test discriminant(E) == QQ(-124848)
@test j_invariant(E) == QQ(0)
@test isvalid(E)

E2 = Weierstrass(QQ(1), QQ(2), QQ(3), QQ(4), QQ(6))

@test a_invariants(E2) == (QQ(1), QQ(2), QQ(3), QQ(4), QQ(6))
@test b_invariants(E2) == (QQ(9), QQ(11), QQ(33), QQ(44))
@test c_invariants(E2) == (QQ(-183), QQ(-4293))
@test discriminant(E2) == QQ(-14212)
@test j_invariant(E2) == QQ(6128487//14212)
@test isvalid(E2)




#Testing changes of variables

E1 = ShortWeierstrass(QQ(0), QQ(17))
E2 = Weierstrass(QQ(0), QQ(0), QQ(0), QQ(0), QQ(17))
P2 = EllipticPoint(QQ(-2), QQ(3), QQ(1), E2)
#E3, phi, psi = tolongWeierstrass(E)

#@test E2 == E3
#@test Eval(phi, P1) == P2
#@test Eval(psi, P2) == P1

#E2 = Weierstrass(QQ(1), QQ(2), QQ(3), QQ(4), QQ(6))
#E3, _, _ = toshortWeierstrass(E2)

#@test j_invariant(E2) == j_invariant(E3)



#Testing isogenies between Weierstrass curves over the rationals

print("Testing isogenies...\n")

E = ShortWeierstrass(QQ(1), QQ(2))

P1 = divisionpolynomial(E, 1)
P7 = divisionpolynomial(E, 7)

@test P1 == 1
@test degree(P7) == (7^2 - 1) // 2

Q = 1//lead(P7) * P7
phi = Isogeny(E, Q)
Eprime = image(phi)
d = degree(phi)


@test d == 7^2
#= Very long
psi = Isogeny(E, Eprime, d)

@test kernel(psi) == Q
=#
id = Isogeny(E, E, 1)

@test kernel(id) == 1

#Over small finite fields

F, x = FiniteField(233, 1, "x")
E = ShortWeierstrass(F(1), F(2))

P1 = divisionpolynomial(E, 1)
P7 = divisionpolynomial(E, 7)

@test P1 == 1
@test degree(P7) == (7^2 - 1) // 2

Q = 1//lead(P7) * P7
phi = Isogeny(E, Q)
Eprime = image(phi)
d = degree(phi)

@test d == 7^2

psi = Isogeny(E, Eprime, d)

#Over large finite fields

@test kernel(psi) == Q^2
@test squarefree_kernel(psi) == Q

p = ZZ(3273390607896141870013189696827599152216642046043064789483291368096133796404674554883270092325904157150886684127560071009217256545885393053328527589431)

F, x = FiniteField(p, 1, "x")
E = ShortWeierstrass(F(1), F(2))

P1 = divisionpolynomial(E, 1)
P7 = divisionpolynomial(E, 7)

@test P1 == 1
@test degree(P7) == (7^2 - 1) // 2

Q = 1//lead(P7) * P7
phi = Isogeny(E, Q)
Eprime = image(phi)
d = degree(phi)

@test d == 7^2

psi = Isogeny(E, Eprime, d)

@test kernel(psi) == Q^2
@test squarefree_kernel(psi) == Q

######################################################################
# Testing points.jl
######################################################################


#Addition laws for points on Weierstrass curves. Example is from Silverman 1 p.59

print("Testing addition laws...\n")


E = ShortWeierstrass(QQ(0), QQ(17))
P1 = EllipticPoint(QQ(-2), QQ(3), QQ(1), E)
P2 = EllipticPoint(QQ(-1), QQ(4), QQ(1), E)
P3 = EllipticPoint(QQ(2), QQ(5), QQ(1), E)
Pinf = infinity(E)

@test P1 + P1 == EllipticPoint(QQ(8), QQ(-23), QQ(1), E)
@test Pinf + Pinf == Pinf
@test P1 + Pinf == Pinf + P1 == P1
@test P2 + P2 == EllipticPoint(QQ(137//64), QQ(-2651//512), QQ(1), E)
@test P2 + P3 == EllipticPoint(QQ(-8//9), QQ(-109//27), QQ(1), E)

mP3 = -P3

@test mP3 == EllipticPoint(QQ(2), QQ(-5), QQ(1), E)
@test P3 + mP3 == Pinf
@test P1 + mP3 == P1 - P3 == EllipticPoint(QQ(4), QQ(9), QQ(1), E)




######################################################################
# Testing montgomery.jl
######################################################################

print("Testing Montgomery curves...\n")

#Testing basic functions


E = Montgomery(QQ(7), QQ(1))
P1 = EllipticPoint(QQ(1), QQ(3), QQ(1), E)
P0 = EllipticPoint(QQ(0), QQ(0), QQ(1), E)

@test base_ring(E) == QQ
@test j_invariant(E) == 24918016//45

@test isvalid(P1) & isvalid(P2)
@test isvalid(E)


#Testing model changes


Z = Nemo.ZZ
#E2, _, _ = tolongWeierstrass(E)

#@test a_invariants(E2) == (0, 7, 0, 1, 0)

#Testing x-only arithmetic


xP1 = xonly(P1)
xP0 = xonly(P0)
xinf = xinfinity(E)

#@test isfixedtorsion(xP0)
#@test !isfixedtorsion(xP1)

#=
No addition for complete projective points on Montgomery curves yet
@test plus(P0, P0) == infinity(E)
@test plus(P1, P1) == P0
=#

@test areequal(xdouble(xP0), xinf)
@test areequal(xdouble(xP1), xP0)
@test areequal(xdouble(xinf), xinf)



@test areequal(xadd(xP1, xP0, xP1), xP1)


@test areequal(0 * xP1, xinfinity(E))
@test areequal(1 * xP1, xP1)
@test areequal(2 * xP1, xP0)
@test areequal(1 * xP0, xP0)
@test areequal(2 * xP0, xinfinity(E))
@test areequal(45 * xinfinity(E), xinfinity(E))
@test areequal(4 * xP1, xinfinity(E))

@test areequal(8 * xP1, xinfinity(E))


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


#=
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


