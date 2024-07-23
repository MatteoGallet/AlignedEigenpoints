## Here we construct examples
## Case (P1|P1) = 0 and (P1|P2) = 0.
print("We define a ring with sufficiently many variables.")

varAn1 = ["s"+str(i)+str(j) for i in range(1,8) for j in range(i, 8)]
varAn2 = ["x", "y", "z", "u1", "u2", "v1", "v2", "w1", "w2", "l1", "l2"]
varAn3 = ["a"+str(i) for i in range(10)]
varAn4 = [str(XX)+str(i) for i in range(1, 8) for XX in ["A", "B", "C"]]
varAn5 = ["l", "m", "n"]
var("xx")
K.<im> = NumberField(xx^2+1)
S = PolynomialRing(K, varAn1+varAn2+varAn3+varAn4+varAn5)
S.inject_variables(verbose=false)


P1, P2, P4  = vector((A1, B1, C1)), vector((A2, B2, C2)), vector((A4, B4, C4))
P3 = u1*P1 + u2*P2
P5 = v1*P1 + v2*P4

mon = [x^3, x^2*y, x*y^2, y^3, x^2*z, x*y*z, y^2*z, x*z^2, y*z^2, z^3]
## the generic cubic
cb = add([mon[i]*S(varAn3[i]) for i in range(10)])

load("auxiliaryProcedures.sage")

# We investigate whether, given three points PP1, PP2, PP4 such that PP2 and PP4 are on the isotropic quadric and PP1 is the polar of the line PP2, PP4, the matrix Phi(PP1), Phi(PP2), Phi(PP4) has rank 5 or 6

# Using the action of SO3, we can suppose that PP2 is (1: i: 0) where i is the imaginary unit
PP2 = vector(S, (1, im, 0))

# PP4 is an arbitraty point on the isotropic conic
PP4 = vector(S, (m^2-n^2, 2*m*n, im*(m^2+n^2)))

# PP1 is the polar with respect to the isotropic conic of the line PP2, PP4
PP1 = vector(S, (-(m^2+n^2), -im*(m^2+n^2), 2*m*n-im*(m^2-n^2)))

# Matrix imposing that PP1, PP2, and PP4 are aigenpoints
M = matrixEigenpoints([PP1, PP2, PP4])


