## HERE WE PROVE THAT CONFIGURATION (4) IS NOT POSSIBLE:
## IT GIVES CONFIGURATION (8).

varAn1 = ["s"+str(i)+str(j) for i in range(1,8) for j in range(i, 8)]
varAn2 = ["x", "y", "z", "u1", "u2", "v1", "v2", "w1", "w2", "l1", "l2", "m1", "m2"]
varAn3 = ["a"+str(i) for i in range(10)]
varAn4 = [str(XX)+str(i) for i in range(1, 8) for XX in ["A", "B", "C"]]
var("xx")
K = QQ
## K.<ii> = NumberField(xx^2+1)
S = PolynomialRing(K, varAn1+varAn2+varAn3+varAn4)
S.inject_variables(verbose=false)

P1 = vector(S, (A1, B1, C1))

load("auxiliaryProcedures.sage")


load("../AlignedEigenpoints/software/ancillary.sage")

test = SymbolicCheck()

## 
P2 = vector(S, (A2, B2, C2))
P3 = u1*P1+u2*P2
P4 = vector(S, (A4, B4, C4))
P5 = v1*P1+v2*P4
P6 = w1*P2+w2*P4

## If we have configuration (4), then we must have:
## delta1(P1, P2, P4) = 0, delta1(P2, P1, P4) = 0, delta1(P4, P1, P2) = 0
J = S.ideal(delta1(P1, P2, P4), delta1(P2, P1, P4), delta1(P4, P1, P2))

## we saturate J w.r.t. some ideals which say that the points P1, P2, P4 
## do not have all the coordinates zero and P1, P2, P4 are not aligned:

J = J.saturation(S.ideal(A1, B1, C1))[0]
J = J.saturation(S.ideal(A2, B2, C2))[0]
J = J.saturation(S.ideal(A4, B4, C4))[0]
J = J.saturation(S.ideal(matrix([P1, P2, P4]).det()))[0]

## hence we have that J = ((P1|P2), (P1|P4), (P2|P4))

assert(J == S.ideal(scalarProd(P1, P2), \
                    scalarProd(P1, P4), \
                    scalarProd(P2, P4)))

## These conditions and the definition of delta2 give that 
## delta2(P1, P2, P3, P4, P5) = 0
## delta2(P2, P1, P3, P4, P6) = 0
## delta2(P4, P1, P5, P2, P6) = 0

## we can also verify directly these equalities:

assert(J.reduce(delta2(P1, P2, P3, P4, P5))==S.ideal(S(0)))
assert(J.reduce(delta2(P2, P1, P3, P4, P6))==S.ideal(S(0)))
assert(J.reduce(delta2(P4, P1, P5, P2, P6))==S.ideal(S(0)))


## CONCLUSION:
## configuration (4) is not possible: it gives configuration (8).