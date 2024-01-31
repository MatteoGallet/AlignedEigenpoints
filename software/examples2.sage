## examples of delta2 = 0


print("We define a ring with sufficiently many variables.")

varAn1 = ["s"+str(i)+str(j) for i in range(1,8) for j in range(i, 8)]
varAn2 = ["x", "y", "z", "u1", "u2", "v1", "v2", "w1", "w2", "l1", "l2"]
varAn3 = ["a"+str(i) for i in range(10)]
varAn4 = [str(XX)+str(i) for i in range(1, 8) for XX in ["A", "B", "C"]]
varAn5 = ["l", "m", "n"]
var("xx")
##K.<im> = NumberField(xx^2+1)
K = QQ
S = PolynomialRing(K, varAn1+varAn2+varAn3+varAn4+varAn5)
S.inject_variables(verbose=false)


P1, P2, P4  = vector((A1, B1, C1)), vector((A2, B2, C2)), vector((A4, B4, C4))
P3 = u1*P1 + u2*P2
P5 = v1*P1 + v2*P4

mon = [x^3, x^2*y, x*y^2, y^3, x^2*z, x*y*z, y^2*z, x*z^2, y*z^2, z^3]
## the generic cubic
cb = add([mon[i]*S(varAn3[i]) for i in range(10)])

load("auxiliaryProcedures.sage")

Jrel = S.ideal(s11-scalarProd(P1, P1), \
               s12-scalarProd(P1, P2), \
               s22-scalarProd(P2, P2), \
               s14-scalarProd(P1, P4), \
               s24-scalarProd(P2, P4), \
               s44-scalarProd(P4, P4), \
               s15-scalarProd(P1, P5), \
               s25-scalarProd(P2, P5), \
               s45-scalarProd(P4, P5), \
               s55-scalarProd(P5, P5))

print(delta2(P1, P2, P3, P4, P5))


threeSLG1 = load("longComput/threeSLG1.sobj")
threeSLG2 = load("longComput/threeSLG2.sobj")
threeSLG3 = load("longComput/threeSLG3.sobj")
SLG1a, SLG1b, SLG1c = tuple(threeSLG1)
SLG2a, SLG2b, SLG2c = tuple(threeSLG2)
SLG3a, SLG3b, SLG3c = tuple(threeSLG3)
