### Here we prove lemma: twoTangentsCiso
### which defines C(r) and says that C(r) has two lines of eigenpoints.
## About 10'' of computations

varAn1 = ["s"+str(i)+str(j) for i in range(1,8) for j in range(i, 8)]
varAn2 = ["x", "y", "z", "u1", "u2", "v1", "v2", "w1", "w2", "l1", "l2", "m1", "m2"]
varAn3 = ["a"+str(i) for i in range(10)]
varAn4 = [str(XX)+str(i) for i in range(1, 8) for XX in ["A", "B", "C"]]
var("xx")
K.<ii> = NumberField(xx^2+1)
S = PolynomialRing(K, varAn1+varAn2+varAn3+varAn4)
S.inject_variables(verbose=false)
##
P1 = vector(S, (A1, B1, C1))

load("auxiliaryProcedures.sage")

P1 = vector((1, ii, 0))

mon = [x^3, x^2*y, x*y^2, y^3, x^2*z, x*y*z, y^2*z, x*z^2, y*z^2, z^3]

load("../AlignedEigenpoints/software/ancillary.sage")

test = SymbolicCheck()

## 
## construction of a generic point on the isotropic conic:
rt1 = l1*(y-ii*x)+l2*z

rt1.subs({x:P1[0], y:P1[1], z:P1[2]})

Ciso = x^2+y^2+z^2
scndP = S.ideal(Ciso, rt1).radical().primary_decomposition()[1]
aux = scndP.gens()[:2]
mm2 = matrix([[aux[0].coefficient(x), aux[0].coefficient(y), \
aux[0].coefficient(z)],\
       [aux[1].coefficient(x), aux[1].coefficient(y), \
       aux[1].coefficient(z)]]).minors(2)

## Generic point on Ciso (depending on two parameters):
PP = vector(S, (mm2[2], -mm2[1], mm2[0]))

assert(scndP.subs({x:PP[0], y:PP[1], z:PP[2]}) == S.ideal(S(0)))
assert(Ciso.subs({x:PP[0], y:PP[1], z:PP[2]}) == S(0))

## we can always assume that l2 != 0, since 
## l2 = 0 gives that PP = P1
assert(matrix([P1, PP.subs(l2=0)]).rank() == 1)

## PP is the point: 
## vector(S, ((-ii)*l1^2 + (-ii)*l2^2, l1^2 - l2^2, 2*l1*l2))

assert(PP == vector(S, ((-ii)*l1^2 + (-ii)*l2^2, l1^2 - l2^2, 2*l1*l2)))

## Now that we know the generic point of Ciso, we 
## define two (distinct) points on Ciso:

PP1 = vector(S, ((-ii)*l1^2 + (-ii)*l2^2, l1^2 - l2^2, 2*l1*l2))
PP2 = vector(S, ((-ii)*m1^2 + (-ii)*m2^2, m1^2 - m2^2, 2*m1*m2))

## and the lines ttg1, ttg2, the first is tangent to Ciso in PP1, 
## the second is tangent ot Ciso in PP2:

## tangent to Ciso in PP1:
ttg1 = (derivative(Ciso, x).subs({x: PP1[0], y:PP1[1], z:PP1[2]}))*(PP1[0]-x)+\
(derivative(Ciso, y).subs({x: PP1[0], y:PP1[1], z:PP1[2]}))*(PP1[1]-y)+\
(derivative(Ciso, z).subs({x: PP1[0], y:PP1[1], z:PP1[2]}))*(PP1[2]-z)

## tangent to Ciso in PP2:
ttg2 = (derivative(Ciso, x).subs({x: PP2[0], y:PP2[1], z:PP2[2]}))*(PP2[0]-x)+\
(derivative(Ciso, y).subs({x: PP2[0], y:PP2[1], z:PP2[2]}))*(PP2[1]-y)+\
(derivative(Ciso, z).subs({x: PP2[0], y:PP2[1], z:PP2[2]}))*(PP2[2]-z)

## We verify that ttg1 is tanget:
assert(S.ideal(ttg1, Ciso).radical() == \
S.ideal(z*l1 + (-ii)*x*l2 - y*l2, x*l1 + ii*y*l1 + ii*z*l2, x^2 + y^2 + z^2))

## we suspect that the same result also holds for ttg2...

## we define the line PP1+PP2:
rr = matrix([PP1, PP2, (x, y, z)]).det().factor()[-1][0]

## Finally, we define the cubic C(r) = r*(r^2-3*(a^2+b^2+c^2)*Ciso)
## and we verify that C(r) has ttg1 and ttg2 as eigenpoints.

Cr = rr*(rr^2-3*(rr.coefficient(x)^2+\
         rr.coefficient(y)^2+rr.coefficient(z)^2)*Ciso)

JJ = S.ideal(matrix([[Cr.derivative(x), Cr.derivative(y), \
                      Cr.derivative(z)], [x, y, z]]).minors(2))

PDJ = JJ.radical().primary_decomposition()

assert(len(PDJ) == 2)
assert(PDJ[0] == S.ideal(ttg1))
assert(PDJ[1] == S.ideal(ttg2))


## CONCLUSION: Given a line r of the plane which intersects Ciso in two 
## distinct points P1 and P2, the eigenpoints of the cubic C(r) given by 
## r*(r^2-3*(a^2+b^2+c^2)*Ciso) are the two tangents to Ciso in P1 and P2.