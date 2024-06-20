## HERE WE PROVE THAT IF 7 EIGENPOINTS ARE SUCH THAT HAVE THE COLLINEARITIES:
## (P1, P2, P3), (P1, P4, P5), (P1, P6, P7), THEN FIVE OF THEM MUST 
## SATISFY A delta2 CONDITION.

varAn1 = ["s"+str(i)+str(j) for i in range(1,8) for j in range(i, 8)]
varAn2 = ["x", "y", "z", "u1", "u2", "v1", "v2", "w1", "w2", "l1", "l2", "m1", "m2"]
varAn3 = ["a"+str(i) for i in range(10)]
varAn4 = [str(XX)+str(i) for i in range(1, 8) for XX in ["A", "B", "C"]]
var("xx")
## K = QQ
K.<ii> = NumberField(xx^2+1)
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
P6 = vector(S, (A6, B6, C6))
P7 = w1*P1+w2*P6
## we study the condition delta1(P1, P2, P4)=0 and 
## delta1(P1, P2, P6) = 0, in order to see if the configuration (3)
## of the figure can be realized in this way.
pl1 = delta1(P1, P2, P4)
pl2 = delta1(P1, P2, P6)

## Here we have two linear equations in the coordinates of P2
## We construct the matrix of the system: 
M = matrix([[pl1.coefficient(aa) for aa in (A2, B2, C2)],\
            [pl2.coefficient(aa) for aa in (A2, B2, C2)]])

## and the solution:
m2 = M.minors(2)
slz = {A2: m2[2], B2: -m2[1], C2: m2[0]}
## we verify that this is the solution of pl1 = pl2 = 0:
assert((pl1.subs(slz), pl2.subs(slz)) == (S(0), S(0)))

## we get that the solution is given PP2:

PP2 = scalarProd(P1, P1)*det(matrix([P1, P4, P6]))*P1

assert(PP2 == P2.subs(slz))

## but with this solution, P1 and PP2 coincide as projective points:

assert(matrix([P1, PP2]).minors(2) == [S(0), S(0), S(0)])

## Finally, we want to consider the case in which M does not have 
## maximal rank. 
J = S.ideal(m2)
pdJ = J.radical().primary_decomposition()

## pdJ has two components: det([P1, P4, P6]) and (P1|P1):

assert(len(pdJ) == 2)

assert(pdJ[1] == S.ideal(det(matrix([P1, P4, P6]))))

assert(pdJ[0] == S.ideal(scalarProd(P1, P1)))

## Since P1, P4, P6 are assumed not collinear, it remain to study 
## the case P1 a point on the isotropic conic.

## We can assume therefore P1 = (1, ii, 0) and we redefine the 
## other points:

P1 = vector(S, (1, ii, 0))
P2 = vector(S, (A2, B2, C2))
P3 = u1*P1+u2*P2
P4 = vector(S, (A4, B4, C4))
P5 = v1*P1+v2*P4
P6 = vector(S, (A6, B6, C6))
P7 = w1*P1+w2*P6

## we study when delta1(P1, P2, P4) and 
## delta1(P1, P2, P6) is zero:

pl1 = delta1(P1, P2, P4)
pl2 = delta1(P1, P2, P6)

## Here we have two linear equations in the coordinates of P2
## We construct the matrix of the system: 
M = matrix([[pl1.coefficient(aa) for aa in (A2, B2, C2)],\
            [pl2.coefficient(aa) for aa in (A2, B2, C2)]])

## and the solution:
slz = {A2: (-ii)*B2*A4 + B2*B4, B2: (A4+ii*B4)*B2, C2: (A4+ii*B4)*C2}
assert((pl1.subs(slz), pl2.subs(slz)) == (S(0), S(0)))

## we redefine the points, according to this condition on A2, B2, C2:
p1 = P1.subs(slz)
p2 = P2.subs(slz)
p3 = P3.subs(slz)
p4 = P4.subs(slz)
p5 = P5.subs(slz)
p6 = P6.subs(slz)
p7 = P7.subs(slz)

## now delta1(p1, p2, p4) and delta1(p1, p2, p6) are zero:
assert((delta1(p1, p2, p4), delta1(p1, p2, p6)) == (S(0), S(0)))

## the matrix Phi(p1, p2, p3, p4, p5, p6, p7) must have rank 9 or less.
## in particular the rank of Phi(p1, p4, p5, p6, p7) must be <=9, 
## hence delta1(p1, p4, p6)*delta2(P1, p4, p5, p6, p7) = 0
## If delta2(P1, p4, p5, p6, p7) = 0, then the configuration is realized
## with a delta2 condition.
## If delta1(p1, p4, p6) = 0 we have:
## (A6 + ii*B6) * (A4 + ii*B4), indeed:
assert(delta1(p1, p4, p6) == -(A6 + ii*B6) * (A4 + ii*B4))
## and we study the two possibilities:
## If A6 + ii*B6 = 0:

slz1 = {A6: -ii*B6}
pp1 = p1.subs(slz1)
pp2 = p2.subs(slz1)
pp3 = p3.subs(slz1)
pp4 = p4.subs(slz1)
pp5 = p5.subs(slz1)
pp6 = p6.subs(slz1)
pp7 = p7.subs(slz1)

## In this case we have delta2(pp1, pp2, pp3, pp6, pp7), 
## so the configuration is realized by a delta2 condition:
assert(delta2(pp1, pp2, pp3, pp6, pp7) == S(0))

## If A4 + ii*B4 = 0:

slz1 = {A4: -ii*B4}
pp1 = p1.subs(slz1)
pp2 = p2.subs(slz1)
pp3 = p3.subs(slz1)
pp4 = p4.subs(slz1)
pp5 = p5.subs(slz1)
pp6 = p6.subs(slz1)
pp7 = p7.subs(slz1)

## In this casewe have delta2(pp1, pp2, pp3, pp4, pp5), 
## so the configuration is realized by a delta2 condition:
assert(delta2(pp1, pp2, pp3, pp4, pp5) == S(0))

### CONCLUSION:
### In order to have the configuration (3), i.e. the alignments:
## (P1, P2, P3), (P1, P4, P5), (P1, P6, P7)
## it is necessary to have a delta2 condition among the points,
## hence the configuration cannot be realized with only delta1 conditions.



