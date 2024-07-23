##
## we know that delta2 = 0 in general gives that 
## P3 = (s14*s15*s22-s12^2*s45)*P1+ s12*(s11*s45-s14*s15)*P2
## hence, for some values of sij, P3 is not defined. Here we want to see 
## when this happens.
## The conclusion is that P3 is always uniquely defined, unless:
## 1) the points of the V-configuration are such that s12=0, s14=0
##    (or s13=0, s15=0, or...). In this case the unique cubic is a 
##    cubic in which the eigenpoints are in configuration (5).
## 2) the points of the V-configuration satisfy the condition 
##    s12 = 0, s23 = 0. In this case the line P1+P2 of the V-configuration
##    is tangent in P2 to Ciso, so the matrix Phi(P1, P2, P3) has 
##    rank 5 and P5 can be chosen in a free way on the line P1+P4 and the 
##    matrix Phi(P1, P2, P3, P4, P5) has rank 9, (unless P1+P4 is tangent
##    to Ciso...)
## 3) The points satisfy the conditions sigma(P1, P2) = 0 and 
##    sigma(P1, P4) = 0. In this case the V configuration is given 
##    by a point P1 not on Ciso and the other points on the two tangents 
##    to Ciso through P1. If no points of the V configurations are on 
##    Ciso, the two tangent lines are lines of eigenpoints, if (say)
##    P2 and P4 are on Ciso, the matrix Phi(P1, P2, P3, P4) has rank 8.
##

varAn1 = ["s"+str(i)+str(j) for i in range(1,8) for j in range(i, 8)]
varAn2 = ["x", "y", "z", "u1", "u2", "v1", "v2", "w1", "w2", "l1", "l2", "m1", "m2"]
varAn3 = ["a"+str(i) for i in range(10)]
varAn4 = [str(XX)+str(i) for i in range(1, 8) for XX in ["A", "B", "C"]]
var("xx")
K = QQ
K.<ii> = NumberField(xx^2+1)
S = PolynomialRing(K, varAn1+varAn2+varAn3+varAn4)
S.inject_variables(verbose=false)


P1 = vector(S, (A1, B1, C1))

load("auxiliaryProcedures.sage")


load("../AlignedEigenpoints/software/ancillary.sage")

test = SymbolicCheck()


P2 = vector(S, (A2, B2, C2))
P4 = vector(S, (A4, B4, C4))

## Non si sa mai:
## Jrel = S.ideal(s11-scalarProd(P1, P1), \
##                s12-scalarProd(P1, P2), \
##                s22-scalarProd(P2, P2), \
##                s14-scalarProd(P1, P4), \
##                s24-scalarProd(P2, P4), \
##                s44-scalarProd(P4, P4), \
##                s13-scalarProd(P1, P3), \
##                s23-scalarProd(P2, P3), \
##                s34-scalarProd(P3, P4), \
##                s33-scalarProd(P3, P3))


## invJrel = {s11: scalarProd(P1, P1), s12: scalarProd(P1, P2), \
##            s13: scalarProd(P1, P3), s14: scalarProd(P1, P4), \
##            s15: scalarProd(P1, P5), s22: scalarProd(P2, P2), \
##            s23: scalarProd(P2, P3), s24: scalarProd(P2, P4), \
##            s25: scalarProd(P2, P5), s33: scalarProd(P3, P3), \
##            s34: scalarProd(P3, P4), s35: scalarProd(P3, P5), \
##            s44: scalarProd(P4, P4), s45: scalarProd(P4, P5), \
##            s55: scalarProd(P5, P5)}

mon = [x^3, x^2*y, x*y^2, y^3, x^2*z, x*y*z, y^2*z, x*z^2, y*z^2, z^3]

print("First computations: about 10 sec...")
sleep(1)

## delta2 in terms of sij
d2s = s12*s13*s45-s14*s15*s23  

## In general, we have that 5 points in a V-configuration are eigenpoints 
## iff P3 = (s14*s15*s22-s12^2*s45)*P1+s12*(s11*s45-s14*s15)*P2, 
## unless the two coefficients of P1 and P2 are 0. Here we study this 
## particular case.
## 
## Hence:
Js = S.ideal(s14*s15*s22-s12^2*s45, s11*s12*s45-s12*s14*s15, d2s)

## We study the ideal Js:
pd = Js.radical().primary_decomposition()

## we get 6 ideals:
assert(len(pd)==6)

## the first ideal is generated by s14 and s12:

assert(pd[0] == S.ideal(s12, s14))

## we construct 5 points in a V configuration such that 
## s12 = 0, s14 = 0:
P2 = vector(S, (A2, B2, C2))
P4 = vector(S, (A4, B4, C4))
P1 = vector(S, list(wedgeProd(P2, P4)))
P3 = u1*P1+u2*P2
P5 = v1*P1+v2*P4

assert(scalarProd(P1, P2) == 0)
assert(scalarProd(P1, P4) == 0)

## we have: delta2(P1, P2, P3, P4, P5) = 0:
assert(delta2(P1, P2, P3, P4, P5) == 0)
## we see therefore that we can define P3 and P5 in an arbitrary way:
## for any choice delta2 is zero.

## from here the computations are exactly the computations in the file 
## config5.sage (up to a permutation of the points) and we see that in 
## this case the eigenpoints of the cubic are in configuration (5) of the 
## figure of the possible configurations.

#### CONCLUSION FOR THE CASE s12=0, s14=0:
## If s12=0, s14=0, then the V-configuration has delta2 = 0 and p3 and P5 
## are free (on the line P1+P2 and P1+P4 respectrively). 

######################## 
## the ideal pd[1] is s12, s22, s13:

assert(pd[1] == S.ideal(s12, s22, s23))

## an easy computation shows that if s22 = 0 and s12 = 0, then s23 = 0, 
## (or if s12 = 0, s23 = 0, then s22 = 0)
## so the ideal is generated by s12 and s22. In particular, the line 
## P1+P2 is tangent to the isotropic conic in P2, hence the matrix 
## Phi(P1, P2, P3) has rank 5, therefore, for all P3 and P5 the rank 
## of the matrix Phi(P1, P2, P3, P4, P5) is <= 9, therefore we have a 
## $V$ configuration which satisfies delta2(P1, P2, P3, P4, P5) = 0 
## for any choice of the points P3 and P5.
###
### CONCLUSION FOR THE CASE s12, s22:
### the V config is of eigenpoints, with every choice of P3 and P5
### (if also P1+P4 is tangent to Ciso, the two tangent lines are 
### lines of eigenpoints, as claimed in the section on eigenpoints 
### of dim > 0.
###### here is an example:
## It shows that in case s12=0 and s22 = 0 we have 7 eigenpoints
## with (apparently) no particular specific properties.

# Q1, Q2 = vector(S, (0, 0, 1)), vector(S, (1, ii, 0))
# Q3 = vector(S, (1, -ii, 0))
# rt1 = det(matrix([Q1, Q2, (x, y, z)]))
# rt2 = det(matrix([Q1, Q3, (x, y, z)]))
# 
# p1 = 1*Q1
# p2 = 1*Q2
# p3 = 2*p1+5*p2
# p4 = vector(S, (3, 1, 5))
# p5 = p1+3*p4
# assert(matrixEigenpoints([p1, p2, p3, p4, p5]).rank() == 9)
# cb = cubic_from_matrix(matrixEigenpoints([p1, p2, p3, p4, p5]))
# eigenpoints(cb)
# assert(cb.is_prime())

#######################
## The ideal pd[2] is S.ideal(s15, s12):

assert(pd[2] == S.ideal(s15, s12))
##
## this case is analog to the case s12, s14 (4<-->5)
##

#######################
## ==> The ideal pd[3] is S.ideal(s45, s14):
assert(pd[3] == S.ideal(s45, s14))

## In this case we have s44 = 0, indeed:
## (P4|v1*P1+v2*P4) = 0
## v1*(P1|P4) + v2*(P4|P4) = 0
## 0 + v2*(P4|P4) = 0
## (P4|P4) = 0

## CONCLUSION: THIS CASE IS THE SAME AS THE CASE pd[1]

## 

########################
## ==> The ideal pd[4] is S.ideal(s45, s15):
assert(pd[4] == S.ideal(s45, s15))

## In this case we have s55 = 0

## CONCLUSION: THIS CASE IS THE SAME AS THE CASE pd[1]

########################
## ==> The ideal pd[5] is:
## (s13*s22 - s12*s23, s14*s15 - s11*s45, s12*s13 - s11*s23, s12^2 - s11*s22)

assert(pd[5] == S.ideal((s13*s22 - s12*s23, s14*s15 - s11*s45, \
                         s12*s13 - s11*s23, s12^2 - s11*s22)))

## In this case we determine the ideal in the variables A1, ..., C4:
## but first we have to re-define the points:
P1 = vector(S, (A1, B1, C1))
P2 = vector(S, (A2, B2, C2))
P4 = vector(S, (A4, B4, C4))
P3 = u1*P1+u2*P2
P5 = v1*P1+v2*P4

J = pd[5].subs({s11: scalarProd(P1, P1), s22: scalarProd(P2, P2), \
                s13: scalarProd(P1, P3), s12: scalarProd(P1, P2), \
                s23: scalarProd(P2, P3), s14: scalarProd(P1, P4), \
                s15: scalarProd(P1, P5), s45: scalarProd(P4, P5)}) 

## here Jrel can help:
Jrel = S.ideal(s11-scalarProd(P1, P1), \
               s12-scalarProd(P1, P2), \
               s22-scalarProd(P2, P2), \
               s14-scalarProd(P1, P4), \
               s24-scalarProd(P2, P4), \
               s44-scalarProd(P4, P4), \
               s13-scalarProd(P1, P3), \
               s23-scalarProd(P2, P3), \
               s34-scalarProd(P3, P4), \
               s33-scalarProd(P3, P3))


## We saturate J w.r.t. the condition that v2 != 0 and 
## w.r.t. the condition that P1, P2, P4 are not aligned.

J = J.saturation(v2)[0]
J = J.saturation(det(matrix([P1, P2, P4])))[0]

## We study the ideal J

## J contains the two polynomials sigma(P1, P2) and sigma(P1, P4):

assert(sigma(P1, P2) in J)
assert(sigma(P1, P4) in J)

## moreover, J and the ideal (sigma(P1, P2), sigma(P1, P4))
## coincide if P1, P2, P4 are not aligned:
assert(J == S.ideal(sigma(P1, P2), \
            sigma(P1, P4)).saturation(det(matrix([P1, P2, P4])))[0])

## hence in this case the V  configuration is such that 
## is obtained by a point P1 and the two lines are the two tangents to 
## Ciso passing through P1. In general, we get that the V-config is 
## such that the two lines are lines of eigenpoints and the cubic splits 
## into the conic (x^2 + y^2 + 2/3*z^2) and a line. This case is discussed
## in the section on eigenpoints of positive dimension.
#### Here is an example:
#### example of P1+P2 tangent to Ciso in P2:

# Q1, Q2, Q3 = vector(S, (0, 0, 1)), vector(S, (1, ii, 0)), vector(S, (1, -ii, 0))
# rt1 = det(matrix([Q1, Q2, (x, y, z)]))
# rt2 = det(matrix([Q1, Q3, (x, y, z)]))
# 
# p1 = 1*Q1
# p2 = Q1+3*Q2
# p3 = 2*Q1+5*Q2
# p4 = Q1+7*Q3
# p5 = Q1-2*Q3
# assert(matrixEigenpoints([p1, p2, p3, p4, p5]).rank() == 9)
# cb = cubic_from_matrix(matrixEigenpoints([p1, p2, p3, p4, p5]))
# cb.factor()  ## is z*(x^2 + y^2 + 2/3*z^2)

##### In the particular case in which p2 and p4 are on Ciso, we 
## get that the matrix Phi(P1, P2, P3, P4, P5) has rank 8 and 
## a generic cubic of the corresponding pencil of cubics is smooth
## and has 7 eigenpoints. The case is studied in the section on 
## V conf of rank 8.

## Here is an example:
# Q1, Q2, Q3 = vector(S, (0, 0, 1)), vector(S, (1, ii, 0)), vector(S, (1, -ii, 0))
# rt1 = det(matrix([Q1, Q2, (x, y, z)]))
# rt2 = det(matrix([Q1, Q3, (x, y, z)]))

# p1 = 1*Q1
# p2 = 1*Q2
# p3 = 2*Q1+5*Q2
# p4 = 1*Q3
# p5 = Q1-2*Q3
# assert(matrixEigenpoints([p1, p2, p3, p4, p5]).rank() == 8)
# M1 = matrixEigenpoints([p1, p2, p3, p4, p5])
# M2 = M1.matrix_from_rows([1, 2, 3, 4, 7, 9, 10, 13])
# assert(M2.rank() == 8)
# M3 = M2.stack(vector(S, (1, 3, 6, -4, 0, 3, 2, 1, -5, 7)))
# M3 = M3.stack(vector(S, mon))
# cb = M3.det()
# eigenpoints(cb)

