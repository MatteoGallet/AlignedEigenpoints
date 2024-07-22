## Here we prove that if the points of a V configuration are such that:
## delta1(P1, P2, P4) = 0
## delta1b(P1, P2, P3) = 0
## delta1b(P1, P4, P5) = 0
## (i.e. this is one of the two cases in which rk Phi(P1, ..., P5) = 8)
## then, in general, we have the only collinearity:
## (1, 2, 3), (1, 4, 5)
## (although delta2(P1, ..., P5) = 0)
## In the one dimensional linear family of all the cubics with 
## P1, ..., P5 eigenpoints, there is one cubic which has the further 
## collinearity (1, 6, 7).
## Moreover, if we have the collinearities (1, 2, 3), (1, 4, 5) and (1, 6, 7)
## we cannot find cubics with eigenpoints with other collinearities.
## This is not in contradiction with the results of three_d_partA.sage, 
## because there the configuration (8) does not involve the the 
## alignment (1, 6, 7).
## This file considers one possible solution of the equation 
## delta1(P1, P2, P4) = 0, 
## The file three_deltas_II.sage considers the same problem, but 
## with the other solution of the above equation.



varAn1 = ["s"+str(i)+str(j) for i in range(1,8) for j in range(i, 8)]
varAn2 = ["x", "y", "z", "u1", "u2", "v1", "v2", "w1", "w2", "l1", "l2"]
varAn3 = ["a"+str(i) for i in range(10)]
varAn4 = [str(XX)+str(i) for i in range(1, 8) for XX in ["A", "B", "C"]]
var("xx")
K = QQ
K.<ii> = NumberField(xx^2+1)
S = PolynomialRing(K, varAn1+varAn2+varAn3+varAn4)
S.inject_variables(verbose=false)

P1 = vector(S, (A1, B1, C1))

load("auxiliaryProcedures.sage")

P1 = vector((1, 0, 0))
P2 = vector(S, (A2, B2, C2))
P4 = vector(S, (A4, B4, C4))
P3 = u1*P1+u2*P2
P5 = v1*P1+v2*P4

load("../AlignedEigenpoints/software/ancillary.sage")

test = SymbolicCheck()

## We want delta1(P1, P2, P4) = 0.
## delta1 is B2*B4+C2*C4:
assert(delta1(P1, P2, P4) == B2*B4+C2*C4)


## It is not possible to have C2 = 0 and C4 = 0
## here we assume C2  != 0. If C2 = 0, necessarily B2 != 0.
## If C2 != 0 and B4 = 0, we get C4 = 0, but P4 != P1, so
## if C2 != 0, we can assume B4 != 0.
## so the other case to consider is C2 = 0, B4 = 0
## in which case the points are 
## p1 = (1, 0, 0), p2 = (A2, B2, 0), p4 = (A4, 0, C4).
## Here we assume C2 != 0
st1 = {A4: A4*C2, B4: B4*C2, C4:-B2*B4}
p1 = P1.subs(st1)
p2 = P2.subs(st1)
p4 = P4.subs(st1)
assert(delta1(p1, p2, p4) == 0)
## we define p3 and p5 according to the known formulas:
 
p3 = (scalarProd(p1, p2)^2+scalarProd(p1, p1)*scalarProd(p2, p2))*p1-\
    2*(scalarProd(p1, p1)*scalarProd(p1, p2))*p2
p5 = (scalarProd(p1, p4)^2+scalarProd(p1, p1)*scalarProd(p4, p4))*p1-\
    2*(scalarProd(p1, p1)*scalarProd(p1, p4))*p4

## the entries of p5 can be divided by B4 (which is, as said, not 0)
assert(gcd(list(p5)) == B4)

p5 = vector(S, [px.quo_rem(B4)[0] for px in list(p5)])

assert(delta1b(p1, p2, p3) == 0)
assert(delta1b(p1, p4, p5) == 0)
## Incidentally, delta2 is also 0:
assert(delta2(p1, p2, p3, p4, p5) == 0)

## an example shows that, in general,  p6 and p7 are not aligned with p1.
## Here is an example. The 7 eigenpoints are such that p6 and p7 
## are not aligned with p1.

#ss3 = {A2:1, B2:-5, C2:-3, A4:7, B4:-5}
#pp1 = p1.subs(ss3)
#pp2 = p2.subs(ss3)
#pp3 = p3.subs(ss3)
#pp4 = p4.subs(ss3)
#pp5 = p5.subs(ss3)
#cb = cubic_from_matrix(matrixEigenpoints([pp1, pp2, pp3, pp4, pp5]).\
#                 stack(matrix([[2, 3, 4, 5, 6, 7, 8, 9, 1, 2]])))
#eigenpoints(cb)


### NOW WE WANT TO SEE WHAT HAPPENS IF WE IMPOSE THE ALIGNMENT:
## p1, p6, p7.
##  We start with p1, p2, p3, p4, p5 as above, such that 
##  delta1(p1, p2, p4), delta1b(p1, p2, p3) and delta1b(p1, p4, p5) are 0

## the matrix Phi(p1, p2, p3, p4, p5) has rank 8. 
M = matrixEigenpoints([p1, p2, p3, p4, p5])
assert(M.rank() == 8)

## We select 8 linearly independent rows: 

mm = M.matrix_from_rows([0, 1, 4, 5, 7, 8, 12, 14])

## in general, mm has rank 8
assert(mm.rank() == 8)  

## let us see when it is not 8:
hj = S.ideal(mm.minors(8))
hj = hj.saturation(S.ideal(matrix([p3, p5]).minors(2)))[0]
hj = hj.saturation(S.ideal(matrix([p1, p3]).minors(2)))[0]
hj = hj.saturation(S.ideal(matrix([p1, p5]).minors(2)))[0]
hj = hj.saturation(S.ideal(matrix([p2, p3]).minors(2)))[0]
## 
## it holds: hj = (1)
assert(hj == S.ideal(1))
## hence the order 8 minor mm has always rank 8.

## The cubics of P^9 which have p1, ..., p5 as eigenpoints are a line L
## of P^9. We want to see if, among the points of this line, we can find
## some which give cubic curves with the eigenpoints p1, p6, p7 collinear.
## Here is the way in which we procede.
## We fix two 10-components vectors like:
## (1, 2, 5, 6, 0, 2, 3, 4, 9, 11) and (-1, 3, 6, 5, 0, 1, 3, 7, 9, -5)
## and we consider the following matrices: mmA and mmB
## (we put an index "1" because later we shall repeat the computation):

mmA_1 = mm.stack(matrix([1, 2, 5, 6, 0, 2, 3, 4, 9, 11]))
mmB_1 = mm.stack(matrix([-1, 3, 6, 5, 0, 1, 3, 7, 9, -5]))

## mmA and mmB have rank 9 (if not, take two other points!)
## From mmA we get one cubic curve cbA, from mmB we get one
## cubic cbB and they are two points of the line L of P^9.
## Both cbA and cbB have, among the eigenpoints, p1, p2, p3, p4, p5.
## Given the generic point P = (x: y: z) of P^2, consider the three 
## matrices obtained from mmA_1 given by
## mmA plus the row phi_p(P)[0], 
## mmA plus the row phi_p(P)[1] and
## mmA plus the row phi_p(P)[2].
## We get three square matrices with determinant GA1, GA2, GA3 such 
## that p1[2]*GA1-p1[1]*GA2+p1[0]*GA3 splits, as a polynomial in 
## x, y, z, into three linear factors, one corresponds to the line 
## p1+p2, the other to the line p1+p4 and the third to the line
## passing through the two remaining eigenpoints of cbA.
## A simplification comes from the fact that p1 is the point 
## (1: 0: 0), hence p1[2]*GA1-p1[1]*GA2+p1[0]*GA3 is GA3. Hence we have
## to factorize it.

## The same construction can be done for mmB and we get a polynomial 
## GB3 which factorizes into three linear factors in x, y, z, the 
## first two correspond again to the lines p1+p2 and p1+p4, the third 
## corresponds to the line through the two remaining eigenpoints of cbB.
## A point of the line L corresponds to the cubic cb = w1*cbA+w2*cbB,
## where w1 and w2 are parameters.
## The cubic cb has p1, p2, p3, p4, p5 as eigenpoints and two other 
## eigenpints, say p6 and p7, that are obtained from the factorization 
## of w1*GA3+w2*GB3. Now the explicit computations (we repeat
## the construction twice):


GA3_1 = det(mmA_1.stack(matrix([phi_p((x, y, z))[2]])))
GB3_1 = det(mmB_1.stack(matrix([phi_p((x, y, z))[2]])))

rr3_1 = list(filter(lambda uu: w1 in uu[0].variables(), \
                                list(factor(w1*GA3_1+w2*GB3_1))))[0][0]

## rr3_1 is a polynomial of degree 1 in x, y, z which gives 
## the line passing through the eigenpoints p6 and p7 of 
## cb (it depends of w1 and w2).

## We want to see if, among the cubics cb = cb(w1, w2), it is possible 
## to find a cubic with p1, p6, p7 aligned. Hence p1 must be a point of rr3_1,
## so the following polynomial must be zero:

hh1 = rr3_1.subs({x:1, y:0, z:0})

## In order to get rid of the choice of the two rows above (i.e. the 
## choice of two points of L), we repeat the construction for two other 
## random rows:

mmA_2 = mm.stack(matrix([1, -5, 1, 2, 0, 1, -2, 1, 3, 7]))
mmB_2 = mm.stack(matrix([-1, -1, 0, 4, 0, 1, 0, 1, 0, -5]))

GA3_2 = det(mmA_2.stack(matrix([phi_p((x, y, z))[2]])))
GB3_2 = det(mmB_2.stack(matrix([phi_p((x, y, z))[2]])))

rr3_2 = list(filter(lambda uu: w1 in uu[0].variables(), \
                                list(factor(w1*GA3_2+w2*GB3_2))))[0][0]

## if p1, p6, p7 are alligned, also the following polynomial must be zero
hh2 = rr3_2.subs({x:1, y:0, z:0})

## hh1 (and hh2) are of the form w1*()+w2*(). Hence hh1 (and hh2) 
## is zero iff w1 and w2 are chosen as solution of the equation 
## w1*()+w2*() = 0 or if the coefficients of w1 and w2 are both zero.
## In the first case, we construct r3_1 and r3_2 as follows:

r3_1 = (w1*GA3_1+w2*GB3_1).subs({w1: hh1.coefficient(w2), \
                           w2: -hh1.coefficient(w1)}).factor()[-2][0]
r3_2 = (w1*GA3_2+w2*GB3_2).subs({w1: hh2.coefficient(w2), \
                           w2: -hh2.coefficient(w1)}).factor()[-2][0]

## (i.e. r3_1 and r3_2 are the line passing through p1, p6, p7.
## They should be equal, because they should not depend of the 
## two points of the line L chosen. Indeed, they are equal:

assert(r3_1 == r3_2)

## Now we have to consider the case in which r3_1 and r3_2 are not defined,
## i.e. when the coefficients of w1 and w2 in hh1 and hh2 are all zero:

HH = [hh1, hh2]
JJ = S.ideal([hh.coefficient(w1) for hh in HH]+[hh.coefficient(w2) for hh in HH])

## If we have some values of the points p1, ..., p5 such that give a zero
## of JJ, we have to study that case.

## But JJ, after saturation, is (1):

JJ = JJ.saturation(S.ideal(matrix([p1, p2]).minors(2)))[0]
JJ = JJ.saturation(S.ideal(matrix([p1, p3]).minors(2)))[0]
JJ = JJ.saturation(S.ideal(matrix([p1, p4]).minors(2)))[0]
JJ = JJ.saturation(S.ideal(matrix([p1, p5]).minors(2)))[0]
JJ = JJ.saturation(S.ideal(matrix([p3, p5]).minors(2)))[0]
assert(JJ == S.ideal(1))

## This computation shows that there are no exceptions to consider. 
## We have that on the line L of P^9 there is a point which 
## corresponds to a cubic which has the alignments ## p1, p2, p3,
## p1, p4, p5 and p1, p6, p7 among the eigenpoints. We can determine 
## the cubic:

MM_1 = (w1*mmA_1+w2*mmB_1).subs({w1: hh1.coefficient(w2), \
                           w2: -hh1.coefficient(w1)})

Mcb1 = MM_1.stack(matrix([[x^3, x^2*y, x*y^2, y^3, x^2*z, x*y*z, y^2*z, \
                         x*z^2, y*z^2, z^3]]))

MM_2 = (w1*mmA_2+w2*mmB_2).subs({w1: hh2.coefficient(w2), \
                           w2: -hh2.coefficient(w1)})

Mcb2 = MM_2.stack(matrix([[x^3, x^2*y, x*y^2, y^3, x^2*z, x*y*z, y^2*z, \
                         x*z^2, y*z^2, z^3]]))

## The cubic is the determinant of Mcb1 (or of Mcb2).

## the following computations require, respectively, about 3000 and 
## 4000 seconds. We omit them but we define below the cubic cb which 
## is obtained:

# ttA = cputime()
# cb1 = Mcb1.det()
# print(cputime()-ttA)
# 
# ttA = cputime()
# cb2 = Mcb2.det()
# print(cputime()-ttA)
# 
# assert(cb1.factor()[-1][0]  == cb2.factor()[-1][0])

cb = z^3*A2*B2^3*C2^2*A4^2 - 3*y*z^2*A2*B2^2*C2^3*A4^2 + 3*y^2*z*A2*B2*C2^4*A4^2 - y^3*A2*C2^5*A4^2 - y^3*A2^2*B2^3*C2*A4*B4 + 3/2*x*y^2*A2*B2^4*C2*A4*B4 + 3/2*x*z^2*A2*B2^4*C2*A4*B4 + 1/2*y^3*B2^5*C2*A4*B4 - 3*y^2*z*A2^2*B2^2*C2^2*A4*B4 + 3/2*y^2*z*B2^4*C2^2*A4*B4 - 3*y*z^2*A2^2*B2*C2^3*A4*B4 + 3*x*y^2*A2*B2^2*C2^3*A4*B4 + 3*x*z^2*A2*B2^2*C2^3*A4*B4 + 1/2*y^3*B2^3*C2^3*A4*B4 + 3/2*y*z^2*B2^3*C2^3*A4*B4 - z^3*A2^2*C2^4*A4*B4 + 3/2*y^2*z*B2^2*C2^4*A4*B4 + 1/2*z^3*B2^2*C2^4*A4*B4 + 3/2*x*y^2*A2*C2^5*A4*B4 + 3/2*x*z^2*A2*C2^5*A4*B4 + 3/2*y*z^2*B2*C2^5*A4*B4 + 1/2*z^3*C2^6*A4*B4 - 1/2*z^3*A2*B2^5*B4^2 + 3/2*y*z^2*A2*B2^4*C2*B4^2 - 3/2*y^2*z*A2*B2^3*C2^2*B4^2 - 1/2*z^3*A2*B2^3*C2^2*B4^2 + 1/2*y^3*A2*B2^2*C2^3*B4^2 + 3/2*y*z^2*A2*B2^2*C2^3*B4^2 - 3/2*y^2*z*A2*B2*C2^4*B4^2 + 1/2*y^3*A2*C2^5*B4^2


## We remark however that cb can be computed in 12 seconds usign the 
## fact that the rows of Mcb1 have big common factors:
Ms = []
for i in range(10):
    gd = gcd([Mcb1[i,j] for j in range(10)])
    Ms.append([Mcb1[i,j].quo_rem(gd)[0] for j in range(10)])

print("12 seconds of computation:")
sleep(1)
ttA = cputime()
cb_alt = matrix(Ms).det().factor()[-1][0]
print(cputime()-ttA)

assert(cb_alt == cb)

## An example shows that cb, in general, has the following aligned eigenpoints:
## p1, p2, p3 and p1, p4, p5 and p1, p6, p7. The cubic is irreducible and 
## singular in p1.

ccb = cb.subs({A2:5, B2:-3, C2:-1, A4:2, B4:-7})

## ## HERE WE CONCLUDE THE FIRST PART OF THE COMPUTATION:
##                         In case C2 != 0, 
## IT IS POSSIBLE TO HAVE THREE ALIGNMENTS:
### HENCE IT IS POSSIBLE TO HAVE THE ALIGNMENTS (1, 2, 3), (1, 4, 5), (1, 6, 7)
### WITH NO OTHER ALIGNMENTS.

## Now we want to see if it is possible to have more then three alignments
## (We continue to assume C2 != 0)

## Recall that cb is our cubic.

## Recall that r3_1 (= r3_2) is the line through p6 and p7.

## Up to a permutation of the indices of the points, if there is another 
## alignment among the eigenpoints, p6 must be on the line p2+p4. Hence
## we can find it, since is the intersection or p2+p4 and r3.

r24 = det(matrix([p2, p4, (x, y, z)]))

E1 = matrix([[r3_1.coefficient(xx) for xx in [x, y, z]], \
             [r24.coefficient(xx) for xx in [x, y, z]]]).minors(2)

## from Cramer (or Rouche' Capelli?) the intersection of the two lines
## is obtained from E1. Since all the entries of E1 can be divided 
## by B2^2*B4 + C2^2*B4 and we know that B2^2*B4 + C2^2*B4=0 implies 
## p3 = p5, we divide with now problems.

p6 = vector(S, (E1[2]/(B2^2*B4 + C2^2*B4), -E1[1]/(B2^2*B4 + C2^2*B4), \
                E1[0]/(B2^2*B4 + C2^2*B4)))

## Now we compute the ideal kJ of the eigenpoints of cb and we 
## saturate it as much as possible (in particular, we saturate it
## w.r.t. the ideals of the points p1, p2, p3, p4, p5):

kJ = S.ideal(matrix([[x, y, z], \
        [cb.derivative(x), cb.derivative(y), cb.derivative(z)]]).minors(2))

kJ = kJ.saturation(B2^2*B4 + C2^2*B4)[0]
kJ = kJ.saturation(S.ideal(matrix([p1, p3]).minors(2)))[0]
kJ = kJ.saturation(S.ideal(matrix([p1, p5]).minors(2)))[0]
kJ = kJ.saturation(S.ideal(z, y))[0]  ##p1
kJ = kJ.saturation(S.ideal(p2[0]*y-p2[1]*x, p2[0]*z-p2[2]*x, \
                           p2[1]*z-p2[2]*y))[0] ## p2
kJ = kJ.saturation(S.ideal(p3[0]*y-p3[1]*x, p3[0]*z-p3[2]*x, \
                           p3[1]*z-p3[2]*y))[0] ## p3
kJ = kJ.saturation(S.ideal(p4[0]*y-p4[1]*x, p4[0]*z-p4[2]*x, \
                           p4[1]*z-p4[2]*y))[0] ## p4
kJ = kJ.saturation(S.ideal(p5[0]*y-p5[1]*x, p5[0]*z-p5[2]*x, \
                           p5[1]*z-p5[2]*y))[0] ## p5

##
## after these computations, kJ is the ideal of the two reminining eigenpoints.
## we want that p1 defined above is an eigenpoint, so the ideal kkJ 
## here defined must be zero:

kkJ = kJ.subs({x:p6[0], y:p6[1], z:p6[2]}).radical()

## we saturate kkJ and we get a primary decomposiiton:
kkJ = kkJ.saturation(S.ideal(matrix([p1, p3]).minors(2)))[0]
kkJ = kkJ.saturation(S.ideal(matrix([p1, p4]).minors(2)))[0]
kkJ = kkJ.saturation(S.ideal(matrix([p1, p5]).minors(2)))[0]
kkJ = kkJ.saturation(S.ideal(matrix([p3, p5]).minors(2)))[0]
kkJ = kkJ.saturation(S.ideal(matrix([p2, p6]).minors(2)))[0]
kkJ = kkJ.saturation(S.ideal(matrix([p4, p6]).minors(2)))[0]


## kkJ is (1), so there are no solutions:

assert(kkJ == S.ideal(1))

#################################################
####                  CONCLUSION 1
#################################################
## CONCLUSION (for the case C2 != 0) when the three deltas are zero, 
## (hence rank of the matrix of the five points is 8), we have 
## also delta2 = 0 and we have the collineariries (1, 2, 3) (1, 4, 5).
## We have a sub-case in which there is the further collinearity (1, 6, 7).
## No other collinearities among the 7 eigenpoints are possible.

## ## The case C2 = 0 is on the file three_deltas_II.sage
