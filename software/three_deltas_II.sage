## This file makes the same computations and obtains the same 
## conclusions of the file thee_deltas_I.sage. The only difference
## is that here we consider a second solution of the equation
## delta(P1, P2, P4) = 0
## for all the explanations, see the file three_deltas_I.sage

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

###### now the case C2 = 0.
## The points are:
## p1 = (1, 0, 0), p2 = (A2, B2, 0), p4 = (A4, 0, C4).
## We are sure that B2 != 0 (since p2 != p1) and C4 != 0 (since p4 != p1)

p1 = vector(S, (1, 0, 0))
p2 = vector(S, (A2, B2, 0))
p4 = vector(S, (A4, 0, C4))

p3 = (scalarProd(p1, p2)^2+scalarProd(p1, p1)*scalarProd(p2, p2))*p1-\
    2*(scalarProd(p1, p1)*scalarProd(p1, p2))*p2
p5 = (scalarProd(p1, p4)^2+scalarProd(p1, p1)*scalarProd(p4, p4))*p1-\
    2*(scalarProd(p1, p1)*scalarProd(p1, p4))*p4

## we redefine the points, since B2 and C4 are not 0.

p3, p5 = p3/B2, p5/C4

assert(delta1b(p1, p2, p3) == 0)
assert(delta1b(p1, p4, p5) == 0)
## Incidentally, delta2 is also 0:
assert(delta2(p1, p2, p3, p4, p5) == 0)

## an example shows that, in general,  p6 and p7 are not aligned with p1.
## Here is an example. The 7 eigenpoints are such that p6 and p7 
## are not aligned with p1

#ss3 = {A2:1, B2:-5, A4:7, C4:-5}
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

mm = M.matrix_from_rows([0, 1, 3, 4, 6, 7, 9, 10])

## in general, mm has rank 8
assert(mm.rank() == 8)  

## let us see when it is not 8:
hj = S.ideal(mm.minors(8))

hj = hj.saturation(S.ideal(matrix([p1, p2]).minors(2)))[0]
hj = hj.saturation(S.ideal(matrix([p1, p3]).minors(2)))[0]
hj = hj.saturation(S.ideal(matrix([p1, p4]).minors(2)))[0]
hj = hj.saturation(S.ideal(matrix([p1, p5]).minors(2)))[0]
hj = hj.saturation(S.ideal(matrix([p2, p3]).minors(2)))[0]

## hj is (1), so mm has always rank 8.

assert(hj == S.ideal(1))

## Hence the order 8 minor mm has always rank 8.
## As above, we construct mmA_1 and mmB_1.
## The construction is the same as above.

mmA_1 = mm.stack(matrix([1, 2, 5, 6, 0, 2, 3, 4, 9, 11]))
mmB_1 = mm.stack(matrix([-1, 3, 6, 5, 0, 1, 3, 7, 9, -5]))

## these two matirces have rank 9
assert(mmA_1.rank() == 9)
assert(mmB_1.rank() == 9)

GA3_1 = det(mmA_1.stack(matrix([phi_p((x, y, z))[2]])))
GB3_1 = det(mmB_1.stack(matrix([phi_p((x, y, z))[2]])))


rr3_1 = list(filter(lambda uu: w1 in uu[0].variables(), \
                                list(factor(w1*GA3_1+w2*GB3_1))))[0][0]


hh1 = rr3_1.subs({x:1, y:0, z:0})



mmA_2 = mm.stack(matrix([1, -5, 1, 2, 0, 1, -2, 1, 3, 7]))
mmB_2 = mm.stack(matrix([-1, -1, 0, 4, 0, 1, 0, 1, 0, -5]))

GA3_2 = det(mmA_2.stack(matrix([phi_p((x, y, z))[2]])))
GB3_2 = det(mmB_2.stack(matrix([phi_p((x, y, z))[2]])))

rr3_2 = list(filter(lambda uu: w1 in uu[0].variables(), \
                                list(factor(w1*GA3_2+w2*GB3_2))))[0][0]

## if p1, p6, p7 are alligned, also the following polynomial must be zero
hh2 = rr3_2.subs({x:1, y:0, z:0})

r3_1 = (w1*GA3_1+w2*GB3_1).subs({w1: hh1.coefficient(w2), \
                           w2: -hh1.coefficient(w1)}).factor()[-1][0]
r3_2 = (w1*GA3_2+w2*GB3_2).subs({w1: hh2.coefficient(w2), \
                           w2: -hh2.coefficient(w1)}).factor()[-1][0]


## (i.e. r3_1 and r3_2 are the line passing through p1, p6, p7.
## They should be equal, because they should not depend of the 
## two points of the line L chosen. Indeed, they are equal:

assert(r3_1 == r3_2)

HH = [hh1, hh2]
JJ = S.ideal([hh.coefficient(w1) for hh in HH]+[hh.coefficient(w2) for hh in HH])

JJ = JJ.saturation(S.ideal(matrix([p1, p2]).minors(2)))[0]
JJ = JJ.saturation(S.ideal(matrix([p1, p3]).minors(2)))[0]
JJ = JJ.saturation(S.ideal(matrix([p1, p4]).minors(2)))[0]

assert(JJ == S.ideal(1))

## This computation shows that there are no exceptions to consider. 


MM_1 = (w1*mmA_1+w2*mmB_1).subs({w1: hh1.coefficient(w2), \
                           w2: -hh1.coefficient(w1)})

Mcb1 = MM_1.stack(matrix([[x^3, x^2*y, x*y^2, y^3, x^2*z, x*y*z, y^2*z, \
                         x*z^2, y*z^2, z^3]]))

MM_2 = (w1*mmA_2+w2*mmB_2).subs({w1: hh2.coefficient(w2), \
                           w2: -hh2.coefficient(w1)})

Mcb2 = MM_2.stack(matrix([[x^3, x^2*y, x*y^2, y^3, x^2*z, x*y*z, y^2*z, \
                         x*z^2, y*z^2, z^3]]))

## The cubic is the determinant of Mcb1 (or of Mcb2).

## one possibility is the following computation (not long)
cb1 = Mcb1.det().factor()[-1][0]
cb2 = Mcb2.det().factor()[-1][0]

## here is an alternative:
Ms = []
for i in range(10):
    gd = gcd([Mcb1[i,j] for j in range(10)])
    Ms.append([Mcb1[i,j].quo_rem(gd)[0] for j in range(10)])

cb_alt = matrix(Ms).det().factor()[-1][0]

## We have:
## cb1 = cb2 and cb2 = cb_alt:
assert(cb1 == cb2)
assert(cb1 == cb_alt)

## The following example shows that in general there are only the 
## collinearities (1, 2, 3), (1, 4, 5), (1, 6, 7)

#ccb = cb_alt.subs({A2:5, B2:-3, A4:2, C4:-7})
#pp1 = p1.subs({A2:5, B2:-3, A4:2, C4:-7})
#pp2 = p2.subs({A2:5, B2:-3, A4:2, C4:-7})
#pp3 = p3.subs({A2:5, B2:-3, A4:2, C4:-7})
#pp4 = p4.subs({A2:5, B2:-3, A4:2, C4:-7})
#pp5 = p5.subs({A2:5, B2:-3, A4:2, C4:-7})

## ## HERE WE CONCLUDE THE FIRST PART OF THE COMPUTATION:
##                         In case C2 = 0, 
## IT IS POSSIBLE TO HAVE THREE ALIGNMENTS:
### HENCE IT IS POSSIBLE TO HAVE THE ALIGNMENTS (1, 2, 3), (1, 4, 5), (1, 6, 7)
### WITH NO OTHER ALIGNMENTS.

## Now we want to see if it is possible to have more then three alignments
## (We continue to assume C2 != 0)

## Recall that cb_alt is our cubic.

## Recall that r3_1 (= r3_2) is the line through p6 and p7.

## Up to a permutation of the indices of the points, if there is another 
## alignment among the eigenpoints, p6 must be on the line p2+p4. Hence
## we can find it, since is the intersection or p2+p4 and r3.

r24 = det(matrix([p2, p4, (x, y, z)]))

E1 = matrix([[r3_1.coefficient(xx) for xx in [x, y, z]], \
             [r24.coefficient(xx) for xx in [x, y, z]]]).minors(2)


p6 = vector(S, (E1[2], -E1[1], E1[0]))

kJ = S.ideal(matrix([[x, y, z], \
        [cb_alt.derivative(x), cb_alt.derivative(y), \
         cb_alt.derivative(z)]]).minors(2))

kJ = kJ.saturation(S.ideal(z, y))[0]  ##p1
kJ = kJ.saturation(S.ideal(p2[0]*y-p2[1]*x, p2[0]*z-p2[2]*x, \
                           p2[1]*z-p2[2]*y))[0] ## p2
kJ = kJ.saturation(S.ideal(p3[0]*y-p3[1]*x, p3[0]*z-p3[2]*x, \
                           p3[1]*z-p3[2]*y))[0] ## p3
kJ = kJ.saturation(S.ideal(p4[0]*y-p4[1]*x, p4[0]*z-p4[2]*x, \
                           p4[1]*z-p4[2]*y))[0] ## p4
kJ = kJ.saturation(S.ideal(p5[0]*y-p5[1]*x, p5[0]*z-p5[2]*x, \
                           p5[1]*z-p5[2]*y))[0] ## p5
kJ = kJ.saturation(S.ideal(matrix([p1, p2]).minors(2)))[0]
kJ = kJ.saturation(S.ideal(matrix([p1, p3]).minors(2)))[0]
kJ = kJ.saturation(S.ideal(matrix([p1, p4]).minors(2)))[0]
kJ = kJ.saturation(S.ideal(matrix([p1, p5]).minors(2)))[0]


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
## CONCLUSION (for the case C2 = 0) when the three deltas are zero, 
## (hence rank of the matrix of the five points is 8), we have 
## also delta2 = 0 and we have the collineariries (1, 2, 3) (1, 4, 5).
## We have a sub-case in which there is the further collinearity (1, 6, 7).
## No other collinearities among the 7 eigenpoints are possible.




