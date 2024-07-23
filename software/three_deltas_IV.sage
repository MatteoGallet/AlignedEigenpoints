## Here we to consider the case in which 
## delta1(P1, P2, P4) = 0,
## delta1b(P1, P2, P3) = 0,
## delta1b(P1, P4, P5) = 0.

## This case is given in proposition prop:d2_6allin

## and corresponds to one of the cases of \Phi(P1,..., P5) with rank 8

## There are other files which complete the solution of the problem.
## They are:
########################################## 
## three_deltas_I.sage
## three_deltas_II.sage
## three_deltas_III.sage
## three_deltas_IV.sage

## We start as in the file rank_8_1.sage.

## This file is parallel to the file three_deltas_III.sage
## and obtain the same results. The only difference is that 
## here we consider a solution to the equation delta1(P1, P2, P4)=0
## which is not considered in the file three_deltas_III.sage

## Hence, again:
## The result of this file is that if in this situation we have that 
## there is a further collinearity between p2, p4, p6, then we have 
## a configuration o type (8).
## at the end of the file we verify the conditions of orthogonality
## of the line containing the eigenpints.

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

## We define the points p1, ..., p5 in such a way that 
## delta1(p1, p2, p4), delta1b(p1, p2, p3), delta1b(p1, p4, p5)
## are zero and are zero in the most general situation.
## In order to have delta1(P1, P2, P4) = 0, since this is:
## B2*B4 + C2*C4
## we have to consider two cases: 
## {B2: 0, C4: 0}
## or 
## {B4:-C2, B2:C4}
########################################
##### Here we consider {B2: 0, C4: 0}
########################################
##
## The case {B4:-C2, B2:C4} is considered in the file three_deltas_III.sage



########################################
##### Here we consider the other case: {B2: 0, C4: 0}
########################################

## we redefine p1, ..., in such a way that the three delta are 0
p1 = P1.subs({B2: 0, C4: 0})
p2 = P2.subs({B2: 0, C4: 0})
p3 = P3.subs({B2: 0, C4: 0})
p4 = P4.subs({B2: 0, C4: 0})
p5 = P5.subs({B2: 0, C4: 0})
ss2 = {u1: delta1b(p1, p2, p3).coefficient(u2), u2:-delta1b(p1, p2, p3).coefficient(u1),\
       v1: delta1b(p1, p4, p5).coefficient(v2), v2:-delta1b(p1, p4, p5).coefficient(v1)}
p1 = p1.subs(ss2)
p2 = p2.subs(ss2)
p3 = p3.subs(ss2)
p4 = p4.subs(ss2)
p5 = p5.subs(ss2)

## Now delta1(p1, p2, p4), delta1b(p1, p2, p3), delta1b(p1, p4, p5) are 0:

assert(delta1(p1, p2, p4)==0)
assert(delta1b(p1, p2, p3)==0)
assert(delta1b(p1, p4, p5)==0)


## It holds: 
## delta2(p1, p2, p3, p4, p5) = 0:

assert(delta2(p1, p2, p3, p4, p5) == 0)


## an example shows that p6 and p7 are not aligned with p1.
## Here is an example. The 7 eigenpoints are such that p6 and p7 
## are not aligned.

## ss3 = {A2:1, A4:-5, B4:-3, C2:7}
## pp1 = p1.subs(ss3)
## pp2 = p2.subs(ss3)
## pp3 = p3.subs(ss3)
## pp4 = p4.subs(ss3)
## pp5 = p5.subs(ss3)
## cb = cubic_from_matrix(matrixEigenpoints([pp1, pp2, pp3, pp4, pp5]).\
##                 stack(matrix([[2, 3, 4, 5, 6, 7, 8, 9, 1, 2]])))
## eigenpoints(cb)


## We want to see if it is possible to have an eigenpoint which is collinear 
## with p2 and p4. We define the point w1*p2+w2*p4 and we want to see
## if the matrix \Phi(p1, ..., p5, w1*p2+w2*p4) can have rank 9.

mm2 = matrixEigenpoints([p1, p2, p3, p4, p5, w1*p2+w2*p4])

## in order to extract the order 10 minors of mm2, we use some strategies:
## all the possible rows to get the order 10 minors of mm2 are:

possible_rows = list(Combinations(18, 10))

## but many of them are unnecessary: 
## for instance if a row contains the rows [0, 1, 2] or [3, 4, 5], ...
## of mm2, the corresponding determinant is 0, moreover 
## \Phi(p1, p2, p3, p4, p5) has rank 8, so if a row contains only one 
## row of position 15, 16, 17, again the determinant is 0.
## Here is a function which says if a row is a "good row".
## 
def is_good_row(rg):
    det_sure_0 = [[0, 1, 2], [3, 4, 5], [6, 7, 8], \
                  [9, 10, 11], [12, 13, 14], [15, 16, 17]]
    for ds in det_sure_0:
        if Set(ds).issubset(Set(rg)) or \
           len(Set(rg).difference(Set(range(15)))) < 2:
            return(false)
    return(true)

## here we extract the "good rows" i.e. the rows that give an order 
## 10 minors whose determinant has a chnce to be not 0.

print("\nTime of computation: 10 seconds\n")
sleep(1)
good_rows = []
for rg in possible_rows:
    if is_good_row(rg):
        good_rows.append(rg)

print("done\n")

## Here we compute the ideal generated by some minors of mm2 of order 10.
## We do not need to compute all the possible minors (which are more then
## 8000, but only some of them. If the informations we get from the ideal 
## they generate are enough, we can conclude.


JJ0 = []
flag = 0
for gr in good_rows:
    JJ0.append(mm2.matrix_from_rows(gr).det())
    flag += 1
    if flag > 100:
        break

## Study of the ideal:


JJ = S.ideal(JJ0)
JJ = JJ.saturation(w1*w2)[0]
JJ = JJ.saturation(matrix([p1, p2, p4]).det())[0]
JJ = JJ.saturation(S.ideal(matrix([p1, p3]).minors(2)))[0].\
        saturation(S.ideal(matrix([p1, p5]).minors(2)))[0]
JJ = JJ.saturation(S.ideal(matrix([p2, p3]).minors(2)))[0]
JJ = JJ.saturation(S.ideal(matrix([p4, p5]).minors(2)))[0]

gJJ = JJ.groebner_basis()

sJJ = [test.get_sqrfree(mm) for mm in gJJ]

pd = S.ideal(sJJ).radical().primary_decomposition()

## we get only one ideal which is principal and is generated by:
## 2*w1*A2^2*A4 + w1*C2^2*A4 + 2*w2*A2*A4^2 + w2*A2*B4^2


assert(len(pd) == 1)

FF = 2*w1*A2^2*A4 + w1*C2^2*A4 + 2*w2*A2*A4^2 + w2*A2*B4^2


assert(pd[0] == S.ideal(FF))



## We define therefore p6 in with w1 and w2 such that give FF = 0

## First of all we verify that it is not possible to have 
## FF.coefficient(w1) = 0 and FF.coefficient(w2) = 0:

JJ = S.ideal(FF.coefficient(w1), FF.coefficient(w2))
JJ = JJ.saturation(S.ideal(matrix([p1, p3]).minors(2)))[0]
JJ = JJ.saturation(S.ideal(matrix([p1, p5]).minors(2)))[0]
JJ = JJ.saturation(S.ideal(matrix([p2, p3]).minors(2)))[0]
## it is not possible to have FF.coefficient(w1) = 0 and 
## FF.coefficient(w2) = 0:
assert(JJ == S.ideal(1))


p6 = (w1*p2+w2*p4).subs({w2: pd[0].gens()[0].coefficient(w1), \
                         w1: -pd[0].gens()[0].coefficient(w2)})

## p6 is aligned with p2 and p4 and also with p3 and p5

assert(det(matrix([p2, p4, p6])) == 0)
assert(det(matrix([p3, p5, p6])) == 0)

## If we define p7 as the point of intersection of p2+p5 and p3+p4, 
## we are going to verify that p7 is an eigenpoint.

p7 = vector(S, intersect_lines(p3, p4, p2, p5))

## here we show that p7 is an eigenpoint.

## We consider the matrix Phi(p1, p2, p3, p4, p5, p6, p7):

mm3 = matrixEigenpoints([p1, p2, p3, p4, p5, p6, p7])

## then we extract from mm3 the rows of position: 
## 0, 3, 4, 6, 7, 9, 12, 13, 15

mm4 = mm3.matrix_from_rows([0, 3, 4, 6, 7, 9, 12, 13, 15])

## and we verify that phi_1(p7) is a linear combination of 
## the rows of mm4:

assert(det(mm4.stack(matrix([phi_p(p7)[0]]))) == 0)

## then we verify that phi_2(p7) is a linear combination of 
## the rows of mm4:

assert(det(mm4.stack(matrix([phi_p(p7)[1]]))) == 0)

## and, finally,  we verify that phi_3(p7) is a linear combination of 
## the rows of mm4:

assert(det(mm4.stack(matrix([phi_p(p7)[2]]))) == 0)

## the above computations show that p7 is an eigenpoint.

## Here we have the alignments of the seven eigenpoints:
##
## [(1, 2, 3), (1, 4, 5), (2, 4, 6), (2, 5, 7), (3, 4, 7), (3, 5, 6)]

assert(allignments([p1, p2, p3, p4, p5, p6, p7]) == \
[(1, 2, 3), (1, 4, 5), (2, 4, 6), (2, 5, 7), (3, 4, 7), (3, 5, 6)])

## so the points are in configuration (8) 
## Note that  p1, p6, p7 are not aligned, 
## but here we have delta2(p1, p2, p3, p4, p5) = 0.


## so the points are in configuration (8) 
## Note that  p1, p6, p7 are not aligned, 
## but here we have delta2(p1, p2, p3, p4, p5) = 0.
##

### we have the following orthogonalities:

## p1+p4 ort to p1+p2
## p2+p4 ort to p3+p5
## p2+p5 ort to p3+p4

assert(scalarProd(wedgeProd(p1, p4), wedgeProd(p1, p2))==0)
assert(scalarProd(wedgeProd(p2, p4), wedgeProd(p3, p5))==0)
assert(scalarProd(wedgeProd(p2, p5), wedgeProd(p3, p4))==0)

## here we verify that p4 can be obtained from p2, p5, p7 via the 
## formula 
## p4 == (p2||p3)(p2|p5)(p3|p5)-(p2|p3)(p2||p5)(p3|p5)+
##       (p2|p3)(p2|p5)(p3||p5)

q4 = wedgeProd(p2, p3)*scalarProd(p2, p5)*scalarProd(p3, p5)\
        -scalarProd(p2, p3)*wedgeProd(p2, p5)*scalarProd(p3, p5)\
        +scalarProd(p2, p3)*scalarProd(p2, p5)*wedgeProd(p3, p5)

assert(matrix([p4, q4]).minors(2) == [0, 0, 0])



## The unique cubic which have p1, ..., p7 as eigenpoints is singular in p1
## as follows from proposition 4.18 (\\ref{proposition:P1_sing})
