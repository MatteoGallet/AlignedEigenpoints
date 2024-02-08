## HERE WE PROVE THAT CONFIGURATION (5) IS POSSIBLE AND MUCH MORE.

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

P1 = vector(S, (A1, B1, C1))
P2 = vector(S, (A2, B2, C2))
P3 = u1*P1+u2*P2
P4 = vector(S, (A4, B4, C4))
P5 = v1*P1+v2*P4
P6 = w1*P2+w2*P4



## we consider the following aligmnents:
## (1, 2, 3), (1, 4, 5), (2, 4, 6), (3, 4, 7)

## In order to have configuration (5), we have to consider these polynomials 
## which must be zero:
## delta1(P1, P2, P4), 
## delta1(P2, P1, P4), 
## delta1(P3, P1, P4).quo_rem(u2)[0], 
## delta2(P4, P1, P5, P6, P2)
## This last condition implies that P7 is collinear with P4 and P3

## moreover, delta1(P3, P1, P4) is divisible by u2:
assert(delta1(P3, P1, P4).quo_rem(u2)[1] == S(0))

## hence we can consider the following ideal: 
J = S.ideal(delta1(P1, P2, P4), delta1(P2, P1, P4), \
            delta1(P3, P1, P4).quo_rem(u2)[0], delta2(P4, P1, P5, P6, P2),\
)

## we forgot the condition delta1(P3, P2, P4), but it is not 
## necessary:
assert(J == J+S.ideal(delta1(P3, P2, P4)))


## We saturate J and we get that J is the ideal generated 
## by (P1|P4) and (P2|P4):

J = J.saturation(matrix([P1, P2, P4]).det())[0]
assert(J == S.ideal(scalarProd(P2, P4), scalarProd(P1, P4)))

## So we define P4 in this way and we re-write the points:

P1 = vector(S, (A1, B1, C1))
P2 = vector(S, (A2, B2, C2))
P3 = u1*P1+u2*P2
P4 = vector(S, list(wedgeProd(P1, P2)))
P5 = v1*P1+v2*P4
P6 = w1*P2+w2*P4


J = S.ideal(delta1(P1, P2, P4), delta1(P2, P1, P4), \
            delta1(P3, P1, P4).quo_rem(u2)[0], delta2(P4, P1, P5, P6, P2))
## J is (0):
assert(J == S.ideal(0))

M = matrixEigenpoints([P1, P2, P3, P4, P5, P6])

## In order to have a cubic with P1, ..., P6 eigenpoints, 
## the matrix M must have rank <= 9.
## We compute one minor of M of order 10:

Nx = M.matrix_from_rows([0, 1, 3, 4, 6, 7, 9, 10, 12, 13])

## the next computation requires about 20'. Here we omit it and we
## give the result:
## ttA = cputime()
## dn = Nx.det()
## print(cputime()-ttA)

## 

## we get the polynomial dn defined below:

dn = expand((27) * A2 * A1 * w1 * v1 * u2^2 * u1^2 * w2^3 * v2^3 * (u1*A1 + u2*A2) * (-w2*C1*B2 + w2*B1*C2 + w1*A2) * (-v2*C1*B2 + v2*B1*C2 + v1*A1) * (B1^2*A2^2 + C1^2*A2^2 - 2*A1*B1*A2*B2 + A1^2*B2^2 + C1^2*B2^2 - 2*A1*C1*A2*C2 - 2*B1*C1*B2*C2 + A1^2*C2^2 + B1^2*C2^2)^7 * (-2*u1*v1*w2*A1^3*A2 - 2*u1*v1*w2*A1*B1^2*A2 - 2*u1*v1*w2*A1*C1^2*A2 + 2*u1*v2*w1*A1^2*A2^2 - 2*u2*v1*w2*A1^2*A2^2 + u1*v2*w1*B1^2*A2^2 - u2*v1*w2*B1^2*A2^2 + u1*v2*w1*C1^2*A2^2 - u2*v1*w2*C1^2*A2^2 + 2*u2*v2*w1*A1*A2^3 - 2*u1*v1*w2*A1^2*B1*B2 - 2*u1*v1*w2*B1^3*B2 - 2*u1*v1*w2*B1*C1^2*B2 + 2*u1*v2*w1*A1*B1*A2*B2 - 2*u2*v1*w2*A1*B1*A2*B2 + 2*u2*v2*w1*B1*A2^2*B2 + u1*v2*w1*A1^2*B2^2 - u2*v1*w2*A1^2*B2^2 + 2*u1*v2*w1*B1^2*B2^2 - 2*u2*v1*w2*B1^2*B2^2 + u1*v2*w1*C1^2*B2^2 - u2*v1*w2*C1^2*B2^2 + 2*u2*v2*w1*A1*A2*B2^2 + 2*u2*v2*w1*B1*B2^3 - 2*u1*v1*w2*A1^2*C1*C2 - 2*u1*v1*w2*B1^2*C1*C2 - 2*u1*v1*w2*C1^3*C2 + 2*u1*v2*w1*A1*C1*A2*C2 - 2*u2*v1*w2*A1*C1*A2*C2 + 2*u2*v2*w1*C1*A2^2*C2 + 2*u1*v2*w1*B1*C1*B2*C2 - 2*u2*v1*w2*B1*C1*B2*C2 + 2*u2*v2*w1*C1*B2^2*C2 + u1*v2*w1*A1^2*C2^2 - u2*v1*w2*A1^2*C2^2 + u1*v2*w1*B1^2*C2^2 - u2*v1*w2*B1^2*C2^2 + 2*u1*v2*w1*C1^2*C2^2 - 2*u2*v1*w2*C1^2*C2^2 + 2*u2*v2*w1*A1*A2*C2^2 + 2*u2*v2*w1*B1*B2*C2^2 + 2*u2*v2*w1*C1*C2^3))

## some factors of dn are specific of the choice of the minor of M
## We consider the last factor:
fdn = dn.factor()
ftC = fdn[-1][0]

## In the computations below we shall show that the rank of M is <=9
## iff the polynomial ftC is zero.
## One way to see this, is to consider the ideal of all the order 10 
## minors of M, but this computation requires too much time.
## Hence we assume that the point A4 is (1,0,0) or (1, ii, 0). In this 
## case the computation of the ideal of all the order 10 minors of the 
## corresponding matrix M is easy to manipulate.

## case A4 = (1, 0, 0)

## Since
## (P1|P4) = 0, (P2|P4) = 0

P1 = vector(S, (0, B1, C1))
P2 = vector(S, (0, B2, C2))
P3 = u1*P1+u2*P2
P4 = vector(S, (1, 0, 0))
P5 = v1*P1+v2*P4
P6 = w1*P2+w2*P4

M1 = matrixEigenpoints([P1, P2, P3, P4, P5, P6])

## The matrix M1 has the 9th row equals to (0,1,0,...,0)
## the 10th row equals to (0,0,0,0,1,0,...,0) and the 11st 
## row given by: (0, 0,..., 0). Hence we can extract from 
## M1 a matrix N1 which do not have the rows 9, 10, 11 and 
## do not have the columns 1 and 4. All the order 10 minors 
## of M1 are 0 iff all the order 8 minors of N1 are 0.

N1 = M1.matrix_from_rows_and_columns([0,1,2,3,4,5,6,7,8,12,13,14,15,16,17], \
     [0,2,3,5,6,7,8,9])

## we can see that the following rows of N1 are linearly dependent:
## row 0 and row 1:
assert(N1.matrix_from_rows([0, 1]).rank() == 1)
## row 3 and row 4:
assert(N1.matrix_from_rows([3, 4]).rank() == 1)
## row 6 and row 7:
assert(N1.matrix_from_rows([6, 7]).rank() == 1)
## row 9 and row 10 and row 11:
assert(N1.matrix_from_rows([9, 10, 11]).rank() == 2)
## row 12 and row 13 and row 14:
assert(N1.matrix_from_rows([12, 13, 14]).rank() == 2)

## hence, in order to compute all the order 8 minors of N1, 
## we can skip several submatrices.

## Here we construct the list of all the rows of 8 elements
## which have to be considered:

L1 = list(Combinations(15,8))

LL1 = []
for lx in L1:
    ll = Set(lx)
    if Set([0,1]).issubset(ll) or Set([3, 4]).issubset(ll) or Set([6, 7]).issubset(ll) or Set([9, 10, 11]).issubset(ll) or Set([12, 13, 14]).issubset(ll):
        continue
    else:
        LL1.append(lx)

## LL1 contains 1362 elements:
assert(len(LL1) == 1362)

## The polynomial ftC constructed above should appear in the order 8
## minors of N1 (when specialized with the condition A1 = 0, A2 = 0)
## We verify this and we collect the order 8 minors of N1 divided by
## the polynomial ftC specialized:
ftCs = ftC.subs({A1:0, A2:0})

## the following computation requires < 20 seconds
print("11 seconds of computations...")
JJ = []
for nr in LL1:
    NN = N1.matrix_from_rows(nr)
    dt = NN.det()
    dvs = dt.quo_rem(ftCs)
    if dvs[1] != 0:
        print("Unexpacted situation. The minor is not a multiple of ftCs")
        print("Do not trust to the next computations!")
        ## al caso lanciare un errore
    else:
        JJ.append(dvs[0])

print("... end of the computation")

## We define the ideal generated by JJ and we saturate it:
JJ = S.ideal(JJ)
JJ = JJ.saturation(u1*u2*v1*v2*w1*w2)[0]
JJ = JJ.saturation(S.ideal(matrix([P1, P2]).minors(2)))[0]

## we get that JJ = (1), so the order 8 minors of N1 
## are 0 iff ftCs is zero.

assert(JJ == S.ideal(S(1)))


## then we have to consider the case in which P4 is (1, i, 0)
## But in this case we have that:
## If Q1, Q2 are points of the plane, if Q4 = wedgeProd(Q1, Q2)
## and if Q4 is on the isotropic conic, i.e. (Q4|Q4) = 0, then
## Q1, Q2, Q4 are aligned:
Q1 = vector(S, (A1, B1, C1))
Q2 = vector(S, (A2, B2, C2))
Q4 = wedgeProd(Q1, Q2)
assert(det(matrix([Q1, Q2, Q4])) == scalarProd(Q4, Q4))

## from this we get that the case P4 = (1, i, 0) does not 
## need to be considered.

### First conclusion:
### Configuration (5) is possible iff ftC is zero

##################
##############
## we go back to the general case:

P1 = vector(S, (A1, B1, C1))
P2 = vector(S, (A2, B2, C2))
P3 = u1*P1+u2*P2
P4 = vector(S, list(wedgeProd(P1, P2)))
P5 = v1*P1+v2*P4
P6 = w1*P2+w2*P4
#
## From the above computation we have that P4 cannot be a point 
## on the isotropic conic (indeed we so that if P4 is on the 
## isotropic conic, then P1, P2, P4 are collinear). Hence 
## (P4|P4) is not zero.
## we have that the condition ftC, which is the condition that 
## implies that P1, P2, P3, P4, P5, P6 in configuration (5) are 
## eigenpoints, can be expressed by:
##
## (
## scalarProd(P1,P3)*(scalarProd(P2,P6)*scalarProd(P4, P5)-\
##                    scalarProd(P2,P5)*scalarProd(P4, P6))+\
## scalarProd(P2,P3)*(scalarProd(P1,P6)*scalarProd(P4, P5)-\
##                    scalarProd(P1,P5)*scalarProd(P4, P6))
## )/scalarProd(P4, P4)

ftC1 = scalarProd(P1,P3)*(scalarProd(P2,P6)*scalarProd(P4, P5)-\
                    scalarProd(P2,P5)*scalarProd(P4, P6))+\
 scalarProd(P2,P3)*(scalarProd(P1,P6)*scalarProd(P4, P5)-\
                    scalarProd(P1,P5)*scalarProd(P4, P6))

assert(ftC1 == scalarProd(P4, P4)*ftC)

## The condition ftC1=0 gives that in order to have 6 points in configuration
## (5) we can choose P1 and P2 in an arbitrary way, P4 = (P1||P2)
## P5 on the line P1+P4, P6 on the line P2+P4. Then P3 is a point on the 
## line P1+P2 determined by a linear equation in u1 and u2 given by ftC1=0.

## Solving ftC1 = 0 w.r.t. u1 and u2, we get that the following U1 
## and U2 are the solution:

U1 = scalarProd(P1,P2)*scalarProd(P2,P6)*scalarProd(P4,P5)-scalarProd(P1,P2)*scalarProd(P2,P5)*scalarProd(P4,P6)+\
scalarProd(P2,P2)*scalarProd(P1,P6)*scalarProd(P4,P5)-scalarProd(P2,P2)*scalarProd(P1,P5)*scalarProd(P4,P6)

U2 = -scalarProd(P1,P1)*scalarProd(P2,P6)*scalarProd(P4,P5)+scalarProd(P1,P1)*scalarProd(P2,P5)*scalarProd(P4,P6)-\
scalarProd(P1,P2)*scalarProd(P1,P6)*scalarProd(P4,P5)+scalarProd(P1,P2)*scalarProd(P1,P5)*scalarProd(P4,P6)

assert(ftC1.subs({u1:U1, u2:U2})==S(0))

## Hence P3 is given by: P3 = U1*P1+U2*P2
## 
## It should be possible to find a formula for P7.
##
## CONCLUSION: CONFIGURATION (5) IS POSSIBLE AND THE CUBICS ARE PARAMETRIZED BY
## P^2xP^2xP^1xP^1, I.E. ARE OF DIMENSION 6

