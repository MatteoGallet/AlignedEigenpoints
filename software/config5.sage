## HERE WE PROVE THAT CONFIGURATION (5) IS POSSIBLE AND MUCH MORE.

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

P1 = vector(S, (A1, B1, C1))
P2 = vector(S, (A2, B2, C2))
P3 = u1*P1+u2*P2
P4 = vector(S, (A4, B4, C4))
P5 = v1*P1+v2*P4
P6 = w1*P2+w2*P4


## we consider the following aligmnents:
## (1, 2, 3), (1, 4, 5), (1, 6, 7), (2, 4, 6)

print("First computation: about 5''")
sleep(1)


## In order to have configuration (5), we have to consider these polynomials 
## which must be zero:
## delta1(P2, P1, P4), 
## delta1(P4, P1, P2), 
## delta1(P6, P1, P2).quo_rem(w2)[0], 
## delta2(P1, P2, P3, P4, P5)
## This last condition implies that P7 is collinear with P1 and P6

## indeed, delta1(P6, P1, P2) is divisible by w2:
assert(delta1(P6, P1, P2).quo_rem(w2)[1] == S(0))

## hence we can consider the following ideal:
J = S.ideal(delta1(P2, P1, P4), delta1(P4, P1, P2), \
delta1(P6, P1, P2).quo_rem(w2)[0], delta2(P1, P2, P3, P4, P5))

## We saturate J and we get that J is the ideal generated 
## by (P1|P2) and (P1|P4):

J = J.saturation(matrix([P1, P2, P4]).det())[0]
assert(J == S.ideal(scalarProd(P1, P4), scalarProd(P1, P2)))

## Moreover, we have that s16 is in J:

assert(scalarProd(P1, P6) in J)

## So we define P1 in this way and we re-write the points:

P2 = vector(S, (A2, B2, C2))
P4 = vector(S, (A4, B4, C4))
P1 = vector(S, list(wedgeProd(P2, P4)))
P3 = u1*P1+u2*P2
P5 = v1*P1+v2*P4
P6 = w1*P2+w2*P4


J = S.ideal(delta1(P2, P1, P4), delta1(P4, P1, P2), \
            delta1(P6, P1, P2), delta2(P1, P2, P3, P4, P5))
## J is (0):
assert(J == S.ideal(0))

## We have several orthogonalities among the lines:

## P1+P2 ort P4+P6
## P1+P6 ort P2+P4
## P1+P4 ort P2+P6

assert(scalarProd(wedgeProd(P1, P2), wedgeProd(P4, P6)) == S(0))
assert(scalarProd(wedgeProd(P1, P6), wedgeProd(P2, P4)) == S(0))
assert(scalarProd(wedgeProd(P1, P4), wedgeProd(P2, P6)) == S(0))

## hence the line P2+P4+P6 is orthogonal to P1+P2, to P1+P4 and to P1+P6


M = matrixEigenpoints([P1, P2, P3, P4, P5, P6])

## In order to have a cubic with P1, ..., P6 eigenpoints, 
## the matrix M must have rank <= 9.
## We compute one minor of M of order 10:

Nx = M.matrix_from_rows([0, 1, 3, 4, 6, 7, 9, 10, 12, 15])

## the next computation requires about 8'. Here we omit it and we
## give the result:
## ttA = cputime()
## dn = Nx.det()
## print(cputime()-ttA)

## 


## we get the polynomial dn defined below:

dn = expand((-27) * A4 * A2 * w2 * w1 * v2 * v1 * u2^2 * u1^2 * (-C2*B4 + B2*C4) * (-B2*A4 + A2*B4) * (A2^2 + B2^2 + C2^2) * (-u1*C2*B4 + u1*B2*C4 + u2*A2) * (-C2*A4^2 - C2*B4^2 + A2*A4*C4 + B2*B4*C4) * (B2^2*A4^2 + C2^2*A4^2 - 2*A2*B2*A4*B4 + A2^2*B4^2 + C2^2*B4^2 - 2*A2*C2*A4*C4 - 2*B2*C2*B4*C4 + A2^2*C4^2 + B2^2*C4^2)^6 * (-2*u2*v1*w1*A2^3*A4 - 2*u2*v1*w1*A2*B2^2*A4 - 2*u2*v1*w1*A2*C2^2*A4 + 2*u1*v2*w1*A2^2*A4^2 - 2*u2*v1*w2*A2^2*A4^2 + u1*v2*w1*B2^2*A4^2 - u2*v1*w2*B2^2*A4^2 + u1*v2*w1*C2^2*A4^2 - u2*v1*w2*C2^2*A4^2 + 2*u1*v2*w2*A2*A4^3 - 2*u2*v1*w1*A2^2*B2*B4 - 2*u2*v1*w1*B2^3*B4 - 2*u2*v1*w1*B2*C2^2*B4 + 2*u1*v2*w1*A2*B2*A4*B4 - 2*u2*v1*w2*A2*B2*A4*B4 + 2*u1*v2*w2*B2*A4^2*B4 + u1*v2*w1*A2^2*B4^2 - u2*v1*w2*A2^2*B4^2 + 2*u1*v2*w1*B2^2*B4^2 - 2*u2*v1*w2*B2^2*B4^2 + u1*v2*w1*C2^2*B4^2 - u2*v1*w2*C2^2*B4^2 + 2*u1*v2*w2*A2*A4*B4^2 + 2*u1*v2*w2*B2*B4^3 - 2*u2*v1*w1*A2^2*C2*C4 - 2*u2*v1*w1*B2^2*C2*C4 - 2*u2*v1*w1*C2^3*C4 + 2*u1*v2*w1*A2*C2*A4*C4 - 2*u2*v1*w2*A2*C2*A4*C4 + 2*u1*v2*w2*C2*A4^2*C4 + 2*u1*v2*w1*B2*C2*B4*C4 - 2*u2*v1*w2*B2*C2*B4*C4 + 2*u1*v2*w2*C2*B4^2*C4 + u1*v2*w1*A2^2*C4^2 - u2*v1*w2*A2^2*C4^2 + u1*v2*w1*B2^2*C4^2 - u2*v1*w2*B2^2*C4^2 + 2*u1*v2*w1*C2^2*C4^2 - 2*u2*v1*w2*C2^2*C4^2 + 2*u1*v2*w2*A2*A4*C4^2 + 2*u1*v2*w2*B2*B4*C4^2 + 2*u1*v2*w2*C2*C4^3))

## some factors of dn are specific of the choice of the minor of M
## We consider the last factor:
fdn = dn.factor()
ftC = fdn[-1][0]
## In the computations below we shall show that the rank of M is <=9
## iff the polynomial ftC is zero.
## One way to see this, is to consider the ideal of all the order 10 
## minors of M, but this computation requires too much time.
## Hence we assume that the point P1 is (1,0,0) or (1, ii, 0). In this 
## case the computation of the ideal of all the order 10 minors of the 
## corresponding matrix M is easy to manipulate.

## case P1 = (1, 0, 0)

## Since
## (P1|P2) = 0, (P1|P4) = 0

P1 = vector(S, (1, 0, 0))
P2 = vector(S, (0, B2, C2))
P3 = u1*P1+u2*P2
P4 = vector(S, (0, B4, C4))
P5 = v1*P1+v2*P4
P6 = w1*P2+w2*P4

M1 = matrixEigenpoints([P1, P2, P3, P4, P5, P6])



## The matrix M1 has the 0th row equals to (0,1,0,...,0)
## the 1st row equals to (0,0,0,0,1,0,...,0) and the 2nd 
## row given by: (0, 0,..., 0). Hence we can extract from 
## M1 a matrix N1 which do not have the rows 0, 1, 2 and 
## do not have the columns 1 and 4. All the order 10 minors 
## of M1 are 0 iff all the order 8 minors of N1 are 0.

N1 = M1.matrix_from_rows_and_columns([3,4,5,6,7,8,9,10,11,12,13,14,15,16,17], \
     [0,2,3,5,6,7,8,9])



## we can see that the following rows of N1 are linearly dependent:
## row 0 and row 1:
assert(N1.matrix_from_rows([0, 1]).rank() == 1)
## row 6 and row 7:
assert(N1.matrix_from_rows([6, 7]).rank() == 1)
## row 12 and row 13:
assert(N1.matrix_from_rows([12, 13]).rank() == 1)
## row 3 and row 4 and row 5:
assert(N1.matrix_from_rows([3, 4, 5]).rank() == 2)
## row 9 and row 10 and row 11:
assert(N1.matrix_from_rows([9, 10, 11]).rank() == 2)
## row 0 and row 12 and row 13:
assert(N1.matrix_from_rows([0, 12, 13]).rank() == 2)

## hence, in order to compute all the order 8 minors of N1, 
## we can skip several submatrices.

## Here we construct the list of all the rows of 8 elements
## which have to be considered:

L1 = list(Combinations(15,8))

LL1 = []
for lx in L1:
    ll = Set(lx)
    if Set([0,1]).issubset(ll) or Set([6, 7]).issubset(ll) or Set([12, 13]).issubset(ll) or Set([3, 4, 5]).issubset(ll) or Set([9, 10, 11]).issubset(ll) or Set([0, 12, 13]).issubset(ll):
        continue
    else:
        LL1.append(lx)

## LL1 contains 1362 elements:
assert(len(LL1) == 1362)

## The polynomial ftC constructed above should appear in the order 8
## minors of N1 (when specialized with the condition A2 = 0, A4 = 0)
## We verify this and we collect the order 8 minors of N1 divided by
## the polynomial ftC specialized:
ftCs = ftC.subs({A2:0, A4:0})

## the following computation requires < 20 seconds
print("16 seconds of computations...")
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
JJ = JJ.saturation(S.ideal(matrix([P2, P4]).minors(2)))[0]

## we get that JJ = (1), so the order 8 minors of N1 
## are 0 iff ftCs is zero.

assert(JJ == S.ideal(S(1)))


## then we have to consider the case in which P1 is (1, i, 0)
## But in this case we have that:
## If Q2, Q4 are points of the plane, if Q1 = wedgeProd(Q2, Q4)
## and if Q1 is on the isotropic conic, i.e. (Q1|Q1) = 0, then
## Q1, Q2, Q4 are aligned:
Q4 = vector(S, (A4, B4, C4))
Q2 = vector(S, (A2, B2, C2))
Q1 = wedgeProd(Q2, Q4)
assert(det(matrix([Q1, Q2, Q4])) == scalarProd(Q1, Q1))

## from this we get that the case P1 = (1, i, 0) does not 
## need to be considered.

### First conclusion:
### Configuration (5) is possible iff ftC is zero

##################
##############
## we go back to the general case:

P2 = vector(S, (A2, B2, C2))
P4 = vector(S, (A4, B4, C4))
P1 = vector(S, list(wedgeProd(P2, P4)))
P3 = u1*P1+u2*P2
P5 = v1*P1+v2*P4
P6 = w1*P2+w2*P4

#
## From the above computations we have that P1 cannot be a point 
## on the isotropic conic (indeed we sow that if P1 is on the 
## isotropic conic, then P1, P2, P4 are collinear). Hence 
## (P1|P1) is not zero.
## we have that the condition ftC, which is the condition that 
## implies that P1, P2, P3, P4, P5, P6 in configuration (5) are 
## eigenpoints, can be expressed by:
## (
## scalarProd(P2,P6)*(scalarProd(P4,P5)*scalarProd(P1, P3)-\
##                    scalarProd(P4,P3)*scalarProd(P1, P5))+\
## scalarProd(P4,P6)*(scalarProd(P2,P5)*scalarProd(P1, P3)-\
##                    scalarProd(P2,P3)*scalarProd(P1, P5))
## )/scalarProd(P1, P1)



ftC1 = scalarProd(P2,P6)*(scalarProd(P4,P5)*scalarProd(P1, P3)-\
                    scalarProd(P4,P3)*scalarProd(P1, P5))+\
 scalarProd(P4,P6)*(scalarProd(P2,P5)*scalarProd(P1, P3)-\
                    scalarProd(P2,P3)*scalarProd(P1, P5))

assert(ftC1 == -2*scalarProd(P1, P1)*ftC)

## since (P1|P1) != 0, we have that ftC = 0 <==> ftC1 = 0.

## Hence P1, ..., P6 in config (5) are eigenpoints iff 
## s26*(s45*s13 - s34*s15) + s46*(s25*s13 - s23*s15)=0

## If we expand this w.r.t. P3 = u1*P1+u2*P2, we get:

## U1 = s12*(s26*s45+s46*s25)-s26*s15*s24-s46*s15*s22
## U2 = s11*(s26*s45+s46*s25)-s26*s15*s14-s46*s15*s12

## and ftC is zero iff u1 = U1, u2 = -U2:


U1 = scalarProd(P1, P2)*(scalarProd(P2, P6)*scalarProd(P4, P5)+\
                         scalarProd(P4, P6)*scalarProd(P2,P5))-\
scalarProd(P2, P6)*scalarProd(P1, P5)*scalarProd(P2, P4)-\
scalarProd(P4, P6)*scalarProd(P1, P5)*scalarProd(P2, P2)

U2 = scalarProd(P1, P1)*(scalarProd(P2, P6)*scalarProd(P4, P5)+\
                         scalarProd(P4, P6)*scalarProd(P2,P5))-\
scalarProd(P2, P6)*scalarProd(P1, P5)*scalarProd(P1, P4)-\
scalarProd(P4, P6)*scalarProd(P1, P5)*scalarProd(P1, P2)

assert(ftC.subs({u1:U1, u2:-U2}) == S(0))




## The condition ftC1=0 gives that in order to have 6 points in configuration
## (5) we can choose P2 and P4 in an arbitrary way, P1 = (P2||P4)
## P3 on the line P1+P2, P5 on the line P1+P4. Then P6 is a point on the 
## line P2+P4 determined by a linear equation in w1 and w2 given by ftC1=0.


## We need the following substitution:
sst_5 = {\
s11:scalarProd(P1, P1), \
s12:scalarProd(P1, P2), \
s13:scalarProd(P1, P3), \
s14:scalarProd(P1, P4), \
s15:scalarProd(P1, P5), \
s22:scalarProd(P2, P2), \
s23:scalarProd(P2, P3), \
s24:scalarProd(P2, P4), \
s25:scalarProd(P2, P5), \
s33:scalarProd(P3, P3), \
s34:scalarProd(P3, P4), \
s35:scalarProd(P3, P5), \
s44:scalarProd(P4, P4), \
s45:scalarProd(P4, P5), \
s55:scalarProd(P5, P5)\
}


## Then we have:
assert(ftC1.coefficient(w1) == \
   (s13*s24*s25 - 2*s15*s22*s34 + s13*s22*s45).subs(sst_5))

assert(ftC1.coefficient(w2) == \
   -(s15*s24*s34 + s15*s23*s44 - s13*s25*s44 - s13*s24*s45).subs(sst_5))


## Hence we can define P6 as follows:
P6 = (s15*s24*s34 + s15*s23*s44 - s13*s25*s44 - s13*s24*s45)*P2+\
     (s13*s24*s25 - 2*s15*s22*s34 + s13*s22*s45)*P4

P6 = P6.subs(sst_5)
## [[[If in the above formula for U1 and U2 we exchange P2 with P6, 
## we get the coordinates of P7:]]]
## Somehow we obtain P7:

L1 = scalarProd(P1, P6)*(scalarProd(P6, P2)*scalarProd(P4, P5)+\
                         scalarProd(P4, P2)*scalarProd(P6,P5))-\
scalarProd(P6, P2)*scalarProd(P1, P5)*scalarProd(P6, P4)-\
scalarProd(P4, P2)*scalarProd(P1, P5)*scalarProd(P6, P6)

L2 = scalarProd(P1, P1)*(scalarProd(P6, P2)*scalarProd(P4, P5)+\
                         scalarProd(P4, P2)*scalarProd(P6,P5))-\
scalarProd(P6, P2)*scalarProd(P1, P5)*scalarProd(P1, P4)-\
scalarProd(P4, P2)*scalarProd(P1, P5)*scalarProd(P1, P6)

P7 = L1*P1 - L2*P6

## hence in order to obtain configuration (5) we take:
## P2, P4 generic, 
## P1 = (P2 || P4)
## P3 on the line P1+P2
## P5 on the line P1+P4
## P6 and P7 given by the above formulas.

## Moreover, as said, we have several orthogonalities among the lines:

## P1+P2 ort P4+P6
## P1+P6 ort P2+P4
## P1+P4 ort P2+P6

## We do not have orthogonality of P1+P2 with P3+P5
## We do not have orthogonality of P1+P6 with P3+P5
## We do not have orthogonality of P1+P4 with P3+P5


assert(scalarProd(wedgeProd(P1, P2), wedgeProd(P3, P5)) != 0)
assert(scalarProd(wedgeProd(P1, P6), wedgeProd(P3, P5)) != 0)
assert(scalarProd(wedgeProd(P1, P4), wedgeProd(P3, P5)) != 0)

## The point P6 is always defined: 

J = S.ideal(list(P6))
J = J.saturation(u1*u2*v1*v2)[0]
J = J.saturation(matrix([P1, P2, P4]).det())[0]

assert(J == S.ideal(1))

## The point P7 is always defined:

print("27 seconds of computations...")
sleep(1)
gg = gcd(list(P7))
P7 = vector(S, [p7.quo_rem(gg)[0] for p7 in P7]) 
J = S.ideal(list(P7))
gJ = J.groebner_basis()
Jss = [0]
for ff in gJ:
    l1 = list(filter(lambda uu: not uu[0] in [u1, u2, v1, v2, det(matrix([P1, P2, P4]))], list(ff.factor())))
    l2 = list(map(lambda uu: uu[0], l1))
    Jss.append(prod(l2))
assert(S.ideal(Jss) == S.ideal(1))



##
## CONCLUSION: CONFIGURATION (5) IS POSSIBLE AND THE CUBICS ARE PARAMETRIZED BY
## (P^2, P^2,P^1,P^1), I.E. ARE OF DIMENSION 6


## a test:
sstA = {A2:3, B2:-2, C2:7, A4:-5, B4:1, C4:8, v1:3, v2:7, u1:9, u2:1}
p1 = P1.subs(sstA)
p2 = P2.subs(sstA)
p3 = P3.subs(sstA)
p4 = P4.subs(sstA)
p5 = P5.subs(sstA)
p6 = P6.subs(sstA)
p7 = P7.subs(sstA)

assert(matrixEigenpoints([p1, p2, p3, p4, p5, p6, p7]).rank() == 9)
cb = cubic_from_matrix(matrixEigenpoints([p1, p2, p3, p4, p5, p6, p7]))
ep = eigenpoints(cb)

Ep = points_from_id(ep)
Ep = list(map(tuple, Ep))

def last_coord_one(pt):
    return((pt[0]/pt[2], pt[1]/pt[2], 1))

assert(len(Ep) == 7)
## the point p1, ..., p7 defined above and the seven points 
## Ep, computed as the eigenpoints of the cubic cb, are the same:

assert(Set(list(map(last_coord_one, [p1, p2, p3, p4, p5, p6, p7]))) == Set(Ep))



#######
#######
#######
####### Here we prove that, if P1 = (P2||P4) and if P1, P2, P3, P4, P5 
## are a V configuration such that rank Phi(P1, ..., P5) = 9, then 
## the unique cubic which has, among the eigenpoints, the points 
## P1, ..., P5, has an eigenpoint P6 which is collinear with P2 and P4.
## 
## It is enough to prove the result for P2 = (1, 0, 0) or P2 = (1, ii, 0).
## In the second case, the line P1 + P2 is tangent in P2 to Ciso, so
## the rank of the matrix Phi(P1, ..., P5) cannot be 9.
## Here we verify this claim:

P2 = vector(S, (1, ii, 0))
P1 = vector(S, (-ii*B1, B1, C1))  ## we know that P1 is orthogonal to P2
assert(S.ideal(matrix([P1, P2, (x, y, z)]).det(), x^2+y^2+z^2).\
saturation(S.ideal(matrix([P1, P2]).minors(2)))[0] == S.ideal(x+ii*y, z^2))

## Hence it remains to consider the case P2 = (1, 0, 0).
## Here we prove that in this case, the point P6 given by the above formula, 
## i.e. 
## P6 = (s15*s24*s34 + s15*s23*s44 - s13*s25*s44 - s13*s24*s45)*P2+\
##      (s13*s24*s25 - 2*s15*s22*s34 + s13*s22*s45)*P4
##
## is such that rank Phi(P1, ..., P6) = 9.

## We redefine the points:

P2 = vector(S, (1, 0, 0))
P4 = vector(S, (A4, B4, C4))
P1 = wedgeProd(P2, P4)
P3 = u1*P1+u2*P2
P5 = v1*P1+v2*P4

P6 = (s15*s24*s34 + s15*s23*s44 - s13*s25*s44 - s13*s24*s45)*P2+\
     (s13*s24*s25 - 2*s15*s22*s34 + s13*s22*s45)*P4

sst_5 = {s11:scalarProd(P1, P1), \
s12:scalarProd(P1, P2), \
s22:scalarProd(P2, P2), \
s14:scalarProd(P1, P4), \
s24:scalarProd(P2, P4), \
s44:scalarProd(P4, P4), \
s13:scalarProd(P1, P3), \
s23:scalarProd(P2, P3), \
s34:scalarProd(P3, P4), \
s33:scalarProd(P3, P3), \
s45:scalarProd(P4, P5), \
s15:scalarProd(P1, P5), \
s25:scalarProd(P2, P5)}

P6 = P6.subs(sst_5)

M = matrixEigenpoints([P1, P2, P3, P4, P5, P6])

## The nine rows 0, 2, 3, 4, 6, 7, 9, 10, 12 of M are linearly independent:
assert(M.matrix_from_rows([0, 2, 3, 4, 6, 7, 9, 10, 12]).rank() == 9)

## now we see that the three new rows given by phi(P6) are lin. dep. w.r.t.
## the nine rows above:

##
print("three seconds of computations...")
sleep(1)

dtA1 = M.matrix_from_rows([0, 2, 3, 4, 6, 7, 9, 10, 12, 15]).det()
dtA2 = M.matrix_from_rows([0, 2, 3, 4, 6, 7, 9, 10, 12, 16]).det()
dtA3 = M.matrix_from_rows([0, 2, 3, 4, 6, 7, 9, 10, 12, 17]).det()

assert((dtA1, dtA2, dtA2) == (0, 0, 0))

## conclusion: P6 is an eigenpoint and P2, P4, P6 are aligned:

assert(det(matrix([P2, P4, P6])) == 0)

## The point P6 is always defined (is not (0,0,0)):

assert(S.ideal(list(P6)).saturation(u1*u2*v1*v2)[0].saturation(scalarProd(P1, P1))[0] == S.ideal(1))

#######
#######
#######
####### Here we prove that again, as above, if P1 = (P2||P4) and if 
## P1, P2, P3, P4, P5 are a V configuration such that rank 
## Phi(P1, ..., P5) = 9, then 
## the unique cubic which has, among the eigenpoints, the points 
## P1, ..., P5, has an eigenpoint P7 which is collinear with P1 and P6
## given by the formula:
## P7 = (s16*(s26*s45+s24*s56)-s26*s15*s46-s24*s15*s66)*P1-\
##      (s11*(s26*s45+s24*s56)-s26*s15*s14-s24*s15*s16)*P6
##

sst_5a = {s11:scalarProd(P1, P1), \
s12:scalarProd(P1, P2), \
s22:scalarProd(P2, P2), \
s14:scalarProd(P1, P4), \
s24:scalarProd(P2, P4), \
s44:scalarProd(P4, P4), \
s13:scalarProd(P1, P3), \
s23:scalarProd(P2, P3), \
s34:scalarProd(P3, P4), \
s33:scalarProd(P3, P3), \
s45:scalarProd(P4, P5), \
s15:scalarProd(P1, P5), \
s25:scalarProd(P2, P5), \
s16:scalarProd(P1, P6), \
s26:scalarProd(P2, P6), \
s36:scalarProd(P3, P6), \
s46:scalarProd(P4, P6), \
s56:scalarProd(P5, P6), \
s66:scalarProd(P6, P6)}

## we define pp7 with the above formula:
## (it is easy to see that this pp7 and the point P7 defined 
## above by the formula P7 = L1*P1 - L2*P6 are the same).

pp7 = (s16*(s26*s45+s24*s56)-s26*s15*s46-s24*s15*s66)*P1-\
      (s11*(s26*s45+s24*s56)-s26*s15*s14-s24*s15*s16)*P6

P7 = pp7.subs(sst_5a)

assert(gcd(list(P7)) == scalarProd(P1, P1)^4*v1)

## Since we know that s11 is not 0, we redefine P7:
P7 = vector(S, [p7.quo_rem(gcd(list(P7)))[0] for p7 in P7])

## now as above (note that now the last point of the matrixEigenpint
## below is P7 and not P6, as above).
M = matrixEigenpoints([P1, P2, P3, P4, P5, P7])

## The nine rows 0, 2, 3, 4, 6, 7, 9, 10, 12 of M are linearly independent:
assert(M.matrix_from_rows([0, 2, 3, 4, 6, 7, 9, 10, 12]).rank() == 9)

## now we see that the three new rows given by phi(P7) are lin. dep. w.r.t.
## the nine rows above:

##

print("180 seconds of computations...")
sleep(1)

dtA1 = M.matrix_from_rows([0, 2, 3, 4, 6, 7, 9, 10, 12, 15]).det()
dtA2 = M.matrix_from_rows([0, 2, 3, 4, 6, 7, 9, 10, 12, 16]).det()
dtA3 = M.matrix_from_rows([0, 2, 3, 4, 6, 7, 9, 10, 12, 17]).det()

assert((dtA1, dtA2, dtA2) == (0, 0, 0))


## conclusion: P7 is an eigenpoint and P1, P6, P7 are aligned:

assert(det(matrix([P1, P6, P7])) == 0)

## The point P7 is always defined (is not (0,0,0)):

assert(S.ideal(list(P7)).saturation(u1*u2*v1*v2)[0].saturation(scalarProd(P1, P1))[0] == S.ideal(1))



### Here we see that if we impose that (P2, P5, P7) are aligned, then 
### we get a (C8) configuration:

## we define a list of triplets, each is given by three points that 
## must not be aligned:

impossible_collin = [[P1, P2, P4], [P1, P2, P6], [P1, P4, P6], [P2, P3, P4], [P2, P4, P5], [P2, P6, P7], [P4, P6, P7], [P1, P2, P7], [P2, P3, P7]]

## We define the condition (P2, P5, P7) aligned and we simplify it, 
## erasing factors that are surely not zero:

J = S.ideal(matrix([P2, P5, P7]).det()).saturation(u1*u2*v1*v2)[0]
for tr in impossible_collin:
    J = J.saturation(matrix(tr).det())[0]


## J is a principal ideal. It has two factors. The fist implies that 
## P3, P4, P7 are aligned, so we have the alignments:
## (1, 2, 3), (1, 4, 5), (2, 4, 6), (2, 5, 7), (3, 4, 7)
## The second factor implies the alignments:
## (1, 2, 3), (1, 4, 5), (2, 4, 6), (2, 5, 7), (3, 5, 6):

pdJ = J.primary_decomposition()
assert(len(pdJ) == 2)
assert(pdJ[0].reduce(matrix([P3, P4, P7]).det()) == 0)
assert(pdJ[1].reduce(matrix([P3, P5, P6]).det()) == 0)


## similar computations can be done for the condition (P3, P4, P7) aligned
## or (P3, P5, P6) aligned. 