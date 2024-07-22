# Here we prove that:
## If 4 distinct points of the isotropic conic Ciso are eigenpoints 
## of a cubic and if among the points we do not have (1: -ii: 0), 
## then all the cubics which have the 4 points as eigenpoints 
## split into the isotropic conic and a line of the plane. 
##
## Furthermore, the same is proved also if one of the points is
## (1: -ii: 0)
## 

## the file is executed in 14''


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

P1 = vector((1, ii, 0))

mon = [x^3, x^2*y, x*y^2, y^3, x^2*z, x*y*z, y^2*z, x*z^2, y*z^2, z^3]

load("../AlignedEigenpoints/software/ancillary.sage")

test = SymbolicCheck()

## construction of a generic point on the isotropic conic:
rt1 = l1*(y-ii*x)+l2*z

rt1.subs({x:P1[0], y:P1[1], z:P1[2]})

Ciso = x^2+y^2+z^2

## Intersection of the line rt1 with Ciso:
scndP = S.ideal(Ciso, rt1).radical().primary_decomposition()[1]
aux = scndP.gens()[:2]
mm2 = matrix([[aux[0].coefficient(x), aux[0].coefficient(y), \
aux[0].coefficient(z)],\
       [aux[1].coefficient(x), aux[1].coefficient(y), \
       aux[1].coefficient(z)]]).minors(2)

## Genric point on Ciso (depending on two parameters):

PP = vector(S, (mm2[2], -mm2[1], mm2[0]))

assert(scndP.subs({x:PP[0], y:PP[1], z:PP[2]}) == S.ideal(S(0)))
assert(Ciso.subs({x:PP[0], y:PP[1], z:PP[2]}) == S(0))

## we can always assume that l2 != 0, since l2=0 gives that 
## PP = P1
assert(matrix([P1, PP.subs(l2=0)]).rank() == 1)

## Definition of three points on Ciso
P2 = PP.subs({l1:u1, l2:1})
P3 = PP.subs({l1:v1, l2:1})
P4 = PP.subs({l1:w1, l2:1})

M = matrixEigenpoints([P1, P2, P3, P4])

### Case 1: we assume that the point (1, -ii, 0) is NOT one of the 
### points P2, P3, P4. Hence we can assume that u1, v1, w1, are not zero.


## The consequence is that the third coordinate of P2, P3 and P4 is 
## not zero.

## since we have: 
## P2[2]*M[3]-P2[1]*M[4]+P2[0]*M[5]
## P3[2]*M[6]-P3[1]*M[7]+P3[0]*M[8]
## P4[2]*M[9]-P4[1]*M[10]+P4[0]*M[11]

## we have that in our hypothesis, M[3], M[6], M[9] are unnecessary
## similarly M[2] is also unnecessary.


M1 = M.matrix_from_rows([0, 1, 4, 5, 7, 8, 10, 11]) 

### M1 has always rank 7:
## It cannot have rank 8, since all the order 8 minors are zero:
assert(Set(M1.minors(8)) == Set([S(0)]))

## If we impose that all the 7-minors of M1 are zero, we get 
## that (if the points P1, P2, P3, P4 are distinct) there are no
## solutions:

mm7 = S.ideal(M1.minors(7))
mm7 = mm7.saturation(S.ideal(u1*v1*w1))[0]
mm7 = mm7.saturation(S.ideal((u1-v1)*(u1-w1)*(v1-w1)))[0]
assert(mm7 == S.ideal(S(1)))

## conclusion 1: If four distinc points of Ciso are eigenpoints of 
## a cubic (and if the points are not (1, -i, 0)), then the cubic 
## splits into Ciso and a line.
## Proof: The matrix M1 has always rank 7, hence all the cubics 
## with P1, P2, P3, P4 eigenpoints are a linear variety in P^9 of 
## dimension 2. But this family contains all the cubics which split 
## into Ciso and a line of the plane, which is of dimension 2, hence
## the two families coincide.

## Case 2 the point  P2 is (1, -ii, 0)

P2 = vector(S, (1, -ii, 0))
P3 = PP.subs({l1:v1, l2:1})
P4 = PP.subs({l1:w1, l2:1})

M = matrixEigenpoints([P1, P2, P3, P4])

## since we have: 
## P3[2]*M[6]-P3[1]*M[7]+P3[0]*M[8]
## P4[2]*M[9]-P4[1]*M[10]+P4[0]*M[11]

## we have that in our hypothesis, M[6], M[9] are unnecessary
## similarly M[2] and M[5] are also unnecessary.


M1 = M.matrix_from_rows([0, 1, 3, 4, 7, 8, 10, 11]) 

### M1 has always rank 7:
## It cannot have rank 8:
assert(Set(M1.minors(8)) == Set([S(0)]))

## If we impose that all the 7-minors of M1 are zero, we get 
## that (if the points P1, P2, P3, P4 are distinct) there are no
## solutions:

mm7 = S.ideal(M1.minors(7))
mm7 = mm7.saturation(S.ideal(v1*w1))[0]
mm7 = mm7.saturation(S.ideal((v1-w1)))[0]
assert(mm7 == S.ideal(S(1)))

## conclusion 2: Also in case one of the point (1, -ii, 0) we 
## obtain the same conclusion than conclusion 1

## Final remark:
## the cubic Ciso*(l1*x+l2*y+l3*z) 
ccs = Ciso*(u1*x+v1*y+w1*z)

Je = S.ideal(matrix([[ccs.derivative(x), ccs.derivative(y), \
         ccs.derivative(z)], [x, y, z]]).minors(2))

