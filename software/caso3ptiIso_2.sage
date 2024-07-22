## Three eigenpoints on Ciso.

### The result obtained by this file is that:
### IT IS NOT POSSIBLE TO HAVE A CONIC OF EIGENPOINTS WHICH INTERSECTS
### THE ISOTROPIC CONIC Ciso IN THREE DISTINCT POINTS (HENCE IS TANGENT 
### IN A POINT TO Ciso AND PASSES THROUGH TWO OTHER PONTS OF Ciso).

### time of computation:
### 12'  if doLongComputations is true, 
### 13'' if doLongComputations is false.


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

## we can always assume that l2 != 0, since l2 = 0 gives that 
## PP = P1
assert(matrix([P1, PP.subs(l2=0)]).rank() == 1)

P2 = PP.subs({l1:u1, l2:1})
P3 = PP.subs({l1:v1, l2:1})

## If doLongComputations is true, we perform all the computations, 
## if it is false, we load the computations done in a previous section.
doLongComputations = false ##true

### Case 1: 
###
### we assume that the point (1, -ii, 0) is NOT one of the 
### points P2, P3. Hence we can assume that u1, v1 are not zero.

## Then we compute the pencil of conics passing through P1, P2, P3 and
## tangent to Ciso in P1:

## conic given by the two lines: P1+P2 and P1+P3
c2 = matrix([P1, P2, [x, y, z]]).det()*matrix([P1, P3, [x, y, z]]).det()

## generic conic tangent to Ciso in P1 and passing through P2 and P3
## (it depends on the parameter l1):
Cg = Ciso+l1*c2

## construction of a generic point (different from (1, ii, 0)) on Cg:
foo = Cg.subs(y=ii*x+w1*z).factor()[-1][0]

## generic point of Cg (depends on the parameter w1):
Pg = vector(S, (foo.coefficient(z), ii*(foo.coefficient(z))+\
                     w1*(-foo.coefficient(x)), -foo.coefficient(x)))

## Pg is the following:
assert(Pg == vector(S, (u1*v1*w1^2*l1 + u1*w1*l1 + v1*w1*l1 + \
                          1/4*w1^2 + l1 + 1/4, \
                        ii*u1*v1*w1^2*l1 + ii*u1*w1*l1 + ii*v1*w1*l1 + \
                          (-1/4*ii)*w1^2 + ii*l1 + (1/4*ii), \
                        (-1/2*ii)*w1)))

Pg1, Pg2 = Pg.subs({w1:w1}), Pg.subs({w1:w2})


## We can assume l1, u1, v1, w1, w2, (w1-w2), (u1-v1),
## (v1*w2+1), (u1*w2+1), (v1*w1+1), (u1*w1+1)
## all different from 0. 
## Indeed, l1 = 0 gives Cg = Ciso (condition already considered elsewere)
## u1 = 0 or v1 = 0 gives that P2 or P3 is the point (1, -ii, 0), and we 
## are in case 1, so this is not possible; 
## w1 = 0 gives Pg1 = P1:
assert(Set(matrix([Pg1.subs(w1=0), P1]).minors(2)) == Set([S(0)]))
## w2 = 0 gives Pg2 = P1:
assert(Set(matrix([Pg2.subs(w2=0), P1]).minors(2)) == Set([S(0)]))
##
## w1-w2 = 0 gives Pg1 = Pg2
## u1-v1 = 0 gives P2 = P3

## v1*w2+1 = 0 gives P3 = Pg2:
assert(S.ideal(matrix([P3, Pg2]).minors(2)).saturation(w2)[0] == \
          S.ideal(v1*w2+1))
## Similarly, (u1*w2+1), (v1*w1+1), (u1*w1+1) give, respectively:
## P2 = Pg2, P3 = Pg1, P2 = Pg1. Hence we define a polynomial dgnCs
## which contains all these degenerate cases and we can saturate our 
## computations w.r.t. this polynomial.

dgnCs = l1*v1*u1*w1*w2*(w1-w2)*(u1-v1)*(v1*w2+1)*(u1*w2+1)*(v1*w1+1)*(u1*w1+1)

## The points P1, P2, P3, Pg1, Pg2 are 5 points on the conic Cg.
## If they are eigenpoints, the following matrix must have
## rank 9 or less:

M = matrixEigenpoints([P1, P2, P3, Pg1, Pg2])

## We have P1[2]*M[0]-P1[1]*M[1]+P1[0]*M[2] = 0 and P1[0] = 1
## so M[2] is lin dep of M[0] and M[1] and can be omitted.

assert(P1[2]*M[0]-P1[1]*M[1]+P1[0]*M[2] == 0)
assert(P1[0] == 1)


## We have P2[2]*M[3]-P2[1]*M[4]+P2[0]*M[5] = 0 and P2[2] = 2*u1
## which, under our hypothesis, is always not zero, 
## so M[3] can be omitted.

assert(P2[2]*M[3]-P2[1]*M[4]+P2[0]*M[5] == 0)
assert(P2[2] == 2*u1)


##
## In a similar way, M[6], M[9], M[12], M[15] can be omitted.
## Hence we construct the square matrix of order 10:

MM1 = M.matrix_from_rows([0, 1, 4, 5, 7, 8, 10, 11, 13, 14])

## The first row of MM1 is ((-3*ii), 3, (3*ii), -3, 0, 0, 0, 0, 0, 0)
## The second row is (0, 0, 0, 0, 1, ii, -1, 0, 0, 0)

assert(MM1[0] == vector(S, ((-3*ii), 3, (3*ii), -3, 0, 0, 0, 0, 0, 0)))
assert(MM1[1] == vector(S, (0, 0, 0, 0, 1, ii, -1, 0, 0, 0)))

## so with elementary row operations we can simplify MM1:
MM1.rescale_row(0, 1/3*ii)
for i in range(2, 10):
    MM1.add_multiple_of_row(i, 0, -MM1[i][0])


for i in range(2, 10):
    MM1.add_multiple_of_row(i, 1, -MM1[i][4])

## Now the 0-th column of MM1 is (1, 0, ..., 0)
## and the 4-th column of MM1 is (0, 1, 0, ..., 0)

assert([MM1[i, 0] for i in range(10)] == [1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
assert([MM1[i, 4] for i in range(10)] == [0, 1, 0, 0, 0, 0, 0, 0, 0, 0])


## We extract from MM1 an order 8 square matrix, 
## extracting the last 8 rows and the columns of position
## 1, 2, 3, 5, 6, 7, 8, 9. 
MM1 = MM1.matrix_from_rows_and_columns([2, 3, 4, 5, 6, 7, 8, 9], \
                                       [1, 2, 3, 5, 6, 7, 8, 9])

## This new matrix MM1 has rank n iff the original MM1 has rank n+2
## (iff M has rank n+2)
## In particular, we want to see if MM1 can have rank <= 6.

if doLongComputations:
   print("A computation of about 8' 40''")
   sleep(1)
   ttA = cputime()
   mm1_7 = MM1.minors(7)
   print("Computation of the 64 order 7 minors:")
   print("time: "+str(cputime()-ttA))
   save(mm1_7, "longComput/mm1_7_3ptiIso.sobj")
else:
   mm1_7 = load("longComput/mm1_7_3ptiIso.sobj")

J7 = S.ideal(mm1_7)
J7 = J7.saturation(dgnCs)[0]
## The above ideal is (1), so it is not possible to have that the 
## matrix MM1 has rank <= 6 (hence it is not possible that M has rank <= 8)
assert(J7 == S.ideal(S(1)))

## In order to have P1, P2, Pg1, Pg2, Pg3 eigenpoints, 
## M must have det zero. But det(M) = det(MM1). So we compute 
## det(MM1) and we saturate it w.r.t. dgnCs.

if doLongComputations:
   print("About  44'' of computation:")
   sleep(1)
   ttA = cputime()
   ddt = MM1.det()
   ddt = S.ideal(ddt).saturation(dgnCs)[0].gens()[0]
   print(cputime()-ttA)
   sleep(1)
else:
   ddt = (16) * (l1 + 1/4)^2 * (u1*w1*w2 + v1*w1*w2 + w1 + w2)

assert(ddt == (16) * (l1 + 1/4)^2 * (u1*w1*w2 + v1*w1*w2 + w1 + w2))


## Hence we have two possibilities: 
## l1 + 1/4 = 0
## or u1*w1*w2 + v1*w1*w2 + w1 + w2 = 0.
## This second condition must be satisfied for every point Pg2 of 
## the conic Cg, i.e. for every w2, and therefore this condition is 
## impossible. 

## Consider the case l1 = -1/4.
## In this case Cg splits into two lines: the line
## (x + ii*y) 
## and the line (x*u1*v1 + ii*y*u1*v1 + ii*z*u1 + ii*z*v1 + x + (-ii)*y)
## The first is the tangent line to Ciso in the point P1
## and the second line is the line passing through the points P2 and P3:

assert(S.ideal(Ciso, x+ii*y).groebner_basis() == [z^2, x + ii*y])

assert(Cg.subs(l1=-1/4) == \
          (x+ii*y)*S(det(matrix([P2, P3, (x, y, z)]))/(2*(u1-v1))))

## The matrix M.matrix_from_rows([0, 1, 4, 5, 7, 8, 10, 11, 13, 14])
## i.e. the original matrix extracted from M,  with the condition 
## l1=-1/4 has rank 9:
MM2 = M.matrix_from_rows([0, 1, 4, 5, 7, 8, 10, 11, 13, 14]).subs(l1 = -1/4)
assert(MM2.rank() == 9)

## and the last row of MM2 is linearly dependent w.r.t. the other 
## rows of MM2. Indeed, the rank of MM2 without the last row is 
## still 9:

MM3 = MM2.matrix_from_rows([0, 1, 2, 3, 4, 5, 6, 7, 8])
assert(MM3.rank() == 9)

## hence the cubic corresponding to the matrix MM2 is:

if doLongComputations:
   print("about 60'' of computations:")
   sleep(1)
   ttA = cputime()
   cbb = det(MM3.stack(matrix([mon])))
   print(cputime()-ttA)
   ttA = cputime()
   sleep(1)
   print("about 17'' of computations")
   sleep(1)
   cbb = S.ideal(cbb).saturation(dgnCs)[0].gens()[0]
   print(cputime()-ttA)
   sleep(1)
else:
   cbb = (x + ii*y) * (u1*v1 - 1) * \
       (x*u1*v1 + ii*y*u1*v1 + ii*z*u1 + ii*z*v1 + x + (-ii)*y)^2

assert(cbb == (x + ii*y) * (u1*v1 - 1) * \
          (x*u1*v1 + ii*y*u1*v1 + ii*z*u1 + ii*z*v1 + x + (-ii)*y)^2)

## and cbb splits into the line (P1+P2)^2 and the tangent to 
## Ciso in the point P1. The eigenpoints are the line (P1+P2) plus 
## other points and is not a conic.

## HENCE IN THIS SITUATION WE DO NOT HAVE CUBICS WHOSE EIGENPOINTS 
## ARE A CONIC.


####################
####################
## Case 2: 
## P3 is the point (1, -ii, 0):
##

P3 = ii*P3.subs(v1=0)

assert(P3 == vector(S, (1, -ii, 0)))

## conic given by two lines: P1+P2 and P1+P3
c2 = matrix([P1, P2, [x, y, z]]).det()*matrix([P1, P3, [x, y, z]]).det()

## generic conic tangent to Ciso in P1 and passing through P2 and P3
## (it depends on the parameter l1):
Cg = Ciso+l1*c2

## construction of a generic point (different from (1, ii, 0)) on Cg:
foo = Cg.subs(y=ii*x+w1*z).factor()[-1][0]

## generic point of Cg (depends on the parameter w1):
Pg = vector(S, (foo.coefficient(z), ii*(foo.coefficient(z))+\
                     w1*(-foo.coefficient(x)), -foo.coefficient(x)))


Pg1, Pg2 = Pg.subs({w1:w1}), Pg.subs({w1:w2})


M = matrixEigenpoints([P1, P2, P3, Pg1, Pg2])


## We have P1[2]*M[0]-P1[1]*M[1]+P1[0]*M[2] = 0 and P1[0] = 1
## so M[2] is lin dep of M[0] and M[1] and can be omitted.

assert(P1[2]*M[0]-P1[1]*M[1]+P1[0]*M[2] == 0)
assert(P1[0] == 1)

## We have P2[2]*M[3]-P2[1]*M[4]+P2[0]*M[5] = 0 and P2[2] = 2*u1
## which, under our hypothesis, is always not zero, 
## so M[3] can be omitted.

assert( P2[2]*M[3]-P2[1]*M[4]+P2[0]*M[5] == 0)
assert(P2[2] == 2*u1)


##
## In a similar way, M[9], M[12], M[15] can be omitted.
##
## P3[2]*M[6]-P3[1]*M[7]+P3[0]*M[8] = 0 and P3[0] = 1, 
## hence M[8] can be omitted.

## Hence we construct the square matrix of order 10:

MM1 = M.matrix_from_rows([0, 1, 4, 5, 6, 7, 10, 11, 13, 14])


## The first row of MM1 is ((-3*ii), 3, (3*ii), -3, 0, 0, 0, 0, 0, 0)
## The second row is (0, 0, 0, 0, 1, ii, -1, 0, 0, 0)

assert(MM1[0] == vector(S, ((-3*ii), 3, (3*ii), -3, 0, 0, 0, 0, 0, 0)))
assert(MM1[1] == vector(S, (0, 0, 0, 0, 1, ii, -1, 0, 0, 0)))


## so with elementary row operations we can simplify MM1:
MM1.rescale_row(0, 1/3*ii)
for i in range(2, 10):
    MM1.add_multiple_of_row(i, 0, -MM1[i][0])


for i in range(2, 10):
    MM1.add_multiple_of_row(i, 1, -MM1[i][4])

## Now the 0-th column of MM1 is (1, 0, ..., 0)
## and the 4-th column of MM1 is (0, 1, 0, ..., 0)

assert([MM1[i, 0] for i in range(10)] == [1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
assert([MM1[i, 4] for i in range(10)] == [0, 1, 0, 0, 0, 0, 0, 0, 0, 0])

## We extract from MM1 an order 8 square matrix, 
## extracting the last 8 rows and the columns of position
## 1, 2, 3, 5, 6, 7, 8, 9. 
MM1 = MM1.matrix_from_rows_and_columns([2, 3, 4, 5, 6, 7, 8, 9], \
                                       [1, 2, 3, 5, 6, 7, 8, 9])

## This new matrix MM1 has rank n iff the original MM1 has rank n+2
## (iff M has rank n+2)
## In particular, we want to see if MM1 can have rank <= 6.

if doLongComputations:
   ttA = cputime()
   mm1_7b = MM1.minors(7)
   print("Computation of the 64 order 7 minors (about 25''):")
   sleep(1)
   print("time: "+str(cputime()-ttA))
   save(mm1_7b, "longComput/mm1_7b_3ptiIso.sobj")
else:
   mm1_7b = load("longComput/mm1_7b_3ptiIso.sobj")

dgnCs = l1*u1*w1*w2*(w1-w2)*(u1-v1)*(v1*w2+1)*(u1*w2+1)*(v1*w1+1)*(u1*w1+1)

J7 = S.ideal(mm1_7b)
J7 = J7.saturation(dgnCs)[0]
## The above ideal is (1), so it is not possible to have that the 
## matrix MM1 has rank <= 6 (hence it is not possible that M has rank <= 8)
assert(J7 == S.ideal(S(1)))

## We want now to compute the determinant of the original MM1 
## (the order 10 matrix) which is equal to the determiant of the order 
## 8 matrix MM1:

dt2 = MM1.det()

dt2 = S.ideal(dt2).saturation(dgnCs)[0].gens()[0]

## we get that dt2 is:
## (16) * (l1 + (-1/4*ii))^2 * (u1*w1*w2 + w1 + w2)

assert(dt2 == (16) * (l1 + (-1/4*ii))^2 * (u1*w1*w2 + w1 + w2))

## as above, u1*w1*w2+w1+w2 = 0, since w2 is generic, does not 
## give solutions.

## The case l1 = 1/4*ii
## In this case Cg splits into the line x+ii*y (tangent line to Cg in P1)
## and the line P2+P3:
assert((-2)*u1*Cg.subs(l1=1/4*ii) == \
      ii*(x+ii*y)*matrix([P2, P3, (x, y, z)]).det())

## We consider the substitution of l1=1/4*ii in the original MM,
## i.e. in MM1 = M.matrix_from_rows([0, 1, 4, 5, 6, 7, 10, 11, 13, 14])

MM2 = M.matrix_from_rows([0, 1, 4, 5, 6, 7, 10, 11, 13, 14]).subs(l1=1/4*ii)

## Also here the last row of MM2 is unnecessary:
## assert(MM2.rank() == MM2.matrix_from_rows(range(9)).rank())

## Hence we compute the cubic:
M3 = (MM2.matrix_from_rows(range(9))).stack(matrix([mon]))
cb2 = M3.det()

## The cubic splits into the line tangent to Ciso in P1 and the 
## reducible conic given by (P2+P3)^2
## 
assert(S.ideal(cb2).saturation(dgnCs)[0].gens()[0] == \
    (x + ii*y) * (z*u1 + (-ii)*x - y)^2)

### CONCLUSION:
### IT IS NOT POSSIBLE TO HAVE A CONIC OF EIGENPOINTS WHICH INTERSECTS
### THE ISOTROPIC CONIC Ciso IN THREE DISTINCT POINTS (HENCE IS TANGENT 
### IN A POINT TO Ciso).



