## Three eigenpoints on Ciso.

### The result obtained by this file is:
### IF A CUBIC HAS AMONG ITS EIGENPOINTS A CONIC WHICH IS TANGENT TO THE
### ISOTROPIC CONIC Ciso IN TWO POINTS P1 AND P2, THEN THE CUBIC 
### IS OF THE FORM:
### r12*(U*r12^2+V*Ciso)
### WHERE r12 IS THE LINE P1+P2 AND U AND V ARE ELEMENTS IN THE FIELD K.
### HERE WE ASSUME P1 = (1, ii, 0) AND P2 != (1, -ii, 0)
### IN THE FILE caso2ptiIso_partII.sage WE CONSIDER THE CASE 
### P2 = (1, -ii, 0) AND WE GET THE SAME RESULT.  



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

doLongComputations = false

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

## Tangent line to Ciso in P1:
rtg1 = scalarProd(P1, vector((x, y, z))-P1)


## Tangent line to Ciso in P2:
rtg2 = scalarProd(P2, vector((x, y, z))-P2)*ii

## Pencil of conics tangent to P1 and P2 to Ciso:
Cg = Ciso + l1*rtg1*rtg2

## If l1 = -1, Cg is the conic given by (P1+P2)^2, so we can assume
## l1+1 != 0
assert(4*Cg.subs(l1=-1) == det(matrix([P1, P2, (x, y, z)]))^2)

## construction of a generic point (different from (1, ii, 0)) on Cg:

foo = Cg.subs(y=ii*x+w1*z).factor()[-1][0]

## generic point of Cg (depends on the parameter w1):
Pg = vector(S, (foo.coefficient(z), ii*(foo.coefficient(z))+\
                     w1*(-foo.coefficient(x)), -foo.coefficient(x)))

## the last coordinate of Pg is ((2*ii)) * (l1 + 1) * w1. 
## If w1 = 0, then Pg = P1, hence we can assume w1 != 0.
assert(matrix([P1, Pg.subs(w1=0)]).rank() == 1)

## Now we define three points on Cg:

Pg1, Pg2 = Pg.subs({w1:w1}), Pg.subs({w1:w2})
Pg3 = Pg.subs({w1:m1})

## and we can assume w1, w2, m1 != 0.

## the following matrix must have rank <= 9:

M = matrixEigenpoints([P1, P2, Pg1, Pg2, Pg3])

## We assume that u1 != 0, so P2 is not the point (1, -ii, 0)

#######################
#####################
## The case P2 = (1, -ii, 0) WILL BE CONSIDERED IN THE FILE:
## caso2ptiIso_partII.sage
#####################
#######################

## Under this hypothesis, we have that we can extract from M
## the rows: 0, 1; 4, 5; 7, 8; 10, 11; 13, 14.
## (Remember that w1, w2, m1 are not zero).

MM1 = M.matrix_from_rows([0, 1, 4, 5, 7, 8, 10, 11, 13, 14])

## we make a copy of MM1
MM2 = M.matrix_from_rows([0, 1, 4, 5, 7, 8, 10, 11, 13, 14])

## Using the fact that the first two rows are good
## (MM2[0,0] is a non zero constant, and MM2[1, 4] is 1)
## we can reduce MM2 with elementary rows and columns operations.

MM2.rescale_row(0, ii/3)
for i in range(2, 10):
    MM2.add_multiple_of_row(i, 0, -MM2[i][0])


for i in range(2, 10):
    MM2.add_multiple_of_row(i, 1, -MM2[i][4])

## We extract from MM2 an order 8 square matrix, 
## extracting the last 8 rows and the columns of position
## 1, 2, 3, 5, 6, 7, 8, 9. 
MM2 = MM2.matrix_from_rows_and_columns([2, 3, 4, 5, 6, 7, 8, 9], \
                                       [1, 2, 3, 5, 6, 7, 8, 9])

## The computation of det(MM2) gives 0 (time of computation: 3' 20'')
if doLongComputations:
   ttA = cputime()
   dtMM2 = MM2.det()
   print("Computation of the determinant of MM2:")
   print("time: "+str(cputime()-ttA))
   sleep(1)
else:
   dtMM2 = S(0)

assert(dtMM2 == S(0))

## Hence M and MM1 have rank <= 9.
## We want to see when M has rank <=8, i.e. when MM2 has rank <=6
## hence we compute the ideal of the order 7-minors of MM2.
## Time of computation: 35'

if doLongComputations:
   ttA = cputime()
   mm2_7 = MM2.minors(7)
   print("Computation of the order 7 minors:")
   print("time: "+str(cputime()-ttA))
   save(mm2_7, "longComput/mm2_7.sobj")
   sleep(1)
else:
   mm2_7 = load("longComput/mm2_7.sobj")
    
## here is a list of factors that cannot be zero:
nonZeroFt = [l1,v1,u1,w1,w2,m1,l1+1, (w1-w2),(u1-v1),\
             (w2-m1),(w1-m1),(v1*w2+1),(u1*w2+1),(v1*w1+1),\
             (u1*w1+1), (u1*m1+1)]

## a procedure:
## input: a polynomial and a list of polynomials:
## output: a list of polynomial obtained deleting from the factors of pol
## those factors which are contained in the list listFt. Each factor is taken 
## with exponent 1.
## example:
## erase_superfluous_factors(x^3*(x+1)^4*(z+y)^3*(x+y+z)^5, [x, z+y, z])
## answer: [x+1, x+y+z]

def erase_superfluous_factors(pol, listFt):
    ftOK = []
    for ft in list(factor(pol)):
        if not ft[0] in listFt:
            ftOK.append(ft[0])
    return(ftOK)

## here from each order 7-minor of MM2 (i.e. every element 
## of mm2_7) we clear off the superfluous factors and 
## we construct the ideal J7 of these polynomials (7'' of computation)

J7 = []
for ff in mm2_7:
    J7.append(prod(erase_superfluous_factors(ff, nonZeroFt)))

## since the ideal J7 (after suitable saturation) is (1), we have 
## that the matrix MM2 cannot have rank 6 or smaller, hence
## M cannot have rank 8 or smaller:
assert(S.ideal(J7).saturation((l1+1)*(w2-m1))[0] == S.ideal(S(1)))

## THE ABOVE COMPUTATIONS SHOW THAT THE MATRIX M HAS RANK 9
## THERE IS ONLY ONE CUBIC GIVEN BY M AND IS OBTAINED 
## FROM THE DETERMINANT OF THE MATRIX Mc GIVEN WHEN WE SUBSTITUTE 
## ONE OF THE ROWS OF MM1 WITH THE LIST mon.

Mc = MM1.matrix_from_rows(range(9))  
Mc = Mc.stack(matrix([mon]))

## the next computation requires 25':

if doLongComputations:
   ttA = cputime()
   dtMc = Mc.det()
   print("Computation of the cubic:")
   print("time: "+str(cputime()-ttA))
   save(dtMc, "longComput/dtMc.sobj")
   sleep(1)
else:
   dtMc = load("longComput/dtMc.sobj")

## from dtMc we erase useless factors:

dt1 = erase_superfluous_factors(dtMc, nonZeroFt)

## We get that dt1 is a list of 4 elements, which are:
Ls = [x*u1 + ii*y*u1 + ii*z,\
 u1*w1*w2 + 1/2*w1 + 1/2*w2,\
 u1^2*l1*m1^2 - l1*m1^2 - 2*u1*m1 - m1^2 - 1,\
 x^2*u1^2*l1 + (2*ii)*x*y*u1^2*l1 - y^2*u1^2*l1 + (2*ii)*x*z*u1*l1 -\
    2*y*z*u1*l1 + 3*x^2*l1 + 3*y^2*l1 + 2*z^2*l1 + 3*x^2 + 3*y^2 + 3*z^2]

assert(dt1 == Ls)

## we select the two factors which are a line and a conic and whose 
## product is the desired cubic.

## The factor dt1[0] is the line r12 given by P1+P2:

r12 = det(matrix([P1, P2, (x, y, z)]))/(2*ii)

assert(dt1[0] == r12)

## The factor Lt[3] is the l1*r12^2+3*(l1+1)*Ciso:

assert(dt1[3] == l1*r12^2+3*(l1+1)*Ciso)

## The cubic is therefore:

Cbc = dt1[0]*dt1[3]

assert(Cbc == r12*(l1*r12^2+3*(l1+1)*Ciso))

### HENCE WE HAVE THAT Cbc IS OF THE FORM r12*(u*r12^2+v*Ciso)
### FOR U, V IN THE FIELD K.

## now we compute the eigenpoints of Cbc:

Je = S.ideal(matrix([[Cbc.derivative(x), \
                      Cbc.derivative(y), \
                      Cbc.derivative(z)], [x, y, z]]).minors(2))

## we get two components: 
PDe = Je.primary_decomposition()
##
assert(PDe[0] == S.ideal(rtg1, rtg2))

assert(PDe[1] == S.ideal(Ciso+l1*rtg1*rtg2))

#################################
#############################
#####################
### Conversely, we assume here that the cubic C is of the 
### form C1 = (l1*Ciso+l2*(u1*x+v1*y+w1*z)^2)*(u1*x+v1*y+w1*z)
### where u1*x+v1*y+w1*z is a generic line of the plane and l1 and l2 
### are parameters.
C1 = (l1*Ciso+l2*(u1*x+v1*y+w1*z)^2)*(u1*x+v1*y+w1*z)

## we compute the eigenpoints of C1:
JJ = S.ideal(matrix([[C1.derivative(x), \
              C1.derivative(y), C1.derivative(z)], [x, y, z]]).minors(2))

## and its primary decomposition:

pdJJ = JJ.primary_decomposition()

## We get two components: 
## The first is:
## Ideal (z*v1 - y*w1, z*u1 - x*w1, y*u1 - x*v1)
## i.e. the point (u1, v1, w1)

assert(pdJJ[0] == S.ideal (z*v1 - y*w1, z*u1 - x*w1, y*u1 - x*v1))

## and the second component is the conic 
## l1*Ciso+3*l2*(u1*x+v1*y+w1*z)^2

assert(pdJJ[1] == S.ideal(l1*Ciso+3*l2*(u1*x+v1*y+w1*z)^2))












