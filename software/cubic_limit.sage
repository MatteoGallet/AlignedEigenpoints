## Here we prove \\ref{prop:limitCubics}
## Hence we see that a cubic which has a line in the eigenpoints is a limit
## of cubics whith three aligned eigenpoints.

varAn1 = ["s"+str(i)+str(j) for i in range(1,8) for j in range(i, 8)]
varAn2 = ["x", "y", "z", "u1", "u2", "v1", "v2", "w1", "w2", "l1", "l2"]
varAn3 = ["a"+str(i) for i in range(10)]
varAn4 = [str(XX)+str(i) for i in range(1, 8) for XX in ["A", "B", "C"]]
var("xx")
K = QQ
K.<ii> = NumberField(xx^2+1)
S = PolynomialRing(K, varAn1+varAn2+varAn3+varAn4)
S.inject_variables(verbose=false)
P = vector(S, (x, y, z))
P1 = vector((A1, B1, C1))

load("auxiliaryProcedures.sage")

mon = [x^3, x^2*y, x*y^2, y^3, x^2*z, x*y*z, y^2*z, x*z^2, y*z^2, z^3]



## CASE 1: CUBICS WITH A LINE OF EIGENPOINTS (NOT TANGENT TO Ciso)

## line z = 0
 
P1 = vector(S, (1, 0, 0))  
P2 = vector(S, (0, 1, 0))  
P3 = P1 + P2

### P1, P2, P3 are aligned points.
## We construct the space of all the cubics which have P1, P2, P3 among
## the eigenpoints.

M = matrixEigenpoints([P1, P2, P3])

monV = vector(S, mon)

## we extract 6 lin. indep. rows:
M1 = M.matrix_from_rows([0, 1, 3, 5, 6, 7])
assert(M1.rank() == 6)

## Now we construct 4 more or less random cubics which are 
## cubics which have P1, P2, P3 among the eigenpoints.

M2 = M1.stack(vector(S, (0, 2, 0, 4, 0, 0, 0, 6, -7, 3)))
M2 = M2.stack(vector(S, (1, 1, 0, 2, 0, 0, 0, 8, 5, 3)))
M2 = M2.stack(vector(S, (4, 7, 0, 1, 0, 0, 0, 11, -4, 2)))
M2 = M2.stack(vector(S, mon))
cb1 = M2.det()

M2 = M1.stack(vector(S, (4, 1, 0, -2, 0, 0, 0, 3, -1, 5)))
M2 = M2.stack(vector(S, (2, -3, 0, 1, 0, 0, 0, 2, 1, 4)))
M2 = M2.stack(vector(S, (1, 2, 0, -1, 0, 0, 0, 9, -1, -2)))
M2 = M2.stack(vector(S, mon))
cb2 = M2.det()

M2 = M1.stack(vector(S, (7, 7, 0, 1, 0, 0, 0, 5, -2, -5)))
M2 = M2.stack(vector(S, (9, 7, 0, 2, 0, 0, 0, 1, -1, 3)))
M2 = M2.stack(vector(S, (5, 3, 0, 8, 0, 0, 0, 1, 11, 3)))
M2 = M2.stack(vector(S, mon))
cb3 = M2.det()

M2 = M1.stack(vector(S, (3, 5, 0, 2, 0, 0, 0, -1, 6, -1)))
M2 = M2.stack(vector(S, (2, 3, 0, 11, 0, 0, 0, -1, 1, 2)))
M2 = M2.stack(vector(S, (1, 2, 0, 1, 0, 0, 0, -1, 1, 4)))
M2 = M2.stack(vector(S, mon))
cb4 = M2.det()

## cb1, cb2, cb3, cb4 are linearly independent:
assert(matrix([[cb1.coefficient(mm) for mm in mon], \
              [cb2.coefficient(mm) for mm in mon], \
              [cb3.coefficient(mm) for mm in mon], \
              [cb4.coefficient(mm) for mm in mon]]).rank() == 4)


## we use cb1, cb2, cb3, cb4 to construct a simpler basis of 
## the 3-dim space of all the cubics with P1, P2, P3 eigenpoints:

Ma = matrix([[cb1.coefficient(mm) for mm in mon], \
             [cb2.coefficient(mm) for mm in mon], \
             [cb3.coefficient(mm) for mm in mon], \
             [cb4.coefficient(mm) for mm in mon]])


Ma = Ma.echelon_form()
## we redefine cb1, cb2, cb3, cb4:
cb1 = add([Ma[0][i]*mon[i] for i in range(10)])
cb2 = add([Ma[1][i]*mon[i] for i in range(10)])
cb3 = add([Ma[2][i]*mon[i] for i in range(10)])
cb4 = add([Ma[3][i]*mon[i] for i in range(10)])

assert(cb1 == x^3+y^3)
assert(cb2 == x*z^2)
assert(cb3 == y*z^2)
assert(cb4 == z^3)

## now cb1, cb2, cb3, cb4 are a good basis. From it we construct
## the generic cubic (which has P1, P2, P3 among the eigenpoints)

cb = u1*cb1+v1*cb2+w1*cb3+l1*cb4


## if u1 = 0, cb is the double line z=0 and a generic line:

assert(cb.subs(u1=0) == z^2*(v1*x+w1*y+l1*z))

## we extract the eigenpoints of cb (and we erase P1, P2, P3 and 
## we assume u1 != 0).

Jc = S.ideal(matrix([[x, y, z], [cb.derivative(x), \
                      cb.derivative(y), cb.derivative(z)]]).minors(2))

Jc = Jc.saturation(S.ideal(matrix([(x, y, z), P1]).minors(2)))[0]
Jc = Jc.saturation(S.ideal(matrix([(x, y, z), P2]).minors(2)))[0]
Jc = Jc.saturation(S.ideal(matrix([(x, y, z), P3]).minors(2)))[0]
Jc = Jc.saturation(u1)[0]

assert(Jc.is_prime())

## Jc is generated by g1, g2, g3 below:
 
g1 = 3*y^2*u1 - 2*x*y*v1 - 2*y^2*w1 + z^2*w1 - 3*y*z*l1
g2 = 3*x^2*u1 - 2*x^2*v1 + z^2*v1 - 2*x*y*w1 - 3*x*z*l1 
g3 = 2*x^3*y*v1 - 2*x^2*y^2*v1 + y^2*z^2*v1 + 2*x^2*y^2*w1 - \
     2*x*y^3*w1 - x^2*z^2*w1 + 3*x^2*y*z*l1 - 3*x*y^2*z*l1


assert(Jc == S.ideal(g1, g2, g3))

## the zeros of Jc should be 4 points which are- (in general) 
## distinct and with no other collinearities).

## FIRST CONCLUSION:
## A cubic of the form r^2*l is the limit of a family of cubics with three 
## aligned points (which depend on one parameter u1).

### CASE 2: CUBICS WITH EIGENPOINTS ON A LINE TANGENT TO Ciso.

### line x+ii*y = 0

P1 = vector(S, (1, ii, 0))
P2 = vector(S, (0, 0, 1))
P3 = P1 + P2

M = matrixEigenpoints([P1, P2, P3])


## we extract 5 lin. indep. rows:
M1 = M.matrix_from_rows([0, 1, 4, 5, 7])
assert(M1.rank() == 5)

## we construct 5 cubics which are a basis for the space of allo
## the cubics with P1, P2, P3 as eigenpoints.

M2 = M1.stack(vector(S, (0, 0, 1, 2, 0, -1, 3, 0, 0, 1)))
M2 = M2.stack(vector(S, (0, 0, -2, 1, 0, 1, 1, 0, 0, 2)))
M2 = M2.stack(vector(S, (0, 0, 4, -3, 0, 1, 1, 0, 0, 3)))
M2 = M2.stack(vector(S, (0, 0, 2, 1, 0, -3, 2, 0, 0, 1)))
M2 = M2.stack(monV)
cb1 = M2.det()

M2 = M1.stack(vector(S, (0, 0, 6, 3, 0, 5, 2, 0, 0, -1)))
M2 = M2.stack(vector(S, (0, 0, -2, 1, 0, -3, -1, 0, 0, 2)))
M2 = M2.stack(vector(S, (0, 0, 4, 6, 0, 2, 5, 0, 0, 3)))
M2 = M2.stack(vector(S, (0, 0, 1, 7, 0, -5, 4, 0, 0, 4)))
M2 = M2.stack(monV)
cb2 = M2.det()

M2 = M1.stack(vector(S, (0, 0, 5, 3, 0, 3, 2, 0, 0, -3)))
M2 = M2.stack(vector(S, (0, 0, -3, 1, 0, -6, -1, 0, 0, 2)))
M2 = M2.stack(vector(S, (0, 0, 5, 6, 0, 5, 5, 0, 0, 6)))
M2 = M2.stack(vector(S, (0, 0, 2, 7, 0, -5, 7, 0, 0, 5)))
M2 = M2.stack(monV)
cb3 = M2.det()

M2 = M1.stack(vector(S, (0, 0, 1, 3, 0, 3, 2, 0, 0, -2)))
M2 = M2.stack(vector(S, (0, 0, -3, 1, 0, -7, -4, 0, 0, 5)))
M2 = M2.stack(vector(S, (0, 0, 5, 3, 0, 5, 8, 0, 0, 2)))
M2 = M2.stack(vector(S, (0, 0, 2, 1, 0, -5, 7, 0, 0, 3)))
M2 = M2.stack(monV)
cb4 = M2.det()

M2 = M1.stack(vector(S, (0, 0, 7, 1, 0, 3, 2, 0, 0, -2)))
M2 = M2.stack(vector(S, (0, 0, -3, 1, 0, -5, -2, 0, 0, 5)))
M2 = M2.stack(vector(S, (0, 0, 1, 3, 0, 5, 1, 0, 0, 7)))
M2 = M2.stack(vector(S, (0, 0, 3, 2, 0, -1, 5, 0, 0, 3)))
M2 = M2.stack(monV)
cb5 = M2.det()

Ma = matrix([[cb1.coefficient(mm) for mm in mon], \
             [cb2.coefficient(mm) for mm in mon], \
             [cb3.coefficient(mm) for mm in mon], \
             [cb4.coefficient(mm) for mm in mon], \
             [cb5.coefficient(mm) for mm in mon]])

assert(Ma.rank() == 5)

## cb1, ..., cb5 are a basis. Now we want a better basis.
## We choose:

ccb1 = (x+ii*y)^2*x
ccb2 = (x+ii*y)^2*y
ccb3 = (x+ii*y)^2*z
ccb4 = z*(x^2+y^2+2/3*z^2)
ccb5 = x^3 + (-ii)*y^3 + z^3

## we verify that also this is a basis:

MMa = matrix([[ccb1.coefficient(mm) for mm in mon], \
             [ccb2.coefficient(mm) for mm in mon], \
             [ccb3.coefficient(mm) for mm in mon], \
             [ccb4.coefficient(mm) for mm in mon], \
             [ccb5.coefficient(mm) for mm in mon]])

assert(Ma.echelon_form() == MMa.echelon_form())

## we define the generic cubic linear combinations of ccb1, ...., ccb5

cb = u1*ccb1+u2*ccb2+v1*ccb3+v2*ccb4+w1*ccb5

## for w1 = 0, we get that cb is of the form t^2*l + L*C(r), so is the 
## generic cubic with a line of eigenpints tangent to Ciso, for w1 != 0
## the ideal of the remaining 4 eigenpoints is given by:

Jc = S.ideal(matrix([[x, y, z], [cb.derivative(x), \
                      cb.derivative(y), cb.derivative(z)]]).minors(2))

Jc = Jc.saturation(S.ideal(matrix([(x, y, z), P1]).minors(2)))[0]
Jc = Jc.saturation(S.ideal(matrix([(x, y, z), P2]).minors(2)))[0]
Jc = Jc.saturation(S.ideal(matrix([(x, y, z), P3]).minors(2)))[0]
Jc = Jc.saturation(w1)[0]

print("Now a computation of 1' 30''")
sleep(1)

pd = Jc.radical().primary_decomposition()

assert(pd[0] == Jc)

## The ideal Jc gives the 4 eigenpoints.



