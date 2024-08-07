## Here we have the proof of 
## \\ref{cubiche_con_retta_autop_tg}
## i.e. the prop. which says when a cubic has a line of eigenpoints
## and of prop. the proposition which says that all the cubic with 
## a line of eigenpoints which is tangent to Ciso has the tancency point
## singular. 


varAn1 = ["s"+str(i)+str(j) for i in range(1,8) for j in range(i, 8)]
varAn2 = ["x", "y", "z", "u1", "u2", "v1", "v2", "w1", "w2", "l1", "l2", "m1", "m2"]
varAn3 = ["a"+str(i) for i in range(10)]
varAn4 = [str(XX)+str(i) for i in range(1, 8) for XX in ["A", "B", "C"]]
var("xx")
K.<ii> = NumberField(xx^2+1)
S = PolynomialRing(K, varAn1+varAn2+varAn3+varAn4)
S.inject_variables(verbose=false)

P1 = vector(S, (A1, B1, C1))

load("auxiliaryProcedures.sage")

P1 = vector((1, ii, 0))

mon = [x^3, x^2*y, x*y^2, y^3, x^2*z, x*y*z, y^2*z, x*z^2, y*z^2, z^3]

load("../AlignedEigenpoints/software/ancillary.sage")

test = SymbolicCheck()

Ciso = x^2+y^2+z^2

F = add([S(varAn3[i])*mon[i] for i in range(10)])

## gradient:
def gdn(F):
    return(vector(S, [F.derivative(xx) for xx in [x, y, z]]))

## a vector which gives the equations of the eigenpoints:
def eig(cub):
    gF = gdn(cub)
    return(vector(S, matrix([gF, [x, y, z]]).minors(2)))


## Tangent line to Ciso in P1:

tg1 = scalarProd(gdn(Ciso).subs(substitution(P1)), vector((x, y, z))-P1)/2

assert(tg1 == x+ii*y)

## the eigenpoints locus of F:
eigF = eig(F)

EE = list(eigF.subs(x = -ii*y)) 

## equations in a0, ..., a9 that must be satisfied
## in order to have tg1 of eigenpoints:

eqz = [[ee.coefficient(mm) for mm in mon] for ee in EE]


eqz1 = []
for ee in eqz:
    eqz1 += ee

F1 = S.ideal(eqz1).reduce(F)

## F1 has the variables x, y, z, a2, a3, a6, a9:
assert(F1.variables() == (x, y, z, a2, a3, a6, a9))

## F1 is a linear combination of the following four cubics:
G1, G2, G3, G4 = F1.coefficient(a2), F1.coefficient(a3), \
                 F1.coefficient(a6), F1.coefficient(a9)

## we define four other cubics which are:
H1, H2, H3 = x*tg1^2, y*tg1^2, z*tg1^2

## and H4 constructed in this way:

## we take a line passing through the point P1:
r = z
## and H4 is the corresponding cubic whose eigenpoints are 
## two tangent lines to Ciso in the points r \cap Ciso.

H4 = r*(r^2-3*(0^2+0^2+1^2)*Ciso)

## we verify that the linear space generated by G1, G2, G3, G4
## coincides with the linear space generated by H1, H2, H3, H4:

for hh in [H1, H2, H3, H4]:
    assert(matrix([[gg.coefficient(mn) for mn in mon] for gg in \
           [G1, G2, G3, G4]+[hh.coefficient(mn) for mn in mon]]).rank()==4)

for gg in [G1, G2, G3, G4]:
    assert(matrix([[hh.coefficient(mn) for mn in mon] for hh in \
           [H1, H2, H3, H4]+[gg.coefficient(mn) for mn in mon]]).rank()==4)


#### CONCLUSION 1:
#### IF A CUBIC HAS A LINE t OF EIGENPOINTS AND t IS TANGENT TO Ciso, 
#### THEN THE EQUATION OF THE CUBIC IS r*(t^2*l+lambda*C(r)), FOR A FIXED r
#### THROUGH P1.
#### 
####
#### NOW THE CONVERSE:

## 
## if rr = u1*x+v1*y+w1*z is the generic cubic of the plane, 
## all the cubics of equation rr*(rr^2-3*(u1^2+v1^2+w1^2)*Ciso)
## are in the linear system generated by H1, H2, H3, H4:

Cb = rr*(rr^2-3*(u1^2+v1^2+w1^2)*Ciso)

assert(matrix([[hh.coefficient(mn) for mn in mon] for hh in \
           [H1, H2, H3, H4]+[Cb.coefficient(mn) for mn in mon]]).rank()==4)

## From this we get that linear combinations of H1, H2, H3, H4 give, 
## in particular, all the cubics C(r), for all possible r.


## Now we take a cubic of the form tg1^2*l+lambda C(r), hence a cubic 
## which is a linear combination of H1, H2, H3, H4.

FF = A1*H1+A2*H2+A3*H3+A4*H4

pdFF = S.ideal(list(eig(FF))).primary_decomposition()

### Among the eigenpoints, we have tg1:
assert(pdFF[0] == S.ideal(tg1))

#### CONCLUSION 2:
#### CUBICS OF THE FORM t^2*l+lambda*C(r) ALWAYS HAVE tg1 OF EIGENPOINTS.

####
#### THIRD COMPUTATION: ALL THE CUBICS FF HAVE THE POINT P1 SINGULAR 
#### AND tg1 TANGENT IN P1.
####
##
## all the cubics FF have P1 singular:
assert(gdn(FF).subs(substitution(P1)) == vector(S, (0,0,0)))

## tg1 is tangent (three intersections of tg1 with FF in P1):
assert(S.ideal(FF, tg1).radical() == S.ideal(x+ii*y, z*A4))
assert(S.ideal(FF, tg1).primary_decomposition()[0] == S.ideal(tg1, z^3))



