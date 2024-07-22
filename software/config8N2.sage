## HERE WE PROVE THAT CONFIGURATION (8) IS POSSIBLE AND MUCH MORE.

## 
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



## print("We define 'mon' and 'Jrel'")
## mon = [x^3, x^2*y, x*y^2, y^3, x^2*z, x*y*z, y^2*z, x*z^2, y*z^2, z^3]

## We define 4 generic points in the plane
P1 = vector((A1, B1, C1))
P2 = vector((A2, B2, C2))
P4 = vector((A4, B4, C4))
P7 = vector((A7, B7, C7))  

## A first property of 4 generic points in the plane:
## It is not possible that every couple of different points are orthogonal, 
## indeed, the following ideal is (1):

JJ = S.ideal(scalarProd(P1, P2), scalarProd(P1, P4), \
        scalarProd(P1, P7), scalarProd(P2, P4),
        scalarProd(P2, P7), scalarProd(P4, P7)).\
     saturation(S.ideal(matrix([P2, P4]).minors(2)))[0].\
     saturation(matrix([P2, P4, P7]).det())[0].\
     saturation(S.ideal(list(P1)))[0]

assert(JJ == S.ideal(1))

### this property will be used later. We call it [AAA] property, 
### for further references in this file.


### another property of three distinct of the plane P1, P2, P4:
### the three vectors (P1||P2), (P1||P4), (P2||P4) are linearly independent.
 
ddt = matrix([wedgeProd(P1, P2), wedgeProd(P1, P4), wedgeProd(P2, P4)]).det()
assert(ddt == matrix([P1, P2, P4]).det()^2)

### this property will be used later. We call it [BBB]


## and we define P5 P6 P7 in order to obtain configuration (8),
## which is given by the alignments:
## (1, 2, 3), (1, 4, 5), (1, 6, 7), (2, 4, 6), (2, 5, 7), (3, 4, 7) 
## 

P3 = intersect_lines(P1, P2, P7, P4)  ## al posto di 5 scritto 3
P5 = intersect_lines(P1, P4, P2, P7) ## al posto di 6 messo 5
P6 = intersect_lines(P1, P7, P2, P4)  ## al posto di 7 messo 6

Jrel = S.ideal(s11-scalarProd(P1, P1), \
               s12-scalarProd(P1, P2), \
               s22-scalarProd(P2, P2), \
               s14-scalarProd(P1, P4), \
               s24-scalarProd(P2, P4), \
               s44-scalarProd(P4, P4), \
               s17-scalarProd(P1, P7), \
               s27-scalarProd(P2, P7), \
               s47-scalarProd(P4, P7), \
               s77-scalarProd(P7, P7))


## If configuration (8) is given by eigenpoints, we must have
## e1 = e2 = e3 = 0, where:
# e1 = delta1(P5, P1, P2), e2 = delta1(P6, P1, P2), e3 = delta1(P3, P1, P4)

e1 = delta1(P5, P1, P2)
e2 = delta1(P6, P1, P2)
e3 = delta1(P3, P1, P4)


## We have:
##e1 == matrix([P1, P2, P7]).det()*matrix([P1, P2, P4]).det()*
##      *(s12*s47-s17*s24)

assert(e1 == matrix([P1, P2, P7]).det()*matrix([P1, P2, P4]).det()*(scalarProd(P1, P2)*scalarProd(P4, P7)-scalarProd(P1, P7)*scalarProd(P2, P4)))

# Similarly,  

assert(e2 == matrix([P1, P2, P7]).det()*matrix([P1, P2, P4]).det()*(scalarProd(P1, P2)*scalarProd(P4, P7)-scalarProd(P1, P4)*scalarProd(P2, P7)))

assert(e3 == -matrix([P1, P4, P7]).det()*matrix([P1, P2, P4]).det()*(scalarProd(P1, P4)*scalarProd(P2, P7)-scalarProd(P1, P7)*scalarProd(P2, P4)))

## we redefine e1, e2, e3 as s12*s47-s17*s24, s12*s47-s14*s27 and 
## s17*s24-s14*s27

e1 = scalarProd(P1, P2)*scalarProd(P4, P7)-scalarProd(P1, P7)*scalarProd(P2, P4)

e2 = scalarProd(P1, P2)*scalarProd(P4, P7)-scalarProd(P1, P4)*scalarProd(P2, P7)

e3 = -scalarProd(P1, P4)*scalarProd(P2, P7)+scalarProd(P1, P7)*scalarProd(P2, P4)

## we have: e1-e2+e3 == 0, hence e3 is not necessary
assert(e1-e2+e3 == 0)

## We solve the system e1 = 0, e2 = 0, linear in A4, B4, C4 and 
## we get the solution that we denote by ss6. We verify that ss6 is the 
## solution of the two equations e1 = 0, e2 = 0

M1 = matrix([[e1.coefficient(A4), e1.coefficient(B4), e1.coefficient(C4)],\
[e2.coefficient(A4), e2.coefficient(B4), e2.coefficient(C4)]])

minM1 = M1.minors(2)

ss6 = {A4: minM1[2], B4: -minM1[1], C4: minM1[0]}

assert(e1.subs({A4:minM1[2], B4: -minM1[1], C4: minM1[0]}) == 0)
assert(e2.subs({A4:minM1[2], B4: -minM1[1], C4: minM1[0]}) == 0)

## The solution ss6 is not unique (degenerate case) iff 
## s12=0, s17=0
## or 
## s27=0, s17=0
## or
## s12=0, s27=0
## as we can obtain from the computations below


pd = S.ideal(minM1).saturation(det(matrix([P1, P2, P7])))[0].radical().primary_decomposition()

assert(pd[0] == S.ideal(scalarProd(P1, P2), scalarProd(P1, P7)))
assert(pd[1] == S.ideal(scalarProd(P2, P7), scalarProd(P1, P7)))
assert(pd[2] == S.ideal(scalarProd(P1, P2), scalarProd(P2, P7)))

#### We will study these conditions later and we will see that 
#### the three conditions are actually only one which says that
#### s12 = 0, s17 = 0, s27 = 0.
#### see (+++)

## Now we redefine the points P4, P3, P5, P6 using the substitution
## ss6 and we obtain the points PP4, PP3, PP5, PP6

PP4 = P4.subs(ss6)
PP3 = P3.subs(ss6)
PP5 = P5.subs(ss6)
PP6 = P6.subs(ss6)


## and we get that 4 deltas2 are zero:

assert(delta2(P1, P2, PP3, PP4, PP5) == 0)
assert(delta2(PP4, P1, PP5, P2, PP6) == 0)
assert(delta2(P7, PP3, PP4, P2, PP5) == 0)
assert(delta2(P2, P1, PP3, PP4, PP6) == 0)


## We have that the lines P1+P2 and PP3+PP4 are orthogonal,
## the lines P1+PP6 and P2+PP4 are orthogonal and
## the lines P1+PP4 and P2+PP5 are orthogonal:
assert(scalarProd(wedgeProd(P1, P2), wedgeProd(PP3, PP4)) == S(0))
assert(scalarProd(wedgeProd(P1, PP6), wedgeProd(P2, PP4)) == S(0))
assert(scalarProd(wedgeProd(P1, PP4), wedgeProd(P2, PP5)) == S(0))

## Moreover,it holds: 
## PP4 == (P1||P2)(P1|P7)(P2|P7)-(P1|P2)(P1||P7)(P2|P7)+(P1|P2)(P1|P7)(P2||P7)

Q4 = wedgeProd(P1, P2)*scalarProd(P1, P7)*scalarProd(P2, P7)\
        -scalarProd(P1, P2)*wedgeProd(P1, P7)*scalarProd(P2, P7)\
        +scalarProd(P1, P2)*scalarProd(P1, P7)*wedgeProd(P2, P7)

assert(matrix([PP4, Q4]).minors(2) == [0, 0, 0])


## We want to see that the rank of the matrix 
## Phi(P1, P2, PP3, PP4, PP5, PP6, P7) is <=9.
## We redefine the points.
p1 = vector(S, (1, 0, 0))

p2 = vector(S, (A2, B2, C2))
p7 = vector(S, (A7, B7, C7))
p4 = wedgeProd(p1, p2)*scalarProd(p1, p7)*scalarProd(p2, p7)\
        -scalarProd(p1, p2)*wedgeProd(p1, p7)*scalarProd(p2, p7)\
        +scalarProd(p1, p2)*scalarProd(p1, p7)*wedgeProd(p2, p7)
p3 = intersect_lines(p1, p2, p4, p7)
p5 = intersect_lines(p1, p4, p2, p7)
p6 = intersect_lines(p1, p7, p2, p4)


## We have: the matrix Phi([p1, p2, p3, p4, p5, p6, p7]) has rank 9:

assert(matrixEigenpoints([p1, p2, p3, p4, p5, p6, p7]).rank() == 9)

## Similarly, we can verify the result for the case in which 
## p1 is (1, ii, 0), but in this case the computation requires about 10 minutes

## We redefine the points.
p1 = vector(S, (1, ii, 0))

p2 = vector(S, (A2, B2, C2))
p7 = vector(S, (A7, B7, C7))
p4 = wedgeProd(p1, p2)*scalarProd(p1, p7)*scalarProd(p2, p7)\
        -scalarProd(p1, p2)*wedgeProd(p1, p7)*scalarProd(p2, p7)\
        +scalarProd(p1, p2)*scalarProd(p1, p7)*wedgeProd(p2, p7)
p3 = intersect_lines(p1, p2, p4, p7)
p5 = intersect_lines(p1, p4, p2, p7)
p6 = intersect_lines(p1, p7, p2, p4)


## We have: the matrix Phi([p1, p2, p3, p4, p5, p6, p7]) has rank 9:

## computation of about 10 minutes:
## assert(matrixEigenpoints([p1, p2, p3, p4, p5, p6, p7]).rank() == 9)


## CONCLUSION: 
## If we fix three points P1, P2, P7 in an arbitrary way and we define
## P4 = (P1||P2)(P1|P7)(P2|P7)-(P1|P2)(P1||P7)(P2|P7)+(P1|P2)(P1|P7)(P2||P7)
## then we define 

## P3 = (P1+P2) *I* (P4+P7)
## P5 = (P1+P4) *I* (P2+P7)
## P6 = (P1+P7) *I* (P2+P4)

## and we get a (C8) configuration. 

## (+++)
## There are some exceptions to study:
## s12=0, s17=0
## or 
## s27=0, s17=0
## or
## s27=0, s12=0

## we define a function to obtain P4 as above:

def fourth_point(Px1, Px2, Px3):
    PP = wedgeProd(Px1, Px2)*scalarProd(Px1, Px3)*scalarProd(Px2, Px3)\
         -scalarProd(Px1, Px2)*wedgeProd(Px1, Px3)*scalarProd(Px2, Px3)\
         +scalarProd(Px1, Px2)*scalarProd(Px1, Px3)*wedgeProd(Px2, Px3)
    return(PP)

## we consider the case 
## (P1 | P2) = 0, (P1 | P7) = 0. The other two are symmetric.
## we redefine the points:

P2 = vector(S, (A2, B2, C2))
P7 = vector(S, (A7, B7, C7))

P1 = wedgeProd(P2, P7)  ## in order to have (P1|P2) = (P1|P7) = 0.

P4 = vector(S, (A4, B4, C4))

## hence P3, P5, P6 are:

P3 = intersect_lines(P1, P2, P4, P7)
P5 = intersect_lines(P1, P4, P2, P7)
P6 = intersect_lines(P1, P7, P2, P4)

## If configuration (8) is given by eigenpoints, we must have
## e1 = e2 = e3 = 0, where:
# e1 = delta1(P3, P1, P4), e2 = delta1(P5, P1, P2), e3 = delta1(P6, P1, P2)

e1 = delta1(P3, P2, P4)
e2 = delta1(P5, P1, P2)
e3 = delta1(P6, P1, P2)


## We have:
## e2 == 0
assert(e2 == S(0))


## we are going to prove that, if s12=0, s17=0, then s27=0.
## e1 can be obtained in different ways:
## as delta1(P3, P2, P4), but also as delta1(P3, P2, P7) or ...
## similarly the others, so we compute three ideals, Je1, Je2 Je3, the
## first is the ideal of all the ways in which delta1(P3, #..) can be computed
## and similarly for the others.

Je1 = S.ideal(delta1(P3, P1, P7), delta1(P3, P1, P4), \
             delta1(P3, P2, P7), delta1(P3, P2, P4)).\
             saturation(matrix([P2, P4, P7]).det())[0]

Je2 = S.ideal(delta1(P5, P1, P2), delta1(P5, P1, P7), \
        delta1(P5, P2, P4), delta1(P5, P4, P7)).\
        saturation(matrix([P2, P4, P7]).det())[0]

Je3 = S.ideal(delta1(P6, P1, P2), delta1(P6, P2, P7), \
             delta1(P6, P1, P4), delta1(P6, P4, P7)).\
             saturation(matrix([P2, P4, P7]).det())[0]

## (Je2 is (0), but we leave it for symmetry)

## Then we see when e1=0, e2=0, e3=0, and precisely, we compute
## the ideal Je1+Je2+Je3.
## It holds: rad(Je1+Je2+Je3) == (s27)

assert((Je1+Je2+Je3).radical() == S.ideal(scalarProd(P2, P7)))

## conclusion: s12 = 0, s17 = 0 ==> s27 = 0.
## for symmetry, it also holds: s12 = 0, s27 = 0 ==> s12 = 0 
## and s27 = 0, s17 = 0 ==> s12 = 0.

## Hence, if we have, among P1, P2, P4, P7, two couples of points of the 
## type (P|Q) = 0, (P|R) = 0, then we have also (Q|R) = 0. As a consequence 
## of [AAA] property, we have that if among P1, P2, P4, P7 we have three
## points P, Q, R such that (P|Q) = 0, (P|R) = 0, (Q|R) = 0, then the 
## remaining point T cannot be orthogonal neither to Q, nor to Q, nor to R.
## Hence fourth_point(P, Q, T) is defined and gives R, 
##  fourth_point(P, R, T) is defined and gives Q, 
##  fourth_point(R, Q, T) is defined and gives P, 
## and the formula fourth_point can always be used.
## Here we also use [BBB], since the three vectors of 
## fourth_point() are linearly independent, and so define a point.



