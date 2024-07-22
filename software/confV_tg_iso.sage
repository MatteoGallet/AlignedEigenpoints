## Here we want to consider the generic case of a V configuration
## where P2 and P4 are on the isotropic conic and the two lines
## rtg1 and rtg2 of the V-conf are tangent to the isotropic conic in P2
## and P4 and rtg1 and rtg2 intersect in P1.
## In this case the matrix Phi(P1,...,P5) has rank 8.
## We construct generic points P1, P2, P3, P4, P5 such that 
## the line P1+P2 is tangent to Ciso in P2, the line P1+P4 
## is tangent to Ciso in P4 and P3 and P5 are on P1+P2
## and P1+P4.
## Successively we prove under which conditions a point on the line
## P2 + P4 is an eigenpoint. 
## The computations done here are also contained in the file 
## "rank_8_twoTangIso.sage" where the computations are done for 
## the specific case P1 = (1, 0, 0). The advantage of this approach 
## is that we have the formula for P6 and P7 in the case in which 
## P2, P4, P6 are aligned.
## The file "rank_8_twoTangIso.sage" also contains the subcase in 
## which there are 6 alignments among the points, which is not considered 
## here, but on the file confV2_tg_iso.sage


## The execution of this file requires about 30'.
print("About 30 minutes of computations")
print("\n\n")

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


## here we start the construction of the generic points P1, ..., P5:
Q1 = vector((1, ii, 0)) ## a known point on Ciso

mon = [x^3, x^2*y, x*y^2, y^3, x^2*z, x*y*z, y^2*z, x*z^2, y*z^2, z^3]

load("../AlignedEigenpoints/software/ancillary.sage")

test = SymbolicCheck()

doLongComputations = false

## construction of a generic point on the isotropic conic:
rt1 = l1*(y-ii*x)+l2*z

rt1.subs({x:Q1[0], y:Q1[1], z:Q1[2]})

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


P2 = PP.subs({l1:u1, l2:u2})

P4 = PP.subs({l1:v1, l2:v2})

## Tangent line to Ciso in P2:
rtg1 = scalarProd(P2, vector((x, y, z))-P2)*(ii)


## Tangent line to Ciso in P4:
rtg2 = scalarProd(P4, vector((x, y, z))-P4)*(ii)

## we compute the intersection point or rtg1 and rtg2 and
## we call it P1:
## (we can assume that (u2*v1-u1*v2) != 0, since P2 and P4 are distinct)

mn2 = matrix([[rtg1.coefficient(xx) for xx in (x, y, z)], \
             [rtg2.coefficient(xx) for xx in (x, y, z)]]).minors(2)


P1 = vector(S, (mn2[2]/(u2*v1 - u1*v2), \
                -mn2[1]/(u2*v1 - u1*v2), \
		mn2[0]/(u2*v1 - u1*v2)))

## P1 is the intersection point of rtg1 and rtg2:

assert((rtg1.subs(substitution(P1)), rtg2.subs(substitution(P1))) \
        == (S(0), S(0)))

## Now we define P3 and P5:

P3 = l2*P2 + l1*P1
P5 = m2*P4 + m1*P1

## ==> NOW WE HAVE CONTRUCTED P1, P2, P3, P4, P5.
## IN A V-CONF S.T. P2 AND P4 ARE ON Ciso AND P1+P2 AND 
## P1+P4 ARE TANGENT TO Ciso.

## We have: P1 = (P2||P4):

assert(matrix([P1, wedgeProd(P2, P4)]).minors(2) == [0, 0, 0])



## In this way we have a V-configuration s.t. Phi(P1, P2, P3, P4, P5)
## has rank 8:

## The following computation requires about 1' 30''

#assert(matrixEigenpoints([P1, P2, P3, P4, P5]).rank() == 8)

## Now we want to see when a point P6 on the line P2+P4 is an eigenpoint:

P6 = w1*P2 + w2*P4  ## in general, this is not an eigenpoint.

## We have to see under which conditions on w1 and w2 P6 is an eigenpoint,
## i.e. the matrix MM = Phi([P1, P2, P3, P4, P5, P6]) has rank <=9.

MM = matrixEigenpoints([P1, P2, P3, P4, P5, P6])

## We select in a "very smart way (!!??)" some order 10-minors of MM in such a 
## way that the ideal they generate gives precisely the conditions 
## for which P6 is an eigenpoint.
## Remember that we can assume (u2*v1-u1*v2) != 0, since P2 and P4 
## are distinct. It turns out that 
## all the determinants of the order 10-minors of MM can be divided by
## (u2*v1-u1*v2)^24.

##
## Here we have a list of 18 lists of 10 rows of MM, which are 
## selected in order to have "good" computations.

Lrw = [[0, 1, 3, 4, 6, 9, 10, 12, 15, 16],\
 [0, 1, 3, 4, 6, 9, 10, 12, 15, 17],\
 [0, 2, 3, 4, 6, 9, 10, 12, 15, 16],\
 [0, 1, 3, 4, 7, 9, 10, 12, 15, 16],\
 [0, 1, 3, 4, 6, 9, 10, 13, 15, 16],\
 [0, 1, 3, 4, 6, 9, 11, 12, 15, 16],\
 [0, 1, 3, 5, 6, 9, 10, 12, 15, 16],\
 [1, 2, 3, 4, 7, 9, 10, 13, 15, 16],\
 [0, 2, 3, 4, 6, 9, 10, 12, 15, 17],\
 [0, 1, 3, 4, 7, 9, 11, 12, 15, 16],\
 [0, 1, 3, 5, 6, 9, 11, 12, 15, 17],\
 [0, 1, 3, 5, 6, 9, 10, 13, 15, 16],\
 [0, 1, 3, 4, 7, 9, 10, 12, 15, 17],\
 [0, 1, 3, 4, 6, 9, 10, 13, 15, 17],\
 [1, 2, 3, 4, 7, 9, 10, 13, 15, 17],\
 [0, 1, 3, 4, 6, 9, 10, 12, 16, 17],\
 [0, 1, 3, 5, 6, 9, 11, 12, 16, 17],\
 [0, 2, 3, 4, 6, 9, 10, 12, 16, 17]]

### the following computation requires about 25 minutes.
## The condition on w1 and w2
Jb = S.ideal(0)
flag = 1
for ll in Lrw:
    print(ll)
    sleep(1)
    ttA = cputime()
    Mx = MM.matrix_from_rows(ll)
    ddt = Mx.det()
    print("computed long determinant")
    sleep(1)
    ddtDiv = ddt.quo_rem((u2*v1-u1*v2)^24) ## it turns out this happens
    if ddtDiv[1] == 0:
        ff = ddtDiv[0]
    else:
        ff = ddt
        print("just in case...") ## In practise, this never happens
    print(cputime()-ttA)
    sleep(1)
    ffId = S.ideal(ff)
    ## we saturate w.r.t. conditions that are surely satisfied.
    ffId = ffId.saturation((u2*v1-u1*v2)*m1*m2*l1*l2*w1*w2)[0]
    Jb = Jb + ffId
    Jb = Jb.saturation((u2*v1-u1*v2)*m1*m2*l1*l2*w1*w2)[0]
    print("computation n. : "+str(flag)+" over "+ str(len(Lrw)))
    print("")
    flag += 1
    sleep(1)


## The only condition we get is: w2*l2*m1 - w1*l1*m2 = 0

assert(Jb == S.ideal(w2*l2*m1 - w1*l1*m2))

## We construct P6 with this condition:
PP6 = P6.subs({w1:l2*m1, w2: l1*m2})

## It holds: 
## PP6 == s11*s15*P3-2*s13*s15*P1+s11*s13*P5

PP6a = scalarProd(P1, P1)*scalarProd(P1, P5)*P3-2*scalarProd(P1, P3)*scalarProd(P1, P5)*P1 + scalarProd(P1, P1)*scalarProd(P1, P3)*P5

assert(S.ideal(matrix([PP6, PP6a]).minors(2)) == S.ideal(0))

## and it also holds:
## PP6 = s15*s34*P2+s13*s25*P4

PP6b = scalarProd(P1, P5)*scalarProd(P3, P4)*P2+scalarProd(P1, P3)*scalarProd(P2, P5)*P4
assert(S.ideal(matrix([PP6, PP6b]).minors(2)) == S.ideal(0))

### ==> We 
## we verify that PP6 (or PP6a or PP6b) is an eigenpoint:

## First we construct a matrix whose determinant is the cubic
## which has P1, P2, P3, P4, P5 as eigenpoints and a row of phi(P6):
## 
MM1 = matrixEigenpoints([P1, P2, P3, P4, P5, PP6]).matrix_from_rows([0, 1, 3, 4, 6, 9, 10, 12, 16])
MM1 = MM1.stack(vector(mon))

## the following computation requires 146 seconds:
print("A computation of 146 seconds")
sleep(1)
cbc = MM1.det()

## Now we verify that PP6 is an eigenpoint:

N1 = matrix([[x, y, z], [cbc.derivative(x), cbc.derivative(y), cbc.derivative(z)]]).subs({x:PP6[0], y:PP6[1], z:PP6[2]})

## PP6 is an eigenpoint:

assert(S.ideal(N1.minors(2)) == S.ideal(0))

## The reminining eigenpoint is P7. We verify that P7 is given by:

## PP7 = s11*s15*P3+s13*s15*P1+s11*s13*P5  ##omogenea!

## and it is also given by:

## PP7b = s15*(s26*s46+s24*s66)*P1+s11*s24*s56*PP6

PP7b = scalarProd(P1, P5)*(scalarProd(P2, PP6)*scalarProd(P4, PP6)+\
       scalarProd(P2, P4)*scalarProd(PP6, PP6))*P1+\
       scalarProd(P1, P1)*scalarProd(P2, P4)*scalarProd(P5, PP6)*PP6



PP7 = scalarProd(P1, P1)*scalarProd(P1, P5)*P3+scalarProd(P1, P3)*scalarProd(P1, P5)*P1 + scalarProd(P1, P1)*scalarProd(P1, P3)*P5

## PP7 and PP7b are the same point:

assert(matrix([PP7, PP7b]).minors(2) == [0, 0, 0])


N2 = matrix([[x, y, z], [cbc.derivative(x), cbc.derivative(y), cbc.derivative(z)]]).subs({x:PP7[0], y:PP7[1], z:PP7[2]})

## PP7 is an eigenpoint:

assert(S.ideal(N2.minors(2)) == S.ideal(0))

## The alignments of the points are:
## [(1, 2, 3), (1, 4, 5), (1, 6, 7), (2, 4, 6)]

assert(allignments([P1, P2, P3, P4, P5, PP6, PP7]) == [(1, 2, 3), (1, 4, 5), (1, 6, 7), (2, 4, 6)])

## First conclusion:
## In the case of a V configuration tangent to Ciso, 
## there is a sub-case in which we have 4 alignments and 
## the eigenpints PP6 and PP7 are given by the above formulas.

## In addition, it is not possible to have further collinearities 
## among the points
## 
## P3, P5, PP7 cannot be collinear:
assert(S.ideal(det(matrix([P3, P5, PP7]))).saturation(m2 * m1 * l2 * l1 * (u2*v1 - u1*v2))[0] == S.ideal(1))

## P3, P4, PP7 cannot be collinear:
assert(S.ideal(det(matrix([P3, P4, PP7]))).saturation(m2 * m1 * l2 * l1 * (u2*v1 - u1*v2))[0] == S.ideal(1))

## P3, P5, PP6 cannot be collinear:
assert(S.ideal(det(matrix([P3, P5, PP6]))).saturation(m2 * m1 * l2 * l1 * (u2*v1 - u1*v2))[0] == S.ideal(1))

## P2, P5, PP7 cannot be collinear:
assert(S.ideal(det(matrix([P2, P5, PP7]))).saturation(m2 * m1 * l2 * l1 * (u2*v1 - u1*v2))[0] == S.ideal(1))


### This concludes the case config V tangent to Ciso.
### There are not sub-cases to consider.

print("End of computations.")