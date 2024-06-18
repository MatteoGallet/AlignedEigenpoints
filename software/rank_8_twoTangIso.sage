## This is the case in which the 5 eigenpoints P1, P2, P3, P4, P5 are
## such that P1+P2 is tangent in the point P2 to the isotropic conic, 
## as well as P1+P4 is tangent to the isotropic conic in P4.
## 
## Here we prove that \\label{lemma:special_case_rank_8}
## and we prove propositions: 
## \Cref{prop:rk8_2B}

## Previous name of the file: contiCasoDegenere2.sage
## connected with the files casoDegenere2.tex, casoDegenere2.pdf

varAn1 = ["s"+str(i)+str(j) for i in range(1,8) for j in range(i, 8)]
varAn2 = ["x", "y", "z", "u1", "u2", "v1", "v2", "w1", "w2", "l1", "l2"]
varAn3 = ["a"+str(i) for i in range(10)]
varAn4 = [str(XX)+str(i) for i in range(1, 8) for XX in ["A", "B", "C"]]
var("xx")
K = QQ
K.<ii> = NumberField(xx^2+1)
S = PolynomialRing(K, varAn1+varAn2+varAn3+varAn4)
S.inject_variables(verbose=false)

P1 = vector(S, (1, 0, 0))
P2 = vector(S, (0, ii, 1))
P4 = vector(S, (0, -ii, 1))
P3 = u1*P1 + u2*P2
P5 = v1*P1 + v2*P4

load("auxiliaryProcedures.sage")

## a remark on delta1 and delta2:
## delta1(P1, P2, P4) is not zero, 
## while delta2(P1, P2, P3, P4, P5) is 0.
assert(delta1(P1, P2, P4) != 0)
assert(delta2(P1, P2, P3, P4, P5) == 0)



## this method is used to find the orthogonality of the lines:
def ortLines(p1, p2, q1, q2):
    rt1 = det(matrix([p1, p2, (x, y, z)]))
    rt2 = det(matrix([q1, q2, (x, y, z)]))
    cf1 = [rt1.coefficient(xx) for xx in (x, y, z)]
    cf2 = [rt2.coefficient(xx) for xx in (x, y, z)]
    return(scalarProd(cf1, cf2) == 0)

mon = [x^3, x^2*y, x*y^2, y^3, x^2*z, x*y*z, y^2*z, x*z^2, y*z^2, z^3]

### In this configuration, the line P1+P2 is tangent to
### the isotropic conic in P2 and the line P1+P4 is tangent
### to the isotropic conic in P4 (for all P3 and P5)

M = matrixEigenpoints([P1, P2, P3, P4, P5])

## 
## M[2] is the zero row
assert(M[2] == vector(S, [0 for i in range(10)]))
## M[3]-ii*M[4] is zero
assert(M[3] - ii*M[4] == vector(S, [0 for i in range(10)]))
## u2*M[6]-ii*u2*M[7]+u1*M[8] is zero
assert(u2*M[6] -ii*u2*M[7] +u1*M[8] == vector(S, [0 for i in range(10)]))
## M[9]+ii*M[10] is zero
assert(M[9] + ii*M[10] == vector(S, [0 for i in range(10)]))
## v2*M[12]+ii*v2*M[13]+v1*M[14] is zero
assert(v2*M[12] +ii*v2*M[13] +v1*M[14] == vector(S, [0 for i in range(10)]))

## therefore in the matrix M we can erase the rows: 2, 3, 6, 9, 12

M = M.matrix_from_rows([0, 1, 4, 5, 7, 8, 10, 11, 13, 14])

## The rank of M cannot be < 8.
m8 = M.minors(8)
assert(S.ideal(m8).saturation(u1*u2*v1*v2)[0] == S.ideal(S.one()))

### Hence we obtain the proof of proposition 4.23(?) (labelled with 
### \label{lemma:special_case_rank_8})

## The matrix M can be manipulated in order to have a simpler matrix.
## Two function for elementary row operations:
def rid(v1, v2, pv):
    return(v2-(v2[pv]/v1[pv])*v1)

def ridM(M, nRiga, pv):
    v1 = M[nRiga]
    M1 = [v1] + [rid(v1, M[i], pv) for i in range(M.nrows()) if i != nRiga]
    return(matrix(M1))

## 


M1 = ridM(M, 0, 1)
M1 = ridM(M1, 1, 4)
M1 = ridM(M1, 2, 2)
M1 = ridM(M1, 3, 3)
M1 = ridM(M1, 6, 5)
M1 = ridM(M1, 7, 6)

### now M1 (which is equiv to M) is better. It has two zero rows:

assert(M1[7]) == vector(S, [0 for _ in range(10)])
assert(M1[9]) == vector(S, [0 for _ in range(10)])

## We erase the two rows from M1
M1 = M1.matrix_from_rows([0, 1, 2, 3, 4, 5, 6, 8])

M1.rescale_row(0, 1/(-6*ii))
M1.rescale_row(2, 1/(-3))
M1.rescale_row(1, 1/(2*ii))
M1.add_multiple_of_row(6, 0, -2*u1*u2^2)
M1.add_multiple_of_row(7, 0, -2*v1*v2^2)
M1 = (M1.matrix_from_rows(range(6)).stack(vector(S, list(map(S, list(M1[6]*1/(u1*u2))))))).\
stack(vector(S, list(map(S, list(M1[7]*1/(v1*v2))))))

## we make a copy of M1 (you never know)
M2 = copy(M1)

## M1 (a 8x10 matrix) is of rank 8 and is the matrix of a linear 
## system whose solutions give all the cubics with the above P1, ..., P5
## as eigenpoints.
## We solve the system with Cramer. We use 
## the information that the determinant of the matrix obtained with
## the columns 1, 3, 4, 5, 6, 7, 8, 9 of M1 cannot be 0.

## Now some computations in order to find the \infty^1 family 
## of cubic curves corresponding to the solution of the system whose 
## coefficients matrix is M1

coeffM = M1.matrix_from_columns([1, 3, 4, 5, 6, 7, 8, 9])
assert(S.ideal(coeffM.det()).saturation(u1*u2*v1*v2)[0] == S.ideal(S.one()))

matB = -M1.matrix_from_columns([0, 2])*matrix(S, [[a0], [a2]])

## solution of the system:
sol = solve_with_Cramer(coeffM, matB)

## let us verify the assert, i.e. that sol is the solution of the system:
assert(vector(coeffM*sol-det(coeffM)*vector(matB.transpose()))== vector(S, [0 for _ in range(8)]))

## the solution of the system depends on the two free parameters a0, a2.
## We rewrite sol, putting matA.det() in position 0 and 2:
i1, i2 = 0, 2

sol1 = list(sol)
sol1 = vector(S, sol1[:i1]+[a0*coeffM.det()]+sol1[i1:i2-1]+[a2*coeffM.det()]+sol1[i2-1:])

## sol1 is the solution of the system of M1
assert(M1*sol1 == vector(S, [0 for _ in range(8)]))

## The family of cubics is:
cb = sum([sol1[i]*mon[i] for i in range(10)])

cb = S.ideal(cb).saturation(u1*u2*v1*v2)[0].gens()[0]

## cb is the linear combination of the following two cubics:
cbb1 = x*(x^2+3/2*y^2+3/2*z^2)
cbb2 = (y + ii*z) * (y + (-ii)*z) * (y*u2*v1 + (-ii)*z*u2*v1 - y*u1*v2 + (-ii)*z*u1*v2 + (2*ii)*x*u2*v2)
## i.e. cb1 = w1*cbb1 + w2*ccb2
## as is obtained by the following computation, which gives that
## w2 = 3*a0 - 2*a2, w1 = (-ii)*4*u2*v2*a0:

assert(S.ideal([(w1*cbb1+w2*cbb2 -cb).coefficient(mm) for mm in mon]).saturation(u1*u2*v1*v2)[0] == S.ideal(w2 - 3*a0 + 2*a2, 4*u2*v2*a0 + (-ii)*w1))
## 
assert(cb == ((-ii)*4*u2*v2*a0)*cbb1 + (3*a0 - 2*a2)*cbb2)

## hence we re-define cb:

cb = l1*cbb1+l2*cbb2

## The above is the equation of ALL the cubics with eigenpoints 
## P1, P2, P3, P4, P5 (it depends on u1, u2, v1, v2 because of P3 and 
## P5) and depends of l1 and l2 (because the matrix Phi(P1, ..., P5)
## has rank 8).
## The equation is inserted in the paper, see the label ***C1_C2***


## cbb1 is a cubic which has as eigenpoints P1 and the line P1+P2 and 
## P1+P4
e_cbb1 = S.ideal(matrix(S, [(x, y, z), [cbb1.derivative(xx) for xx in (x, y, z)]]).minors(2)).saturation(u1*u2*v1*v2)[0]

assert(e_cbb1 == S.ideal(y^2*z + z^3, y^3 + y*z^2))
## hence we can assume l2 != 0 

## We summarize what we have obtained so far:
## We have that all the cubics that come from a rank 8 matrix 
## are given by a linear combination of cbb1 and cbb2
## this family of cubics is called cb and depends on two parameters:
## l1 and l2 (and on u1, u2, v1, v2).

## We compute the eigenpoints of cb, using the definition:

ej = S.ideal(matrix(S, [(x, y, z), [cb.derivative(xx) for xx in (x, y, z)]]).minors(2))

ej = ej.saturation(u1*u2*v1*v2*l2)[0]

## we erase from ej the known eigenpoints:
for pp in [P1, P2, P3, P4, P5]:
    ej = ej.saturation(S.ideal(matrix([pp, (x, y, z)]).minors(2)))[0]

## we get precisely a new ideal:
ep = ej.radical().primary_decomposition()

assert(len(ep) == 1)

## ep[0] is the intersection of a line and a conic 
## and gives the eigenpoints P6 and P7

rt = y*u2*v1 + (-ii)*z*u2*v1 + y*u1*v2 + ii*z*u1*v2 

conic = 12*x*y*u1*v2*l2 + (12*ii)*x*z*u1*v2*l2 + (-8*ii)*x^2*u2*v2*l2 + (4*ii)*y^2*u2*v2*l2 + (4*ii)*z^2*u2*v2*l2 + 3*y^2*l1 + 3*z^2*l1


assert(ep[0] == S.ideal(conic, rt).saturation(u1*u2*v1*v2)[0])
## ep[2] is an ideal of 2 points which are on the line rt and rt passes 
## through P1

## rt, which contains P6 and P7,  passes through P1:
assert(rt.subs(substitution(P1)).is_zero())
## hence we always have the alignment: P1, P6, P7.

## The line  rt = y*u2*v1 + (-ii)*z*u2*v1 + y*u1*v2 + (ii)*z*u1*v2
## is orthogonal to the line  P3+P5, as follows from these computations:

r35 = matrix([P3, P5, (x, y, z)]).det()

## r35 and rt are orthogonal:
assert(scalarProd([r35.coefficient(xx) for xx in (x, y, z)], \
                  [rt.coefficient(xx) for xx in (x, y, z)]).is_zero())

         
### Now we study the possible configurations of eigenpoints of 
### the family cb
  
### In general, the conic which generate the ideal ep[0] is irreducible:

assert(conic.is_prime())

## This means that the coordinates of P6 and P7 are not rational in the 
## parameters.

## Now we want to see if there are other alignments among the points 
## We have to consider three cases:
##  (1) P2, P4, P6 aligned
##  (2) P2, P5, P6 aligned
##  (3) P3, P5, P6 (or P3, P5, P7) aligned 


###   ===> case 1: P2, P4, P6 collinear.

## equation line P2+P4: x = 0
assert(matrix([P2, P4, (x, y, z)]).det() == 2*ii*x)

## construction of the point P6 intersection of the lines rt and P2+P4:
P6 = vector(S, (0, rt.coefficient(z), -rt.coefficient(y)))

## P6 always exists:
assert(S.ideal(list(P6)).saturation(u1*u2*v1*v2)[0] == S.ideal(S.one()))

## if P6 is an eigenpoint, it must be a point of the conic. 
## Hence u2*v2*l2 + (-3/4*ii)*l1 = 0, i.e. 
## l1 = u2*v2, l2 = 3/4*ii
assert(conic.subs(substitution(P6)) == ((16*ii))*v2*v1*u2*u1*(u2*v2*l2 + (-3/4*ii)*l1))

## the family of  cubics in which P6 is aligned with P2 and P4:
cb1 = S(cb.subs({l1: u2*v2, l2: 3/4*ii}))
## cb1 is smooth
assert(sing_loc(cb1).saturation(u1*u2*v1*v2)[0] == S.ideal(x, y, z))

## construction of P7 for cb1:

e_cb1 = S.ideal(matrix([(x, y, z), [cb1.derivative(xx) for xx in (x, y, z)]]).minors(2)).\
saturation(u1*u2*v1*v2)[0]
for pp in [P1, P2, P3, P4, P5, P6]:
    e_cb1 = e_cb1.saturation(S.ideal(matrix([pp, (x, y, z)]).minors(2)))[0]

pl1, pl2 = tuple(e_cb1.gens()[:2])


slz = matrix([[pl1.coefficient(xx) for xx in (x, y, z)], \
              [pl2.coefficient(xx) for xx in (x, y, z)]]).minors(2)

P7 = 1/6*vector(S, [slz[2], -slz[1], slz[0]])

## P7 always exists:
assert(S.ideal(list(P7)).saturation(u1*u2*v1*v2)[0] == S.ideal(S.one()))

## In this situation, when P2, P4, P6 are aligned, we have:

## P6 = (0, (-ii)*u2*v1 + ii*u1*v2, -u2*v1 - u1*v2)
## P7 = ((3*ii)*u1*v1, -u2*v1 + u1*v2, ii*u2*v1 + ii*u1*v2)
assert(P6 == vector((0, (-ii)*u2*v1 + ii*u1*v2, -u2*v1 - u1*v2)))
assert(P7 == vector(S, ((3*ii)*u1*v1, -u2*v1 + u1*v2, ii*u2*v1 + ii*u1*v2)))

## 
print("Points P6 and P7 when P2, P4, P6 are aligned:")
print(P6)
print(P7)

print("")

## P6 and P7 are obtained from the formulas:
## P6 = s11*s15*P3-2*s13*s15*P1+s11*s15*P5

P6a = scalarProd(P1, P1)*scalarProd(P1, P5)*P3-\
       2*scalarProd(P1, P3)*scalarProd(P1, P5)*P1 +\
       scalarProd(P1, P1)*scalarProd(P1, P3)*P5

assert(S.ideal(matrix([P6, P6a]).minors(2)) == S.ideal(0))


## P7 = s11*s15*P3+s13*s15*P1+s11*s15*P5

P7a = scalarProd(P1, P1)*scalarProd(P1, P5)*P3+\
      scalarProd(P1, P3)*scalarProd(P1, P5)*P1 + \
      scalarProd(P1, P1)*scalarProd(P1, P3)*P5

assert(S.ideal(matrix([P7, P7a]).minors(2)) == S.ideal(0))

## These formulas are proved, for the general case, in the file:
## confV_tg_iso.sage


## The seven eigenpoints of cb1 have 4 alignments:
lp = [P1, P2, P3, P4, P5, P6, P7]

assert(allignments(lp) == [(1, 2, 3), (1, 4, 5), (1, 6, 7), (2, 4, 6)])

## the seven eigenpoints cannot have more then the 4 alignments above:
for i in range(5):
    for j in range(i+1, 6):
        for k in range(j+1, 7):
            if not ((i+1, j+1, k+1) in [(1, 2, 3), (1, 4, 5), (1, 6, 7), (2, 4, 6)]):
                assert(S.ideal(matrix([lp[i], lp[j], lp[k]]).det()).saturation(u1*u2*v1*v2)[0] == S.ideal(S.one()))



## We have the following orthogonalities:
## P1+P2 ort P2+P4
## P1+P6 ort P2+P4
## P1+P4 ort P2+P4
## P1+P6 ort P3+P5

assert(scalarProd(wedgeProd(P1, P2), wedgeProd(P2, P4)) == 0)
assert(scalarProd(wedgeProd(P1, P6), wedgeProd(P2, P4)) == 0)
assert(scalarProd(wedgeProd(P1, P4), wedgeProd(P2, P4)) == 0)
assert(scalarProd(wedgeProd(P1, P6), wedgeProd(P3, P5)) == 0)

## hence we have:

## P1 = (P2 || P4)

assert(matrix([P1, wedgeProd(P2, P4)]).minors(2) == [0, 0, 0])

## conclusion:
## In the case P2, P4, P6 aligned, 
## we have P1 = (P2 || P4)

## we have a configuration of type (5) and we cannot have other configurations
## of the 7 eigenpoints (but there is a configuration with lines of eigepoints).


### 
###   ===> case 2: case in which P2, P5, P6 are collinear

## equation line P2+P5: r25 = y*v1 + (-ii)*z*v1 + (2*ii)*x*v2
r25 = matrix([P2, P5, (x, y, z)]).det()
assert(r25 == y*v1 + (-ii)*z*v1 + (2*ii)*x*v2)

## construction of the point P6 intersection of the lines rt and P2+P5:

slz6 = matrix([[r25.coefficient(xx) for xx in (x, y, z)], [rt.coefficient(xx) for xx in (x, y, z)]]).minors(2)

P6 = vector(S, (slz6[2], -slz6[1], slz6[0]))
## P6 is the intersection point of rt and r25
assert(r25.subs(substitution(P6)).is_zero() and  rt.subs(substitution(P6)).is_zero())

## P6 always exists:
assert(S.ideal(list(P6)).saturation(u1*u2*v1*v2)[0] == S.ideal(S.one()))

## if P6 is an eigenpoint, it must be a point of the conic. Hence 
# u1*v1*l2 + u2*v2*l2 + (-3/4*ii)*l1)
assert(conic.subs(substitution(P6)) == \
((-64*ii)) * v1 * u2 * u1 * v2^3 * (u1*v1*l2 + u2*v2*l2 + (-3/4*ii)*l1))


## Hence the substitution is:
## {l1: u1*v1+u2*v2: l2: 3/4*ii}

cb1 = S(cb.subs({l1: u1*v1+u2*v2, l2: 3/4*ii}))

## if we assume (2*u1*v1+3*u2*v2)*(u1*v1+u2*v2)*(u1*v1+3*u2*v2) != 0, 
## then  cb1 is smooth:
irr = sing_loc(cb1).saturation(u1*u2*v1*v2)[0].saturation(2*u1*v1+3*u2*v2)[0].\
saturation(u1*v1+u2*v2)[0].saturation(u1*v1+3*u2*v2)[0]
assert(irr == S.ideal(x, y, z))

## construction of P7:

e_cb1 = S.ideal(matrix([(x, y, z), [cb1.derivative(xx) for xx in (x, y, z)]]).minors(2)).\
saturation(u1*u2*v1*v2)[0]
for pp in [P1, P2, P3, P4, P5, P6]:
    e_cb1 = e_cb1.saturation(S.ideal(matrix([pp, (x, y, z)]).minors(2)))[0]

pl1, pl2 = tuple(e_cb1.gens()[:2])


slz = matrix([[pl1.coefficient(xx) for xx in (x, y, z)], \
              [pl2.coefficient(xx) for xx in (x, y, z)]]).minors(2)

P7 = vector(S, [slz[2], -slz[1], slz[0]])

## P7 is an eigenpoint of cb1:
assert(S.ideal(matrix([(x, y, z), [cb1.derivative(xx) for xx in (x, y, z)]]).minors(2)).subs(substitution(P7))==S.ideal(0))

## P7 always exists:
assert(S.ideal(list(P7)).saturation(u1*u2*v1*v2)[0] == S.ideal(S.one()))

## here are P6 and P7:

## P6 = ((2*ii)*u1*v1*v2, -2*u2*v1*v2 + 2*u1*v2^2, (2*ii)*u2*v1*v2 + (2*ii)*u1*v2^2),
##P7 = ((2*ii)*u1*v1, -u2*v1 + u1*v2, ii*u2*v1 + ii*u1*v2)

assert(P6 == vector(S, ((2*ii)*u1*v1*v2, -2*u2*v1*v2 + 2*u1*v2^2, (2*ii)*u2*v1*v2 + (2*ii)*u1*v2^2)))
assert(P7 == vector(S,  ((2*ii)*u1*v1, -u2*v1 + u1*v2, ii*u2*v1 + ii*u1*v2)))

## The seven eigenpoints of cb1 have the 6 alignments [(1, 2, 3), (1, 4, 5), (1, 6, 7), (2, 5, 6), (3, 4, 6), (3, 5, 7)]
lp = [P1, P2, P3, P4, P5, P6, P7]
assert(allignments(lp) == [(1, 2, 3), (1, 4, 5), (1, 6, 7), (2, 5, 6), (3, 4, 6), (3, 5, 7)])



### end of case P2, P5, P6 collinear. We have configuration (8)

## we have the following orthogonalities:
assert(scalarProd(wedgeProd(P1, P4), wedgeProd(P3, P4))==0)
assert(scalarProd(wedgeProd(P1, P2), wedgeProd(P2, P5))==0)
assert(scalarProd(wedgeProd(P3, P5), wedgeProd(P1, P6))==0)


## We also have the relation P3 == QQ3, where:
QQ3 = wedgeProd(P1, P5)*scalarProd(P1, P6)*scalarProd(P5, P6)-\
      scalarProd(P1, P5)*wedgeProd(P1, P6)*scalarProd(P5, P6)+\
      scalarProd(P1, P5)*scalarProd(P1, P6)*wedgeProd(P5, P6)

assert(QQ3 != vector(S, (0, 0, 0)))
assert(matrix([P3, QQ3]).minors(2) == [0, 0, 0])

###############################
###############################

## case 3: P3, P5, P7 aligned.

## case P3, P5 aligned with a point of the line P1+P6+P7.
## Here we call P7 the point of rt which is aligned with 
## P3 and P5

## equation line P3+P5: r35 = y*u2*v1 + (-ii)*z*u2*v1 - y*u1*v2 + (-ii)*z*u1*v2 + (2*ii)*x*u2*v2
r35 = matrix([P3, P5, (x, y, z)]).det()
assert(r35 == y*u2*v1 + (-ii)*z*u2*v1 - y*u1*v2 + (-ii)*z*u1*v2 + (2*ii)*x*u2*v2)

## construction of the point P7 intersection of the lines rt and P3+P5:

slz7 = matrix([[r35.coefficient(xx) for xx in (x, y, z)], [rt.coefficient(xx) for xx in (x, y, z)]]).minors(2)
P7 = vector(S, (slz7[2], -slz7[1], slz7[0]))
## P7 is the intersection point of rt and r35
assert(r35.subs(substitution(P7)).is_zero() and  rt.subs(substitution(P7)).is_zero())

## P7 always exists:
assert(S.ideal(list(P7)).saturation(u1*u2*v1*v2)[0] == S.ideal(S.one()))

## if P7 is an eigenpoint, it must be a point of the conic. Hence
## u1*v1*l2 + u2*v2*l2 + (-3/4*ii)*l1
## indeed:

assert(conic.subs(substitution(P7)) == \
((-64*ii)) * v1 * u1 * v2^3 * u2^3 * (u1*v1*l2 + u2*v2*l2 + (-3/4*ii)*l1))

## Hence the substitution is:
## {l1: u1*v1+u2*v2, l2: 3/4*ii}

cb1bis = S(cb.subs({l1: u1*v1+u2*v2, l2: 3/4*ii}))

## we obtain that the new cubic is the same the cubic cb1 of the previous case, so this case is already studied.
assert(cb1bis == cb1)

##
## In conclusion, we have two possible configurations: (5) and (8)
## (and a case of two lines of eigenpoints)
## the line P3+P5 and the line P1+P6+P7 are always orthogonal



## Nota:
## Un esempio di cubica singolare che salta fuori e':
## x^3 + 9/2*x*y^2 + (-21/4*ii)*y^3 + 15/4*y^2*z + 9/2*x*z^2 + (-21/4*ii)*y*z^2 + 15/4*z^3
## che e' una cubica che ha 7 autopunti in conf (8) e uno di questi e' il
## punto singolare.
