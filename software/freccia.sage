print("We define a ring with sufficientely many variables.")

varAn1 = ["s"+str(i)+str(j) for i in range(1,8) for j in range(i, 8)]
varAn2 = ["x", "y", "z", "u1", "u2", "v1", "v2", "w1", "w2", "l1", "l2"]
varAn3 = ["a"+str(i) for i in range(10)]
varAn4 = [str(XX)+str(i) for i in range(1, 8) for XX in ["A", "B", "C"]]

K = QQ
S = PolynomialRing(K, varAn1+varAn2+varAn3+varAn4)
S.inject_variables(verbose=false)

print("")
print("Here we want to compute condition delta1 and condition delta2, i.e. ")
print("   ==> WE WANT TO SEE WHEN FIVE POINTS P1, P2, P3, P4, P5 WITH THE CONDITIONS:")
print("   ==> P3 COLLINEAR WITH P1 AND P2, P5 COLLINEAR WITH P1 AND P4, ARE EIGENPOINTS.")
print("")

print("We define 3 generic points: P1, P2, P4")
print("")
P1 = vector((A1, B1, C1))


#load("auxiliaryProcedures1.sage")

P2 = vector((A2, B2, C2))
P4 = vector((A4, B4, C4))

print("")
print("")
print("We define P3 and P5 collinear with P1, P2 and with P1, P4")

P3 = u1*P1+u2*P2
P5 = v1*P1+v2*P4

## an auxiliary procedure:
def substitution(P):
    return {x: P[0], y:P[1], z:P[2]}


def phi_p(Pt):
    F = a0*x^3 + a1*x^2*y + a2*x*y^2 + a3*y^3 + a4*x^2*z \
    + a5*x*y*z + a6*y^2*z + a7*x*z^2 + a8*y*z^2 + a9*z^3
    ## the ideal of the minors:
    R = F.parent()
    Maux = matrix([[x, y, z], [diff(F, x), diff(F, y), diff(F, z)]])
    Jaux = R.ideal(Maux.minors(2))
    Laux = Jaux.subs(substitution(Pt)).gens()
    vv = [a0,a1,a2,a3,a4,a5,a6,a7,a8,a9]
    terna = []
    for aa in Laux:
        terna.append(vector([aa.coefficient(cf) for cf in vv]))
    return(terna)


print("We assume that the vectors:")
print("phi_p(P1)[0], phi_p(P1)[1], phi_p(P2)[0], phi_p(P2)[1],")
print("phi_p(P3)[0], phi_p(P3)[1], phi_p(P4)[0], phi_p(P4)[1], phi_p(P5)[0]\n\
are linearly independent and we contruct a 9 X 10 matrix  M \n\
with these nine rows.")

M = matrix([phi_p(P1)[0], phi_p(P1)[1], phi_p(P2)[0], phi_p(P2)[1], \
phi_p(P3)[0], phi_p(P3)[1], phi_p(P4)[0], phi_p(P4)[1], phi_p(P5)[0]])


################################
## We assume M is a 9x10 matrix. The following method computes 
## all the 9x9 minors of M.
## Same result as M.minors(9), but here we print the computation time 
## of each of the 10 minors.
## Times of computations on my pc:
## 1/10: 31, 2/10: 83, 3/10: 65, 4/10: 141, 5/10: 80, 
## 6/10: 75, 7/10: 124, 8/10: 251, 9/10: 79, 10/10: 407. 

def compute_9_minors(M):
    initialTime = cputime()
    print("We compute the 9x10 minors of the matrix M")
    coffeeBreak = 1
    columnsM = range(10)
    minorsM = []
    for i in range(10):
        cpCol = list(columnsM[:])
        cpCol.remove(i)
        ttA = cputime()
        print("    We start the computation number "+str(i+1)+"/10 ...")
        sleep(coffeeBreak)
        dMx = (M.matrix_from_columns(cpCol)).det()
        minorsM.append(dMx)
        print("    ... end of this partial computation. Time: "\
               +str(cputime()-ttA))
        sleep(coffeeBreak)
    print("global time: "+ str(cputime()-initialTime))
    return(minorsM)

print("We compute the 10 minors of the matrix M")
minM = compute_9_minors(M)


print("We compute the gcd of the ten minors minM")
ttA = cputime()
commFct = gcd(minM)
print("Time of computation: "+str(cputime()-ttA))

print("The common factor is of the form \n\
       A1*A2*A3*A4*u1^2*u2^2*v1*v2*[P1,P2,P4]^2\n\
       (where [P1,P2,P4] is the determinant of the matrix whose rows are\n\
       P1, P2, P4")

print("Since u1, u2, and [P1,P2,P4] are surely not 0, we can erase then \n\
      from the minors and compute the determinant of the matrix MA\n\
      with this semplification (where MA is M plus the line phi(P5)[1]")

d = commFct/(A1*A2*A4*(u1*A1+u2*A2))
minMreduced = list(map(lambda uu: S(uu/d), minM))

detMA = 0
for i in range(10):
    detMA = detMA + (-1)^i*phi_p(P5)[1][i]*minMreduced[i]

print(detMA.factor())

print("We consider only relevant factors of detMA:")

dd = A4 * A2 * A1 * v2 * v1 * (v1*A1 + v2*A4) * (u1*A1 + u2*A2) * (matrix([P1, P2, P4]).det())^3
detMAreduced = S(detMA/dd)


##  The scalar product (of vectors of three components):
def scalarProd(P, Q):
    return(P[0]*Q[0]+P[1]*Q[1]+P[2]*Q[2])

## we define an ideal 
Jrel = S.ideal(s11-scalarProd(P1, P1), s12-scalarProd(P1, P2), 
       s14-scalarProd(P1, P4), s22-scalarProd(P2, P2), 
       s24-scalarProd(P2, P4), s44-scalarProd(P4, P4))


print(list(map(lambda vv: Jrel.reduce(vv), list(map(lambda uu: uu[0], list(detMAreduced.factor()))))))
