## Here we assume that it is defined the ring:
## K = QQ
## R.<a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,x,y,z,\
## u1,u2,v1,v2,w1,w2,q1,q2,A1,B1,C1,A2,B2,C2,\
## A3,B3,C3,A4,B4,C4,A5,B5,C5, A6, B6, C6> = PolynomialRing(K)

## we assume that a point P1 is defined in a ring.
R = P1.base_ring()

##  The scalar product (of vectors of three components):
def scalarProd(P, Q):
    return(P[0]*Q[0]+P[1]*Q[1]+P[2]*Q[2])
    
## the wedge product of two vectors of K^3:
def wedgeProd(pa, pb):
    aux = matrix([pa, pb]).minors(2)
    return(vector((aux[2], -aux[1], aux[0])))

## the function delta1:
def delta1(pa, pb, pc):
    return(scalarProd(pa, pa)*scalarProd(pb, pc)-\
    scalarProd(pa, pb)*scalarProd(pa, pc))

## the function delta1b:
def delta1b(pa, pb, pc):
    return(scalarProd(pa, pa)*scalarProd(pb, pc)+\
    scalarProd(pa, pb)*scalarProd(pa, pc))
 
## the function delta2:
def delta2(pt1, pt2, pt3, pt4, pt5):
    return(scalarProd(pt1, pt2)*scalarProd(pt1, pt3)*scalarProd(pt4, pt5)-\
scalarProd(pt1, pt4)*scalarProd(pt1, pt5)*scalarProd(pt2, pt3))


## the function (pt1|pt1)(pt2|pt2)-(pt1|pt2)^2
def sigma(pt1, pt2):
    return(scalarProd(pt1, pt1)*scalarProd(pt2, pt2)-scalarProd(pt1, pt2)^2)




## an auxiliary procedure:
def substitution(P):
    return {x: P[0], y:P[1], z:P[2]}

## given a list Lp of points P1, P2, ... gives the matrix
## of the linear system in the coefficient of the cubic
## which gives the conditions that P1, P2, ... are eigenpoints
## we assume R is defined.
def matrixEigenpoints(Lp):
    F = a0*x^3 + a1*x^2*y + a2*x*y^2 + a3*y^3 + a4*x^2*z \
    + a5*x*y*z + a6*y^2*z + a7*x*z^2 + a8*y*z^2 + a9*z^3
    ## the ideal of the minors:
    Maux = matrix([[x, y, z], [diff(F, x), diff(F, y), diff(F, z)]])
    Jaux = R.ideal(Maux.minors(2))
    Laux = [Jaux.subs(substitution(Pt)).gens() for Pt in Lp]
    H = []
    for aa in Laux:
        H += aa
    vv = [a0,a1,a2,a3,a4,a5,a6,a7,a8,a9]
    matr = []
    for ii in range(len(H)):
        row = []
        for cf in vv:
            row.append(H[ii].coefficient(cf))
        matr.append(row)
    return(matrix(matr))


## given a rank 9 matrix of elements of K, computes the cubic curve.
def cubic_from_matrix(Ms):
    Ms = matrix(K, Ms)
    if rank(Ms) != 9:
        return("wrong rank!")
    basis = Ms.right_kernel().basis()
    ## basis dovrebbe avere un solo vettore:
    if len(basis) != 1:
        return("wrong basis!")
    mon = [x^3, x^2*y, x*y^2, y^3, x^2*z, x*y*z, y^2*z, x*z^2, y*z^2, z^3]
    Fs = add([basis[0][i]*mon[i] for i in range(10)])
    return(Fs)

## given a cubic (with coefficients in K) computes its singular locus.
def sing_loc(Fs):
    PtSing = R.ideal([Fs, Fs.derivative(x), Fs.derivative(y), Fs.derivative(z)])
    return(PtSing.radical())

## given a cubic (with coefficients in K) computes the eigenpoints
def eigenpoints(Fs, withoutRadical=false):
    Ns = matrix([[x, y, z], [diff(Fs, x), diff(Fs, y), diff(Fs, z)]])
    Js = R.ideal(Ns.minors(2))
    if withoutRadical:
        Pd = Js.primary_decomposition()
    else:
        Pd = Js.radical().primary_decomposition()
    return(Pd)

## given a list of ideals of points, returns the points in coordinates.
def points_from_id(Lid):
    n = len(Lid)
    newEig = [vector((Lid[i].reduce(x), Lid[i].reduce(y), \
        Lid[i].reduce(z))).subs({x:1, y:1, z:1}) for i in range(n)]
    return(newEig)

## given a list of ideals of points, returns the points in coordinates.
def points_from_id6(Lid):
    n = len(Lid)
    newEig = [vector((Lid[i].reduce(A6), Lid[i].reduce(B6), \
        Lid[i].reduce(C6))).subs({A6:1, B6:1, C6:1}) for i in range(n)]
    return(newEig)


## given a list of points, returns the triplets of collinear points.
def allignments(Lpt):
    n = len(Lpt)
    triplet = list(Combinations(n,3))
    allignments = []
    for tr in triplet:
        threePoints = (Lpt[tr[0]], Lpt[tr[1]], Lpt[tr[2]])
        if det(matrix(threePoints)) == 0:
            allignments.append((tr[0]+1, tr[1]+1, tr[2]+1))
    return(allignments)


## given a point ll and a list of points LL2, gives the index
## of the point ll in the list LL2.
def find_index(ll, LL2):
    for i in range(len(LL2)):
        mm = matrix([ll, LL2[i]])
        if mm.rank() < 2:
            return(i)
    return("something went wrong")

## given two lists of points L1 and L2 (the points of L1 are
## contained in the points of L2), return the list L2 ordered
## in such a way that the first points are - in order - the
## points of L1 (and the remaining points are in arbitrary position).
def order_points(L1, L2):
    aux1 = [find_index(ll, L2) for ll in L1]
    aux2 = list(Set(range(len(L2))).difference(Set(aux1)))
    return([L2[i] for i in aux1]+[L2[i] for i in aux2])


def my_factor(f):
    if f == 0:
        return(0)
    else:
        return(f.factor())


## given 4 points Q1, Q2, Q3, Q4 such that (Q1+Q2) intersects (Q3+Q4)
## in a single point, determines the coordinates of this point.
def intersect_lines(Q1, Q2, Q3, Q4):
    twoLines = [det(matrix([[x, y, z], Q1, Q2])), \
                det(matrix([[x, y, z], Q3, Q4]))]
    Msist = matrix([[twoLines[0].coefficient(vr) for vr in [x, y, z]],\
                    [twoLines[1].coefficient(vr) for vr in [x, y, z]]])
    LL = Msist.minors(2)
    dd = gcd(LL)
    UU = list(map(lambda uu: R(uu/dd), LL))
    return(vector((UU[2], -UU[1], UU[0])))



# Given a square non singular matrix matA of order n and an column 
# matrix matB of n rows, the method gives a vector X of n components 
# which is the solution of the system
# matA*(d*X) = matB
# where d = det(matA)

# Example: Consider the system:
#     3*w1 +4*w2 -1*l1 -1*l2 = 0
#     1*w1 -2*w2 -3*l1 -5*l2 = 0

#i.e. the system to solve is:
#     3*w1 +4*w2 = 1*l1 +1*l2
#     1*w1 -2*w2 = 3*l1 +5*l2
#
# we define:
#  matA = matrix(S,[[3, 4], [1, -2]])
#  matB = matrix(S, [[l1+l2], [3*l1+5*l2]])

# Then solve_with_Cramer(matA, matB) gives:
# (-14*l1 - 22*l2, 8*l1 + 14*l2)
# which means that:
# det(matA)*(w1, w2) = (-14*l1 - 22*l2, 8*l1 + 14*l2)
# i.e.
# (-10)*(w1, w2) = (-14*l1 - 22*l2, 8*l1 + 14*l2)



def solve_with_Cramer(matA, matB): 
    n = matA.ncols()
    dt_matA = matA.det()
    sol = []
    for i in range(n):
        Ml = matA.matrix_from_columns(range(0,i))
        Mr = matA.matrix_from_columns(range(i+1, n))
        Maa = (Ml.augment(dt_matA*matB)).augment(Mr)
        sol.append(Maa.det()/dt_matA)
    return(vector(sol))

## here we compute the determinant of a 10 x 10 matrix
## with some tricks

def find_det(MM, withMessages = false):
    if withMessages:
        initialTime = cputime()
        print("We compute the determinant with some tricks")
        pausa = 1
    M2 = MM.matrix_from_rows(range(9))
    columnsM2 = range(10)
    minorsM2 = []
    for i in range(10):
        cpCol = list(columnsM2[:])
        cpCol.remove(i)
        if withMessages:
            ttA = cputime()
            print("    We start the computation number "+str(i+1)+"/10 ...")
            sleep(pausa)
        dMx = (M2.matrix_from_columns(cpCol)).det()
        minorsM2.append(dMx)
        if withMessages:
            print("    ... end of this partial computation. Time: "\
	    +str(cputime()-ttA))
            sleep(pausa)
    if withMessages:
        print("Now we compute a gcd of 10 polynomials")
        ttA = cputime()
        sleep(pausa)
    cmDiv = gcd(minorsM2)
    if withMessages:
        print(cmDiv.factor())
        print("gcd computed. Time: "+str(cputime()-ttA))
        print("Now we compute the 10 x 10 determinant:")
        ttA = cputime()
        sleep(pausa)
    minorsM2 = list(map(lambda uu: R(uu/cmDiv), minorsM2))
    DT = 0
    for i in range(10):
        DT = DT + (-1)^i*MM[9, i]*minorsM2[i]
    if withMessages:    
        print("Determinant computed. Time: "+str(cputime()-ttA))
        print("Total time of computations: "+str(cputime()-initialTime))
    return(DT)


    

## input: a 9 X 10 matrix M2
## output: a list (m1, m2, ..., m9) 
## where m_i are the max minors of M2. More precisely,
## m1 is the minor obtained from M2 after erasing the first
## column, m2 is the minor obtained from M2 after erasing the
## second column, and so on.

def find_9_minors(M2, withMessages = false, withCommonFactors = false):
    if withMessages:
        initialTime = cputime()
        print("We compute the determinant with some tricks")
        pausa = 1
    columnsM2 = range(10)
    minorsM2 = []
    for i in range(10):
        cpCol = list(columnsM2[:])
        cpCol.remove(i)  ##### parte rimuovendo la prima colonna e va avanti
        if withMessages:
            ttA = cputime()
            print("    We start the computation number "+str(i+1)+"/10 ...")
            sleep(pausa)
        dMx = (M2.matrix_from_columns(cpCol)).det()
        minorsM2.append(dMx)
        if withMessages:
            print("    ... end of this partial computation. Time: "\
	    +str(cputime()-ttA))
            sleep(pausa)
    #### This can be useful if we need to consider all possible cases"	    
    if withCommonFactors:
        print("There can be some common factors among the output of")
        print("the procedure find_9_minors()")
        return(minorsM2)
    #### 
    if withMessages:
        print("Now we compute a gcd of 10 polynomials")
        ttA = cputime()
        sleep(pausa)
    cmDiv = gcd(minorsM2)
    if withMessages:
        print("gcd computed. Time: "+str(cputime()-ttA))
    minorsM2 = list(map(lambda uu: S(uu/cmDiv), minorsM2))
    return(minorsM2)








### calcolo del determinante furbo.
def detFurbo(M1s):
    R = M1s[0,0].parent() ## cerco anello dove vive tutto
    initialTime = cputime()
    ## parto da M1s, matrice 10x10. Voglio il suo determinante furbo
    pausa = 1
    M2 = M1s.matrix_from_rows(range(9))
    columnsM2 = range(10)
    minorsM2 = []
    for i in columnsM2: ##range(10):
        cpCol = list(range(10))
        cpCol.remove(i)
        ttA = cputime()
        print("Inizio conto "+str(i)+"...")
        sleep(pausa)
        dMx = (M2.matrix_from_columns(cpCol)).det()
        minorsM2.append(dMx)
        print("... fine. Tempo: "+str(cputime()-ttA))
        sleep(pausa)
    print("Inizio calcolo gcd")
    ttA = cputime()
    sleep(pausa)
    cmDiv = gcd(minorsM2)
    print("Calcolato gcd. Tempo: "+str(cputime()-ttA))
    print("Inizio calcolo determinante 10 x 10:")
    ttA = cputime()
    sleep(pausa)
    minorsM2 = list(map(lambda uu: R(uu/cmDiv), minorsM2))
    DT = 0
    for i in range(10):
        DT = DT + (-1)^i*M1s[9, i]*minorsM2[i]
    print("Calcolato. Tempo: "+str(cputime()-ttA))
    print("Tempo totale: "+str(cputime()-initialTime))
    return(DT)


## computes phi_p(Pt), according to the definition of \Phi given in the
## paper (i.e. given a point Pt, if we require that it is an eigenpoint
## of the cubic a0*x^3+ ... +a9*z^3, we get three linear equations in
## a0, ..., a9. This method gives the coefficients of these lin. eq.

def phi_p(Pt):
    F = a0*x^3 + a1*x^2*y + a2*x*y^2 + a3*y^3 + a4*x^2*z \
    + a5*x*y*z + a6*y^2*z + a7*x*z^2 + a8*y*z^2 + a9*z^3
    ## the ideal of the minors:
    Maux = matrix([[x, y, z], [diff(F, x), diff(F, y), diff(F, z)]])
    Jaux = R.ideal(Maux.minors(2))
    Laux = Jaux.subs(substitution(Pt)).gens()
    vv = [a0,a1,a2,a3,a4,a5,a6,a7,a8,a9]
    terna = []
    for aa in Laux:
        terna.append(vector([aa.coefficient(cf) for cf in vv]))
    return(terna)

## transforme (if possible) a projective point in an affine point
def affinize(pp):
    if pp[2] == 0:
        return("forget it!")
    return([pp[0]/pp[2], pp[1]/pp[2]])