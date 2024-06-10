##
## Proof of prop. 4.7:
## if P1, P2, P4 are three points, then 
## 5 <= rank(Phi(P1, P2, P4)) <= 6
## and, if rank(Phi(P1, P2, P4)) = 5, then P1, P2, P4 are aligned.
##

print("We define a ring with sufficiently many variables.")

varAn1 = ["s"+str(i)+str(j) for i in range(1,8) for j in range(i, 8)]
varAn2 = ["x", "y", "z", "u1", "u2", "v1", "v2", "w1", "w2", "l1", "l2"]
varAn3 = ["a"+str(i) for i in range(10)]
varAn4 = [str(XX)+str(i) for i in range(1, 8) for XX in ["A", "B", "C"]]
varAn5 = ["l", "m", "n"]
var_ii = ["ii"]
K = QQ
S = PolynomialRing(K, varAn1+varAn2+varAn3+varAn4+varAn5+var_ii)
S.inject_variables(verbose=false)

P1 = vector((A1, B1, C1))

load("auxiliaryProcedures.sage")

## case 1: 
##                  P4 = (1, 0, 0)


P1, P2 = vector((A1, B1, C1)), vector((A2, B2, C2))
P4 = vector((1, 0, 0))

load("../AlignedEigenpoints/software/ancillary.sage")
test = SymbolicCheck()
M = test.condition_matrix([P1, P2, P4], R, standard="all")

J5 = R.ideal(M.minors(5))
J6 = R.ideal(M.minors(6))

J5 = J5.saturation(R.ideal(matrix([P1, P2]).minors(2)))[0].saturation(R.ideal(matrix([P1, P4]).minors(2)))[0].saturation(R.ideal(matrix([P2, P4]).minors(2)))[0].radical()

# J5 is the ideal (1), so the matrix cannot have rank < 5:
assert(J5 == R.ideal(R.one()))

# When J6 is satisfied, it means that M has rank 5. 

J6 = J6.saturation(R.ideal(matrix([P1, P2]).minors(2)))[0].saturation(R.ideal(matrix([P1, P4]).minors(2)))[0].saturation(R.ideal(matrix([P2, P4]).minors(2)))[0]

pd = J6.radical().primary_decomposition()

## we have two possibilities: 
## either P1, P2, P4 are aligned and the line is tangent to the isotropic
## conic in P2 (hence P2 orthogonal to P4, P2 orthogonal to P2 and 
## P1, P2, P4 aligned):

PD0 = R.ideal(scalarProd(P2, P4), scalarProd(P2, P2), det(matrix([P1, P2, P4]))).saturation(R.ideal(list(P2)))[0].radical().primary_decomposition()

assert(len(PD0)==1)

assert(pd[0] == PD0[0])

## or
## P1, P2, P4 are aligned and the line is tangent to the isotropic
## conic in P1 (hence P1 orthogonal to P4, P1 orthogonal to P1 and 
## P1, P2, P4 aligned):

PD1 = R.ideal(scalarProd(P1, P4), scalarProd(P1, P1), det(matrix([P1, P2, P4]))).saturation(R.ideal(list(P1)))[0].radical().primary_decomposition()

assert(len(PD1) == 1)
assert(pd[1] == PD1[0])


## case 2:
##              P4 = (1, i, 0)


P1, P2 = vector((A1, B1, C1)), vector((A2, B2, C2))
P4 = vector((1, ii, 0))

M = test.condition_matrix([P1, P2, P4], R, standard="all")

J5 = R.ideal(M.minors(5)) + R.ideal(ii^2+1)
J6 = R.ideal(M.minors(6)) + R.ideal(ii^2+1)

J5 = J5.saturation(R.ideal(matrix([P1, P2]).minors(2)))[0].saturation(R.ideal(matrix([P1, P4]).minors(2)))[0].saturation(R.ideal(matrix([P2, P4]).minors(2)))[0].radical()


# J5 is the ideal (1), so the matrix cannot have rank < 5:
assert(J5 == R.ideal(R.one()))

# When J6 is satisfied, it means that M has rank 5. 

print("a computation of 6 seconds:")
sleep(1)

J6 = J6.saturation(R.ideal(matrix([P1, P2]).minors(2)))[0].saturation(R.ideal(matrix([P1, P4]).minors(2)))[0].saturation(R.ideal(matrix([P2, P4]).minors(2)))[0]


## when J6 is satisfied, we have that P1, P2, P4 are aligned
## and the line P1+P2+P4 is tangent to the isotropic conic in P4:

PD = R.ideal(ii^2+1, scalarProd(P2, P4), scalarProd(P1, P4), matrix([P1, P2, P4]).det()).radical().primary_decomposition()

J6 = J6.radical()

assert(len(PD) == 1)


assert(J6 == PD[0])
