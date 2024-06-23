## Here we give a direct proof that if P1, P2 are arbitrary, if P4 is 
## such that delta1(P1, P2, P4) = 0 P3 is chosen in such a way that 
## delta1b(P1, P2, P3) = 0 and P5 is chosen in such a way that 
## delta1b(P1, P4, P5) = 0, then the matrix Phi(P1, ..., P5) 
## has rank 8. 

varAn1 = ["s"+str(i)+str(j) for i in range(1,8) for j in range(i, 8)]
varAn2 = ["x", "y", "z", "u1", "u2", "v1", "v2", "w1", "w2", "l1", "l2"]
varAn3 = ["a"+str(i) for i in range(10)]
varAn4 = [str(XX)+str(i) for i in range(1, 8) for XX in ["A", "B", "C"]]
var("xx")
K = QQ
#K.<ii> = NumberField(xx^2+1)
S = PolynomialRing(K, varAn1+varAn2+varAn3+varAn4)
S.inject_variables(verbose=false)

P1 = vector(S, (A1, B1, C1))

load("auxiliaryProcedures.sage")

P1 = vector(S, (1, 0, 0))  ## possible, up to rotations
P2 = vector(S, (A2, B2, C2))

## we want delta1(P1, P2, P4) = 0.
## Hence: (P4| s11*P1-s12*P2) = 0, hence:

Ptmp4 = vector(S, (A4, B4, C4))
Q4 = scalarProd(P1,P1)*P2-scalarProd(P1,P2)*P1
aa4, bb4, cc4 = scalarProd(Ptmp4, Q4).coefficient(A4), \
                scalarProd(Ptmp4, Q4).coefficient(B4), \
                scalarProd(Ptmp4, Q4).coefficient(C4)

## two alternative definitions of P4 (solving delta1()=0 w.r.t. 
## A4 or C4)

P4 = vector(S, (-bb4*B4-cc4*C4, aa4*B4, aa4*C4)) ## this gives P5 = (0, 0, 0)
P4 = vector(S, (cc4*A4, cc4*B4, -aa4*A4-bb4*B4))
P4 = vector(S, (bb4*A4, -aa4*A4-cc4*C4, bb4*C4))

assert(delta1(P1, P2, P4) == 0)

P3 = (scalarProd(P1, P2)^2+scalarProd(P1, P1)*scalarProd(P2, P2))*P1-2*scalarProd(P1, P1)*scalarProd(P1, P2)*P2
assert(delta1b(P1, P2, P3) == 0)

P5 = (scalarProd(P1, P4)^2+scalarProd(P1, P1)*scalarProd(P4, P4))*P1-2*scalarProd(P1, P1)*scalarProd(P1, P4)*P4
assert(delta1b(P1, P4, P5) == 0)

M = matrixEigenpoints([P1, P2, P3, P4, P5])

M1 = M.matrix_from_rows([0, 1, 3, 4, 6, 7, 9, 11, 12, 13])
assert(M1.rank() == 8)

m9 = M1.minors(9)
assert(Set(m9) == Set([0]))