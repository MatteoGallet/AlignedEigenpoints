## verifica che se vale delta1(P1, P2, P4)=0,  (P1, P2)(P1,P3+(P1,P1)*(P2,P3) = 0
## e analogo per P1, P4, P5, allora tutti i minoridi ordine 9 sono zero
## nel caso generale.


varAn1 = ["s"+str(i)+str(j) for i in range(1,8) for j in range(i, 8)]
varAn2 = ["x", "y", "z", "u1", "u2", "v1", "v2", "w1", "w2", "l1", "l2"]
varAn3 = ["a"+str(i) for i in range(10)]
varAn4 = [str(XX)+str(i) for i in range(1, 8) for XX in ["A", "B", "C"]]
var("xx")
K = QQ
#### K.<ii> = NumberField(xx^2+1)
S = PolynomialRing(K, varAn1+varAn2+varAn3+varAn4)
S.inject_variables(verbose=false)

P1 = vector(S, (A1, B1, C1))

load("auxiliaryProcedures.sage")
load("../AlignedEigenpoints/software/ancillary.sage")

test = SymbolicCheck()

P1 = vector((1, 0, 0))
P1 = vector((0, 1, 0))
P1 = vector((0, 0, 1))

P1 = vector((A1, B1, C1))
P2 = vector(S, (A2, B2, C2))
P4 = vector(S, (A4, B4, C4))
P3 = u1*P1+u2*P2
P5 = v1*P1+v2*P4

def qr_gener(pol, divisor):
    pol1 = pol
    qr = pol1.quo_rem(divisor)
    while qr[1] == 0:
        pol1 = qr[0]
        qr = pol1.quo_rem(divisor)
    return(pol1)

J1 = S.ideal(scalarProd(P1, P2)*scalarProd(P1, P3)+scalarProd(P1, P1)*scalarProd(P2, P3), 
scalarProd(P1, P4)*scalarProd(P1, P5)+scalarProd(P1, P1)*scalarProd(P4, P5),
delta1(P1, P2, P4)).saturation(matrix([P1, P2, P4]).det())[0]

J1 = J1.saturation(u1*u2*v1*v2)[0]

M = test.condition_matrix([P1, P2, P3, P4, P5], S)

#matrix([phi_p(P1)[0], phi_p(P1)[1],phi_p(P2)[0], phi_p(P2)[1],phi_p(P3)[0], phi_p(P3)[1],\
#           phi_p(P4)[0], phi_p(P4)[1],phi_p(P5)[0], phi_p(P5)[1]])


##rg = list(Combinations(10, 9))
##cl = list(Combinations(10, 9))

print("Long computation (about 53'):")
sleep(1)

ttA = cputime()
min9 = M.minors(9)
##min9 = [(M.matrix_from_rows_and_columns(rrg, ccl)).det() for rrg in rg for ccl in cl]
print(cputime()-ttA)

print()

print("Second long computation (about 30''): ")
sleep(1)
ttA = cputime()
dt = matrix([P1, P2, P4]).det()
min9d = [qr_gener(mm, dt) for mm in min9]
print(cputime()-ttA)

print()

print("Third long computation (about 30''): ")
sleep(1)
ttA = cputime()
rid = [J1.reduce(mm) for mm in min9d]

print(cputime()-ttA)


print(Set(rid) == Set([0]))