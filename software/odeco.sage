## Here we consider the case of ODECO.

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


Q1 = vector(S, (u1, 0, 0))
Q2 = vector(S, (0, v1, 0))
Q3 = vector(S, (0, 0, w1))
Q4 = Q1+Q2+Q3
Q5 = Q1+Q2
Q6 = Q1+Q3
Q7 = Q2+Q3

mon = [x^3, x^2*y, x*y^2, y^3, x^2*z, x*y*z, y^2*z, x*z^2, y*z^2, z^3]

ME = matrixEigenpoints([Q1, Q2, Q3, Q4, Q5, Q6, Q7])

ME1 = ME.matrix_from_rows([1, 5, 7, 8, 9, 10, 12, 15, 16])
assert(ME.rank() == 9)
cb = ME1.stack(vector(S, mon)).det()
cb1 = cb.factor()[-1][0]
assert(cb1 == z^3*u1*v1 + y^3*u1*w1 + x^3*v1*w1)


assert(allignments([Q1, Q2, Q3, Q4, Q5, Q6, Q7])==\
       [(1, 2, 5), (1, 3, 6), (1, 4, 7), (2, 3, 7), (2, 4, 6), (3, 4, 5)])

assert(scalarProd(wedgeProd(Q1, Q2), wedgeProd(Q3, Q5)) == S(0))
assert(scalarProd(wedgeProd(Q2, Q6), wedgeProd(Q1, Q3)) == S(0))
assert(scalarProd(wedgeProd(Q2, Q3), wedgeProd(Q1, Q7)) == S(0))

## Moreover,it holds: 
## PP4 == (P1||P2)(P1|P7)(P2|P7)-(P1|P2)(P1||P7)(P2|P7)+(P1|P2)(P1|P7)(P2||P7)

QQ3 = wedgeProd(Q1, Q2)*scalarProd(Q2, Q4)*scalarProd(Q1, Q4)\
    -scalarProd(Q1, Q2)*wedgeProd(Q2, Q4)*scalarProd(Q1, Q4)\
    +scalarProd(Q1, Q2)*scalarProd(Q2, Q4)*wedgeProd(Q1, Q4)

assert(QQ3 != vector(S, (0, 0, 0)))
assert(matrix([QQ3, Q3]).minors(2) == [0, 0, 0])

ort = []
PT = [Q1, Q2, Q3, Q4, Q5, Q6, Q7]
for i in range(6):
    for j in range(i, 7):
        if scalarProd(PT[i], PT[j]) == 0:
            ort.append((i, j))

assert(ort == [(0, 1), (0, 2), (0, 6), (1, 2), (1, 5), (2, 4)])