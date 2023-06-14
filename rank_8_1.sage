varAn1 = ["s"+str(i)+str(j) for i in range(1,8) for j in range(i, 8)]
varAn2 = ["x", "y", "z", "u1", "u2", "v1", "v2", "w1", "w2", "l1", "l2"]
varAn3 = ["a"+str(i) for i in range(10)]
varAn4 = [str(XX)+str(i) for i in range(1, 8) for XX in ["A", "B", "C"]]
var("xx")
K = QQ
K.<ii> = NumberField(xx^2+1)
S = PolynomialRing(K, varAn1+varAn2+varAn3+varAn4)
S.inject_variables(verbose=false)

P1 = vector(S, (A1, B1, C1))

load("auxiliaryProcedures.sage")

P1 = vector((1, 0, 0))
P2 = vector(S, (A2, B2, C2))
P4 = vector(S, (A4, B4, C4))
P3 = u1*P1+u2*P2
P5 = v1*P1+v2*P4

load("../AlignedEigenpoints/software/ancillary.sage")

test = SymbolicCheck()

# Construction of the matrix of all the linear conditions.

M1 = test.condition_matrix([P1, P2, P3, P4, P5], S, standard="all")

## Since the first three rows of M1 are, respectively,
## (0, 1, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 1, 0, 0, 0, 0, 0), and
## (0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
## in order to compute the minors of order 9 of M1, we can compute
## the minors of order 7 of the matrix obtained from the rows 3, 4, ..., 14
## of M1 and all the columns of M1 except columns 1 and 4.
MM = M1.matrix_from_rows_and_columns(range(3, 15), [0, 2, 3, 5, 6, 7, 8, 9])

## Since we have:
assert(tuple(P2[2]*MM[0]-P2[1]*MM[1]+P2[0]*MM[2]) == (0, 0, 0, 0, 0, 0, 0, 0))
assert(tuple(P3[2]*MM[3]-P3[1]*MM[4]+P3[0]*MM[5]) == (0, 0, 0, 0, 0, 0, 0, 0))
assert(tuple(P4[2]*MM[6]-P4[1]*MM[7]+P4[0]*MM[8]) == (0, 0, 0, 0, 0, 0, 0, 0))
assert(tuple(P5[2]*MM[9]-P5[1]*MM[10]+P5[0]*MM[11]) == (0, 0, 0, 0, 0, 0, 0, 0))

## a square submatrix of order 7 of MM has surely determinant zero if
## it contains the three rows 0, 1, 2 or the three rows 3, 4, 5 or the
## three rows 6, 7, 8 or the three rows 9, 10, 11. Hence we construct
## all the submatrix of MM or order seven which do not contain these
## these triplets of rows.

rg = Combinations(12, 7).list()

## given a list of (seven) rows st, the method checks if it
## contains the triplet 0, 1, 2 or ...
def is_min_sure_zero(st):
    return(Set([0, 1, 2]).issubset(Set(st)) or\
           Set([3, 4, 5]).issubset(Set(st)) or\
           Set([6, 7, 8]).issubset(Set(st)) or\
           Set([9, 10, 11]).issubset(Set(st)))

## select the "good" rows
rg1 = list(filter(lambda uu: not is_min_sure_zero(uu), rg))

## select the "good" columns
cl1 = Combinations(8, 7).list()

print("First 'long' computation: computation of minors of order 7:")
ttA = cputime()
min7 = []
for rr in rg1:
    for cc in cl1:
        min7.append(MM.matrix_from_rows_and_columns(rr, cc).det())
print(cputime()-ttA)
print()


sleep(1)

min7 = list(filter(lambda uu: uu != 0, min7))

print("Possible 'long' computation: computation of squarefree polynomials")
ttA = cputime()
min7s = [test.get_sqrfree(mm) for mm in min7]
print(cputime()-ttA)
print()


sleep(1)

print("Possible 'long' computation: saturation w.r.t. u and v")
ttA = cputime()
min7ss = [test.clear_uv(mm) for mm in min7s]
print(cputime()-ttA)
print()


sleep(1)

print("")
print("                      Case 1: A2 = 0")
print("Possible 'long' computation: saturation of the ideal w.r.t. u and v")

ttA = cputime()
J = (S.ideal(min7ss).subs(A2=0))
J1 = J.saturation(u1*v1*u2*v2)[0]
print(cputime()-ttA)
print()

J1 = J1.saturation((matrix([P1, P2, P4]).subs(A2=0)).det())[0]

PD1 = J1.radical().primary_decomposition()


print("Computed PD1, primary dec of the ideal when A2 = 0")

print("")

print("                      Case 2: A2 != 0")
print("")

## given a polynomial, the method computes the saturation w.r.t. A2
def clear_A2(poly):
    return S.ideal(poly).saturation(A2)[0].gens()[0]

print("Possible 'long' computation: saturation w.r.t. A2")
ttA = cputime()
mn7 = [clear_A2(mm) for mm in min7ss]
print(cputime()-ttA)
sleep(1)
print("")

print("Possible 'long' computation: saturation of the ideal w.r.t. u, v, A2")
ttA = cputime()
JJ = S.ideal(mn7)
J2 = JJ.saturation(u1*v1*u2*v2*A2)[0]
cputime()-ttA

sleep(1)

print("Possible 'long' computation: saturation of the ideal w.r.t. P1, P2, P4 aligned")
ttA = cputime()
J2 = J2.saturation(matrix([P1, P2, P4]).det())[0]
cputime()-ttA

sleep(1)

print("Possible 'long' computation: radical and pr. dec.")
ttA = cputime()
PD2 = J2.radical().primary_decomposition()
cputime()-ttA

sleep(1)

print("Computed PD2, primary dec of the ideal when A2 != 0")

