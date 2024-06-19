## 
## We prove that, if P1, ..., P5 is a V-configuration and P1 is on 
## the isotropic conic, then Phi(P1, ..., P5) cannot have rakn <9.
## (this result is used in theorem:rank_V)

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

P1 = vector((1, ii, 0))
P2 = vector(S, (A2, B2, C2))
P4 = vector(S, (A4, B4, C4))
P3 = u1*P1+u2*P2
P5 = v1*P1+v2*P4



load("../AlignedEigenpoints/software/ancillary.sage")

test = SymbolicCheck()

# Construction of the matrix of all the linear conditions.

M1 = test.condition_matrix([P1, P2, P3, P4, P5], S, standard="all")

## Manipulation of M1 with elementary rows and columns operations,
## in order to have a simpler matrix.
## We use the fact that M1 has the following three rows:
## ((-3*ii), 3, (3*ii), -3, 0, 0, 0, 0, 0, 0)
## (0, 0, 0, 0, 1, ii, -1, 0, 0, 0)
## (0, 0, 0, 0, ii, -1, (-ii), 0, 0, 0)
## The first and the second are linearly independent.


M1.rescale_row(0, 1/3)

for j in range(3, 15):
    M1.add_multiple_of_row(j, 0, -M1[j, 1])

assert([M1[j, 1] for j in range(3, 15)] == [0 for j in range(3, 15)])

for j in range(3, 15):
    M1.add_multiple_of_row(j, 1, -M1[j, 4])

assert([M1[j, 4] for j in range(3, 15)] == [0 for j in range(3, 15)])

## In order to compute the minors of order 9 of M1, now, we can compute
## the minors of order 7 of the matrix obtained from the rows 3, 4, ..., 14
## of M1 and all of its columns, except columns 0 and 4.


MM = M1.matrix_from_rows_and_columns(range(3, 15), [0, 2, 3, 5, 6, 7, 8, 9])

## Since we have:
assert(tuple(P2[2]*MM[0]-P2[1]*MM[1]+P2[0]*MM[2]) == (0, 0, 0, 0, 0, 0, 0, 0))
assert(tuple(P3[2]*MM[3]-P3[1]*MM[4]+P3[0]*MM[5]) == (0, 0, 0, 0, 0, 0, 0, 0))
assert(tuple(P4[2]*MM[6]-P4[1]*MM[7]+P4[0]*MM[8]) == (0, 0, 0, 0, 0, 0, 0, 0))
assert(tuple(P5[2]*MM[9]-P5[1]*MM[10]+P5[0]*MM[11]) == (0, 0, 0, 0, 0, 0, 0, 0))

## a square submatrix of order 7 of MM has surely determinant zero if
## it contains the three rows 0, 1, 2 or the three rows 3, 4, 5 or the
## three rows 6, 7, 8 or the three rows 9, 10, 11. Hence we construct
## all the submatrices of MM of order seven which do not contain these
## triplets of rows.

## given a list of (seven) rows st, the method checks if it
## contains the triplet 0, 1, 2 or ...
def is_min_sure_zero(st):
    return(Set([0, 1, 2]).issubset(Set(st)) or\
           Set([3, 4, 5]).issubset(Set(st)) or\
           Set([6, 7, 8]).issubset(Set(st)) or\
           Set([9, 10, 11]).issubset(Set(st)))

## select the "good" rows
rg = Combinations(12, 7).list()
rg1 = list(filter(lambda uu: not is_min_sure_zero(uu), rg))

## select the "good" columns
cl1 = Combinations(8, 7).list()

print("First 'long' computation: computation of minors of order 7")
print("About 6 min")
sleep(1)

ttA = cputime()
min7 = [MM.matrix_from_rows_and_columns(rr, cc).det() for rr in rg1 \
for cc in cl1]
print(cputime()-ttA)
print()

min7 = list(filter(lambda uu: uu != 0, min7))

## Some preprocessing. We divide each element of min7
## by u1, u2, v1, v2, and dt as much as possible (dt is the condition
## P1, P2, P4 aligned).
## In this way the factorization of the elements of min7 becomes simpler.


dt = matrix([P1, P2, P4]).det()

## Given a polynomial "pol" and another polynomial "divisor", this
## method divides pol by divisor as much as possible
## The method returns a couple. The first component is the divided
## polynomial, the second component is true/false and says if the given
## polynomial was divided at least one time.

def qr_gener(pol, divisor):
    pol1 = pol
    qr = pol1.quo_rem(divisor)
    while qr[1] == 0:
        pol1 = qr[0]
        qr = pol1.quo_rem(divisor)
    return(pol1)

## division by u1 (1 minute)

print("some divisions. Can take some time (~4 min).")
sleep(1)
ttA = cputime()

min7 = [qr_gener(mm, u1) for mm in min7]
min7 = [qr_gener(mm, u2) for mm in min7]
min7 = [qr_gener(mm, v1) for mm in min7]
min7 = [qr_gener(mm, v2) for mm in min7]
min7 = [qr_gener(mm, dt) for mm in min7]

print(cputime()-ttA)
print()

print("Possible 'long' computation: computation of squarefree polynomials")
print("About 55 min")
sleep(1)
ttA = cputime()
min7s = [test.get_sqrfree(mm) for mm in min7]
print(cputime()-ttA)
print()

print("Possible 'long' computation: saturation of ideal w.r.t. u and v")
print("About 60 min")
sleep(1)

Ja = S.ideal(min7s)

ttA = cputime()
JJs = Ja.saturation(u1*u2*v1*v2)[0]
print(cputime()-ttA)
print()

print("Possible 'long' computation: saturation of ideal w.r.t. det(mat())")
print("About 1 min 20 sec")
sleep(1)

ttA = cputime()
JJsa = JJs.saturation(dt)
print(cputime()-ttA)

print()
print("Final result:")
print("We do not have solutions: the ideal is (1)")
assert(JJsa[0] == R.ideal(R.one()))


