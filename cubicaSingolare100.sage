#
# Qui si vuol vedere quali sono le condizioni affinche' nei 5 punti a V
# il punto P1 sia singolare. 
# Qui si trova il caso P1 = (1, 0, 0)

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
## If we assume that P1 is singular, then we have the further condition:
## (1, 0, 0, 0, 0, 0, 0, 0, 0, 0) to add to M1, 
## hence we can extract from M1 the matrix without the first three rows and
## without the columns 0, 1, 4
MM = M1.matrix_from_rows_and_columns(range(3, 15), [2, 3, 5, 6, 7, 8, 9])


## Since we have:
assert(tuple(P2[2]*MM[0]-P2[1]*MM[1]+P2[0]*MM[2]) == (0, 0, 0, 0, 0, 0, 0))
assert(tuple(P3[2]*MM[3]-P3[1]*MM[4]+P3[0]*MM[5]) == (0, 0, 0, 0, 0, 0, 0))
assert(tuple(P4[2]*MM[6]-P4[1]*MM[7]+P4[0]*MM[8]) == (0, 0, 0, 0, 0, 0, 0))
assert(tuple(P5[2]*MM[9]-P5[1]*MM[10]+P5[0]*MM[11]) == (0, 0, 0, 0, 0, 0, 0))
## in the computation of the order 10 minors of M1 and hence of order 7 minors
## of MM, we can erase many matrices:

rg = Combinations(12, 7)

## given a list of (seven) rows st, the method checks if it
## contains the triplet 0, 1, 2 or ...
def is_min_sure_zero(st):
    return(Set([0, 1, 2]).issubset(Set(st)) or\
           Set([3, 4, 5]).issubset(Set(st)) or\
           Set([6, 7, 8]).issubset(Set(st)) or\
           Set([9, 10, 11]).issubset(Set(st)))

## select the "good" rows

rg1 = filter(lambda uu: not is_min_sure_zero(uu), rg)

print("First 'long' computation: computation of minors of order 7:")
sleep(1)
ttA = cputime()
min7 = [MM.matrix_from_rows(rr).det() for rr in rg1]

print(cputime()-ttA)
print()

def qr_gener(pol, divisor):
    pol1 = pol
    qr = pol1.quo_rem(divisor)
    while qr[1] == 0:
        pol1 = qr[0]
        qr = pol1.quo_rem(divisor)
    return(pol1)

## division by u1 (1 minute)
dt = matrix([P1, P2, P4]).det()

print("some divisions. Can take some time.")
sleep(1)
ttA = cputime()

min7 = filter(lambda uu: uu != 0, min7)

min7 = [qr_gener(mm, u1) for mm in min7]
min7 = [qr_gener(mm, u2) for mm in min7]
min7 = [qr_gener(mm, v1) for mm in min7]
min7 = [qr_gener(mm, v2) for mm in min7]
min7 = [qr_gener(mm, dt) for mm in min7]
print(cputime()-ttA)

print()

print("Calcolo base gr. Lunga")

sleep(1)
J7 = S.ideal(min7)
ttA = cputime()
gJ7 = J7.groebner_basis()
print(cputime()-ttA)

print("Calcolata base groebner")

## division by u1 (1 minute)
ttA = cputime()

gJ7 = [qr_gener(mm, u1) for mm in gJ7]
gJ7 = [qr_gener(mm, u2) for mm in gJ7]
gJ7 = [qr_gener(mm, v1) for mm in gJ7]
gJ7 = [qr_gener(mm, v2) for mm in gJ7]
gJ7 = [qr_gener(mm, dt) for mm in gJ7]
print("Calcolate saturazioni parziali")
print(cputime()-ttA)
sleep(1)
print()
sgJ7 = S.ideal(gJ7).saturation(u1*u2*v1*v2)[0]
sgJ7 = S.ideal(gJ7).saturation(dt)[0]
print("Calcolo dec primaria")
sleep(1)

PD = sgJ7.radical().primary_decomposition()
print("Calcolata dec. primaria")

print(len(PD) == 3)
print(PD[0] == S.ideal(delta1(P1, P2, P4), delta1b(P1, P4, P5)))
print(PD[1] == S.ideal(delta1(P1, P2, P4), delta1b(P1, P2, P3)))
print(PD[2] == S.ideal(delta1b(P1, P2, P3), delta1b(P1, P4, P5)))