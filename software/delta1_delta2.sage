## conto per delta1 e delta2... 24 apr 23

## Qui si trovano tre polinomi di grado 3 in x, y, z
## in modo che i loro zeri comuni siano gli autopunti.
## Il conto e' fatto sia per delta1, sia per delta2.
## In entrambi i casi il conto e' ripetuto per 3 possibilita'
## (che nascono dal fatto di prendere la prima, la seconda e la
## terza componente di Phi(P5)).
## In particolare, se Ga, Gb, Gc sono i tre polinomi relativi al
## caso delta2 = 0 (che in fondo a questo file sono chiamati
## Gg1a, Gg1b, Gg1c) si trova subito che C1*Gg1a-B1*Gg1b+A1*Gg1c
## si fattorizza in tre fattori lineari in x, y, z che danno tre rette.
## La prima contiene i punti P1, P2, P3, la seconda i punti P1, P4, P5
## e la terza i punti P1, P6, P7.



print("We define a ring with sufficiently many variables.")


varAn1 = ["s"+str(i)+str(j) for i in range(1,8) for j in range(i, 8)]
varAn2 = ["x", "y", "z", "u1", "u2", "v1", "v2", "w1", "w2", "l1", "l2"]
varAn3 = ["a"+str(i) for i in range(10)]
varAn4 = [str(XX)+str(i) for i in range(1, 8) for XX in ["A", "B", "C"]]

K = QQ
S = PolynomialRing(K, varAn1+varAn2+varAn3+varAn4)
S.inject_variables(verbose=false)
P1 = vector(S, (A1, B1, C1))

P1, P2, P4  = vector((A1, B1, C1)), vector((A2, B2, C2)), vector((A4, B4, C4))
P3 = u1*P1 + u2*P2
P5 = v1*P1 + v2*P4
P6 = vector(S, (x, y, z))

load("auxiliaryProcedures.sage")

print("We define 'mon' and 'Jrel'")
mon = [x^3, x^2*y, x*y^2, y^3, x^2*z, x*y*z, y^2*z, x*z^2, y*z^2, z^3]
Jrel = S.ideal(s11-scalarProd(P1, P1), \
               s12-scalarProd(P1, P2), \
               s22-scalarProd(P2, P2), \
               s14-scalarProd(P1, P4), \
               s24-scalarProd(P2, P4), \
               s44-scalarProd(P4, P4), \
               s15-scalarProd(P1, P5), \
               s25-scalarProd(P2, P5), \
               s45-scalarProd(P4, P5), \
               s55-scalarProd(P5, P5))


print("""We assume (we know) that the vectors:
phi_p(P1)[0], phi_p(P1)[1], phi_p(P2)[0], phi_p(P2)[1], 
phi_p(P3)[0], phi_p(P3)[1], phi_p(P4)[0], phi_p(P4)[1]
are linearly independent and we contruct three  
9 X 10 matrices M1, M2, M3 with these eight rows plus, respectively, 
the row phi_p(P5)[0], phi_p(P5)[1], phi_p(P5)[2]""")

M1 = matrix([phi_p(P1)[0], phi_p(P1)[1], phi_p(P2)[0], phi_p(P2)[1], \
phi_p(P3)[0], phi_p(P3)[1], phi_p(P4)[0], phi_p(P4)[1], phi_p(P5)[0]])


M2 = matrix([phi_p(P1)[0], phi_p(P1)[1], phi_p(P2)[0], phi_p(P2)[1], \
phi_p(P3)[0], phi_p(P3)[1], phi_p(P4)[0], phi_p(P4)[1], phi_p(P5)[1]])


M3 = matrix([phi_p(P1)[0], phi_p(P1)[1], phi_p(P2)[0], phi_p(P2)[1], \
phi_p(P3)[0], phi_p(P3)[1], phi_p(P4)[0], phi_p(P4)[1], phi_p(P5)[2]])

print("\nM1, M2, M3 constructed\n")

print("""We construct 3 vectors minM1, minM2, minM3, each with 10 
components, which are the maximal minors of M1, M2, M3, respectively""")


doCompleteComputation = false


if doCompleteComputation: 
    sleep(1)
    ttA = cputime()
    minM1 = M1.minors(9)
    print(cputime()-ttA)
    minM1 = [numerator(mm/(A1*A2*A4*u1^2*u2^2*v1*v2*(u1*A1+u2*A2)*det(matrix([P1, P2, P4])^2))) for mm in minM1]
    minM1 = [minM1[9-i] for i in range(10)]

    sleep(1)
    ttA = cputime()
    minM2 = M2.minors(9)
    print(cputime()-ttA)
    minM2 = [numerator(mm/(A1*A2*A4*u1^2*u2^2*v1*v2*(u1*A1+u2*A2)*det(matrix([P1, P2, P4])^2))) for mm in minM2]
    minM2 = [minM2[9-i] for i in range(10)]

    sleep(1)
    ttA = cputime()
    minM3 = M3.minors(9)
    print(cputime()-ttA)
    minM3 = [numerator(mm/(A1*A2*A4*u1^2*u2^2*v1*v2*(u1*A1+u2*A2)*det(matrix([P1, P2, P4])^2))) for mm in minM3]
    minM3 = [minM3[9-i] for i in range(10)]

    threeMin = [minM1, minM2, minM3]
    save(threeMin, "longComput/threeMin.sobj")
else:
    threeMin = load("longComput/threeMin.sobj")
    minM1, minM2, minM3 = tuple(threeMin)

print("""Computation of the determinant of the matrices obtained 
from M1 adding the row phi_p(P6)[0], phi_p(P6)[1], phi_p(P6)[2].
These three determinants (which are polynomials of degree 3 in x, y, z)
are denoted by G1a, G1b, G1c. The same computation is repeated 
using the matrix M2 and M3 in place of M1. In this way the polynomials
G2a, G2b, G2c and G3a, G3b, G3c are obtained""")


row1, row2, row3 = tuple(phi_p(P6))
G1a = add([(-1)^i*row1[i]*minM1[i] for i in range(10)])
G1b = add([(-1)^i*row2[i]*minM1[i] for i in range(10)])
G1c = add([(-1)^i*row3[i]*minM1[i] for i in range(10)])

G2a = add([(-1)^i*row1[i]*minM2[i] for i in range(10)])
G2b = add([(-1)^i*row2[i]*minM2[i] for i in range(10)])
G2c = add([(-1)^i*row3[i]*minM2[i] for i in range(10)])

G3a = add([(-1)^i*row1[i]*minM3[i] for i in range(10)])
G3b = add([(-1)^i*row2[i]*minM3[i] for i in range(10)])
G3c = add([(-1)^i*row3[i]*minM3[i] for i in range(10)])

print("\nG1a, G1b, G1c, ..., G3c computed\n")

print("""It holds: 
z*G1a-y*G1b+x*G1c = 0 
z*G2a-y*G2b+x*G2c = 0 
z*G3a-y*G3b+x*G3c = 0 
(since we know that z*phi_p(P6)[0]-y*phi_p(P6)[1]+x*phi_p(P6)[2] = 0).
Analogously, 
P5[2]*G1a-P5[1]*G2a+P5[0]*G3a = 0
P5[2]*G1b-P5[1]*G2b+P5[0]*G3b = 0
P5[2]*G1c-P5[1]*G2c+P5[0]*G3c = 0.""")

print(z*G1a-y*G1b+x*G1c==0, z*G2a-y*G2b+x*G2c==0, z*G3a-y*G3b+x*G3c==0)

print(P5[2]*G1a-P5[1]*G2a+P5[0]*G3a==0, \
      P5[2]*G1b-P5[1]*G2b+P5[0]*G3b==0, \
      P5[2]*G1c-P5[1]*G2c+P5[0]*G3c==0)


print("""In order to improve the computations, we convert the 9 polynomials
G1a,..., G3c into 9 vectors, each of 10 components given by the coefficients
of these polynomials as polynomials in x, y, z. We construct therefore the 
lists listG1a, listG1b, ..., listG3c""")


listG1a = [G1a.coefficient(mn) for mn in mon]
listG1b = [G1b.coefficient(mn) for mn in mon]
listG1c = [G1c.coefficient(mn) for mn in mon]

listG2a = [G2a.coefficient(mn) for mn in mon]
listG2b = [G2b.coefficient(mn) for mn in mon]
listG2c = [G2c.coefficient(mn) for mn in mon]

listG3a = [G3a.coefficient(mn) for mn in mon]
listG3b = [G3b.coefficient(mn) for mn in mon]
listG3c = [G3c.coefficient(mn) for mn in mon]

print("\nConstructed the lists listG1a, listG1b, ..., listG3c\n")



### 		     END OF THE COMMON PART



print("""The lists listG1a, ..., listG3c are the 10 coefficient of nine 
polynomials in x, y, z of degree 3.
We want to see the relations between the first three polynomials defined 
in this way, the second three and the third three in the two cases delta1 = 0
and delta2 = 0.

Case delta1 = 0""")

print("""\nCondition for delta1 = 0:
sst= {A4: -D1.coefficient(B4)*B4-D1.coefficient(C4)*C4, 
      B4: D1.coefficient(A4)*B4, 
      C4: D1.coefficient(A4)*C4,
      v1: D1.coefficient(A4)*v1}""")

D1 = delta1(P1, P2, P4)
sst= {A4: -D1.coefficient(B4)*B4-D1.coefficient(C4)*C4, \
      B4: D1.coefficient(A4)*B4, \
      C4: D1.coefficient(A4)*C4,\
      v1: D1.coefficient(A4)*v1}   ##### aggiunta per non modificare P5
print("\nIs the substitution correct?\n")
print(D1.subs(sst) == 0)

print("""\nWe substitute sst into listG1a, ..., listG3c. 
Time of computation (if required): 6-7 minutes""")

####
if doCompleteComputation:
    ttA = cputime()
    print("start computation for G1...")
    pLG1a = [cf.subs(sst) for cf in listG1a]
    pLG1b = [cf.subs(sst) for cf in listG1b]
    pLG1c = [cf.subs(sst) for cf in listG1c]
    print("... end computation for G1")
    sleep(1)
    
    print("start computation for G2...")
    pLG2a = [cf.subs(sst) for cf in listG2a]
    pLG2b = [cf.subs(sst) for cf in listG2b]
    pLG2c = [cf.subs(sst) for cf in listG2c]
    print("... end computation for G2")
    sleep(1)
    
    print("start computation for G3...")
    pLG3a = [cf.subs(sst) for cf in listG3a]
    pLG3b = [cf.subs(sst) for cf in listG3b]
    pLG3c = [cf.subs(sst) for cf in listG3c]
    print("... end computation for G3")
    
    print(cputime()-ttA)
    
    threepLG1 = [pLG1a, pLG1b, pLG1c]
    threepLG2 = [pLG2a, pLG2b, pLG2c]
    threepLG3 = [pLG3a, pLG3b, pLG3c]
    save(threepLG1, "longComput/threepLG1.sobj")
    save(threepLG2, "longComput/threepLG2.sobj")
    save(threepLG3, "longComput/threepLG3.sobj")
else:
    threepLG1 = load("longComput/threepLG1.sobj")
    threepLG2 = load("longComput/threepLG2.sobj")
    threepLG3 = load("longComput/threepLG3.sobj")
    pLG1a, pLG1b, pLG1c = tuple(threepLG1)
    pLG2a, pLG2b, pLG2c = tuple(threepLG2)
    pLG3a, pLG3b, pLG3c = tuple(threepLG3)



print("""\nNow we have 9 lists: pLG1a, ..., pLG3c. Each list 
is the list of the 10 coefficient of a polynomial of degree 3 in x, y, z
in the parameters A1, A2, ... with the condition delta1 = 0.
We eliminate a common factor from the 10 elements of each list.
It turns out that pLG1a, pLG1b and pLG1c have the same common factor
(similarly the other triplets of lists). We call UG1a, ..., UG3c the 
nine lists and dd1, dd2, dd3 the common factors of the lists pLG*""")

dd1 = gcd(pLG1a)

UG1a = [S(pLG1a[i]/dd1) for i in range(10)]
UG1b = [S(pLG1b[i]/dd1) for i in range(10)]
UG1c = [S(pLG1c[i]/dd1) for i in range(10)]

dd2 = gcd(pLG2a)

UG2a = [S(pLG2a[i]/dd2) for i in range(10)]
UG2b = [S(pLG2b[i]/dd2) for i in range(10)]
UG2c = [S(pLG2c[i]/dd2) for i in range(10)]

dd3 = gcd(pLG3a)

UG3a = [S(pLG3a[i]/dd3) for i in range(10)]
UG3b = [S(pLG3b[i]/dd3) for i in range(10)]
UG3c = [S(pLG3c[i]/dd3) for i in range(10)]

print("""Now we reconstruct from UG1a, ... the corresponding polynomials
and we call them:
UGg1a, UGg1b, UGg1c;   
UGg2a, UGg2b, UGg2c;   
UGg3a, UGg3b, UGg3c:""")   
 

UGg1a = add([UG1a[i]*mon[i] for i in range(10)])
UGg1b = add([UG1b[i]*mon[i] for i in range(10)])
UGg1c = add([UG1c[i]*mon[i] for i in range(10)])

UGg2a = add([UG2a[i]*mon[i] for i in range(10)])
UGg2b = add([UG2b[i]*mon[i] for i in range(10)])
UGg2c = add([UG2c[i]*mon[i] for i in range(10)])

UGg3a = add([UG3a[i]*mon[i] for i in range(10)])
UGg3b = add([UG3b[i]*mon[i] for i in range(10)])
UGg3c = add([UG3c[i]*mon[i] for i in range(10)])

print("""UG1a, UG1b, ... reconstructed.
We verify that 
z*UGg1a-y*UGg1b+x*UGg1c = 0
z*UGg2a-y*UGg2b+x*UGg2c = 0
z*UGg3a-y*UGg3b+x*UGg3c = 0
""")

print(z*UGg1a-y*UGg1b+x*UGg1c == 0, z*UGg2a-y*UGg2b+x*UGg2c == 0, \
z*UGg3a-y*UGg3b+x*UGg3c == 0)

print("""\nThe three triplets of polynomials:
(UGg1a, UGg1b, UGg1c), (UGg2a, UGg2b, UGg2c), (UGg3a, UGg3b, UGg3c)
are the same (up to sign)""")
print((UGg1a, UGg1b, UGg1c) == (-UGg2a, -UGg2b, -UGg2c))
print((UGg1a, UGg1b, UGg1c) == (-UGg3a, -UGg3b, -UGg3c))

print("""\nCONCLUSION FOR THE CASE delta1 = 0:

If (x, y, z) are a common zero to (UGg1a, UGg1b, UGg1c), 
then it is also a common zero to (UGg2a, UGg2b, UGg2c) and to
(UGg3a, UGg3b, UGg3c). Therefore the 18x10 matrix whose rows are:
phi_p(P1), ..., phi_p(P6) has rank <=9

There are some coefficients (dd1, dd2, dd3) which have to be studied.""")

print(dd1.factor())
print("")
print(dd2.factor())
print("")
print(dd3.factor())
print("")



################ CASE delta2 = 0 ##################

print("\nCase delta2 = 0")



print("""The points P1, ..., P5 have to satisfy the condition
delta2(P1,..., P5) = 0, which is of the form U1*u1+U2*u2 = 0,
for suitable U1 and U2. We compute U1 and U2 and we make the
substitution {u1:U2, u2:-U1} into listG1a, ..., listG3b (we follow
this way, since the direct substitution of u1 and u2 into
G1a, ..., G3c seems too difficult).
We obtain 9 lists of 10 elements, which are denoted by
SLG1a, SLG1b, ..., SLG3c""")

U1 = delta2(P1, P2, P3, P4, P5).coefficient(u1)
U2 = delta2(P1, P2, P3, P4, P5).coefficient(u2)

print("""The polynomials U1 and U2 are respectively:
U1 = (P1|P2)*((P1|P1)*(P4,P5)-(P1,P4)*(P1,P5))
U2 = (P1|P2)^2*(P4|P5)-(P1|P4)*(P1|P5)*(P2|P2)
(we obtain this with the instructions: Jrel.reduce(U1) and 
Jrel.reduce(U2))""")
print(U1 == scalarProd(P1, P2)*(scalarProd(P1, P1)*scalarProd(P4,P5)-\
            scalarProd(P1, P4)*scalarProd(P1, P5)))
print(U2 == scalarProd(P1, P2)^2*scalarProd(P4, P5)-\
            scalarProd(P1, P4)*scalarProd(P1, P5)*scalarProd(P2, P2))


if doCompleteComputation:
    ttA = cputime()
    print("start computation for G1...")
    SLG1a = [cf.subs({u1: U2, u2: -U1}) for cf in listG1a]
    SLG1b = [cf.subs({u1: U2, u2: -U1}) for cf in listG1b]
    SLG1c = [cf.subs({u1: U2, u2: -U1}) for cf in listG1c]
    print("... end computation for G1")
    sleep(1)
    
    print("start computation for G2...")
    SLG2a = [cf.subs({u1: U2, u2: -U1}) for cf in listG2a]
    SLG2b = [cf.subs({u1: U2, u2: -U1}) for cf in listG2b]
    SLG2c = [cf.subs({u1: U2, u2: -U1}) for cf in listG2c]
    print("... end computation for G2")
    sleep(1)
    
    print("start computation for G3...")
    SLG3a = [cf.subs({u1: U2, u2: -U1}) for cf in listG3a]
    SLG3b = [cf.subs({u1: U2, u2: -U1}) for cf in listG3b]
    SLG3c = [cf.subs({u1: U2, u2: -U1}) for cf in listG3c]
    print("... end computation for G3")
    
    print(cputime()-ttA)
    
    threeSLG1 = [SLG1a, SLG1b, SLG1c]
    threeSLG2 = [SLG2a, SLG2b, SLG2c]
    threeSLG3 = [SLG3a, SLG3b, SLG3c]
    save(threeSLG1, "longComput/threeSLG1.sobj")
    save(threeSLG2, "longComput/threeSLG2.sobj")
    save(threeSLG3, "longComput/threeSLG3.sobj")
else:
    threeSLG1 = load("longComput/threeSLG1.sobj")
    threeSLG2 = load("longComput/threeSLG2.sobj")
    threeSLG3 = load("longComput/threeSLG3.sobj")
    SLG1a, SLG1b, SLG1c = tuple(threeSLG1)
    SLG2a, SLG2b, SLG2c = tuple(threeSLG2)
    SLG3a, SLG3b, SLG3c = tuple(threeSLG3)

print("\nThe 9 lists SLG1a, ...SLG3c are now computed or loaded\n")

print("""The meaning of SLG1a, ... is the following:
SLG1a, SLG1b, SLG1c, SLG2a, ... are vectors (lists) of 10 components.
Each of them represents a polynomial of degree 3 in x, y, z (which is
obtained computing the scalar product of SLG1a, SLG1b, ... with the 
10-components vector 'mon'. These 9 polynomials are the determinant
of the matrix whose rows are: 
phi_p(P1)[0], phi_p(P1)[1], phi_p(P2)[0], phi_p(P2)[1], 
phi_p(P3)[0], phi_p(P3)[1], phi_p(P4)[0], phi_p(P4)[1],
phi_p(P5)[i], phi_p(P6)[j]
for i = 0, 1, 2 and j = 0, 1, 2, in which the substitution u1 = U2, 
u2 = -U1 is done.
For i = 0, 1, 2 and j = 0 we obtain the polynomials coming from SLG1a, b, c
For i = 0, 1, 2 and j = 1 we obtain the polynomials coming from SLG2a, b, c
For i = 0, 1, 2 and j = 2 we obtain the polynomials coming from SLG3a, b, c""")

print("""Now we reconstruct from SLG1a, ... the corresponding polynomials
and we call them:
Gg1a, Gg1b, Gg1c;   
Gg2a, Gg2b, Gg2c;   
Gg3a, Gg3b, Gg3c:""")   
 

Gg1a = add([SLG1a[i]*mon[i] for i in range(10)])
Gg1b = add([SLG1b[i]*mon[i] for i in range(10)])
Gg1c = add([SLG1c[i]*mon[i] for i in range(10)])

Gg2a = add([SLG2a[i]*mon[i] for i in range(10)])
Gg2b = add([SLG2b[i]*mon[i] for i in range(10)])
Gg2c = add([SLG2c[i]*mon[i] for i in range(10)])

Gg3a = add([SLG3a[i]*mon[i] for i in range(10)])
Gg3b = add([SLG3b[i]*mon[i] for i in range(10)])
Gg3c = add([SLG3c[i]*mon[i] for i in range(10)])

print("""Gg1a, Gg1b, ... reconstructed
We verify that 
z*Gg1a-y*Gg1b+x*Gg1c = 0
z*Gg2a-y*Gg2b+x*Gg2c = 0
z*Gg3a-y*Gg3b+x*Gg3c = 0
""")

print(z*Gg1a-y*Gg1b+x*Gg1c == 0, z*Gg2a-y*Gg2b+x*Gg2c == 0, \
z*Gg3a-y*Gg3b+x*Gg3c == 0)

print("""\nWe look for a relation between the triplets (vectors) 
(Gg1a, Gg1b, Gg1c), (Gg2a, Gg2b, Gg2c) and (Gg3a, Gg3b, Gg3c)
""")

lc1a, lc2a, lc3a = Gg1a.coefficient(x^3), Gg2a.coefficient(x^3), \
                   Gg3a.coefficient(x^3)
d12 = gcd(lc1a, lc2a)
d13 = gcd(lc1a, lc3a)
d23 = gcd(lc2a, lc3a)

alpha = -S(lc3a/d23) ##  = -S(lc3a/d13)
beta = S(lc2a/d23) ## = S(lc2a/d12)
gamma = -S(lc1a/d13) ## = -S(lc1a/d12)
alphaR = Jrel.reduce(alpha)
betaR = Jrel.reduce(beta)
gammaR = Jrel.reduce(gamma)
print("""After defining some coefficients, we have:
beta*vector((Gg1a, Gg1b, Gg1c))+gamma*vector((Gg2a, Gg2b, Gg2c)) = 0
alpha*vector((Gg1a, Gg1b, Gg1c))-gamma*vector((Gg3a, Gg3b, Gg3c)) = 0
alpha*vector((Gg2a, Gg2b, Gg2c))+beta*vector((Gg3a, Gg3b, Gg3c)) = 0
where alpha, beta, gamma, expressed in terms of the scalar product 
of the points, are:
""")
print(alphaR)
print(betaR)
print(gammaR)
print("\nWe verify the above relations:\n")
print(beta*vector((Gg1a, Gg1b, Gg1c))+gamma*vector((Gg2a, Gg2b, Gg2c)) == 0)
print(alpha*vector((Gg1a, Gg1b, Gg1c))-gamma*vector((Gg3a, Gg3b, Gg3c)) == 0)
print(alpha*vector((Gg2a, Gg2b, Gg2c))+beta*vector((Gg3a, Gg3b, Gg3c)) == 0)

print("""\nWe define the following vector:
WW = vector((s14*s25+s12*s45, -s14*s15-s11*s45, -s12*s15+s11*s25))
and we verify that:

vector((alphaR, betaR, gammaR)) == WW * matrix([P1, P2, P4])
""")

WW = vector((s14*s25+s12*s45, -s14*s15-s11*s45, -s12*s15+s11*s25))

print(vector((alphaR, betaR, gammaR)) == 27 * WW * matrix([P1, P2, P4]))

print("""\nHence the coefficients alpha, beta, gamma can be defined only
in terms of the matrix whose rows are the points P1, P2, P4 and their 
scalar products.""")

print("""So far we have proved that the eigenpoints are the common zero 
of three polynomials of degree 3 in x, y, z, which are: 
{Gg1a, Gg1b, Gg1c} (or {Gg2a, Gg2b, Gg2c} or {Gg3a, Gg3b, Gg3c}), at least
for the general case, there are some particular cases which should be
considered""")












