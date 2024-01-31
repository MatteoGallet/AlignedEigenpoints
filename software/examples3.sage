## Here we show that it is possible to have seven rational eigenpoints
## p1, ..., p7 such that the condition delta2(p1, p2, p3, p4, p5) = 0
## is satisfied. 

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



U1 = delta2(P1, P2, P3, P4, P5).coefficient(u1)
U2 = delta2(P1, P2, P3, P4, P5).coefficient(u2)

PP3 = P3.subs({u1:U2, u2:-U1})

### P1, P2, PP3, P4, P5 are generic such that
### delta2(P1, P2, PP3, P4, P5) = 0.

## Now we make the following substitution (more or less random numbers)
## v1 and v2 are not jet assigned.

sst = {A1:2, B1:1, C1:-3, A2:4, B2:2, C2:5, A4:-7, B4:5, C4:2}
sst = {A1:2, B1:1, C1:-8, A2:4, B2:2, C2:5, A4:-11, B4:5, C4:2}

p1 = P1.subs(sst)
p2 = P2.subs(sst)
p3 = PP3.subs(sst)
p4 = P4.subs(sst)
p5 = P5.subs(sst)
p6 = vector(R, (x, y, z))

## we construct three polynomials of degree 3 in x, y, z whose
## common zeros are the eigenpoints:

ga = det(matrix([phi_p(p1)[0], phi_p(p1)[1], phi_p(p2)[0], phi_p(p2)[1], phi_p(p3)[0], phi_p(p3)[1], \
           phi_p(p4)[0], phi_p(p4)[1], phi_p(p5)[0], phi_p(p6)[0]]))
gb = det(matrix([phi_p(p1)[0], phi_p(p1)[1], phi_p(p2)[0], phi_p(p2)[1], phi_p(p3)[0], phi_p(p3)[1], \
           phi_p(p4)[0], phi_p(p4)[1], phi_p(p5)[0], phi_p(p6)[1]]))
gc = det(matrix([phi_p(p1)[0], phi_p(p1)[1], phi_p(p2)[0], phi_p(p2)[1], phi_p(p3)[0], phi_p(p3)[1], \
           phi_p(p4)[0], phi_p(p4)[1], phi_p(p5)[0], phi_p(p6)[2]]))

## this factorize into three lines. One is the line p1+p2+p3, another is
## the line p1+p4+p5. We are interested in the last one which is p1+p6+p7
rt1 = p1[2]*ga-p1[1]*gb+p1[0]*gc

## this is the line p1+p6+p7:
f1 = rt1.factor()[-1][0]

## The following is a conic through the points p4, p5, p6, p7 plus the line
## p1+p2+p3:
cn1 = p2[2]*ga-p2[1]*gb+p2[0]*gc

## here is the conic thorough p6 and p7
f2 = cn1.factor()[-1][0]


## we want to eliminate the variable x from f1 and f2
ff = f2*f1.coefficient(x)^2

elimX = S(ff.subs(x = (-y*f1.coefficient(y)-z*f1.coefficient(z))/f1.coefficient(x)))

FF = elimX.discriminant(y).factor()

## FF splits into a product of several squares plus the last factor,
## which, in general, is not a square.

psQ = FF[-1][0]

## we give random values to v1 and v2 and we select the cases in which
## psQ is a square. If this is the case, f1 and f2 meet in 7 rational
## points.

possibiliCasiOK = []
for j in range(1): ## use, for instance, range(1000):
    rnd = [floor(-500+1000*random()) for i in range(2)]
    if rnd[0] == 0 or rnd[1] == 0:
        continue
    sst1 = {v1:rnd[0], v2:rnd[1]}
    pp1 = p1.subs(sst1)
    pp2 = p2.subs(sst1)
    pp3 = p3.subs(sst1)
    pp4 = p4.subs(sst1)
    pp5 = p5.subs(sst1)
    if det(matrix([pp1, pp2, pp4])) == 0:
        continue
    if ZZ(psQ.subs(sst1)).is_square():
        possibiliCasiOK.append(sst1)
        print("maybe there is an example")

## for instance we get the following values for v1 and v2:
## {v1:1, v2:1}  ## ci sono 6 allineamenti
## {v1: 4, v2:3}  ### in this case p4 = p6 and p5 = p7
## {v1: 2, v2:1}  ## ci sono 6 allineamenti
## {v1: -2, v2: 13}  ## esempio con solo 3 allineamenti: quelli giusti
## {v1: 53, v2: 39}
## {v1: 28, v2:19}
## {v1:8, v2:11}
## {v1: -7, v2: 9}
## {v1: 29, v2: 25}

## For instance, the first case (v1 = 1, v2 = 1) gives the cubic:
## 35*x^3 + 78*x^2*y + 54*x*y^2 - 45*y^3 + 6*x^2*z + 72*x*y*z - 48*y^2*z - 72*x*z^2 - 36*y*z^2 + 80*z^3

## whose eigenpoints are:
## [-2, -1, 8], [4, 2, 5], [12, 6, 1], [-11, 5, 2], [3, -2, 2], [-8, 38, 11], [0, -2, 1]

print("In this file you get for instance the cubic:\n")
print(35*x^3 + 78*x^2*y + 54*x*y^2 - 45*y^3 + 6*x^2*z + 72*x*y*z - 48*y^2*z - 72*x*z^2 - 36*y*z^2 + 80*z^3)

print("""\nwhich has 7 rational eigenpoints and delta2 = 0.
Moreover, in the file it is explained how to find other cubics with the 
same properties""")
