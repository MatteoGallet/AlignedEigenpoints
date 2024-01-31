## Here we construct examples
## Case (P1|P1) = 0 and (P1|P2) = 0.
print("We define a ring with sufficiently many variables.")

varAn1 = ["s"+str(i)+str(j) for i in range(1,8) for j in range(i, 8)]
varAn2 = ["x", "y", "z", "u1", "u2", "v1", "v2", "w1", "w2", "l1", "l2"]
varAn3 = ["a"+str(i) for i in range(10)]
varAn4 = [str(XX)+str(i) for i in range(1, 8) for XX in ["A", "B", "C"]]
varAn5 = ["l", "m", "n"]
var("xx")
K.<im> = NumberField(xx^2+1)
S = PolynomialRing(K, varAn1+varAn2+varAn3+varAn4+varAn5)
S.inject_variables(verbose=false)


P1, P2, P4  = vector((A1, B1, C1)), vector((A2, B2, C2)), vector((A4, B4, C4))
P3 = u1*P1 + u2*P2
P5 = v1*P1 + v2*P4

mon = [x^3, x^2*y, x*y^2, y^3, x^2*z, x*y*z, y^2*z, x*z^2, y*z^2, z^3]
## the generic cubic
cb = add([mon[i]*S(varAn3[i]) for i in range(10)])

load("auxiliaryProcedures.sage")

## we define PP1 and PP2 such that (PP1|PP1) = 0 and (PP1|PP2) = 0,
## so the rank of the matrix [Phi(PP1), Phi(PP2), Phi(PP3)] is 5.

PP1 = vector(S, (m^2-n^2, 2*m*n, im*(m^2+n^2)))
s2 = scalarProd(PP1, P2)
PP2 = P2.subs({A2: -(s2.coefficient(B2)*B2+s2.coefficient(C2)*C2), \
               B2:s2.coefficient(A2)*B2, C2:s2.coefficient(A2)*C2})
PP3 = u1*PP1+u2*PP2

print("We have two generic points such that (PP1|PP1) = 0 and (PP1|PP2) = 0\n")
print((scalarProd(PP1, PP1), scalarProd(PP1, PP2))==(0,0))
print("")
M = matrix([phi_p(PP1)[0], phi_p(PP1)[1], phi_p(PP2)[0], phi_p(PP2)[1], \
            phi_p(PP3)[0], phi_p(PP3)[1]])
## the following instruction shows that the rank is less
## than 6 (hence 5):
## print(M.minors(6))

## specific examples:
sst = {m:2, n:3, B2:5, C2:-7, u1:3, u2:7}
p1, p2, p3 = PP1.subs(sst), PP2.subs(sst), PP3.subs(sst)
p4 = vector(S, (1, 5, 11))
p5 = vector(S, (-5, 1, 7))

## with these points we get a matrix of the Phi of the points of
## rank 9:
M = matrix([phi_p(p1)[0], phi_p(p1)[1], phi_p(p2)[0], phi_p(p2)[1], \
            phi_p(p3)[0], phi_p(p3)[1], phi_p(p4)[0], phi_p(p4)[1], \
	    phi_p(p5)[0], phi_p(p5)[1]])

print("""We define specific points p1, ..., p5 such that the rank of 
the matrix of the Phi of the points have rank 9 (p4, p5 are 
random points of the plane):\n""")

print(M.rank()==9)

## We construct the cubic with p1, ..., p5 as eigenpoints:

nEqz = M.nrows()
sistLin = [add([M[j][i]*S(varAn3[i]) for i in range(10)]) for j in range(nEqz)]

cb1 = S.ideal(sistLin).reduce(cb)

freeVar = Set(cb1.variables()).difference(Set([x, y, z]))
randomValues = [3, 7, 11, -12, 23]
sstA = {freeVar[i]:randomValues[i] for i in range(len(freeVar))}
cb1 = cb1.subs(sstA)

EP = eigenpoints(cb1)

print("\n The eigenpoints:\n")
for ep in EP:
    print(ep)
    print("")


print("NUOVO")
## specific examples:
sst = {m:1, n:2, B2:-1, C2:1, u1:-1, u2:1}
p1, p2, p3 = PP1.subs(sst), PP2.subs(sst), PP3.subs(sst)
p4 = vector(S, (1, 1, 1))
p5 = vector(S, (1, -1, 1))

## with these points we get a matrix of the Phi of the points of
## rank 9:
M = matrix([phi_p(p1)[0], phi_p(p1)[1], phi_p(p2)[0], phi_p(p2)[1], \
            phi_p(p3)[0], phi_p(p3)[1], phi_p(p4)[0], phi_p(p4)[1], \
	    phi_p(p5)[0], phi_p(p5)[1]])

print("""We define specific points p1, ..., p5 such that the rank of 
the matrix of the Phi of the points have rank 9 (p4, p5 are 
random points of the plane):\n""")

print(M.rank()==9)

## We construct the cubic with p1, ..., p5 as eigenpoints:

nEqz = M.nrows()
sistLin = [add([M[j][i]*S(varAn3[i]) for i in range(10)]) for j in range(nEqz)]

cb1 = S.ideal(sistLin).reduce(cb)

freeVar = Set(cb1.variables()).difference(Set([x, y, z]))
randomValues = [3, 7, 11, -12, 23]
sstA = {freeVar[i]:randomValues[i] for i in range(len(freeVar))}
cb1 = cb1.subs(sstA)

EP = eigenpoints(cb1)

print("\n The eigenpoints:\n")
for ep in EP:
    print(ep)
    print("")


print("")
print("Second example: the line p1 p2 is a line of eigenpoints.")
sst = {m:2, n:3, B2:5, C2:-7, u1:3, u2:7}
p1, p2, p3 = PP1.subs(sst), PP2.subs(sst), PP3.subs(sst)
p4 = 2*p1+5*p2
p5 = vector(S, (-5, 1, 7))
## with these points we get a matrix of the Phi of the points of
## rank 8:
M = matrix([phi_p(p1)[0], phi_p(p1)[1], phi_p(p2)[0], phi_p(p2)[1], \
            phi_p(p3)[0], phi_p(p3)[1], phi_p(p4)[0], phi_p(p4)[1],\
	    phi_p(p5)[0], phi_p(p5)[1]])

print("""We define specific points p1, ..., p5 such that the rank of 
the matrix of the Phi of the points have rank 8 (p3 and p4 are random points
on the line p1+p2, p5 is a random point of the plane:\n""")

print(M.rank()==8)

## We construct the cubic with p1, ..., p4 as eigenpoints:

nEqz = M.nrows()
sistLin = [add([M[j][i]*S(varAn3[i]) for i in range(10)]) for j in range(nEqz)]

cb1 = S.ideal(sistLin).reduce(cb)

freeVar = Set(cb1.variables()).difference(Set([x, y, z]))
randomValues = [3, 7, 11, -12, 23]
sstA = {freeVar[i]:randomValues[i] for i in range(len(freeVar))}
cb1 = cb1.subs(sstA)

EP = eigenpoints(cb1)

print("""\nThe result is a line of eigenpoints (which contains p1, p2, p3, p4)
and the point p5""")
print("\n The eigenpoints:\n")
for ep in EP:
    print(ep)
    print("")


print("")
## specific examples:
sst = {m:2, n:3, B2:5, C2:-7, u1:3, u2:7}
p1, p2, p3 = PP1.subs(sst), PP2.subs(sst), PP3.subs(sst)
p4 = vector(S, (1, 5, 11))
p5 = 3*p1+8*p4

## with these points we get a matrix of the Phi of the points of
## rank 9:
M = matrix([phi_p(p1)[0], phi_p(p1)[1], phi_p(p2)[0], phi_p(p2)[1], \
            phi_p(p3)[0], phi_p(p3)[1], phi_p(p4)[0], phi_p(p4)[1], \
	    phi_p(p5)[0], phi_p(p5)[1]])

print("""We define specific points p1, ..., p5 such that the rank of 
the matrix of the Phi of the points have rank 9 (p4 is a random point
of the plane, p5 is a random point on the line p1+p4):\n""")

print(M.rank()==9)

## We construct the cubic with p1, ..., p5 as eigenpoints:

nEqz = M.nrows()
sistLin = [add([M[j][i]*S(varAn3[i]) for i in range(10)]) for j in range(nEqz)]

cb1 = S.ideal(sistLin).reduce(cb)

freeVar = Set(cb1.variables()).difference(Set([x, y, z]))
randomValues = [3, 7, 11, -12, 23]
sstA = {freeVar[i]:randomValues[i] for i in range(len(freeVar))}
cb1 = cb1.subs(sstA)

EP = eigenpoints(cb1)

print("\n The eigenpoints:\n")
for ep in EP:
    print(ep)
    print("")

########################


##### altro tenmtativo: p1, p2, p3 allineati, p2, p4, p5 allineati
##### in questo esempio dovrebbe succedere che delta1(p2, p1, p4) = 0,
##### ma invece cio' non succede.
print("")
## specific examples:
sst = {m:2, n:3, B2:5, C2:-7, u1:3, u2:7}
p1, p2, p3 = PP1.subs(sst), PP2.subs(sst), PP3.subs(sst)
p4 = vector(S, (1, 5, 11))
p5 = 3*p2+8*p4

## with these points we get a matrix of the Phi of the points of
## rank 9:
M = matrix([phi_p(p1)[0], phi_p(p1)[1], phi_p(p2)[0], phi_p(p2)[1], \
            phi_p(p3)[0], phi_p(p3)[1], phi_p(p4)[0], phi_p(p4)[1], \
	    phi_p(p5)[0], phi_p(p5)[1]])

print("""We define specific points p1, ..., p5 such that the rank of 
the matrix of the Phi of the points have rank 9 (p4 is a random point
of the plane, p5 is a random point on the line p1+p4):\n""")

print(M.rank()==9)

## We construct the cubic with p1, ..., p5 as eigenpoints:

nEqz = M.nrows()
sistLin = [add([M[j][i]*S(varAn3[i]) for i in range(10)]) for j in range(nEqz)]

cb1 = S.ideal(sistLin).reduce(cb)

freeVar = Set(cb1.variables()).difference(Set([x, y, z]))
randomValues = [3, 7, 11, -12, 23]
sstA = {freeVar[i]:randomValues[i] for i in range(len(freeVar))}
cb1 = cb1.subs(sstA)

EP = eigenpoints(cb1)

print("\n The eigenpoints:\n")
for ep in EP:
    print(ep)
    print("")

print("""Here we would like to get zero, but we get something different
from zero:""")
print(delta1(p2, p1, p4))
print("""However, in this case, delta2(p2, p1, p3, p4, p5) is zero:""")
print(delta2(p2, p1, p3, p4, p5) == 0)