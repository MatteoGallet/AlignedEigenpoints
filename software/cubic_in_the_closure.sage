## In this file we construct a one dimensional family of cubic surfaces
## (depending on a parameter u1) which has four eigenpoints, three (p1, p2, p3)
## are fixed and alligned, the remaining one (called p4) is chosen on a
## fixed line (passing through two fixed points q1 and q2).
## We want to see what happens in the case in which p4 --> q,
## where q is the point of intersection of the lines p1+p2 and q1+q2.


print("We define a ring whith suficientely many variables.")

varAn1 = ["s"+str(i)+str(j) for i in range(1,8) for j in range(i, 8)]
varAn2 = ["x", "y", "z", "u1", "u2", "v1", "v2", "w1", "w2", "l1", "l2"]
varAn3 = ["a"+str(i) for i in range(10)]
varAn4 = [str(XX)+str(i) for i in range(1, 8) for XX in ["A", "B", "C"]]

K = QQ
S = PolynomialRing(K, varAn1+varAn2+varAn3+varAn4)
S.inject_variables(verbose=false)


P1 = vector(S, (A1, B1, C1))

load("auxiliaryProcedures.sage")

## here we start the computations for the  specific case.
## We define three collinear points in a random way:

p1, p2 = vector((1, 2, -3)), vector((4, 1, 5))
p3 = 3*p1+p2

## we define two random points not alligned with p1 and p2

q1 = vector((3, 1, 1))
q2 = vector((-2, 4, 7))

print("We define three collinear points p1, p2, p3:\n")

for pp in [p1, p2, p3]:
    print(pp)

print("\nwe define two other points q1 and q2:\n")
print(q1)
print(q2)

print("\nq1 is not on the line p1+p2 and q2 is not on the line p1+2:\n")

print(det(matrix([p1, p2, q1])) != 0)
print(det(matrix([p1, p2, q2])) != 0)
print("")

## and now we define p4, a point on the line q1+q2
p4 = q1+u1*q2

print("""We define a point p4 on the line q1+q2, which depends on one
parameter u1:\n""")
print(p4)
## p4 depends on the parameter u1 and is the intersection point of
## (p1+p2) and (q1+q2) when u1 =  15/143

print("\nIt is on the line p1+p2 when the following expression is zero:\n")
print(det(matrix([p1, p2, p4])))

## we define the matrix 8x10 of the conditions that p1, p2, p3, p4 are
## eigenpoints:

print("\nNow we construct a cubic which has p1, p2, p3, p4 as eigenpoints.\n")

M = matrix([phi_p(p1)[0], phi_p(p1)[1], phi_p(p2)[0], phi_p(p2)[1], \
           phi_p(p3)[0], phi_p(p3)[1], phi_p(p4)[0], phi_p(p4)[1]])

## we define the 10 monomials of degree 3 in x, y, z:
mon = [x^3, x^2*y, x*y^2, y^3, x^2*z, x*y*z, y^2*z, x*z^2, y*z^2, z^3]


## M is a rank 8 matrix with 8 rows. We complete it into a 9x10 matrix M1
## adding a random row, hence in this way we get a rank 9 matrix:


M1 = matrix(S, 9, 10)
for i in range(8):
    for j in range(10):
        M1[i, j] = M[i, j]
	
row9 = [2, 4, -1, 5, 6, 9, 0, 7, 2, -3]  ## a random row.
for j in range(10):
    M1[8, j] = row9[j]

## we compute the unique cubic cb1 whose coefficients are the solution of 
## the 9 linear equations in a0, ..., a9 whose asociated matrix is M1:

m9 = M1.minors(9)
cb1 = add([(-1)^i*mon[i]*m9[9-i] for i in range(10)])

## we verify that the cubic cb1 has p1, p2, p3, p4 as eigenpoints.

V1 = matrix([(cb1.derivative(x), cb1.derivative(y), cb1.derivative(z)), (x, y, z)])

w1 = vector(V1.minors(2))

print("""\nWe have constructed a cubic which has p1, p2, p3, p4 as eigenpoints.
If this is correct, here we should get 4 times the zero vector:\n""")
for pp in [p1, p2, p3, p4]:
    print(w1.subs({x:pp[0], y:pp[1], z:pp[2]}))
    
print("")

## cb1 is a polynomial of degree 3 in x, y, z multiplied by
## (2*u1 - 3)  and (143*u1 - 15).
## we select the polynomial without the factors:
cb1 = cb1.factor()[-1][0]


print("""Here is the cubic (which depends on the parameter u1) which has
p1, p2, p3, p4 as eigenvectors\n""")

print(cb1)

## we verify that again this cubic has p1, p2, p3, p4 as eigenpoints:

V1 = matrix([(cb1.derivative(x), cb1.derivative(y), cb1.derivative(z)), (x, y, z)])
w1 = vector(V1.minors(2))

print("Here we should get again 4 times the zero vector:\n")
for pp in [p1, p2, p3, p4]:
    print(w1.subs({x:pp[0], y:pp[1], z:pp[2]}))
    
print("")

## we particularize the cubic cb1 on the value u1 = 15/143 of the
## parameter, in such a way that p4 goes on the line p1+p2

cb2 = cb1.subs({u1:15/143})

## cb2 splits into (p1+p2)^2*l, where l is a line of the plane.

print("\nHere we have a cubic of the form (p1+p2)^2*l:")
print(cb2.factor())
print("")

## the cubic cb1 for u1=0 has 7 eigenpoints:

print("the cubic cb1 for u1=0 has 7 eigenpoints (three in the first ideal):\n")
ep = eigenpoints(cb1.subs({u1:0}))

for ee in ep:
    print(ee)
    print("")
    