print("""Here we prove that if P1, P2 are such that 
(P1|P1) = 0, (P1|P2) = 0 (i.e. P1 is orthogonal to r = P1+P2),
then the line r is tangent to the isotropic conic in P1""")
##The ideal given by the intersection of the line r = P1+P2 
## and the isotropic conic.
cis = x^2+y^2+z^2  ## isotropic conic 
J1 = R.ideal(det(matrix([P1, P2, (x, y, z)])), scalarProd(P1, P1), \
        scalarProd(P1, P2), cis)
## We saturate w.r.t. the condition that P1 and P2 are distinct:
J1 = J1.saturation(R.ideal(matrix([P1, P2]).minors(2)))[0]
## We saturate w.r.t. the condition that P1 and (x, y, z) are distinct:
J1 = J1.saturation(R.ideal(matrix([P1, (x, y, z)]).minors(2)))[0]
print(J1 == R.ideal(1))

print("\n\n")

print("""Here we prove that if P1 and P2 are such that 
sigma(P1, P2) = 0, then the line r = P1 + P2 is tangent to the 
isotropic conic and conversely""")

## Proof:
Q3 = l1*P1+l2*P2  ## generic point of the line r = P1+P2

intP = cis.subs({x: Q3[0], y: Q3[1], z: Q3[2]}) ##intersection points r and cis
dscr = intP.discriminant(l1)
print(dscr == -4*l2^2*sigma(P1, P2))

print("\n\n")

