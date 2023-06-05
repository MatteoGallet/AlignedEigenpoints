### conto per provare che conf 6 non esiste.

print("We define a ring with sufficiently many variables.")

varAn1 = ["s"+str(i)+str(j) for i in range(1,8) for j in range(i, 8)]
varAn2 = ["x", "y", "z", "u1", "u2", "v1", "v2", "w1", "w2", "l1", "l2"]
varAn3 = ["a"+str(i) for i in range(10)]
varAn4 = [str(XX)+str(i) for i in range(1, 8) for XX in ["A", "B", "C"]]
K = QQ
S = PolynomialRing(K, varAn1+varAn2+varAn3+varAn4)
S.inject_variables(verbose=false)
P = vector(S, (x, y, z))



load("auxiliaryProcedures.sage")

### we redefine the point for this situation.
print("""Is the configuration (6) possible?
We define the points P1, P2, P3, P4 generic and P5 as the intersection 
of P1+P4 and P2+P3 and P6 the intersection of P1+P2 and P3+P4.""")
P1 = vector((A1, B1, C1))
P2 = vector((A2, B2, C2))
P3 = vector((A3, B3, C3))
P4 = vector((A4, B4, C4))

P5 = intersect_lines(P1, P4, P2, P3)
P6 = intersect_lines(P1, P2, P3, P4)

print("""If (6) is possible, we must have:
delta1(P1, P2, P5)=0, 
delta1(P2, P1, P5)=0, 
delta1(P6, P2, P4)=0, 
delta1(P4, P1, P6)=0, 
delta1(P5, P1, P3)=0, 
delta1(P3, P5, P4)=0
We define the ideal generated by these 6 polynomials and we saturate it 
w.r.t. the polynomials 
det(matrix([P1, P2, P3]))
det(matrix([P1, P2, P4]))
det(matrix([P2, P3, P4]))
det(matrix([P1, P3, P4]))
and w.r.t. matrix([P1, P4]).minors(2) and matrix([P3, P4]).minors(2)
""")
L1 = list(map(factor,[delta1(P1, P2, P5), delta1(P2, P1, P5), \
                      delta1(P6, P2, P4), delta1(P4, P1, P6), \
		      delta1(P5, P1, P3), delta1(P3, P5, P4)]))

L2 = list(map(lambda uu: uu[-1][0], L1)) ## first rough saturation
J1 = S.ideal(L2)
J2 = J1.saturation(det(matrix([P1, P2, P3])))[0]
J3 = J2.saturation(S.ideal(matrix([P1, P4]).minors(2)))[0].saturation(S.ideal(matrix([P3, P4]).minors(2)))[0]

print(""" If we get Ideal(1), configuration (6) is not possible\n""")
print(J3)

print("""\nWe conclude that configuration (6) is not possible""")

