## data 24 giugno 22
## dati tre punti allineati, si calcola il caso in cui la matrice M
## ha rango minore di 6

print("Definiamo l'anello opportuno e i punti")

varAn1 = ["s"+str(i)+str(j) for i in range(1,8) for j in range(i, 8)]
varAn2 = ["x", "y", "z", "u1", "u2", "v1", "v2", "w1", "w2", "l1", "l2"]
varAn3 = ["a"+str(i) for i in range(10)]
varAn4 = [str(XX)+str(i) for i in range(1, 8) for XX in ["A", "B", "C"]]
K = QQ
S = PolynomialRing(K, varAn1+varAn2+varAn3+varAn4)
S.inject_variables(verbose=false)

P1, P2 = vector((A1, B1, C1)), vector((A2, B2, C2))
P3 = u1*P1 + u2*P2

mon = [x^3, x^2*y, x*y^2, y^3, x^2*z, x*y*z, y^2*z, x*z^2, y*z^2, z^3]
load("auxiliaryProcedures1.sage")

print("""Definiamo la matrice del sistema lineare che da' le cubiche 
con tre punti allineati e i suoi minori di ordine 6. La matrice 
usa un risultato di teoria che dice che basta prendere due righe 
per Phi(P1), due per Phi(P2) e due per Phi(P3)""")

M = matrix(phi_p(P1)[:2]+phi_p(P2)[:2]+phi_p(P3)[:2])
M6 = M.minors(6)
print("")
print("""Semplifichiamo l'ideale dei minori, dividendo per A1*A2*A3*u1^2*u2^2
e costriamo l'ideale corrispondente""")
print("")
MM6 = list(map(lambda uu: uu/(A1*A2*(u1*A1+u2*A2)*u1^2*u2^2), M6))
J = S.ideal(MM6)
print("Saturiamo rispetto alla condizione P1 = P2")
print("Qualche secondo di attesa...")
sleep(1)
ttA = cputime()
print("")
Js = J.saturation(S.ideal(matrix([P1, P2]).minors(2)))
print("Tempo calcolo saturazione: "+str(cputime()-ttA))
print("")
print("Calcolo della dec primaria")
PD = Js[0].radical().primary_decomposition()
print("")
print("""Si ottengono tre ideali che sono della forma 
      ((Pi|P1), (Pi|P2), (Pi|P3))
con i = 1, 2, 3""")
print("")
print(PD[0] == S.ideal(scalarProd(P2, P2), scalarProd(P1, P2)))
print(PD[1] == S.ideal(scalarProd(P1, P1), scalarProd(P1, P2)))
print(PD[2] == S.ideal(scalarProd(P1, P3), scalarProd(P2, P3), scalarProd(P3, P3)).primary_decomposition()[1])

print("\n\nFINE FILE")