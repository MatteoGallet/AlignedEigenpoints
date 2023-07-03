varAn1 = ["s"+str(i)+str(j) for i in range(1,8) for j in range(i, 8)]
varAn2 = ["x", "y", "z", "u1", "u2", "v1", "v2", "w1", "w2", "l1", "l2"]
varAn3 = ["a"+str(i) for i in range(10)]
varAn4 = [str(XX)+str(i) for i in range(1, 8) for XX in ["A", "B", "C"]]
var("xx")
K = QQ
K.<ii> = NumberField(xx^2+1)
S = PolynomialRing(K, varAn1+varAn2+varAn3+varAn4)
S.inject_variables(verbose=false)

P1 = vector(S, (1, 0, 0))
P2 = vector(S, (0, ii, 1))
P4 = vector(S, (0, -ii, 1))
P3 = u1*P1 + u2*P2
P5 = v1*P1 + v2*P4

load("auxiliaryProcedures.sage")

## this method is used to find the orthogonality of the lines:
def ortLines(p1, p2, q1, q2):
    rt1 = det(matrix([p1, p2, (x, y, z)]))
    rt2 = det(matrix([q1, q2, (x, y, z)]))
    cf1 = [rt1.coefficient(xx) for xx in (x, y, z)]
    cf2 = [rt2.coefficient(xx) for xx in (x, y, z)]
    return(scalarProd(cf1, cf2) == 0)

mon = [x^3, x^2*y, x*y^2, y^3, x^2*z, x*y*z, y^2*z, x*z^2, y*z^2, z^3]

### In this configuration, the line P1+P2 is tangent to
### the isotropic conic in P2 and the line P1+P4 is tangent
### to the isotropic conic in P4 (for all P3 and P5)

M = matrixEigenpoints([P1, P2, P3, P4, P5])

## 
## M[2] is the zero row
assert(M[2] == vector(S, [0 for i in range(10)]))
## M[3]-ii*M[4] is zero
assert(M[3] - ii*M[4] == vector(S, [0 for i in range(10)]))
## u2*M[6]-ii*u2*M[7]+u1*M[8] is zero
assert(u2*M[6] -ii*u2*M[7] +u1*M[8] == vector(S, [0 for i in range(10)]))
## M[9]+ii*M[10] is zero
assert(M[9] + ii*M[10] == vector(S, [0 for i in range(10)]))
## v2*M[12]+ii*v2*M[13]+v1*M[14] is zero
assert(v2*M[12] +ii*v2*M[13] +v1*M[14] == vector(S, [0 for i in range(10)]))
## therefore in the matrix M we can erase the rows: 2, 3, 6, 9, 12

M = M.matrix_from_rows([0, 1, 4, 5, 7, 8, 10, 11, 13, 14])

## M non puo' avere rango piu' basso di 8:
m8 = M.minors(8)
print("La matrice non ha rango minore di 8:")
print(S.ideal(m8).saturation(u1*u2*v1*v2)[0] == S.ideal(S.one()))

print("""\n\nQui finisce dimostazione proposizione \\ref{frecciaFissata}
nel file casoDegenere2.tex""")

print("")




## manipolation of the matrix M:
def rid(v1, v2, pv):
    return(v2-(v2[pv]/v1[pv])*v1)

def ridM(M, nRiga, pv):
    v1 = M[nRiga]
    M1 = [v1] + [rid(v1, M[i], pv) for i in range(M.nrows()) if i != nRiga]
    return(matrix(M1))

## operazioni per ottenere la \ref{matriceBella}
print("Qui si comincia a manipolare la matrice per avere \\ref{matriceBella}")

M1 = ridM(M, 0, 1)
M1 = ridM(M1, 1, 4)
M1 = ridM(M1, 2, 2)
M1 = ridM(M1, 3, 3)
M1 = ridM(M1, 6, 5)
M1 = ridM(M1, 7, 6)

### adesso M1 (equiv ad M) e' piu' bella. Ha due righe tutte nulle

show(M1)

## tolgo le due righe nulle di M1
M1 = M1.matrix_from_rows([0, 1, 2, 3, 4, 5, 6, 8])
M2 = copy(M1)

## M1 (matrice 8x10) di rango 8 e' la matrice di un sistema lineare
## le cui soluzioni danno tutte le cubiche che hanno la configurazione
## P1, P2, P3, P4, P5 di autopunti.
## risolviamo il sistema con Kramer. Usiamo il fatto che il determinante di
## M1 fatto con le colonne 1, 2, 3, 4, 5, 6, 7, 8 (togliendo quindi
## ad M1 la prima e l'ultima colonna) non e' zero.

## matrice dei coefficienti e suo determinante:
coeffMat = [[M1[i, j] for i in range(8)] for j in range(1, 9)]
dCf = matrix(coeffMat).det()

## colonna dei termini noti:
tn = -dCf*(l1*vector(S, [M1[i, 0] for i in range(8)])+l2*vector(S, [M1[i, 9] for i in range(8)]))

## soluzione con kramer:
krm = []
for k in range(8):
    dt = matrix([coeffMat[i] for i in range(k)] + [tn] + [coeffMat[i] for i in range(k+1, 8)]).det()
    krm.append(S(dt/dCf))

print("Verifica soluzione con Kramer: Tutto zero?")
assert(M1*vector([l1*dCf]+krm+[l2*dCf])== vector(S, [0 for _ in range(8)]))

## la soluzione del sistema (dipende dai parametri liberi l1 e l2):
sol = [l1*dCf]+krm+[l2*dCf]

## la famiglia di cubiche soluzione:
cb = sum([sol[i]*mon[i] for i in range(10)])
cb = factor(cb)[-1][0]

print("\nCalcolata la famiglia cb di cubiche con configurazione data\n")
print("Famiglia denominata con \\ref{famigliaCubDim1}")

### autopunti di cb:
ej = S.ideal(matrix(S, [(x, y, z), [cb.derivative(xx) for xx in (x, y, z)]]).minors(2))

ep = ej.saturation(u1*u2*v1*v2)[0].saturation(l2)[0].radical().primary_decomposition()


## si trovano i 5 autopunti P1, P2, P3, P4, P5 e inoltre
## due autopunti che sono dati dall'ideale:
dueP = ep[-1]

## e sono due punti su una retta di equazione:
## y*u2*v1 + (-ii)*z*u2*v1 + y*u1*v2 + (ii)*z*u1*v2
## e su una conica di equazione:
## 3*y*z*v1*l1 + (-3*ii)*z^2*v1*l1 - 6*x*y*v1*l2 + (6*ii)*x*z*v1*l2 + (-4*ii)*x^2*v2*l2 + (2*ii)*y^2*v2*l2 + (2*ii)*z^2*v2*l2

assert(dueP == S.ideal(y*u2*v1 + (-ii)*z*u2*v1 + y*u1*v2 + (ii)*z*u1*v2,\
3*y*z*v1*l1 + (-3*ii)*z^2*v1*l1 - 6*x*y*v1*l2 + (6*ii)*x*z*v1*l2 + (-4*ii)*x^2*v2*l2 + (2*ii)*y^2*v2*l2 + (2*ii)*z^2*v2*l2).saturation(v1*v2)[0])

## y*u2*v1 + (-ii)*z*u2*v1 + y*u1*v2 + (ii)*z*u1*v2
## e' la retta passante per P1 e ortogonale alla retta
## P3+P5, come segue da questi conti:
r35 = matrix([P3, P5, (x, y, z)]).det()
rt = dueP.gens()[0]  ### retta per P1, P6, P7

## la retta rt passa per P1:
assert(rt.subs({x:P1[0], y:P1[1], z:P1[2]}).is_zero())

print("\n\n\n")
print(scalarProd([r35.coefficient(xx) for xx in (x, y, z)], \
                  [rt.coefficient(xx) for xx in (x, y, z)]) == 0)

assert(scalarProd([r35.coefficient(xx) for xx in (x, y, z)], \
                  [rt.coefficient(xx) for xx in (x, y, z)]).is_zero())
		  
### studio delle possibili configurazioni dei punti nelle cubiche
### della famiglia cb.
###
### In generale, la conica che compare in dueP e' irriducibile
assert(dueP.gens()[1].is_prime())

### punto 12341231234 (riferimento sul file casoDegenere2.tex)


###   ===> caso P2, P4, P6 allineati

##
## costruzione di Prt: generico punto su retta rt = P1+P6+P7
ss6 = {x:w1*rt.coefficient(z), y:w2*rt.coefficient(z), \
          z:-rt.coefficient(x)*w1-rt.coefficient(y)*w2}

assert(rt.subs(ss6).is_zero())

Prt = vector(S, (x, y, z)).subs(ss6)

## condition of alignment between P2, P4, Prt. gives w1 = 0 or
## u1*v2-u2*v1 = 0
assert(S.ideal(det(matrix([P2, P4, Prt]))) == S.ideal(w1*(u1*v2-u2*v1)))

### Caso 1A  w1 = 0
## hence we define PP6 on the line rt aligned with P2 and P4
## and we see when (for which values of u1 and u2) it is an eigenpoint

PP6 = Prt.subs({w1:0, w2:1})
J2 = dueP.subs({x:PP6[0], y:PP6[1], z:PP6[2]}).saturation(u1*u2*v1*v2)[0]

ff = J2.gens()[0]
print(ff.subs({u1: ff.coefficient(u2), u2:-ff.coefficient(u1)})==0)

## here is the corresponding cubic curve:
cb2 = cb.subs({u1: ff.coefficient(u2), u2:-ff.coefficient(u1)})
cb2 = cb2.factor()[-1][0]
PP3 = P3.subs({u1: ff.coefficient(u2), u2:-ff.coefficient(u1)})

## and here we have the eigenpoints of cb2:
EJ = S.ideal(matrix(S, [(x, y, z), [cb2.derivative(xx) for xx in (x, y, z)]]).minors(2))
## P67 are the 2 ideals corresponding to the two new eigenpoints
P67 = EJ.saturation(u1*u2*v1*v2)[0].saturation(S.ideal(l1, l2))[0].\
radical().primary_decomposition()[-2:]

## here we obtain true, hence we define the point PP6
print(P67[0].subs({x:0, y: 3*v1*l1 - 2*v2*l2, z: -(2*ii)*v2*l2})== S.ideal(0))

PP6 = vector(S, (0, 3*v1*l1 - 2*v2*l2, -(2*ii)*v2*l2))

print(latex(PP6))

## and now we define PP7:
dueG = P67[1].gens()[:2]
mn2 = matrix([[dueG[0].coefficient(cc) for cc in (x, y, z)], \
        [dueG[1].coefficient(cc) for cc in (x, y, z)]]).minors(2)
sst7 = {x:mn2[2], y:-mn2[1], z: mn2[0]}
print(P67[1].subs(sst7)==S.ideal(0))
PP7 = vector(S, (x, y, z)).subs(sst7)

print(latex(PP7))

print("\n\nHere we have the alignments of the 7 points:")
print(allignments([P1, P2, PP3, P4, P5, PP6, PP7]))
print("\nHence we have a cubic whose eigenpoints are in configuration (5)")

## in the output we give the orthogonality of the points:
print("")
PNT = [P1, P2, PP3, P4, P5, PP6, PP7]
for i in range(7):
    for j in range(i, 7):
        if scalarProd(PNT[i], PNT[j]) == 0:
            print("point "+str(i+1)+" orthogonal to point "+str(j+1))

## in the output we show the kind of orthogonality among
## the lines:

print("\nOrthogonalities among the lines in case P2, P4, P6 aligned:\n ")
RT = [[P1, P2], [P1, P4], [P1, PP6], [P2, P5], [PP3, P4], [PP3, P5]]
lines = ["line123", "line145", "line167", "line25", "line34", "line35"]
for i in range(6):
    for j in range(i, 6):
        if ortLines(*(RT[i]+RT[j])):
            print("line "+lines[i]+" ortogonal to line "+lines[j])



print("\nEnd of case P2 P4 P6 aligned. conf (5)\n")

####################
####################
###
### ==> caso P2, P5, P6 allineati

## costruzione di Prt: generico punto su retta rt = P1+P6+P7
ss6 = {x:w1*rt.coefficient(z), y:w2*rt.coefficient(z), \
          z:-rt.coefficient(x)*w1-rt.coefficient(y)*w2}
rt.subs(ss6)
Prt = vector(S, (x, y, z)).subs(ss6)

## condition of alignment between P2, P5, Prt. gives 
## {w1: -(ii)*u1*v1, w2: u2*v1 - u1*v2}:
det(matrix([P2, P5, Prt])).factor()

## hence we define PP6 on the line rt aligned with P2 and P5
## and we see when (for which values of u1 and u2) it is an eigenpoint

PP6 = Prt.subs({w1: -(ii)*u1*v1, w2: u2*v1 - u1*v2})
J2 = dueP.subs({x:PP6[0], y:PP6[1], z:PP6[2]}).saturation(u1*u2*v1*v2)[0]
ff = J2.gens()[0]
ff = ff.factor()[-1][0]
print(ff.subs({u1: ff.coefficient(u2), u2:-ff.coefficient(u1)})==0)
PP3 = P3.subs({u1: ff.coefficient(u2), u2:-ff.coefficient(u1)})

cb2 = cb.subs({u1: ff.coefficient(u2), u2:-ff.coefficient(u1)})
cb2 = cb2.factor()[-1][0]
## and here we have the eigenpoints of cb2:
EJ = S.ideal(matrix(S, [(x, y, z), [cb2.derivative(xx) for xx in (x, y, z)]]).minors(2))
## We saturate EJ with some conditions that are not compatible and we get the 7 eigenpoints:
P67 = EJ.saturation(u1*u2*v1*v2)[0].saturation(S.ideal(l1, l2))[0].\
saturation(S.ideal(3*l1 + 4*l2, v1 + v2))[0].\
saturation(S.ideal(3*l1 - 4*l2, v1 - v2))[0].\
saturation(S.ideal(3*v2*l1 - 4*v1*l2))[0].radical().primary_decomposition()[:2]


## computation of PP6 and PP7:
dueG = P67[1].gens()[:2]
mn2 = matrix([[dueG[0].coefficient(cc) for cc in (x, y, z)], \
        [dueG[1].coefficient(cc) for cc in (x, y, z)]]).minors(2)
sst6 = {x:mn2[2], y:-mn2[1], z: mn2[0]}
print(P67[1].subs(sst6)==S.ideal(0))
PP6 = vector(S, (x, y, z)).subs(sst6)


dueG = P67[0].gens()[:2]
mn2 = matrix([[dueG[0].coefficient(cc) for cc in (x, y, z)], \
        [dueG[1].coefficient(cc) for cc in (x, y, z)]]).minors(2)
sst7 = {x:mn2[2], y:-mn2[1], z: mn2[0]}
print(P67[0].subs(sst7)==S.ideal(0))
PP7 = vector(S, (x, y, z)).subs(sst7)


print("\n\nHere we have the alignments of the 7 points  (configuration (8)):")
print(allignments([P1, P2, PP3, P4, P5, PP6, PP7]))


## in the output we give the orthogonality of the points:
print("")
PNT = [P1, P2, PP3, P4, P5, PP6, PP7]
for i in range(7):
    for j in range(i, 7):
        if scalarProd(PNT[i], PNT[j]) == 0:
            print("point "+str(i+1)+" orthogonal to point "+str(j+1))

## in the output we show the kind of orthogonality among
## the lines:

print("\nOrthogonalities among the lines in case P2, P4, P6 aligned:\n ")
RT = [[P1, P2], [P1, P4], [P1, PP6], [P2, P5], [PP3, P4], [PP3, P5]]
lines = ["line123", "line145", "line167", "line25", "line34", "line35"]
for i in range(6):
    for j in range(i, 6):
        if ortLines(*(RT[i]+RT[j])):
            print("line "+lines[i]+" ortogonal to line "+lines[j])



print("\nEnd of case of alignment P2, P4, P6 conf (8)\n")




####################
###################
######
###### ==> caso P3, P5, P6 allineati
######
####################

## costruzione di Prt: generico punto su retta rt = P1+P6+P7
ss6 = {x:w1*rt.coefficient(z), y:w2*rt.coefficient(z), \
          z:-rt.coefficient(x)*w1-rt.coefficient(y)*w2}
rt.subs(ss6)
Prt = vector(S, (x, y, z)).subs(ss6)

## condition of alignment between P3, P5, Prt. gives 
## sAux:
ff1 = det(matrix([P3, P5, Prt])).factor()[-1][0]

sAux = {w1: ff1.coefficient(w2), w2: -ff1.coefficient(w1)}

print("\nsostituz ok?")
print(ff1.subs(sAux) == S(0))
print("")
## hence we define PP6 on the line rt aligned with P3 and P5
## and we see when (for which values of u1 and u2) it is an eigenpoint

PP6 = Prt.subs(sAux)

J2 = dueP.subs({x:PP6[0], y:PP6[1], z:PP6[2]}).saturation(u1*u2*v1*v2)[0]


ff = J2.gens()[0]
ff = ff.factor()[-1][0]
print(ff.subs({u1: ff.coefficient(u2), u2:-ff.coefficient(u1)})==0)
PP3 = P3.subs({u1: ff.coefficient(u2), u2:-ff.coefficient(u1)})

cb2 = cb.subs({u1: ff.coefficient(u2), u2:-ff.coefficient(u1)})
cb2 = cb2.factor()[-1][0]
## and here we have the eigenpoints of cb2:
EJ = S.ideal(matrix(S, [(x, y, z), [cb2.derivative(xx) for xx in (x, y, z)]]).minors(2))
## We saturate EJ with some conditions that are not compatible and we get the 7 eigenpoints:

P67 = EJ.saturation(u1*u2*v1*v2)[0].saturation(S.ideal(l1, l2))[0].\
saturation(S.ideal(3*l1 + 4*l2, v1 + v2))[0].\
saturation(S.ideal(3*l1 - 4*l2, v1 - v2))[0].\
saturation(S.ideal(3*v2*l1 - 4*v1*l2))[0].radical().primary_decomposition()[:2]

## computation of PP6 and PP7:
dueG = P67[1].gens()[:2]
mn2 = matrix([[dueG[0].coefficient(cc) for cc in (x, y, z)], \
        [dueG[1].coefficient(cc) for cc in (x, y, z)]]).minors(2)
sst6 = {x:mn2[2], y:-mn2[1], z: mn2[0]}
print(P67[1].subs(sst6)==S.ideal(0))
PP6 = vector(S, (x, y, z)).subs(sst6)


dueG = P67[0].gens()[:2]
mn2 = matrix([[dueG[0].coefficient(cc) for cc in (x, y, z)], \
        [dueG[1].coefficient(cc) for cc in (x, y, z)]]).minors(2)
sst7 = {x:mn2[2], y:-mn2[1], z: mn2[0]}
print(P67[0].subs(sst7)==S.ideal(0))
PP7 = vector(S, (x, y, z)).subs(sst7)


print("\n\nHere we have the alignments of the 7 points (configuration (8)):")
print(allignments([P1, P2, PP3, P4, P5, PP6, PP7]))



## in the output we give the orthogonality of the points:
print("")
PNT = [P1, P2, PP3, P4, P5, PP6, PP7]
for i in range(7):
    for j in range(i, 7):
        if scalarProd(PNT[i], PNT[j]) == 0:
            print("point "+str(i+1)+" orthogonal to point "+str(j+1))

## in the output we show the kind of orthogonality among
## the lines:

print("\nOrthogonalities among the lines in case P2, P4, P6 aligned:\n ")
RT = [[P1, P2], [P1, P4], [P1, PP6], [P2, P5], [PP3, P4], [PP3, P5]]
lines = ["line123", "line145", "line167", "line25", "line34", "line35"]
for i in range(6):
    for j in range(i, 6):
        if ortLines(*(RT[i]+RT[j])):
            print("line "+lines[i]+" ortogonal to line "+lines[j])



print("\nEnd of case of alignment P3, P5, P6 conf (8)\n")
