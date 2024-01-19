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

P1 = vector(S, (A1, B1, C1))
P2 = vector(S, (A2, B2, C2))
P4 = vector(S, (A4, B4, C4))
P3 = u1*P1+u2*P2
P5 = v1*P1+v2*P4

load("ancillary.sage")

test = SymbolicCheck()

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

## We verify formula (11), namely
## delta1(P1, P2, P4) = < P4, s11*P2 - s12*P1 >

S11 = scalarProd(P1, P1)
S12 = scalarProd(P1, P2)
assert(delta1(P1, P2, P4) == scalarProd(P4, S11*P2 - S12*P1))

## We verify formula (12), namely
## delta1(P1, P2, P4) = < P2, s11*P4 - s14*P1 >

S11 = scalarProd(P1, P1)
S14 = scalarProd(P1, P4)
assert(delta1(P1, P2, P4) == scalarProd(P2, S11*P4 - S14*P1))

## We verify formula (13), namely

## if delta1b(P1, P2, P3) = 0 then
## P1 is on C_iso and P1+P2+P3 is tangent to C_iso at P1, or
## P2 = (s13^2 + s11*s33)*P1 - 2*s11*s13*P2, or
## P3 = (s12^2 + s11*s22)*P1 - 2*s11*s12*P2

S11 = scalarProd(P1, P1)
S12 = scalarProd(P1, P2)
S13 = scalarProd(P1, P3)
S22 = scalarProd(P2, P2)
S33 = scalarProd(P3, P3)

## if delta1b(P1, P2, P3) = 0 then
## P1 is on C_iso and P1+P2+P3 is tangent to C_iso at P1, or
## P2 = (s13^2 + s11*s33)*P1 - 2*s11*s13*P2, or
## P3 = (s12^2 + s11*s22)*P1 - 2*s11*s12*P2

## We have
assert(delta1b(P1, P2, P3) == 2*S11*S12*u1 + (S12^2 + S11*S22)*u2)
## If s11*s12 != 0 and s12^2 + s11*s22 != 0, then
## delta1b(P1, P2, P3) = 0 implies P3 = (s12^2 + s11*s22)*P1 - 2*s11*s12*P2

## Since the situation is symmetric with respect to swapping P2 and P3, we get
## If s11*s13 != 0 and s13^2 + s11*s33 != 0, then
## delta1b(P1, P2, P3) = 0 implies P2 = (s13^2 + s11*s33)*P1 - 2*s11*s13*P2, or

## We are left with analyzing the case
## s11*s12 = 0 and s12^2 + s11*s22 = 0 and s11*s13 = 0 and s13^2 + s11*s33 = 0
I = S.ideal(S11*S12, S12^2 + S11*S22, S11*S13, S13^2 + S11*S33)
## We saturate with respect to the condition u1 = 0, since P3 != P2
I = I.saturation(u1)[0]
## We get that the only relevant situation after saturation is s12 = 0 and s11 = 0
assert(I.radical() == S.ideal(S12, S11))
## namely, the point P1 is on the isotropic conic and P1 is orthogonal to P2
## (and so to all points of the line P1 + P2 + P3)
## This in turn implies that the line P1 + P2 + P2 is tangent to the isotropic conic
## since we see that s12^2 - s11*s22 = 0 in this case, namely the condition
## \\sigma(P1, P2) = 0 holds

## We now show that each of the conditions
## a. P1 is on C_iso and P1+P2+P3 is tangent to C_iso at P1
## b. P2 = (s13^2 + s11*s33)*P1 - 2*s11*s13*P2
## c. P3 = (s12^2 + s11*s22)*P1 - 2*s11*s12*P2
## implies delta1b(P1, P2, P3) = 0

## c. can be readily checked
PP3 = (S12^2 + S11*S22)*P1 - 2*S11*S12*P2
assert(delta1b(P1, P2, PP3) == 0)
## the symmetry of the situation with respect to swapping P2 and P3 implies b.
## We are left to check case a.
assert(delta1(P1, P2, P3) in S.ideal(S12, S12^2 - S11*S22))

## We verify formula (14), namely
## delta2(P1, P2, P3, P4, P5) = 0 if and only if
## MISSING CONDITION
## P3 = (s14*s15*s22 - s12^2*s45)*P1 + (s11*s12*s45 - s12*s14*s15)*P2 or
## P5 = -(s12*s13*s44 - s14^2*s23)*P1 + (s11*s14*s23 - s14*s12*s13)*P4
## It is enough to show that delta2(P1, P2, P3, P4, P5) = (s11*s12*s45 - s12*s14*s15)*u1 - (s14*s15*s22 - s12^2*s45)*u2

S11 = scalarProd(P1, P1)
S12 = scalarProd(P1, P2)
S14 = scalarProd(P1, P4)
S15 = scalarProd(P1, P5)
S22 = scalarProd(P2, P2)
S45 = scalarProd(P4, P5)

assert(delta2(P1, P2, P3, P4, P5) == (S11*S12*S45 - S12*S14*S15)*u1 - (S14*S15*S22 - S12^2*S45)*u2)
