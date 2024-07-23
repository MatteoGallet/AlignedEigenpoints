#!/usr/bin/env python
# coding: utf-8

# # Configuration $(C_3)$

# In[3]:


load("basic_functions.sage")


# We define 7 generic points in configuration $(C_3)$ and we study the case in which $\delta_1(P_1, P_2, P_4)=0$ 
# and $\delta_1(P_1, P_2, P_6)=0$. We get two linear equations in the coordinates of $P_2$. If the matrix has 
# maximal rank, we construct the unique solution of the system and we call it `PP2`. But we verify that the point 
# `PP2` coincides with $P_1$, which is impossible.

# In[4]:


P1 = vector(S, (A1, B1, C1))
P2 = vector(S, (A2, B2, C2))
P3 = u1*P1+u2*P2
P4 = vector(S, (A4, B4, C4))
P5 = v1*P1+v2*P4
P6 = vector(S, (A6, B6, C6))
P7 = w1*P1+w2*P6
## we study the condition delta1(P1, P2, P4)=0 and 
## delta1(P1, P2, P6) = 0, in order to see if the configuration (3)
## of the figure can be realized in this way.
pl1 = delta1(P1, P2, P4)
pl2 = delta1(P1, P2, P6)

## Here we have two linear equations in the coordinates of P2
## We construct the matrix of the system: 
M = matrix([[pl1.coefficient(aa) for aa in (A2, B2, C2)],\
            [pl2.coefficient(aa) for aa in (A2, B2, C2)]])

## and the solution:
m2 = M.minors(2)
slz = {A2: m2[2], B2: -m2[1], C2: m2[0]}
## we verify that this is the solution of pl1 = pl2 = 0:
assert((pl1.subs(slz), pl2.subs(slz)) == (S(0), S(0)))

## we get that the solution is given PP2:

PP2 = scalar_product(P1, P1)*det(matrix([P1, P4, P6]))*P1

assert(PP2 == P2.subs(slz))

## but with this solution, P1 and PP2 coincide as projective points:

assert(matrix([P1, PP2]).minors(2) == [S(0), S(0), S(0)])


# Hence we consider the case in which the above matrix $M$ does not have maximal rank.
# In this block we are going to see that in this case $P_1$ is on the isotropic conic, 
# so we redefine the points adding the condition that $P_1 = (1: i: 0)$

# In[5]:


## Finally, we want to consider the case in which M does not have 
## maximal rank. 
J = S.ideal(m2)
pdJ = J.radical().primary_decomposition()

## pdJ has two components: det([P1, P4, P6]) and (P1|P1):

assert(len(pdJ) == 2)

assert(pdJ[1] == S.ideal(det(matrix([P1, P4, P6]))))

assert(pdJ[0] == S.ideal(scalar_product(P1, P1)))

## Since P1, P4, P6 are assumed not collinear, it remain to study 
## the case P1 a point on the isotropic conic.

## We can assume therefore P1 = (1, ii, 0) and we redefine the 
## other points:

P1 = vector(S, (1, ii, 0))
P2 = vector(S, (A2, B2, C2))
P3 = u1*P1+u2*P2
P4 = vector(S, (A4, B4, C4))
P5 = v1*P1+v2*P4
P6 = vector(S, (A6, B6, C6))
P7 = w1*P1+w2*P6


# We consider again the case $\delta_1(P_1, P_2, P_4)=0$ and 
# $\delta_1(P_1, P_2, P_6)=0$, we solve the linear system in the 
# coordinates of $P_2$ and we define one more time the seven points (that 
# now will be called $p_1, \dotsc, p_7$), 
# using the coordinates of $P_2$ obtained in this way.

# In[6]:


## we study when delta1(P1, P2, P4) and 
## delta1(P1, P2, P6) is zero:

pl1 = delta1(P1, P2, P4)
pl2 = delta1(P1, P2, P6)

## Here we have two linear equations in the coordinates of P2
## We construct the matrix of the system: 
M = matrix([[pl1.coefficient(aa) for aa in (A2, B2, C2)],\
            [pl2.coefficient(aa) for aa in (A2, B2, C2)]])

## and the solution:
slz = {A2: (-ii)*B2*A4 + B2*B4, B2: (A4+ii*B4)*B2, C2: (A4+ii*B4)*C2}
assert((pl1.subs(slz), pl2.subs(slz)) == (S(0), S(0)))

## we redefine the points, according to this condition on A2, B2, C2:
p1 = P1.subs(slz)
p2 = P2.subs(slz)
p3 = P3.subs(slz)
p4 = P4.subs(slz)
p5 = P5.subs(slz)
p6 = P6.subs(slz)
p7 = P7.subs(slz)

## now delta1(p1, p2, p4) and delta1(p1, p2, p6) are zero:
assert((delta1(p1, p2, p4), delta1(p1, p2, p6)) == (S(0), S(0)))


# We have that $\Phi(p_1, p_4, p_5, p_6, p_7)$ must have rank $\leq 9$ hence 
# $\delta_1(p_1, p_4, p_6) \delta_2(p_1, p_4, p_5, p_6, p_7) = 0$. If the second 
# factor is zero, we have a $\delta_2=0$ and we are done, so we assume 
# $\delta_1(p_1, p_4, p_6) = 0$. But it holds:
# $\delta_1(p_1, p_4, p_6) = (A_6+iB_6)(A_4+iB_4)$, so we have to study two cases.

# In[7]:


assert(delta1(p1, p4, p6) == -(A6 + ii*B6) * (A4 + ii*B4))


# ## Case $A_6+iB_6 = 0$

# In this case we redefine the points and we see that we get a $\delta_2=0$

# In[7]:


slz1 = {A6: -ii*B6}
pp1 = p1.subs(slz1)
pp2 = p2.subs(slz1)
pp3 = p3.subs(slz1)
pp4 = p4.subs(slz1)
pp5 = p5.subs(slz1)
pp6 = p6.subs(slz1)
pp7 = p7.subs(slz1)

## In this case we have delta2(pp1, pp2, pp3, pp6, pp7), 
## so the configuration is realized by a delta2 condition:
assert(delta2(pp1, pp2, pp3, pp6, pp7) == S(0))


# ## Case $A_4+iB_4=0$

# In[8]:


slz1 = {A4: -ii*B4}
pp1 = p1.subs(slz1)
pp2 = p2.subs(slz1)
pp3 = p3.subs(slz1)
pp4 = p4.subs(slz1)
pp5 = p5.subs(slz1)

## In this case we have delta2(pp1, pp2, pp3, pp4, pp5)=0, 
## so the configuration is realized by a delta2 condition:
assert(delta2(pp1, pp2, pp3, pp4, pp5) == S(0))


# Hence we conclude that among the seven points in configuration (C3) we always have a $\delta_2=0$ condition.
