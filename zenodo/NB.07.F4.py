#!/usr/bin/env python
# coding: utf-8

# # Configuration $(C_6)$

# In[2]:


load("basic_functions.sage")


# ## Six points in configuration $(C_6)$

# We assume the collinearities given by:
# 
# $(1, 2, 3), (1, 4, 5), (2, 4, 6), (3, 5, 6)$
# 
# We construct $P_2, \dots, P_5$, four generic points and we 
# define $P_1$ and $P_6$ in such a way that $P_1, \dots, P_6$ are 
# in configuration $(C_6)$.
# 
# We verify that $P_1$ and $P_6$ are always be defined.
# 
# If config $(C_6)$ is realizable by eigenpoints, then 
# $\delta_1(P_2, P_3, P_4), \delta_1(P_3, P_2, P_5)$, 
# $\delta_1(P_5, P_3, P_4), \delta_1(P_4, P_5, P_2)$, 
# $\delta_1(P_1, P_3, P_5), \delta_1(P_6, P_4, P_5)$
# must be zero.
# 
# We compute these polynomials and we simplify them, since two of 
# them have some factors which are surely not zero.
# 
#     

# In[41]:


## we construct 6 points in configuration (6)
## P2, P3, P4, P5 are generic, so there are no collinearities among them.

P2 = vector(S, (A2, B2, C2))
P3 = vector(S, (A3, B3, C3))
P4 = vector(S, (A4, B4, C4))
P5 = vector(S, (A5, B5, C5))

## P1 and P6 are as follows:
P1 = vector(S, list(intersection_lines(P2, P3, P4, P5)))
P6 = vector(S, list(intersection_lines(P2, P4, P3, P5)))

# P1 and P6 are always defined:
J1 = S.ideal(list(P1))
J1 = J1.saturation(S.ideal(matrix([P2, P3]).minors(2)))[0]
J1 = J1.saturation(S.ideal(matrix([P4, P5]).minors(2)))[0]
J1 = J1.saturation(matrix([P2, P4, P5]).det())[0]
assert(J1 == S.ideal(1))

J1 = S.ideal(list(P6))
J1 = J1.saturation(S.ideal(matrix([P2, P3]).minors(2)))[0]
J1 = J1.saturation(S.ideal(matrix([P3, P5]).minors(2)))[0]
J1 = J1.saturation(matrix([P2, P4, P5]).det())[0]
assert(J1 == S.ideal(1))

d1 = delta1(P2, P3, P4)
d2 = delta1(P3, P2, P5)
d3 = delta1(P5, P3, P4)
d4 = delta1(P4, P5, P2)
d5 = delta1(P1, P3, P5)
d6 = delta1(P6, P4, P5)

## we simplify d1, ..., d6:
## we note indeed that delta1(P6, P4, P5) is divisible by 
## det(matrix([P3, P4, P5]))*det(matrix([P2, P4, P5]))
## and  delta1(P1, P3, P5) is divisible by: 
## det(matrix([P3, P4, P5]))*det(matrix([P2, P3, P5]))

assert(d6.quo_rem(det(matrix([P3, P4, P5]))*det(matrix([P2, P4, P5])))[1] == 0)
d6 = d6.quo_rem(det(matrix([P3, P4, P5]))*det(matrix([P2, P4, P5])))[0]

assert(d5.quo_rem(det(matrix([P3, P4, P5]))*det(matrix([P2, P3, P5])))[1] == 0)
d5 = d5.quo_rem(det(matrix([P3, P4, P5]))*det(matrix([P2, P3, P5])))[0]


# Now we define the ideal $J$ generated by $d_1, \dotsc, d_6$ 
# and we saturate it with polynomials which cannot be zero, and we get 
# the ideal $(1)$, so $(C_6)$ cannot be realized.

# In[42]:


J = S.ideal(d1, d2, d3, d4, d5, d6)
J = J.saturation(det(matrix([P2, P3, P4])))[0]
J = J.saturation(det(matrix([P2, P3, P5])))[0]
J = J.saturation(det(matrix([P3, P4, P5])))[0]

## we get that now J = (1):
assert(J == S.ideal(1))
