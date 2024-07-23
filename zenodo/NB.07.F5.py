#!/usr/bin/env python
# coding: utf-8

# # Configuration $(C_7)$

# In[1]:


load("basic_functions.sage")


# We consider the configuration $(C_7)$ given by:
# 
# $(1, 2, 3), (1, 4, 5), (1, 6, 7), (2, 4, 6), (2, 5, 7)$
# 
# and we construct 7 generic points in this configuration. So $P_1, P_2, P5, P6$
# are generic, $P_3$ is on the line $P_1, P_2$, $P_4$ and $P_7$ are 
# intersection points.
# 
# We verify that $P_4$ and $P_7$ are always defined.

# In[6]:


P1 = vector(S, (A1, B1, C1))
P2 = vector(S, (A2, B2, C2))
P5 = vector(S, (A5, B5, C5))
P6 = vector(S, (A6, B6, C6))

P4 = vector(S, list(intersection_lines(P1, P5, P2, P6)))
P7 = vector(S, list(intersection_lines(P1, P6, P2, P5)))

P3 = u1*P1+u2*P2

# P4 and P7 are always defined:
J1 = S.ideal(list(P4))
J1 = J1.saturation(matrix([P2, P5, P6]).det())[0]
J1 = J1.saturation(S.ideal(matrix([P1, P5]).minors(2)))[0]
J1 = J1.saturation(S.ideal(matrix([P2, P6]).minors(2)))[0]
assert(J1 == S.ideal(1))

J1 = S.ideal(list(P7))
J1 = J1.saturation(matrix([P2, P5, P6]).det())[0]
J1 = J1.saturation(S.ideal(matrix([P2, P5]).minors(2)))[0]
J1 = J1.saturation(S.ideal(matrix([P1, P6]).minors(2)))[0]
assert(J1 == S.ideal(1))


# In the considered configuration we have:
# * $\delta_1(P_5, P_1, P_2)=0$,
# * $\delta_1(P_6, P_1, P_2)=0$,
# * $\delta_1(P_4, P_1, P_2)=0$,
# * $\delta_1(P_7, P_1, P_2)=0$,
# * $\delta_2(P_1, P_2, P_3, P_4, P_5)=0$,
# * $\delta_2(P_2, P_1, P_3, P_5, P_7)=0$.
# 
# $\delta_1(P_7, P_1, P_2)$ and $\delta_1(P_4, P_1, P_2)$ can be simplified

# In[ ]:


d1 = delta1(P5, P1, P2)
d2 = delta1(P6, P1, P2)
d3 = delta1(P4, P1, P2)
d4 = delta1(P7, P1, P2)
d5 = delta2(P1, P2, P3, P4, P5)
d6 = delta2(P2, P1, P3, P5, P7)

assert(d4.quo_rem(det(matrix([P1, P2, P6]))*det(matrix([P1, P2, P5])))[1] == 0)
assert(d3.quo_rem(det(matrix([P1, P2, P6]))*det(matrix([P1, P2, P5])))[1] == 0)


# We should consider the ideal $(d_1, \dots, d_6)$ but it is too big, 
# so we split the computations and first we define the ideal 
# $J = (d_1, d_2, d_3, d_4)$

# In[ ]:


d4 = d4.quo_rem(det(matrix([P1, P2, P6]))*det(matrix([P1, P2, P5])))[0]
d3 = d3.quo_rem(det(matrix([P1, P2, P6]))*det(matrix([P1, P2, P5])))[0]

J = S.ideal(d1, d2, d3, d4)


# We saturate $J$:

# In[8]:


J = J.saturation(det(matrix([P1, P2, P5])))[0]
J = J.saturation(det(matrix([P2, P5, P6])))[0]


# and now, that it is simpler, we add $d_5$ and $d_6$. We get the ideal $(1)$ so the configuration is not possible.
# 

# In[10]:


J1 = J + S.ideal(d5, d6)
J1 = J1.saturation(S.ideal(det(matrix([P1, P2, P5]))))[0]
J1 = J1.saturation(u2)[0]
J1 = J1.saturation(S.ideal(det(matrix([P1, P5, P6]))))[0]

assert(J1 == S.ideal(1))

