#!/usr/bin/env python
# coding: utf-8

# # Configuration $(C_8)$

# In[1]:


load("basic_functions.sage")


# We assume all the $V$-configurations of points obtained from $P_1, \dotsc, P_7$ have the matrix of condition of rank 9

# ## First property of 4 points
# 
# Given 4 generic points of the plane, it is not possible that every
# couple of different points are orthogonal.

# In[2]:


P1 = vector((A1, B1, C1))
P2 = vector((A2, B2, C2))
P4 = vector((A4, B4, C4))
P7 = vector((A7, B7, C7))  
 
## The following ideal is (1):

JJ = S.ideal(
    scalar_product(P1, P2), scalar_product(P1, P4),
    scalar_product(P1, P7), scalar_product(P2, P4),
    scalar_product(P2, P7), scalar_product(P4, P7)
).saturation(
    S.ideal(
        matrix([P2, P4]).minors(2)
    )
)[0].saturation(
    matrix([P2, P4, P7]).det()
)[0].saturation(
    S.ideal(list(P1))
)[0]

assert(JJ == S.ideal(1))


# ## Property of 3 points

# ### Given three distinct not collinear points of the plane $P_1, P_2, P_4$:
# ### the three vectors $P_1 \times P_2$, $P_1 \times P_4$, $P_2 \times P_4$ 
# ### are linearly independent.

# In[3]:


ddt = matrix(
    [
        wedge_product(P1, P2),
        wedge_product(P1, P4), 
        wedge_product(P2, P4)
    ]
).det()

assert(ddt == matrix([P1, P2, P4]).det()^2)


# ## A property of 7 eigenpoints in conf. $(C_8)$

# ### If $s_{12} = 0$ and $s_{17} = 0$, then also $s_{27}=0$.

# We define 7 points in fonfiguration $(C_8)$, which is assumed the following:
# 
# $(1, 2, 3), (1, 4, 5), (1, 6, 7), (2, 4, 6), (2, 5, 7), (3, 4, 7)$
# 
# We take $P_2, P_7, P_4$ generic, while $P_1 = P_2 \times P_7$ 
# (since $s_{12}=0$, $s_{17}=0$). $P_3, P_5, P_6$ as intersection of 
# suitable lines

# In[4]:


P2 = vector(S, (A2, B2, C2))
P7 = vector(S, (A7, B7, C7))
P1 = wedge_product(P2, P7)  
P4 = vector(S, (A4, B4, C4))

## hence P3, P5, P6 are forced:

P3 = intersection_lines(P1, P2, P4, P7)
P5 = intersection_lines(P1, P4, P2, P7)
P6 = intersection_lines(P1, P7, P2, P4)

## P1, ..., P7 are in config (C8)
assert(
    alignments([P1, P2, P3, P4, P5, P6, P7]) == 
    [(1, 2, 3), (1, 4, 5), (1, 6, 7), (2, 4, 6), (2, 5, 7), (3, 4, 7)]
)


# It turns out that $P_3$ is not defined
# precisely when $s_{22}=0$ and $s_{27}=0$, which gives $P_1=P_2$
# 
# It turns out that $P_5$ is not defined
# precisely when $s_{24}=0$ and $s_{47}=0$, which gives $P_1=P_4$
# 
# It turns out that $P_6$ is not defined
# precisely when $s_{27}=0$ and $s_{77}=0$, which gives $P_1=P_7$

# In[28]:


J1 = S.ideal(list(P3))
J1 = J1.saturation(S.ideal(matrix([P2, P7]).minors(2)))[0]
J1 = J1.saturation(S.ideal(matrix([P4, P7]).minors(2)))[0]
J1 = J1.saturation(matrix([P2, P4, P7]).det())[0]
assert(J1 == S.ideal(scalar_product(P2, P7), scalar_product(P2, P2)))
assert([J1.reduce(mm) for mm in matrix([P1, P2]).minors(2)] == [S(0), S(0), S(0)])

J1 = S.ideal(list(P5))
J1 = J1.saturation(S.ideal(matrix([P2, P4]).minors(2)))[0]
J1 = J1.saturation(S.ideal(matrix([P4, P7]).minors(2)))[0]
J1 = J1.saturation(matrix([P2, P4, P7]).det())[0]
assert(J1 == S.ideal(scalar_product(P4, P7), scalar_product(P2, P4)))
assert([J1.reduce(mm) for mm in matrix([P1, P4]).minors(2)] == [S(0), S(0), S(0)])

J1 = S.ideal(list(P6))
J1 = J1.saturation(S.ideal(matrix([P2, P7]).minors(2)))[0]
J1 = J1.saturation(S.ideal(matrix([P4, P7]).minors(2)))[0]
J1 = J1.saturation(matrix([P2, P4, P7]).det())[0]
assert(J1 == S.ideal(scalar_product(P2, P7), scalar_product(P7, P7)))
assert([J1.reduce(mm) for mm in matrix([P1, P7]).minors(2)] == [S(0), S(0), S(0)])


# In our hypotheses (matrix of the conditions of the $V$-conf. always of rank 9), we have that
# if configuration $(C_8)$ is given by eigenpoints, we must have
# $e_1 = e_2 = e_3 = 0$, where:
# $e_1 = \delta_1(P_3, P_1, P_4)$, $e_2 = \delta_1(P_5, P_1, P_2)$, 
# $e3 = \delta_1(P_6, P_1, P_2)$

# In[6]:


e1 = delta1(P3, P2, P4)
e2 = delta1(P5, P1, P2)
e3 = delta1(P6, P1, P2)


# We have: $e_2 = 0$

# In[29]:


assert(e2 == S(0))


# We are going to prove that, if $s_{12}=0, s_{17}=0$, then $s_{27}=0$.
# 
# $e_1$ can be obtained in different ways:
# as $\delta_1(P_3, P_2, P_4)$, but also as $\delta_1(P_3, P_2, P_7)$ or \dots
# similarly the others, so we compute three ideals, Je1, Je2 Je3, the
# first is the ideal of all the ways in which $\delta_1(P_3, \dotsc)$ can be computed and similarly for the others.

# In[30]:


Je1 = S.ideal(
    delta1(P3, P1, P7), delta1(P3, P1, P4),
    delta1(P3, P2, P7), delta1(P3, P2, P4)
).saturation(matrix([P2, P4, P7]).det())[0]

Je2 = S.ideal(
    delta1(P5, P1, P2), delta1(P5, P1, P7),
    delta1(P5, P2, P4), delta1(P5, P4, P7)
).saturation(matrix([P2, P4, P7]).det())[0]

Je3 = S.ideal(
    delta1(P6, P1, P2), delta1(P6, P2, P7),
    delta1(P6, P1, P4), delta1(P6, P4, P7)
).saturation(matrix([P2, P4, P7]).det())[0]

## (Je2 is (0), but we leave it for symmetry)


# Then we see when $e_1=0, e_2=0, e_3=0$, and precisely, we compute
# the ideal Je1+Je2+Je3 and we see that it is the ideal $s_{27}$ 
# (up to radical)

# In[31]:


assert((Je1+Je2+Je3).radical() == S.ideal(scalar_product(P2, P7)))


# Conclusion: 
# * $s_{12} = 0, s_{17} = 0$ implies $s_{27} = 0$.
#   
# By symmetry, it also holds: 
# 
# * $s_{12} = 0, s_{27} = 0$ implies $s_{12} = 0$
# * $s_{27} = 0, s_{17} = 0$ implies $s_{12} = 0$
