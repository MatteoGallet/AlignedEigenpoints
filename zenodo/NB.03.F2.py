#!/usr/bin/env python
# coding: utf-8

# # Proposition

# Let $P_1, P_2, P_3, P_4$ be four distinct points of the plane such that
# $P_1, P_2, P_3$ are aligned and let $r$ be the line joining them.
# 
# If $\mathrm{rk} \,\Phi(P_1, P_2, P_3, P_4) \leq 7$ then $r$ is tangent to the isotropic conic in one of the three points $P_1$,  $P_2$, and $P_3$.

# In[2]:


load("basic_functions.sage")


# We distinguish two cases: $P_1 = (1: 0: 0)$ and $P_1 = (1: i: 0)$.

# ## Case $P_1 = (1: 0: 0)$

# We define four points, so that $P_1$, $P_2$, and $P_3$ are aligned.

# In[3]:


P1 = vector((1, 0, 0))
P2 = vector((A2, B2, C2))
P3 = u1*P1 + u2*P2
P4 = vector((A4, B4, C4))


# We define the matrix of conditions of $P_1$, $P_2$, $P_3$ and $P_4$.

# In[4]:


M = condition_matrix([P1, P2, P3, P4], S, standard="all")


# We compute the ideal of minors of order $8$ of $M$.

# In[5]:


J8 = S.ideal(M.minors(8))


# We saturate $J_8$ with respect to the conditions that the points $P_1$, $P_2$, $P_3$, and $P_4$ are distinct.

# In[8]:


J8 = J8.saturation(
    S.ideal(matrix(S, [P1, P2]).minors(2))
)[0].saturation(
    S.ideal(matrix(S, [P1, P3]).minors(2))
)[0].saturation(
    S.ideal(matrix(S, [P2, P3]).minors(2))
)[0].saturation(
    S.ideal(matrix(S, [P1, P4]).minors(2))
)[0].saturation(
    S.ideal(matrix(S, [P2, P4]).minors(2))
)[0].saturation(
    S.ideal(matrix(S, [P3, P4]).minors(2))
)[0]


# We saturate $J_8$ with respect to the conditions that the points $P_1$, $P_2$, and $P_4$ are not aligned.

# In[11]:


J8 = J8.saturation(
    S.ideal(matrix(S, [P1, P2, P4]).det())
)[0]


# The condition imposed by $J_8$ is equivalent to the one 
# that the line joining $P_1$, $P_2$, and $P_3$ is tangent to the isotropic conic in one of the three points,
# namely, $\sigma(P_1, P_2) = 0$ and $\left\langle P_1, P_2 \right\rangle = \left\langle P_1, P_3 \right\rangle = 0$.

# In[20]:


assert(J8 == S.ideal(
    [
        sigma(P1, P2),
        scalar_product(P1, P3)*scalar_product(P1, P2)
    ]
))


# ## Case $P_1 = (1: i: 0)$

# We define four points, so that $P_1$, $P_2$, and $P_3$ are aligned.

# In[21]:


P1 = vector((1, ii, 0))
P2 = vector((A2, B2, C2))
P3 = u1*P1 + u2*P2
P4 = vector((A4, B4, C4))


# We define the matrix of conditions of $P_1$, $P_2$, and $P_4$.

# In[22]:


M = condition_matrix([P1, P2, P3, P4], S, standard="all")


# We compute the ideal of minors of order $8$ of $M$.

# In[23]:


J8 = S.ideal(M.minors(8))


# We saturate $J_8$ with respect to the conditions that the points $P_1$, $P_2$, $P_3$, and $P_4$ are distinct.

# In[24]:


J8 = J8.saturation(
    S.ideal(matrix(S, [P1, P2]).minors(2))
)[0].saturation(
    S.ideal(matrix(S, [P1, P3]).minors(2))
)[0].saturation(
    S.ideal(matrix(S, [P2, P3]).minors(2))
)[0].saturation(
    S.ideal(matrix(S, [P1, P4]).minors(2))
)[0].saturation(
    S.ideal(matrix(S, [P2, P4]).minors(2))
)[0].saturation(
    S.ideal(matrix(S, [P3, P4]).minors(2))
)[0]


# We saturate $J_8$ with respect to the conditions that the points $P_1$, $P_2$, and $P_4$ are not aligned.

# In[25]:


J8 = J8.saturation(
    S.ideal(matrix(S, [P1, P2, P4]).det())
)[0]


# The condition imposed by $J_8$ is equivalent to the one 
# that the line joining $P_1$, $P_2$, and $P_3$ is tangent to the isotropic conic in one of the three points,
# namely, $\left\langle P_1, P_2 \right\rangle = 0$.

# In[28]:


assert(J8 == S.ideal(
    [
        scalar_product(P1, P2)
    ]
))

