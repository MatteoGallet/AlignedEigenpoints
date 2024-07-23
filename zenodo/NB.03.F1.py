#!/usr/bin/env python
# coding: utf-8

# # Proposition

# Let $P_1, P_2, P_4$ be three distinct points of the plane. Then:
# * $5 \leq \mathrm{rk} \,\Phi(P_1, P_2, P_4) \leq 6$;
# * if $\mathrm{rk} \, \Phi(P_1, P_2, P_4) = 5$, then $P_1, P_2, P_4$
#   are aligned and the line joining them is tangent to the isotropic conic
#   in one of the three points.

# In[5]:


load("basic_functions.sage")


# ## Proof of $\mathrm{rk} \,\Phi(P_1, P_2, P_4) \geq 5$

# We distinguish two cases: $P_1 = (1: 0: 0)$ and $P_1 = (1: i: 0)$.

# ### Case $P_1 = (1: 0: 0)$

# We define the three points.

# In[31]:


P1 = vector((1, 0, 0))
P2 = vector((A2, B2, C2))
P4 = vector((A4, B4, C4))


# We define the matrix of conditions of $P_1$, $P_2$, and $P_4$.

# In[32]:


M = condition_matrix([P1, P2, P4], S, standard="all")


# We compute the ideal of minors of order $5$ of $M$.

# In[33]:


J5 = S.ideal(M.minors(5))


# In[34]:


J5 = J5.saturation(
    S.ideal(matrix([P1, P2]).minors(2))
)[0].saturation(
    S.ideal(matrix([P1, P4]).minors(2))
)[0].saturation(
    S.ideal(matrix([P2, P4]).minors(2))
)[0].radical()


# $J_5$ is the ideal $(1)$, so the matrix $M$ cannot have rank $< 5$.

# In[35]:


assert(J5 == S.ideal(S.one()))


# ### Case $P_1 = (1: i: 0)$

# We define the three points.

# In[6]:


P1 = vector((1, ii, 0))
P2 = vector((A2, B2, C2))
P4 = vector((A4, B4, C4))


# We define the matrix of conditions of $P_1$, $P_2$, and $P_4$.

# In[51]:


M = condition_matrix([P1, P2, P4], S, standard="all")


# We compute the ideal of minors of order $5$ of $M$.

# In[52]:


J5 = S.ideal(M.minors(5))


# In[53]:


J5 = J5.saturation(
    S.ideal(matrix([P1, P2]).minors(2))
)[0].saturation(
    S.ideal(matrix([P1, P4]).minors(2))
)[0].saturation(
    S.ideal(matrix([P2, P4]).minors(2))
)[0].radical()


# $J_5$ is the ideal $(1)$, so the matrix $M$ cannot have rank $< 5$.

# In[54]:


assert(J5 == S.ideal(S.one()))


# ## Proof of $\mathrm{rk} \,\Phi(P_1, P_2, P_4) = 5$ if and only if $P_1$, $P_2$, and $P_4$ are aligned and the line joining them is tangent to the isotropic conic in one of the three points.

# We distinguish two cases: $P_1 = (1: 0: 0)$ and $P_1 = (1: i: 0)$.

# ### Case $P_1 = (1: 0: 0)$

# We define the three points.

# In[7]:


P1 = vector((1, 0, 0))
P2 = vector((A2, B2, C2))
P4 = vector((A4, B4, C4))


# We define the matrix of conditions of $P_1$, $P_2$, and $P_4$.

# In[8]:


M = condition_matrix([P1, P2, P4], S, standard="all")


# We compute the ideal of minors of order $6$ of $M$.

# In[9]:


J6 = S.ideal(M.minors(6))


# In[12]:


J6 = J6.saturation(
    S.ideal(matrix([P1, P2]).minors(2))
)[0].saturation(
    S.ideal(matrix([P1, P4]).minors(2))
)[0].saturation(
    S.ideal(matrix([P2, P4]).minors(2))
)[0]


# We compute the primary decomposition of $J_6$

# In[13]:


pd = J6.radical().primary_decomposition()


# We claim we have only two possibilities: 
# * either $P_1$, $P_2$, $P_4$ are aligned and the line is tangent to the isotropic
# conic in $P_2$ (hence $P_2$ orthogonal to $P_4$, $P_2$ orthogonal to $P_2$ and 
# $P_1$, $P_2$, $P_4$ aligned):

# In[44]:


H6 = S.ideal(
    scalar_product(P1, P2), 
    scalar_product(P2, P2), 
    det(matrix([P1, P2, P4]))
).saturation(
    S.ideal(list(P2))
)[0].radical()
PD1 = H6.primary_decomposition()


# * or $P_1$, $P_2$, $P_4$ are aligned and the line is tangent to the isotropic
# conic in $P_4$ (hence $P_1$ orthogonal to $P_4$, $P_4$ orthogonal to $P_4$ and 
# $P_1$, $P_2$, $P_4$ aligned):

# In[48]:


K6 = S.ideal(
    scalar_product(P1, P4), 
    scalar_product(P4, P4), 
    det(matrix([P1, P2, P4]))
).saturation(
    S.ideal(list(P4))
)[0].radical()
PD2 = K6.primary_decomposition()


# We check that indeed these are the only two possibilities.

# In[49]:


assert(Set(pd) == Set(PD1 + PD2))


# ### Case $P_1 = (1: i: 0)$

# We define the three points.

# In[6]:


P1 = vector((1, ii, 0))
P2 = vector((A2, B2, C2))
P4 = vector((A4, B4, C4))


# We define the matrix of conditions of $P_1$, $P_2$, and $P_4$.

# In[51]:


M = condition_matrix([P1, P2, P4], S, standard="all")


# We compute the ideal of minors of order $5$ of $M$.

# In[ ]:


J6 = S.ideal(M.minors(6))


# In[ ]:


J6 = J6.saturation(
    S.ideal(matrix([P1, P2]).minors(2))
)[0].saturation(
    S.ideal(matrix([P1, P4]).minors(2))
)[0].saturation(
    S.ideal(matrix([P2, P4]).minors(2))
)[0].radical()


# When $J_6$ is satisfied, we have that $P_1$, $P_2$, $P_4$ are aligned
# and the line $P_1 \vee P_2 \vee P_4$ is tangent to the isotropic conic in $P_4$.

# In[67]:


K6 = S.ideal(
    scalar_product(P1, P2), 
    scalar_product(P1, P4), 
    matrix([P1, P2, P4]).det()
).saturation(
    S.ideal(list(P4))
)[0].radical()


# In[68]:


assert(J6 == K6)

