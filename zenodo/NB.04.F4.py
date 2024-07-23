#!/usr/bin/env python
# coding: utf-8

# # Remark

# Example of five points for which the matrix of conditions has rank $8$

# In[3]:


load("basic_functions.sage")


# In[4]:


P1 = vector(S, (1, 0, 0))
P2 = vector(S, (A2, B2, C2))


# We want $\delta_1(P_1, P_2, P_4) = 0$.
# Hence $\left\langle P_4,  s_{11} P_1 - s_{12} P_2 \right\rangle = 0$, so:

# In[5]:


Ptmp4 = vector(S, (A4, B4, C4))
Q4 = scalar_product(P1,P1)*P2 - scalar_product(P1,P2)*P1
aa4, bb4, cc4 = (
    scalar_product(Ptmp4, Q4).coefficient(A4), 
    scalar_product(Ptmp4, Q4).coefficient(B4),
    scalar_product(Ptmp4, Q4).coefficient(C4)
)


# Two alternative definitions of $P_4$ (solving $\delta_1 = 0$ w.r.t. $A_4$ or $C_4$)

# In[6]:


P4 = vector(S, (cc4*A4, cc4*B4, -aa4*A4-bb4*B4))
P4 = vector(S, (bb4*A4, -aa4*A4-cc4*C4, bb4*C4))


# In[7]:


assert(delta1(P1, P2, P4) == 0)


# In[8]:


P3 = (
    (scalar_product(P1, P2)^2+scalar_product(P1, P1)*scalar_product(P2, P2))*P1
    - 2*scalar_product(P1, P1)*scalar_product(P1, P2)*P2
)
assert(delta1b(P1, P2, P3) == 0)


# In[9]:


P5 = (
    (scalar_product(P1, P4)^2+scalar_product(P1, P1)*scalar_product(P4, P4))*P1
    - 2*scalar_product(P1, P1)*scalar_product(P1, P4)*P4
)
assert(delta1b(P1, P4, P5) == 0)


# In[10]:


M = condition_matrix([P1, P2, P3, P4, P5], S, standard="all")


# In[12]:


M1 = M.matrix_from_rows([0, 1, 3, 4, 6, 7, 9, 11, 12, 13])
assert(M1.rank() == 8)


# In[13]:


m9 = M1.minors(9)
assert(Set(m9) == Set([0]))

