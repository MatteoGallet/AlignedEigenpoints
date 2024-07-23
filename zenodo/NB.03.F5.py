#!/usr/bin/env python
# coding: utf-8

# #  Lemma

# Let $P_1, \dots, P_5$ be a $V$-configuration of points and assume that
# $$
#   \left\langle P_1, P_2 \right\rangle = 0 \,, \quad
#   \left\langle P_2, P_2 \right\rangle = 0 \,, \quad
#   \left\langle P_1, P_4 \right\rangle = 0 \,, \quad
#   \left\langle P_4, P_4 \right\rangle = 0 \,.
# $$
# Then the matrix $\Phi(P_1, \dots, P_5)$ has rank $8$.

# This is the case in which the $5$ eigenpoints $P_1$, $P_2$, $P_3$, $P_4$, $P_5$ are
# such that $P_1 \vee P_2$ is tangent at the point $P_2$ to the isotropic conic, 
# and $P_1 \vee P_4$ is tangent to the isotropic conic at $P_4$.

# In[4]:


load("basic_functions.sage")


# The argument in the paper shows that we can choose the $5$ points as follows

# In[16]:


P1 = vector(S, (1, 0, 0))
P2 = vector(S, (0, ii, 1))
P4 = vector(S, (0, -ii, 1))
P3 = u1*P1 + u2*P2
P5 = v1*P1 + v2*P4


# A remark on $\delta_1$ and $\delta_2$:
# $\delta_1(P_1, P_2, P_4$ is not zero, 
# while $\delta_2(P_1, P_2, P_3, P_4, P_5)$ is zero.

# In[17]:


assert(delta1(P1, P2, P4) != 0)
assert(delta2(P1, P2, P3, P4, P5) == 0)


# We define the matrix of conditions of $P_1, \dotsc, P_5$.

# In[18]:


M = condition_matrix([P1, P2, P3, P4, P5], S, standard="all")


# Dependencies between the rows of $M$

# In[19]:


# M[2] is the zero row
assert(M[2] == vector(S, [0 for i in range(10)]))
# M[3]-ii*M[4] is zero
assert(M[3] - ii*M[4] == vector(S, [0 for i in range(10)]))
# u2*M[6]-ii*u2*M[7]+u1*M[8] is zero
assert(u2*M[6] -ii*u2*M[7] +u1*M[8] == vector(S, [0 for i in range(10)]))
# M[9]+ii*M[10] is zero
assert(M[9] + ii*M[10] == vector(S, [0 for i in range(10)]))
# v2*M[12]+ii*v2*M[13]+v1*M[14] is zero
assert(v2*M[12] +ii*v2*M[13] +v1*M[14] == vector(S, [0 for i in range(10)]))


# Therefore in the matrix $M$ we can erase the rows $2$, $3$, $6$, $9$, $12$.

# In[20]:


M = M.matrix_from_rows([0, 1, 4, 5, 7, 8, 10, 11, 13, 14])


# We compute the ideal of order $8$ minors of $M$

# In[21]:


m8 = M.minors(8)
J8 = S.ideal(m8)


# The ideal of order $8$ minors of $M$, once saturated by the conditions that the points are distinct, is the whole ring, so $M$ cannot have lower than $8$

# In[22]:


assert(S.ideal(m8).saturation(u1*u2*v1*v2)[0] == S.ideal(S.one()))


# On the other hand, the rank of $M$ is also $\leq 8$

# In[29]:


assert(M.det() == S.zero())
assert(S.ideal(M.minors(9)) == S.ideal(S.zero()))


# Hence the rank of $M$ is precisely $8$.
