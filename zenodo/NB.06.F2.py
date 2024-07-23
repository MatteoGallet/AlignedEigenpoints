#!/usr/bin/env python
# coding: utf-8

# # Lemma

# Suppose that $P_1, P_2, P_3, P_4$ are four distinct points belonging to a line $t$.
# A cubic $C$ has $P_1, \dotsc, P_4$ among its eigenpoints if and only if
# all the points of $t$ are eigenpoints of $C$.
# Moreover,
# $$  6 \leq \mathrm{rk} \,\Phi(P_1, P_2, P_3, P_4) \leq 7 $$
# The rank is $6$ if and only if $\sigma(P_1, P_2) = 0$, i.e.\ if
# and only if $t$ is tangent to the isotropic conic.

# In[1]:


load("basic_functions.sage")


# We can assume that $P_1 = (1: 0: 0)$ since at least one of the four points is not on the isotropic conic.
# 
# Then we define $P_2$, $P_3$ and $P_4$ such that are all collinear.
# Note that $u_1$, $u_2$, $v_1$, $v_2$ can be assumed not zero, 
# since we want distinct points.

# In[2]:


P1 = vector((1, 0, 0))
P2 = vector(S, (A2, B2, C2))
P3 = u1*P1 + u2*P2
P4 = v1*P1 + v2*P2


# Construction of the matrix of the linear conditions:

# In[3]:


M = condition_matrix([P1, P2, P3, P4], S, standard="all")


# In order to get the rank of M, we can erase the rows 0, 1, 2 and the 
# columns 1 and 4 (with no conditions on the parameters) 
# and we get a new matrix $MM$ which has rank $r$ iff M has 
# rank $r+2$.

# In[4]:


MM = M.matrix_from_rows_and_columns(
    [3, 4, 5, 6, 7, 8, 9, 10, 11],
    [0, 2, 3, 5, 6, 7, 8, 9]
)


# MM has rank $\leq 5$:

# In[5]:


J6 = S.ideal(MM.minors(6))
assert(J6 == S.ideal(S.zero()))


# Now we want to see when $MM$ has rank $< 5$

# In[6]:


J5 = S.ideal(MM.minors(5))
J5 = J5.saturation(v1*v2*u1*u2*(u2*v1-u1*v2))[0]
J5r = J5.radical()


# $J5r$ has only one ideal which is the ideal generated by $\sigma(P_1, P_2)$:

# In[7]:


assert(J5r == S.ideal(sigma(P1, P2)))


# $MM$ cannot have rank $\leq 4$:

# In[8]:


J4 = S.ideal(MM.minors(4))
J4 = J4.saturation(v1*v2*u1*u2*(u2*v1-u1*v2))[0]
J4 = J4.saturation(S.ideal(matrix([P1, P3]).minors(2)))[0]


# In[9]:


assert(J4 == S.ideal(S.one()))


# A further remark is that if $P_5$ is another point collinear with $P_1$, $P_2$, then the 
# rank of the matrix $\Phi(P_1, P_2, P_3, P_4, P_5)$ is again $\leq 7$, therefore all the points of $t = P_1 \vee P_2$ are eigenpoints for the cubics defined by $\Lambda(M)$.

# In[10]:


P5 = w1*P1+w2*P2
M1 = condition_matrix([P1, P2, P3, P4, P5], S, standard="all")
MM1 = M1.matrix_from_rows_and_columns(
    [3, 4, 5, 6, 7, 8, 9, 10, 11],
    [0, 2, 3, 5, 6, 7, 8, 9]
)


# The rank of MM1 is $\leq 5$, so the rank of M1 is $\leq 7$:

# In[11]:


JJ6 = S.ideal(MM1.minors(6))
assert(JJ6 == S.ideal(S.zero()))


# Also $M_1$ has rank $\leq 6$ iff $\sigma(P_1, P_2) = 0$:

# In[12]:


JJ5 = S.ideal(MM1.minors(5))
JJ5 = JJ5.saturation(v1*v2*u1*u2*(u2*v1-u1*v2))[0]
JJ5r = JJ5.radical()


# JJ5r is the ideal generated by $\sigma(P_1, P_2)$:

# In[13]:


assert(JJ5r == S.ideal(sigma(P1, P2)))


# Here we can conclude that the matrix $M$ can have rank $6$ or $7$ and if $M$ has rank $6$, then the line $t = P_1 \vee P_2$ is tangent to the isotropic conic and the line $t$ is a line of eigenpoints for the cubics obtained by $M$.
# 
# Now we want to see that if a line $t$ is tangent to the isotropic conic in a point $P_1$ and if $P_1, P_2, P_3, P_4$ are four distinct points on $t$, then the rank of $\Phi(P_1, \dotsc, P_4)$ is $6$.

# In[14]:


P1 = vector((1, ii, 0))


# The tangent line to the isotropic conic in $P_1$ is $x+iy$. We define 3 other points on it.

# In[15]:


tg = x+ii*y
P2 = vector(S, (A2*ii, -A2, C2))
P3 = u1*P1+u2*P2
P4 = v1*P1+v2*P2


# Then we define the condition matrix and we verify that it has rank $6$:

# In[16]:


M = condition_matrix([P1, P2, P3, P4], S, standard="all")
assert(M.rank() == 6)
