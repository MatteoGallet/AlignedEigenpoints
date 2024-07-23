#!/usr/bin/env python
# coding: utf-8

# # Theorem

# The variety $\mathcal{L}$ is an irreducible hypersurface.

# Here we provide some computations that aid the proof of the result.

# In[1]:


load("basic_functions.sage")


# In[2]:


P1 = vector(S, (A1, B1, C1))
P2 = vector(S, (A2, B2, C2))
P3 = u1*P1+u2*P2
M = condition_matrix([P1, P2, P3], S, standard="all")


# The following columns of $M$ are linearly dependent:
# 0, 1, 2, 4, 5, 7 or 1, 2, 3, 5, 6, 8 or 4, 5, 6, 7, 8, 9.
# We can verify this directly, as follows:

# In[3]:


for cols in [
    [0, 1, 2, 4, 5, 7],
    [1, 2, 3, 5, 6, 8],
    [4, 5, 6, 7, 8, 9]
]:
    Ma = M.matrix_from_columns(cols)
    assert(
        (
            S.ideal(Ma.minors(6))
        ).is_zero()
    )


# or we can see the dependencies of the columns in a more explicit way. We select the 10 columns of $M$:

# In[4]:


c = {}
for i in range(10):
    c[i] = M.matrix_from_columns([i])


# We call $\alpha$, $\beta$, $\gamma$ the entries of $P_1 \times P_2$

# In[5]:


alpha, beta, gamma = tuple(wedge_product(P1, P2))


# We call $N_1$ and $N_2$ the following matrices: 
# $$
# \left(\begin{array}{ccc}
#     \alpha & 0 & 0 \\
#     0 & \beta & 0 \\
#     0 & 0 & \gamma
# \end{array} \right)
# \quad \text{and} \quad
# \left(\begin{array}{ccc}
#     0 & \alpha & 0 \\
#     \gamma & 0 & 0 \\
#     0 & 0 & \beta
# \end{array} \right)
# $$

# In[6]:


N1 = matrix([[alpha, 0, 0], [0, beta, 0], [0, 0, gamma]])
N2 = matrix([[0, alpha, 0], [gamma, 0, 0], [0, 0, beta]])


# Then we see that 
# * c0, c1, c2, c4, c5, c7
# * c1, c2, c3, c5, c6, c8
# * c4, c5, c6, c7, c8, c9
# 
# are linearly dependend:

# In[7]:


L1 = c[0].augment(c[2]).augment(c[7])
L2 = c[1].augment(c[4]).augment(c[5])
assert(
    (
        (L1*N1 + 2*L2*N2)*wedge_product(P1, P2)
    ).is_zero()
)


# In[8]:


L1 = c[1].augment(c[3]).augment(c[8])
L2 = c[2].augment(c[5]).augment(c[6])
assert(
    (
        (L1*N1 + 2*L2*N2)*wedge_product(P1, P2)
    ).is_zero()
)


# In[9]:


L1 = c[4].augment(c[6]).augment(c[9])
L2 = c[5].augment(c[7]).augment(c[8])
assert(
    (
        (L1*N1 + 2*L2*N2)*wedge_product(P1, P2)
    ).is_zero()
)

