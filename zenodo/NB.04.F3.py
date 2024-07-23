#!/usr/bin/env python
# coding: utf-8

# # Proposition

# Let $P_1, \dotsc, P_5$ be a $V$-configuration such that it holds
# $$
#    \delta_1(P_1, P_2, P_4)=\overline{\delta}_1(P_1, P_2, P_3) =
#   \overline{\delta}_1(P_1, P_4, P_5) = 0
# $$
# Then $P_4$ is orthogonal to $s_{11} \, P_2 - s_{12} \, P_1$ and one of the four conditions obtained by considering
# $$
#    P_3 = (s_{12}^2+s_{11}s_{22}) \, P_1 - 2s_{11}s_{12} \, P_2 \,, \quad
#    P_5 = (s_{14}^2+s_{11}s_{44}) \, P_1 - 2s_{11}s_{14} \, P_4 \,.
# $$
#   and swapping in the latter formulas $2 \leftrightarrow 3$ and $4 \leftrightarrow 5$ holds.

# To prove this result, we need to show that none of the following situations may happen:
# * $\left\langle P_1, P_1 \right\rangle = \sigma(P_1, P_2) = 0$ and $P_5 = (s_{14}^2+s_{11}s_{44}) \, P_1 - 2s_{11}s_{14} \, P_4$
# * $\left\langle P_1, P_1 \right\rangle = \sigma(P_1, P_2) = 0$ and $P_4 = (s_{15}^2+s_{11}s_{55}) \, P_1 - 2s_{11}s_{15} \, P_5$
# * $\left\langle P_1, P_1 \right\rangle = \sigma(P_1, P_4) = 0$ and $P_3 = (s_{12}^2+s_{11}s_{22}) \, P_1 - 2s_{11}s_{12} \, P_2$
# * $\left\langle P_1, P_1 \right\rangle = \sigma(P_1, P_4) = 0$ and $P_2 = (s_{13}^2+s_{11}s_{33}) \, P_1 - 2s_{11}s_{13} \, P_3$

# In[5]:


load("basic_functions.sage")


# We define five points, so that $P_1$, $P_2$, and $P_3$ are aligned and $P_1$, $P_4$, and $P_5$ are aligned.

# In[6]:


P1 = vector(S, (A1, B1, C1))
P2 = vector(S, (A2, B2, C2))
P3 = u1*P1 + u2*P2
P4 = vector(S, (A4, B4, C4))
P5 = v1*P1 + v2*P4


# In[7]:


S11 = scalar_product(P1, P1)
S12 = scalar_product(P1, P2)
S13 = scalar_product(P1, P3)
S14 = scalar_product(P1, P4)
S15 = scalar_product(P1, P5)
S22 = scalar_product(P2, P2)
S33 = scalar_product(P3, P3)
S44 = scalar_product(P4, P4)
S55 = scalar_product(P5, P5)


# ## Proof that 
# ## $$ \left\langle P_1, P_1 \right\rangle = \sigma(P_1, P_2) = 0 \text{ and } P_5 = (s_{14}^2+s_{11}s_{44}) \, P_1 - 2s_{11}s_{14} \, P_4 $$
# ## cannot happen.

# In[13]:


I = S.ideal(scalar_product(P1, P1), sigma(P1, P2))
m2 = matrix([(S14^2 + S11*S44)*P1 - 2*S11*S14*P4, P5]).minors(2)
J = S.ideal(m2)


# In[14]:


assert(
    (I + J).radical().saturation(
        S.ideal(list(P1))
    )[0].saturation(
        u1*u2*v1*v2
    )[0].saturation(
        matrix([P1, P2, P4]).det()
    )[0] == S.ideal(S.one())
)


# ## Proof that 
# ## $$ \left\langle P_1, P_1 \right\rangle = \sigma(P_1, P_2) = 0 \text{ and } P_4 = (s_{15}^2+s_{11}s_{55}) \, P_1 - 2s_{11}s_{15} \, P_5 $$
# ## cannot happen.

# In[15]:


I = S.ideal(scalar_product(P1, P1), sigma(P1, P2))
m2 = matrix([(S15^2 + S11*S55)*P1 - 2*S11*S15*P5, P4]).minors(2)
J = S.ideal(m2)


# In[16]:


assert(
    (I + J).radical().saturation(
        S.ideal(list(P1))
    )[0].saturation(
        u1*u2*v1*v2
    )[0].saturation(
        matrix([P1, P2, P4]).det()
    )[0] == S.ideal(S.one())
)


# ## Proof that 
# ## $$ \left\langle P_1, P_1 \right\rangle = \sigma(P_1, P_4) = 0 \text{ and } P_3 = (s_{12}^2+s_{11}s_{22}) \, P_1 - 2s_{11}s_{12} \, P_2 $$
# ## cannot happen.

# In[17]:


I = S.ideal(scalar_product(P1, P1), sigma(P1, P4))
m2 = matrix([(S12^2 + S11*S22)*P1 - 2*S11*S12*P2, P3]).minors(2)
J = S.ideal(m2)


# In[18]:


assert(
    (I + J).radical().saturation(
        S.ideal(list(P1))
    )[0].saturation(
        u1*u2*v1*v2
    )[0].saturation(
        matrix([P1, P2, P4]).det()
    )[0] == S.ideal(S.one())
)


# ## Proof that 
# ## $$ \left\langle P_1, P_1 \right\rangle = \sigma(P_1, P_4) = 0 \text{ and } P_2 = (s_{13}^2+s_{11}s_{33}) \, P_1 - 2s_{11}s_{13} \, P_3 $$
# ## cannot happen.

# In[19]:


I = S.ideal(scalar_product(P1, P1), sigma(P1, P4))
m2 = matrix([(S13^2 + S11*S33)*P1 - 2*S11*S13*P3, P2]).minors(2)
J = S.ideal(m2)


# In[20]:


assert(
    (I + J).radical().saturation(
        S.ideal(list(P1))
    )[0].saturation(
        u1*u2*v1*v2
    )[0].saturation(
        matrix([P1, P2, P4]).det()
    )[0] == S.ideal(S.one())
)

