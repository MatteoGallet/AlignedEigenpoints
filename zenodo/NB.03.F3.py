#!/usr/bin/env python
# coding: utf-8

# # Lemma

# * $\delta_1(P_1, P_2, P_4) = 0$ iff $\left\langle P_4, s_{11}P_2-s_{12}P_1 \right\rangle = 0$ iff $\left\langle P_2, s_{11}P_4-s_{14}P_1 \right\rangle = 0$
# * $\bar{\delta}_1(P_1, P_2, P_3) = 0$ iff ($\left\langle P_1, P_1 \right\rangle$ and $P_1 \vee P_2 \vee P_3$ is tangent to the isotropic conic in $P_1$) or $P_3 = (s_{12}^2 + s_{11} s_{22}) P_1 - 2 s_{11} s_{12}P_2$ or $P_2 = (s_{13}^2 + s_{11} s_{33}) P_1 - 2 s_{11} s_{13}P_3$

# In[1]:


load("basic_functions.sage")


# ## Proof of $\delta_1(P_1, P_2, P_4) = 0$ iff $\left\langle P_4, s_{11}P_2-s_{12}P_1 \right\rangle = 0$ iff $\left\langle P_2, s_{11}P_4-s_{14}P_1 \right\rangle = 0$ 

# We define three points $P_1$, $P_2$, and $P_4$.

# In[27]:


P1 = vector([A1, B1, C1])
P2 = vector([A2, B2, C2])
P4 = vector([A4, B4, C4])


# In[28]:


S11 = scalar_product(P1, P1)
S12 = scalar_product(P1, P2)
S14 = scalar_product(P1, P4)


# In[29]:


assert(delta1(P1, P2, P4) == scalar_product(P4, S11*P2 - S12*P1))
assert(delta1(P1, P2, P4) == scalar_product(P2, S11*P4 - S14*P1))


# ## Proof of $\bar{\delta}_1(P_1, P_2, P_3) = 0$ iff ($\left\langle P_1, P_1 \right\rangle$ and $P_1 \vee P_2 \vee P_3$ is tangent to the isotropic conic in $P_1$) or $P_3 = (s_{12}^2 + s_{11} s_{22}) P_1 - 2 s_{11} s_{12}P_2$ or $P_2 = (s_{13}^2 + s_{11} s_{33}) P_1 - 2 s_{11} s_{13}P_3$

# We define three aligned points $P_1$, $P_2$, and $P_3$.

# In[30]:


P1 = vector([A1, B1, C1])
P2 = vector([A2, B2, C2])
P3 = u1*P1 + u2*P2


# In[34]:


S11 = scalar_product(P1, P1)
S12 = scalar_product(P1, P2)
S13 = scalar_product(P1, P3)
S22 = scalar_product(P2, P2)
S33 = scalar_product(P3, P3)


# ### We prove that if $s_{12}^2 + s_{11} s_{22} \neq 0$ and $s_{11} s_{12} \neq 0$, then $\bar{\delta}_1(P_1, P_2, P_3) = 0$ iff $P_3 = (s_{12}^2 + s_{11} s_{22}) P_1 - 2 s_{11} s_{12}P_2$.

# In[32]:


m2 = matrix([(S12^2 + S11*S22)*P1 - 2*S11*S12*P2, P3]).minors(2)
J2 = S.ideal(m2).saturation(
    S.ideal(matrix([P1, P2]).minors(2))
)[0].saturation(
    S.ideal(matrix([P1, P3]).minors(2))
)[0].saturation(
    S.ideal(matrix([P2, P3]).minors(2))
)[0].saturation(
    S.ideal([(S12^2 + S11*S22), 2*S11*S12])
)[0].saturation(
    S.ideal(u1, u2)
)[0]


# In[33]:


assert(J2 == S.ideal(delta1b(P1, P2, P3)))


# ### We prove that if $s_{13}^2 + s_{11} s_{33} \neq 0$ and $s_{11} s_{13} \neq 0$, then $\bar{\delta}_1(P_1, P_2, P_3) = 0$ iff $P_2 = (s_{13}^2 + s_{11} s_{33}) P_1 - 2 s_{11} s_{13}P_3$.

# In[35]:


m2 = matrix([(S13^2 + S11*S33)*P1 - 2*S11*S13*P3, P2]).minors(2)
J2 = S.ideal(m2).saturation(
    S.ideal(matrix([P1, P2]).minors(2))
)[0].saturation(
    S.ideal(matrix([P1, P3]).minors(2))
)[0].saturation(
    S.ideal(matrix([P2, P3]).minors(2))
)[0].saturation(
    S.ideal([(S13^2 + S11*S33), 2*S11*S13])
)[0].saturation(
    S.ideal(u1, u2)
)[0]


# In[36]:


assert(J2 == S.ideal(delta1b(P1, P2, P3)))


# ### We prove that if $s_{12}^2 + s_{11} s_{22} = s_{11} s_{12} = s_{13}^2 + s_{11} s_{33} = s_{11} s_{13} = 0$, then $\bar{\delta}_1(P_1, P_2, P_3) = 0$ iff $\left\langle P_1, P_1 \right\rangle$ and $P_1 \vee P_2 \vee P_3$ is tangent to the isotropic conic in $P_1$

# In[88]:


I = S.ideal(delta1b(P1, P2, P3), S12^2 + S11*S22, S11*S12, S13^2 + S11*S33)
I = I.saturation(
    S.ideal(matrix([P1, P2]).minors(2))
)[0].saturation(
    S.ideal(matrix([P1, P3]).minors(2))
)[0].saturation(
    S.ideal(matrix([P2, P3]).minors(2))
)[0].saturation(
    S.ideal(u1, u2)
)[0].radical()


# In[89]:


J = S.ideal(scalar_product(P1, P1), sigma(P1, P2)).radical()


# In[90]:


assert(I == J)

