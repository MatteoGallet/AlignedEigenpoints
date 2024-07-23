#!/usr/bin/env python
# coding: utf-8

# # Lemma

# Let $r = ax+by+cz$ be a line of the plane 
# and suppose it intersects $\mathcal{Q}_{\mathrm{iso}}$ in two
# distinct points $P_1$ and $P_2$. Consider the cubic:
# $$
#   C(r) = \bigl( r^2-3\left(a^2+b^2+c^2\right) \mathcal{Q}_{\mathrm{iso}} \bigr) \, r \,
# $$
# then, the two tangent lines to $\mathcal{Q}_{\mathrm{iso}}$ in $P_1$ and $P_2$
# are contained in the eigenscheme of~$C(r)$.

# About 10'' of computations

# In[1]:


load("basic_functions.sage")


# Construction of a generic point on the isotropic conic (which will be discovered of 
# the form $(l_1^2+l_2^2, i(l_1^2-l_2^2), 2il_1l_2)$:

# In[2]:


P1 = vector(S, (1, ii, 0))


# In[3]:


rt1 = l1*(y-ii*x)+l2*z
rt1.subs(substitution(P1))

scndP = S.ideal(Ciso, rt1).radical().primary_decomposition()[1]
aux = scndP.gens()[:2]
mm2 = matrix(
    [
        [aux[0].coefficient(x), aux[0].coefficient(y), aux[0].coefficient(z)],
        [aux[1].coefficient(x), aux[1].coefficient(y), aux[1].coefficient(z)]
    ]
).minors(2)

## Generic point on the isotropic conic (depending on two parameters):
PP = vector(S, (mm2[2], -mm2[1], mm2[0]))


# In[4]:


assert(scndP.subs(substitution(PP)) == S.ideal(S(0)))
assert(Ciso.subs(substitution(PP)) == S(0))
assert(PP*ii == vector(S, (l1^2 + l2^2, ii*l1^2 + (-ii)*l2^2, (2*ii)*l1*l2)))


# We can always assume that $l_2 \neq 0$, since $l_2 = 0$ gives that $PP = P_1$:

# In[5]:


assert(matrix([P1, PP.subs(l2=0)]).rank() == 1)


# Now that we know the generic point of $\mathcal{Q}_{\mathrm{iso}}$, we define two (distinct) points on 
# the isotropic conic $\mathcal{Q}_{\mathrm{iso}}$:

# In[6]:


PP1 = vector(S, ((-ii)*l1^2 + (-ii)*l2^2, l1^2 - l2^2, 2*l1*l2))
PP2 = vector(S, ((-ii)*m1^2 + (-ii)*m2^2, m1^2 - m2^2, 2*m1*m2))


# And we defines the lines ttg1 and ttg2, the first is tangent to $\mathcal{Q}_{\mathrm{iso}}$ in PP1, 
# the second is tangent ot $\mathcal{Q}_{\mathrm{iso}}$ in PP2.

# Tangent to isotropic conic in PP1:

# In[7]:


subP1 = substitution(PP1)
ttg1 = (
    (derivative(Ciso, x).subs(subP1))*(PP1[0]-x) 
    + (derivative(Ciso, y).subs(subP1))*(PP1[1]-y)
    + (derivative(Ciso, z).subs(subP1))*(PP1[2]-z)
)


# Tangent to isotropic conic in PP2:

# In[8]:


subP2 = substitution(PP2)
ttg2 = (
    (derivative(Ciso, x).subs(subP2))*(PP2[0]-x)
    + (derivative(Ciso, y).subs(subP2))*(PP2[1]-y)
    + (derivative(Ciso, z).subs(subP2))*(PP2[2]-z)
)


# Just to be sure, we verify that ttg1 is tangent:

# In[9]:


assert(
    S.ideal(ttg1, Ciso).radical() == S.ideal(z*l1 + (-ii)*x*l2 - y*l2, x*l1 + ii*y*l1 + ii*z*l2, x^2 + y^2 + z^2)
)


# We define the line $PP1 \vee PP2$:

# In[10]:


rr = matrix([PP1, PP2, (x, y, z)]).det().factor()[-1][0]


# Finally, we define the cubic 
# $$C(r) = r(r^2-3(a^2+b^2+c^2)\mathcal{Q}_{iso})$$
# and we verify that $C(r)$ has ttg1 and ttg2 as eigenpoints (10 seconds of computation)

# In[11]:


Cr = rr*(rr^2-3*(rr.coefficient(x)^2 + rr.coefficient(y)^2 + rr.coefficient(z)^2)*Ciso)

JJ = S.ideal(matrix(
    [
        [Cr.derivative(x), Cr.derivative(y), Cr.derivative(z)],
        [x, y, z]
    ]
).minors(2))

PDJ = JJ.radical().primary_decomposition()


# In[12]:


assert(len(PDJ) == 2)
assert(Set(PDJ) == Set([S.ideal(ttg1), S.ideal(ttg2)]))


# This concludes the proof of the lemma.
