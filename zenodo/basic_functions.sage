###############
# Auxiliary file for the paper
# "Eigenpoint collinearities of plane cubics"
# by
# Valentina Beorchia, Matteo Gallet, Alessandro Logar
# Last modified: 23/07/2024
# Code written for SageMath >= 10.1
###############

import itertools
import operator

varAn1 = ["s"+str(i)+str(j) for i in range(1,8) for j in range(i, 8)]
varAn2 = ["x", "y", "z", "u1", "u2", "v1", "v2", "w1", "w2", "l1", "l2", "m1", "m2"]
varAn3 = ["a"+str(i) for i in range(10)]
varAn4 = [str(XX)+str(i) for i in range(1, 8) for XX in ["A", "B", "C"]]
var("xx")
K = QQ
K.<ii> = NumberField(xx^2+1)
S = PolynomialRing(K, varAn1+varAn2+varAn3+varAn4)
S.inject_variables(verbose=false)

# The ten degree 3 monomials in x, y , z
mon = [x^3, x^2*y, x*y^2, y^3, x^2*z, x*y*z, y^2*z, x*z^2, y*z^2, z^3]

# The generic cubic of P^2
F = add([S(varAn3[i])*mon[i] for i in range(10)])

# Definition of the isotropic conic Ciso
Ciso = x^2 + y^2 + z^2


def scalar_product(P, Q):
    '''
    Scalar product between two vectors with three coordinates.

    INPUT:
    - P : vector with three coordinates
    - Q : vector with three coordinates

    OUTPUT:
    - the scalar product of P and Q
    '''

    return P[0]*Q[0] + P[1]*Q[1] + P[2]*Q[2]


def wedge_product(P, Q):
    '''
    The cross product of two vectors with three coordinates.

    INPUT:
    - P : vector with three coordinates
    - Q : vector with three coordinates

    OUTPUT:
    - the cross product of P and Q
    '''

    aux = matrix([P, Q]).minors(2)
    return vector((aux[2], -aux[1], aux[0]))


def delta1(P, Q, T):
    '''
    The function delta_1.

    INPUT:
    - P : vector with three coordinates
    - Q : vector with three coordinates
    - T : vector with three coordinates

    OUTPUT: delta_1(P, Q, T)
    '''

    return scalar_product(P, P)*scalar_product(Q, T) - scalar_product(P, Q)*scalar_product(P, T)


def delta1b(P, Q, T):
    return scalar_product(P, P)*scalar_product(Q, T) + scalar_product(P, Q)*scalar_product(P, T)


def delta2(P1, P2, P3, P4, P5):
    '''
    The function delta_2.

    INPUT:
    - P1 : vector with three coordinates.
    - P2 : vector with three coordinates.
    - P3 : vector with three coordinates.
    - P4 : vector with three coordinates.
    - P5 : vector with three coordinates.

    OUTPUT: delta_2(P1, P2, P3, P4, P5)
    '''

    return scalar_product(P1, P2)*scalar_product(P1, P3)*scalar_product(P4, P5) - scalar_product(P1, P4)*scalar_product(P1, P5)*scalar_product(P2, P3)

def sigma(P, Q):
    '''
    The function sigma.

    INPUT:
    - P : vector with three coordinates
    - Q : vector with three coordinates

    OUTPUT: sigma(P, Q)
    '''

    return scalar_product(P, P)*scalar_product(Q, Q) - scalar_product(P, Q)^2

def substitution(P):
    '''
    Return a substitution of variables for a given point.

    INPUT:
    - P : vector with three coordinates

    OUTPUT:
    - a dictionary defining the substitution
    '''

    return {x:P[0], y:P[1], z:P[2]}

def phi(P, R):
    '''
    Compute the eigenscheme conditions determined by P.

    INPUT:
    - P : vector with three coordinates
    - R : the ring of polynomials

    OUTPUT:
    - a list of three vectors, each of which is a condition determined by P
    '''

    # A generic ternary cubic
    F = a0*x^3 + a1*x^2*y + a2*x*y^2 + a3*y^3 + a4*x^2*z + a5*x*y*z + a6*y^2*z + a7*x*z^2 + a8*y*z^2 + a9*z^3

    # The ideal of the eigenscheme of the cubic
    matrix_eigenscheme = matrix([[x, y, z], [diff(F, x), diff(F, y), diff(F, z)]])
    ideal_eigenscheme = R.ideal(matrix_eigenscheme.minors(2))

    # instantiate the eigenscheme ideal at P
    instance_ideal_eigenscheme = ideal_eigenscheme.subs(substitution(P))

    vars_a = [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9]
    return [vector([poly.coefficient(a_variable) for a_variable in vars_a]) for poly in instance_ideal_eigenscheme.gens()]

def condition_matrix(list_points, R, standard="first_two", list_rows=[]):
    '''
    Return the condition matrix Phi(list_points).

    INPUT:
    - list_points : a list of points
    - R : a ring of polynomials
    - standard : a string
    - list_rows : a list

    OUTPUT:
    - the matrix of conditions determined by list_points
    '''

    if standard == "all":
        return matrix(R, list(itertools.chain(*[phi(P, R) for P in list_points])))

    if standard == "first_two":
        return matrix(R, list(itertools.chain(*[phi(P, R)[:2] for P in list_points])))

    if standard == "selected":
        def get_rows(lst_rows, mat):
            if len(lst_rows) == 1:
                return (operator.itemgetter(*lst_rows)(mat),)
            else:
                return operator.itemgetter(*lst_rows)(mat)

        return matrix(R, list(itertools.chain(*[get_rows(list_rows[i], phi(list_points[i], R)) for i in range(len(list_points))])))
    else:
        raise Error("Value of standard not defined")


def get_sqrfree(poly):
    '''
    Compute the squarefree factorization of a polynomial.

    INPUT:
    - poly : a polynomial

    OUTPUT:
    - the squarefree factorization of poly
    '''

    if poly.is_zero():
        return S(0)
    return prod(map(lambda pair: pair[0], list(poly.factor())))


def poly_saturate(pol, divisor):
    '''
    Given a polynomial "pol" and another polynomial "divisor",
    divide pol by divisor as much as possible.

    INPUT:
    - pol : a polynomial
    - divisor : a polynomial

    OUTPUT:
    - pol1, where pol = pol1 * divisor^r and r is the maximal possible power
    '''

    pol1 = pol
    qr = pol1.quo_rem(divisor)
    while qr[1] == 0:
        pol1 = qr[0]
        qr = pol1.quo_rem(divisor)
    return(pol1)


def cubic_from_matrix(M):
    '''
    Computes the cubic from a matrix of conditions of rank 9.

    INPUT:
    - M : a matrix

    OUTPUT:
    - a cubic determined by the one-dimensional kernel of M
    '''

    M = matrix(K, M)
    if rank(M) != 9:
        raiseError("Wrong rank!")
    basis = M.right_kernel().basis()
    # basis should have a single vector
    if len(basis) != 1:
        raiseError("Wrong basis!")
    F = sum([basis[0][i]*mon[i] for i in range(10)])
    return(F)


def orthogonal_lines(p1, p2, q1, q2):
    '''
    Check whether two lines are orthogonal.

    INPUT:
    - p1, p2, q1, q2 : vectors representing points

    OUTPUT:
    - True if the lines p1,p2 and q1,q2 are orthogonal, False otherwise
    '''

    rt1 = matrix([p1, p2, (x, y, z)]).det()
    rt2 = matrix([q1, q2, (x, y, z)]).det()
    cf1 = [rt1.coefficient(xx) for xx in (x, y, z)]
    cf2 = [rt2.coefficient(xx) for xx in (x, y, z)]
    return(scalar_product(cf1, cf2) == 0)


def clear_uv(poly):
    '''
    Remove all the powers of the variables u1, u2, v1, v2 in a polynomial.

    INPUT:
    - poly : a polynomial

    OUTPUT
    - pol1 such that poly = pol1 * u1^r1 * u2^r2 * v1^s1 * v2^s2 where r1, r2, s1, s2 are maximal
    '''

    return S.ideal(poly).saturation(u1*v1*u2*v2)[0].gens()[0]


def alignments(Lpt):
    '''
    Compute the triples of aligned points among a list of points

    INPUT:
    - Lpt : a list of points

    OUTPUT:
    - the list of triples among Lpt that are aligned
    '''

    n = len(Lpt)
    triplet = list(Combinations(n,3))
    alignments = []
    for tr in triplet:
        threePoints = (Lpt[tr[0]], Lpt[tr[1]], Lpt[tr[2]])
        if det(matrix(threePoints)) == 0:
            alignments.append((tr[0]+1, tr[1]+1, tr[2]+1))
    return(alignments)


def gdn(F):
    '''
    Compute the gradient of a given polynomial.

    INPUT:
    - F : a polynomial in the variables x, y, z

    OUTPUT:
    - the gradient of F
    '''

    return(vector(S, [F.derivative(xx) for xx in [x, y, z]]))


def eig(F):
    '''
    Compute the equations of the eigenscheme of a polynomial.

    INPUT:
    - F : a polynomial in the variables x, y, z

    OUTPUT:
    - the list of equations of the eigenscheme of F
    '''

    gF = gdn(F)
    return(vector(S, matrix([gF, [x, y, z]]).minors(2)))


def intersection_lines(Q1, Q2, Q3, Q4):
    '''
    Compute the intersection of two lines.

    INPUT:
    - Q1 : vector with three coordinates
    - Q2 : vector with three coordinates
    - Q3 : vector with three coordinates
    - Q4 : vector with three coordinates

    OUTPUT:
    - the intersection of the lines Q1,Q2 and Q3,Q4
    '''

    twoLines = [
        matrix([[x, y, z], Q1, Q2]).det(),
        matrix([[x, y, z], Q3, Q4]).det()
        ]
    Msist = matrix(
        [
            [twoLines[0].coefficient(vr) for vr in [x, y, z]],
            [twoLines[1].coefficient(vr) for vr in [x, y, z]]
            ]
        )
    LL = Msist.minors(2)
    return(vector((LL[2], -LL[1], LL[0])))


def small_random():
    '''
    Return a random integer in the interval [1, 10].
    '''
    return(floor(random()*11+1))
