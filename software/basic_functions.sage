import itertools


def scalar_product(P, Q):
    '''
    Scalar product between two vectors with three coordinates.

    INPUT:
    - P : vector with three coordinates
    - Q : vector with three coordinates

    OUTPUT: the scalar product of P and Q
    '''

    return P[0]*Q[0] + P[1]*Q[1] + P[2]*Q[2]

def wedge_product(P, Q):
    '''
    The wedge product of two vectors with three coordinates.

    INPUT:
    - P : vector with three coordinates
    - Q : vector with three coordinates

    OUTPUT: the wedge product of P and Q
    '''

    aux = matrix([P, Q]).minors(2)
    return vector((aux[2], -aux[1], aux[0]))

def delta1(P, Q, T):
    '''
    The function delta_1. See Definition {...}

    INPUT:
    - P : vector with three coordinates
    - Q : vector with three coordinates
    - T : vector with three coordinates

    OUTPUT: delta_1(P, Q, T)
    '''

    return scalar_product(P, P)*scalar_product(Q, T) - scalar_product(P, Q)*scalar_product(P, T)

def delta2(P1, P2, P3, P4, P5):
    '''
    The function delta_2. See Definition {...}.

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
    The function sigma. See Definition {...}.

    INPUT:
    - P : vector with three coordinates
    - Q : vector with three coordinates

    OUTPUT: sigma(P, Q)
    '''

    return scalar_product(P, P)*scalar_product(Q, Q) - scalar_product(P, Q)^2

def point_substitution(P):
    '''
    Point substitution.
    '''

    return {x:P[0], y:P[1], z:P[2]}

def phi(P, R):
    '''
    Computes the eigenscheme conditions determined by P.

    INPUT:
    - P : a vector with three coordinates
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
    instance_ideal_eigenscheme = ideal_eigenscheme.subs(point_substitution(P))

    vars_a = [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9]
    return [vector([poly.coefficient(a_variable) for a_variable in vars_a]) for poly in instance_ideal_eigenscheme.gens()]

def condition_matrix(list_points, R, standard="first_two", list_rows=[]):
    '''
    Returns the condition matrix Phi(list_points).

    INPUT:
    - list_points : a list of ...
    ...
    '''

    import operator
    import itertools

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
