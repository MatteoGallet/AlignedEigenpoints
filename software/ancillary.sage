###############
# Auxiliary file for the paper
# "Ternary symmetric tensors with an aligned triple of eigenpoints"
# by
# Valentina Beorchia, Matteo Gallet, Alessandro Logar
# Last modified: 25/05/2023
###############

class SymbolicCheck:

    total_checks = 7

    def __init__(self):
        self.passed_checks = 0

    def scalar_product(self, P, Q):
        '''
        Scalar product between two vectors with three coordinates.

        INPUT:
        - P : vector with three coordinates
        - Q : vector with three coordinates

        OUTPUT: the scalar product of P and Q
        '''

        return P[0]*Q[0] + P[1]*Q[1] + P[2]*Q[2]

    def wedge_product(self, P, Q):
        '''
        The wedge product of two vectors with three coordinates.

        INPUT:
        - P : vector with three coordinates
        - Q : vector with three coordinates

        OUTPUT: the wedge product of P and Q
        '''

        aux = matrix([P, Q]).minors(2)
        return vector((aux[2], -aux[1], aux[0]))

    def delta1(self, P, Q, T):
        '''
        The function delta_1. See Definition {...}

        INPUT:
        - P : vector with three coordinates
        - Q : vector with three coordinates
        - T : vector with three coordinates

        OUTPUT: delta_1(P, Q, T)
        '''

        return self.scalar_product(P, P)*self.scalar_product(Q, T) - self.scalar_product(P, Q)*self.scalar_product(P, T)

    def delta2(self, P1, P2, P3, P4, P5):
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

        return self.scalar_product(P1, P2)*self.scalar_product(P1, P3)*self.scalar_product(P4, P5) - self.scalar_product(P1, P4)*self.scalar_product(P1, P5)*self.scalar_product(P2, P3)

    def sigma(self, P, Q):
        '''
        The function sigma. See Definition {...}.

        INPUT:
        - P : vector with three coordinates
        - Q : vector with three coordinates

        OUTPUT: sigma(P, Q)
        '''

        return self.scalar_product(P, P)*self.scalar_product(Q, Q) - self.scalar_product(P, Q)^2

    def point_substitution(self, P):
        '''
        Point substitution.
        '''

        return {x:P[0], y:P[1], z:P[2]}

    def all_checks(self):
        '''
        Performs all the checks.
        '''

        self.passed_checks = 0

        if self.check_scalar_product():
            self.passed_checks += 1

        if self.check_wedge_product():
            self.passed_checks += 1

        if self.check_delta_1():
            self.passed_checks += 1

        if self.check_delta_2():
            self.passed_checks += 1

        if self.check_sigma():
            self.passed_checks += 1

        if self.check_orthogonal_tangent():
            self.passed_checks += 1

        if self.check_conditions_delta():
            self.passed_checks += 1

    def everything_alright(self):
        '''
        Determines whether all checks have passed.
        '''

        return self.passed_checks == self.total_checks

    def check_scalar_product(self):
        '''
        Check for scalar product function.
        '''

        return self.scalar_product([1,2,3], [4,5,6]) == 32

    def check_wedge_product(self):
        '''
        Check for wedge product function.
        '''

        return self.wedge_product([1,2,3], [4,5,6]) == self.wedge_product([1,2,3], [4,5,6])

    def check_delta_1(self):
        '''
        Check for delta_1 function.
        '''

        P = [1,2,3]
        Q = [4,5,6]
        T = [7,8,9]

        return self.delta1(P, Q, T) == self.delta1(P, Q, T)

    def check_delta_2(self):
        '''
        Check for delta_1 function.
        '''

        P1 = [1,2,3]
        P2 = [3,4,5]
        P3 = [5,6,7]
        P4 = [7,8,9]
        P5 = [9,10,11]

        return self.delta2(P1, P2, P3, P4, P5) == self.delta2(P1, P2, P3, P4, P5)

    def check_sigma(self):
        '''
        Check for sigma function.
        '''

        P = [1,2,3]
        Q = [4,5,6]

        return self.sigma(P, Q) == self.sigma(P, Q)

    def check_orthogonal_tangent(self):
        '''
        Proof of Proposition {prop:orthogonal_tangent}.
        '''

        R = PolynomialRing(QQ, "x, y, z, A1, B1, C1, A2, B2, C2")
        R. inject_variables(verbose=False)
        # Definition of the isotropic conic
        isotropic_conic = x^2 + y^2 + z^2
        # Two generic points in the plane
        P = vector(R, [A1, B1 ,C1])
        Q = vector(R, [A2, B2 ,C2])
        # The ideal given by the intersection of the line r through P and Q and the isotropic conic,
        # together with the condititions that P is orthogonal to itself and to Q.
        J = R.ideal(det(matrix([P, Q, (x, y, z)])), self.scalar_product(P, P), self.scalar_product(P, Q), isotropic_conic)
        # We saturate w.r.t. the condition that P and Q are distinct:
        J = J.saturation(R.ideal(matrix([P, Q]).minors(2)))[0]
        # We saturate w.r.t. the condition that P1 and (x, y, z) are distinct:
        J = J.saturation(R.ideal(matrix([P, (x, y, z)]).minors(2)))[0]
        # If what is left after the saturations is the the whole ring,
        # then r is tangent to the isotropic conic at P.
        return J == R.ideal(R.one())

    def check_conditions_delta(self):
        '''
        Proof of ...
        '''

        var_xyz = ["x", "y", "z", "u1", "u2", "v1", "v2"]
        var_a = ["a"+str(i) for i in range(10)]
        var_ABC = [str(XX)+str(i) for i in range(1, 5) for XX in ["A", "B", "C"]]

        K = QQ
        R = PolynomialRing(K, var_xyz + var_a + var_ABC)
        R.inject_variables(verbose=false)

        # Generic plane cubic
        F = a0*x^3 + a1*x^2*y + a2*x*y^2 + a3*y^3 + a4*x^2*z + a5*x*y*z + a6*y^2*z + a7*x*z^2 + a8*y*z^2 + a9*z^3

        # Three generic points in IP^2
        P1 = vector((A1, B1, C1))
        P2 = vector((A2, B2, C2))
        P4 = vector((A4, B4, C4))

        # We define P3 to be collinear with P1, P2 and P4 to be collinear with P1, P4")
        P3 = u1*P1 + u2*P2
        P5 = v1*P1 + v2*P4

        # Definition of the eigenscheme
        gradients = matrix([[x, y, z], [F.derivative(x), F.derivative(y), F.derivative(z)]])
        eigenscheme = gradients.minors(2)

        # Conditions on the cubics to have P1, P2, P3, P4, and P5 as eigenpoints
        conditionsEigenpoints = [eq.subs(self.point_substitution(pt)) for pt in [P1, P2, P3, P4, P5] for eq in eigenscheme]

        # conditionsEigenpoints are linear in the coefficients of the cubics,
        # so we encode them in a 15x10 matrix

        vv = [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9]
        M = matrix([
                    [condition.coefficient(cf) for cf in vv]
                    for condition in conditionsEigenpoints
                ])

        assert(M.nrows() == 15)
        assert(M.ncols() == 10)

        # Lemmas ... from the paper ensure that if we compute the determinant of a single 10x10 submatrix of M
        # obtained by picking exactly two out of the three rows determined by each of the five points,
        # then we know how any other determinant of this kind looks like

        # This computation takes two minutes
        m1 = M.matrix_from_rows_and_columns([0,1,3,4,6,7,9,10,12,13], range(10)).det()

        # The considerations from the paper then ensure that we can remove from m1 the factors
        # A4 * A2 * A1 * v2^2 * v1^2 * u2^2 * u1^2 * (v1*A1 + v2*A4) * (u1*A1 + u2*A2)
        # Morever, since we suppose that P1, P2, and P4 to be non-collinear,
        # we can factor from m1 also the determinant of the matrix [ P1, P2, P4 ]

        assert(m1.quo_rem(A4 * A2 * A1 * v2^2 * v1^2 * u2^2 * u1^2 * (v1*A1 + v2*A4) * (u1*A1 + u2*A2))[1].is_zero())
        m1 = m1.quo_rem(A4 * A2 * A1 * v2^2 * v1^2 * u2^2 * u1^2 * (v1*A1 + v2*A4) * (u1*A1 + u2*A2))[0]
        assert(m1.quo_rem(matrix([P1, P2, P4]).det()^5)[1].is_zero())
        m1 = m1.quo_rem(matrix([P1, P2, P4]).det()^5)[0]

        return m1 == -54*self.delta1(P1, P2, P4)*self.delta2(P1, P2, P3, P4, P5)

# Running the tests
my_test = SymbolicCheck()
my_test.all_checks()
print(my_test.everything_alright())
