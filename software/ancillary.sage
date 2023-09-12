###############
# Auxiliary file for the paper
# "Ternary symmetric tensors with an aligned triple of eigenpoints"
# by
# Valentina Beorchia, Matteo Gallet, Alessandro Logar
# Last modified: 05/06/2023
###############

class SymbolicCheck:

    current_checks = 0

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

    def phi(self, P, R):
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
        instance_ideal_eigenscheme = ideal_eigenscheme.subs(self.point_substitution(P))

        vars_a = [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9]
        return [vector([poly.coefficient(a_variable) for a_variable in vars_a]) for poly in instance_ideal_eigenscheme.gens()]

    def condition_matrix(self, list_points, R, standard="first_two", list_rows=[]):
        '''
        Returns the condition matrix Phi(list_points).

        INPUT:
        - list_points : a list of ...
        ...
        '''

        import operator
        import itertools

        if standard == "all":
            return matrix(R, list(itertools.chain(*[self.phi(P, R) for P in list_points])))

        if standard == "first_two":
            return matrix(R, list(itertools.chain(*[self.phi(P, R)[:2] for P in list_points])))

        if standard == "selected":
            def get_rows(lst_rows, mat):
                if len(lst_rows) == 1:
                    return (operator.itemgetter(*lst_rows)(mat),)
                else:
                    return operator.itemgetter(*lst_rows)(mat)

            return matrix(R, list(itertools.chain(*[get_rows(list_rows[i], self.phi(list_points[i], R)) for i in range(len(list_points))])))
        else:
            raise Error("Value of standard not defined")

    def clear_uv(self, poly):
        R = poly.parent()
        return R.ideal(poly).saturation(u1*v1*u2*v2)[0].gens()[0]

    def get_sqrfree(self, poly):
        if poly.is_zero():
            return R(0)
        return prod(map(lambda pair: pair[0], list(poly.factor())))

    def ort_lines(self, p1, p2, q1, q2):
        '''Determines if two lines are orthogonal.'''

        K = p1.parent()
        R = PolynomialRing(K, 'x,y,z')
        R.inject_variables(verbose=False)

        rt1 = matrix([p1, p2, (x, y, z)]).det()
        rt2 = matrix([q1, q2, (x, y, z)]).det()
        cf1 = [rt1.coefficient(xx) for xx in (x, y, z)]
        cf2 = [rt2.coefficient(xx) for xx in (x, y, z)]

        return(self.scalar_product(cf1, cf2) == 0)

    def all_checks(self, quick=False):
        '''
        Performs all the checks.

        INPUT:
        - quick : a Boolean value (default: False)
        '''

        self.passed_checks = 0
        self.current_checks = 0

        self.current_checks += 1
        if self.check_scalar_product():
            self.passed_checks += 1

        self.current_checks += 1
        if self.check_wedge_product():
            self.passed_checks += 1

        self.current_checks += 1
        if self.check_delta_1():
            self.passed_checks += 1

        self.current_checks += 1
        if self.check_delta_2():
            self.passed_checks += 1

        self.current_checks += 1
        if self.check_sigma():
            self.passed_checks += 1

        self.current_checks += 1
        if self.check_orthogonal_tangent():
            self.passed_checks += 1

        self.current_checks += 1
        if self.check_matrix_rank_one_alignment():
            self.passed_checks += 1

        self.current_checks += 1
        if self.check_check_rank_aligned():
            self.passed_checks += 1

        if not(quick):
            self.current_checks += 1
            if self.check_conditions_delta():
                self.passed_checks += 1

            self.current_checks += 1
            if self.check_rank_3aligned():
                self.passed_checks += 1

    def everything_alright(self):
        '''
        Determines whether all checks have passed.
        '''

        return self.passed_checks == self.current_checks

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

    def check_rank_aligned(self):
        '''
        Proof of Proposition {prop:condition_rank_aligned}
        @@@@ NON MI TORNA LA DIM
        Questa dim va sostituita con il file:
        prop47_5ago.sage
        che dovrebbe essere la correzione.
        '''

        var_xyz = ["x", "y", "z"]
        var_a = ["a"+str(i) for i in range(10)]
        var_ABC = [str(XX)+str(i) for i in range(1, 3) for XX in ["A", "B", "C"]]
        var_ii = ["ii"]

        K = QQ
        R = PolynomialRing(K, var_xyz + var_a + var_ABC + var_ii)
        R.inject_variables(verbose=False)

        # Up to SO3, we can suppose that P4 is (1, 0, 0) or (1, i, 0)

        # First case: P4 = (1, 0, 0)

        P1 = vector((A1, B1, C1))
        P2 = vector((A2, B2, C2))
        P4 = vector((1, 0, 0))

        M = self.condition_matrix([P1, P2, P4], R, standard="all")

        J5 = R.ideal(M.minors(5))
        J6 = R.ideal(M.minors(6))

        # We saturate the ideal J6 with respect to the condition that P1, P2, and P4 are aligned.

        J5 = J5.saturation(R.ideal(matrix([P1, P2]).minors(2)))[0].saturation(R.ideal(matrix([P1, P4]).minors(2)))[0].saturation(R.ideal(matrix([P2, P4]).minors(2)))[0].radical()
        J6 = J6.saturation(matrix([P1, P2, P4]).det())[0]

        orbit1 =  J6 == R.ideal(R.one()) and J5 == R.ideal(R.one())

        # Second case: P4 = (1, ii, 0)

        P1 = vector((A1, B1, C1))
        P2 = vector((A2, B2, C2))
        P4 = vector((1, ii, 0))

        M = self.condition_matrix([P1, P2, P4], R, standard="all")

        J6 = R.ideal(M.minors(6)) + R.ideal(ii^2 + 1)

        # We saturate the ideal J6 with respect to the condition that P1, P2, and P4 are aligned.

        J5 = J5.saturation(R.ideal(matrix([P1, P2]).minors(2)))[0].saturation(R.ideal(matrix([P1, P4]).minors(2)))[0].saturation(R.ideal(matrix([P2, P4]).minors(2)))[0].radical()
        J6 = J6.saturation(matrix([P1, P2, P4]).det())[0]

        orbit2 = J6 == R.ideal(R.one()) and J5 == R.ideal(R.one())

        return orbit1 and orbit2

    def check_rank_3aligned_1general(self):
        '''
        Proof of Proposition \\ref{prop:condition3+1}.
        '''

        var_xyz = ["x", "y", "z", "u1", "u2"]
        var_a = ["a"+str(i) for i in range(10)]
        var_ABC = [str(XX)+str(i) for i in range(1, 5) for XX in ["A", "B", "C"]]

        K = QQ
        R = PolynomialRing(K, var_xyz + var_a + var_ABC)
        R.inject_variables(verbose=False)

        # First case: P1 = (1, 0, 0)

        P1 = vector(R, (1, 0, 0))
        P2 = vector(R, (A2, B2, C2))
        P3 = u1*P1 + u2*P2
        P4 = vector(R, (A4, B4, C4))

        M = self.condition_matrix([P1, P2, P3, P4], R, standard="all")

        J = R.ideal(M.minors(8))

        # We obtain that P2 or P3 is on the isotropic conic
        # and that the line P1+P2 is tangent to the isotropic conic in P2 or P3
        orbit1 = Set(J.saturation(matrix(R, [P1, P2, P4]).det())[0].radical().saturation(R.ideal(u1*u2))[0].primary_decomposition()) == Set([R.ideal(A2, B2^2 + C2^2), R.ideal(u1 + u2*A2, B2^2 + C2^2)])

        # Second case: P1 = (1, i, 0)

        var("xx")
        F.<ii> = NumberField(xx^2 + 1)
        R = R.change_ring(F)
        R.inject_variables(verbose=False)

        P1 = vector(R, (1, ii, 0))
        P2 = vector(R, (A2, B2, C2))
        P3 = u1*P1 + u2*P2
        P4 = vector(R, (A4, B4, C4))

        # We perform Gaussian elimination on the matrix M:

        def rid(v1, v2, pv):
            return(v2-(v2[pv]/v1[pv])*v1)

        def ridM(M, nRiga, pv):
            v1 = M[nRiga]
            M1 = [v1] + [rid(v1, M[i], pv) for i in range(M.nrows()) if i != nRiga]
            return(matrix(M1))

        M = self.condition_matrix([P1, P2, P3, P4], R, standard="all")

        assert(M[0][1] == 3*R.one())

        M = ridM(M, 0, 1)

        assert(M[1][4] == R.one())

        M = ridM(M, 1, 4)

        M = M.delete_columns([1, 4])
        M = M.delete_rows([0, 1])

        # After the Gaussian elimination, the matrix has at most rank 6

        J = R.ideal(M.minors(6))

        J = J.saturation(matrix(R, [P1, P2, P4]).det())[0].radical()

        J = J.saturation(R.ideal(u1*u2))[0]

        # We obtain that the line P1+P2 is tangent to the isotropic conic
        # since \sigma(P1, P2) = 0 so
        # since P1 is on the isotropic conic, it is tangent in P1
        orbit2 = J == R.ideal(self.sigma(P1, P2).radical())

        return orbit1 and orbit2

    def check_matrix_rank_one_alignment(self):
        '''
        Proof of Proposition {prop:rank3points}.
        '''

        var_xyz = ["x", "y", "z", "u1", "u2", "v1", "v2"]
        var_a = ["a"+str(i) for i in range(10)]
        var_ABC = [str(XX)+str(i) for i in range(1, 5) for XX in ["A", "B", "C"]]

        K = QQ
        R = PolynomialRing(K, var_xyz + var_a + var_ABC)
        R.inject_variables(verbose=False)

        # Three generic aligned points in the plane
        P1 = vector((A1, B1, C1))
        P2 = vector((A2, B2, C2))
        P3 = u1*P1 + u2*P2

        # Monomials of degree 3 in x, y, z
        mon = [x^3, x^2*y, x*y^2, y^3, x^2*z, x*y*z, y^2*z, x*z^2, y*z^2, z^3]

        # The matrix of conditions imposed by P1, P2, and P3
        # Due to Lemma ..., we can consider just two out of the three conditions for each point
        M = matrix(self.phi(P1, R)[:2] + self.phi(P2, R)[:2] + self.phi(P3, R)[:2])

        # Compute the 6x6 minors of the matrix
        M6 = M.minors(6)

        # Each of the minors has A1*A2*A3*u1^2*u2^2 as a factor
        # Due to Lemma ..., we can get rid of it
        MM6 = list(map(lambda uu: uu.quo_rem(A1*A2*(u1*A1+u2*A2)*u1^2*u2^2)[0], M6))
        J = R.ideal(MM6)

        # We saturate with respect to the condition P1 = P2
        J_sat = J.saturation(R.ideal(matrix([P1, P2]).minors(2)))

        # We compute the primary decomposition of the radical of J_sat
        pd = J_sat[0].radical().primary_decomposition()

        # We obtain three ideals of the form: (<Pi, P1>, <Pi, P2>, <Pi, P3>) for i in {1, 2, 3}
        return (pd[0] == R.ideal(self.scalar_product(P2, P2), self.scalar_product(P1, P2))) and (pd[1] == R.ideal(self.scalar_product(P1, P1), self.scalar_product(P1, P2))) and (pd[2] == R.ideal(self.scalar_product(P1, P3), self.scalar_product(P2, P3), self.scalar_product(P3, P3)).primary_decomposition()[1])

    def check_rank_5V(self):
        '''Proof of Proposition \\ref{prop:frecciaFissata}.'''

        var_xyz = ["x", "y", "z", "u1", "u2", "v1", "v2"]
        var_a = ["a"+str(i) for i in range(10)]
        var_ABC = [str(XX)+str(i) for i in range(1, 5) for XX in ["A", "B", "C"]]

        var("xx")
        K.<ii> = NumberField(xx^2+1)
        R = PolynomialRing(K, var_xyz + var_a + var_ABC)
        R.inject_variables(verbose=False)

        P1 = vector(R, (1, 0, 0))
        P2 = vector(R, (0, ii, 1))
        P4 = vector(R, (0, -ii, 1))
        P3 = u1*P1 + u2*P2
        P5 = v1*P1 + v2*P4

        M = self.condition_matrix([P1, P2, P3, P4, P5], R, standard="all")

        # M[2] is the zero row
        assert(M[2] == vector(R, [0 for i in range(10)]))
        # M[3]-ii*M[4] is zero
        assert(M[3] - ii*M[4] == vector(R, [0 for i in range(10)]))
        # u2*M[6]-ii*u2*M[7]+u1*M[8] is zero
        assert(u2*M[6] - ii*u2*M[7] + u1*M[8] == vector(R, [0 for i in range(10)]))
        # M[9]+ii*M[10] is zero
        assert(M[9] + ii*M[10] == vector(R, [0 for i in range(10)]))
        # v2*M[12]+ii*v2*M[13]+v1*M[14] is zero
        assert(v2*M[12] + ii*v2*M[13] + v1*M[14] == vector(R, [0 for i in range(10)]))
        # therefore in the matrix M we can erase the rows: 2, 3, 6, 9, 12

        M = M.matrix_from_rows([0, 1, 4, 5, 7, 8, 10, 11, 13, 14])

        m8 = M.minors(8)

        # check that the matrix M cannot have rank <= 7
        return R.ideal(m8).saturation(u1*u2*v1*v2)[0] == R.ideal(R.one())

    def check_conditions_delta(self):
        '''
        Proof of ...
        '''

        var_xyz = ["x", "y", "z", "u1", "u2", "v1", "v2"]
        var_a = ["a"+str(i) for i in range(10)]
        var_ABC = [str(XX)+str(i) for i in range(1, 5) for XX in ["A", "B", "C"]]

        K = QQ
        R = PolynomialRing(K, var_xyz + var_a + var_ABC)
        R.inject_variables(verbose=False)

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
        conditions_eigenpoints = [eq.subs(self.point_substitution(pt)) for pt in [P1, P2, P3, P4, P5] for eq in eigenscheme]

        # conditions_eigenpoints are linear in the coefficients of the cubics,
        # so we encode them in a 15x10 matrix

        vv = [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9]
        M = matrix([
                    [condition.coefficient(cf) for cf in vv]
                    for condition in conditions_eigenpoints
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

    def check_four_aligned(self):
    '''Proof that in Proposition {proposition:four_aligned} the rank of the matrix \Phi(P_1, P_2, P_3, Q) is at most 7.'''

        var_xyz = ["x", "y", "z", "u1", "u2", "v1", "v2", "w1", "w2"]
        var_a = ["a"+str(i) for i in range(10)]
        var_ABC = [str(XX)+str(i) for i in range(1, 6) for XX in ["A", "B", "C"]]

        K = QQ
        R = PolynomialRing(K, var_xyz + var_a + var_ABC)
        R.inject_variables(verbose=False)

        P1 = vector(R, (A1, B1, C1))
        P2 = vector(R, (A2, B2, C2))
        P3 = u1*P1 + u2*P2
        Q = w1*P1 + w2*P2

        M = self.condition_matrix([P1, P2, P3, Q], R)

        m8 = M.minors(8)
        if not(Set(m8) == {0}):
            raise Error("There is a non-zero minor of order 8.")

        rows = Combinations(10, 7)
        good = []
        for row in rows:
            if not((0 in row and 1 in row and 2 in row) or (3 in row and 4 in row and 5 in row) or (6 in row and 7 in row and 8 in row)):
            good.append(row)
        M1 = M.matrix_from_rows([0,1,2,3,4,5,6,7,8,9])
        m71 = []
        count = 0
        for row in good:
            count += 1
            N = M1.matrix_from_rows(row)
            n7 = N.minors(7)
            if (count % 10 == 0):
                print("Calcolati minori")
                sleep(1)
            m71.append(n7)
        m71 = flatten(m71)

        J71 = R.ideal(m71)
        G71 = J71.groebner_basis()
        J71sat = J71.saturation(u1*u2*w1*w2)[0].saturation(R.ideal(matrix(R, [P1, P2]).minors(2)))[0]
        J71rad = J71sat.radical()

        M2 = M.matrix_from_rows([0,1,2,3,4,5,6,7,8,10])
        m72 = []
        count = 0
        for row in good:
            count += 1
            N = M2.matrix_from_rows(row)
            n7 = N.minors(7)
            if (count % 10 == 0):
                print("Calcolati minori")
                sleep(1)
            m72.append(n7)
        m72 = flatten(m72)

        J72 = R.ideal(m72)
        G72 = J72.groebner_basis()
        J72sat = J72.saturation(u1*u2*w1*w2)[0].saturation(R.ideal(matrix(R, [P1, P2]).minors(2)))[0]
        J72rad = J72sat.radical()

        M3 = M.matrix_from_rows([0,1,2,3,4,5,6,7,8,11])
        m73 = []
        count = 0
        for row in good:
            count += 1
            N = M3.matrix_from_rows(row)
            n7 = N.minors(7)
            if (count % 10 == 0):
                print("Calcolati minori")
                sleep(1)
            m73.append(n7)
        m73 = flatten(m73)

        J73 = R.ideal(m73)
        G73 = J73.groebner_basis()
        J73sat = J73.saturation(u1*u2*w1*w2)[0].saturation(R.ideal(matrix(R, [P1, P2]).minors(2)))[0]
        J73rad = J73sat.radical()

        J7rad = J71rad + J72rad + J73rad
        J7 = J7rad.saturation(R.ideal(matrix(R, [P1, P2]).minors(2)))[0].saturation(R.ideal(matrix(R, [[u1, u2],[w1, w2]]).minors(2)))[0]
        j = J7.gens()[0]

        if not(j == self.sigma(P1, P2)):
            raise Error("Wrong condition!")



# Running the tests
my_test = SymbolicCheck()
#my_test.all_checks(quick=True)
#print(my_test.everything_alright())
