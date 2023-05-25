class SymbolicCheck:

    total_checks = 1

    def __init__(self):
        self.passed_checks = 0

    def scalar_product(self, P1, P2):
        return P1[0]*P2[0] + P1[1]*P2[1] + P1[2]*P2[2]

    def all_checks(self):
        '''
        Performs all the checks
        '''

        self.passed_checks = 0

        if self.check_orthogonal_tangent():
            self.passed_checks += 1

    def everything_alright(self):
        return self.passed_checks == self.total_checks

    def check_orthogonal_tangent(self):
        '''
        Proof of Proposition {prop:orthogonal_tangent}
        '''

        R = PolynomialRing(QQ, 'x, y, z, A1, B1, C1, A2, B2, C2')
        R. inject_variables(verbose=False)
        isotropic_conic = x^2 + y^2 + z^2
        P = vector(R, [A1, B1 ,C1])
        Q = vector(R, [A2, B2 ,C2])
        # The ideal given by the intersection of the line r through P and Q and the isotropic conic.
        J = R.ideal(det(matrix([P, Q, (x, y, z)])), self.scalar_product(P, P), self.scalar_product(P, Q), isotropic_conic)
        # We saturate w.r.t. the condition that P and Q are distinct:
        J = J.saturation(R.ideal(matrix([P, Q]).minors(2)))[0]
        # We saturate w.r.t. the condition that P1 and (x, y, z) are distinct:
        J = J.saturation(R.ideal(matrix([P, (x, y, z)]).minors(2)))[0]
        return J == R.ideal(R.one())
