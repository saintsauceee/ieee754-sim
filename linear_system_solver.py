class LinearSystemSolver:
    """ Solve Ax = b via different methods. """

    def solve(self, A, b, method="gepp"):
        if method == "genp":
            return self._genp(A, b)
        elif method == "gepp":
            return self._gepp(A, b)
        else:
            raise ValueError("No a valid method.")

    
    def _genp(self, A, b):
        """ 
            Cost: 2/3 n^3 flops

            L = n x n unit lower triangular matrix (elements are the multipliers)
            U = n x n upper triangular matrix obtained after forward elimination
        """
        pass

    def _gepp(self, A, b):
        """ 
            GENP can lead to unncessary loss of accuracy due to large multipliers (as eps increases).
            GEPP leverages partial pivoting in GENP s.t. |multiplier| <= 1 (thus eps is normalized).
        
            Cost: 2/3 n^3 flops + 1/2 n^2 comparisons

            P = P_n-1 ... P_2 P_1 (P_k is I_n with rows k, q swapped)
            L = n x n unit lower triangular matrix (elements are the permuted multipliers)
            U = n x n upper triangular matrix obtained after forward elimination
        """
        pass

    def _lu_factor(self, A, pivot=True):
        return self._lupp(A) if pivot else self._lunp(A)

    def _lunp(self, A):
        """
            We want to solve Ax = b for an arbitrary b.
            => Ax = b
            => LUx = b (since A = LU)

            Let y = Ux, then solve Ly = b for y, then solve Ux = y for x.
            This way, we can solve Ax = b for an arbitrary b using LUx = b.
        """
        
        pass

    def _lupp(self, A):
        """
            We want to solve PAx = Pb for an arbitrary b.
            => PAx = Pb
            => LUx = b (since PA = LU)

            Let y = Ux, then solve Ly = Pb for y, then solve Ux = y for x.
            This way, we can solve PAx = Pb for an arbitrary b using LUx = b.
        """
        
        pass