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

            GEPP is basically the numerically stable version of GENP.

            We can find a (n x n) matrix E that satisfies (A + E)x_c = b and usually norm(E) ~ eps * norm(A).
            Thus, x_c (computed x) is basically close to x in their vector space (exact sln of a nearby problem).
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

    def thomas(self, a, d, c, b, inplace=True):
        """ Solve a tridiagonal linear system. """

        n = len(d)

        # forward
        for k in range(n - 1):
            mult = a[k] / d[k]
            d[k + 1] -= a[k] * c[k]
            b[k + 1] -= a[k] * b[k]

        # back
        x = [0.0] * n
        x[-1] = b[-1] / d[-1]
        for k in range(n - 2, -1, -1):
            x[k] = (b[k] - c[k]* x[k + 1]) / d[k]
        
        return x
    
    def diag_dominant(self, A, by="col", tol=0.0):
        """ Check if (n x n) matrix A is diagonally dominant by column. """

        n = len(A)
        for i in range(n):
            diag = A[i][i]
            if by == "row":
                off = sum(abs(A[i][j]) for j in range(n) if j != i)
            elif by == "col":
                off = sum(abs(A[j][i]) for j in range(n) if j != i)
            if not (diag + tol) >= off:
                return False
        
        return True