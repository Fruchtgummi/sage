import matrix_generic_dense

from sage.misc.cachefunc import cached_method

class Matrix_local_generic_dense(matrix_generic_dense.Matrix_generic_dense):
    def echelonize(self):
        raise NotImplementedError("echelonization not supported for matrices over inexact local rings")

    def _solve_right_nonsingular_square(self, B, check_rank=True):
        if not self.is_square():
            raise ValueError

        D,_,pivot_cols = self.augment(B)._complete_pivoting_echelonization(protected_columns=1)
        assert len(pivot_cols) == self.nrows()
        return D.matrix_from_rows_and_columns(pivot_cols, range(self.ncols(),D.ncols()))

    @cached_method
    def _complete_pivoting_echelonization(self, protected_columns=0):
        """
        Return this matrix in reduced row echelon form using full pivoting,
        i.e., swapping rows and columns.

        INPUT:

        - ``protected_columns`` -- an integer (default: ``0``), the last
          ``protected_columns`` may not be used as pivot columns; used in
          :meth:`solve_right`.

        OUTPUT:

        Returns a triple ``A, pivot_rows, pivot_cols``, where ``A`` is the
        matrix in reduced row echelon form, ``pivot_rows`` and ``pivot_cols``
        are the rows and columns used as pivots in this process, respectively.

        ALGORITHM:

        TODO: talk about losing precision when eliminating rows with a pivot of
        non-minimal valuation in its row.

        EXAMPLES::

            sage: k = Qp(2,3)
            sage: A = matrix(k, [[2,4],[8,1]]); A
            [  2 + O(2^4) 2^2 + O(2^5)]
            [2^3 + O(2^6)   1 + O(2^3)]
            sage: A._complete_pivoting_echelonization()
            (
            [1 + O(2^3)     O(2^6)]
            [    O(2^4) 1 + O(2^3)], [1, 0], [1, 0]
            )
            sage: A._complete_pivoting_echelonization(protected_columns=1)
            (
            [1 + O(2^3) 2 + O(2^4)]
            [    O(2^6) 1 + O(2^3)], [0], [0]
            )
            sage: A = matrix(k, [[2,4,6],[8,1,3]]); A
            [      2 + O(2^4)     2^2 + O(2^5) 2 + 2^2 + O(2^4)]
            [    2^3 + O(2^6)       1 + O(2^3)   1 + 2 + O(2^3)]
            sage: A._complete_pivoting_echelonization()
            (
            [      1 + O(2^3)           O(2^6)   1 + 2 + O(2^3)]
            [          O(2^4)       1 + O(2^3) 1 + 2^2 + O(2^3)], [1, 0], [1, 0]
            )
            sage: A = matrix(k, [[2,4],[6,8],[1,3]]); A
            [      2 + O(2^4)     2^2 + O(2^5)]
            [2 + 2^2 + O(2^4)     2^3 + O(2^6)]
            [      1 + O(2^3)   1 + 2 + O(2^3)]
            sage: A._complete_pivoting_echelonization()
            (
            [1 + O(2^3)     O(2^3)]
            [    O(2^3) 1 + O(2^3)]
            [    O(2^4)     O(2^4)], [2, 1], [0, 1]
            )

        """
        if not self.base_ring().is_field():
            raise NotImplementedError("only implemented over fields")

        from copy import copy
        A = copy(self)

        # used to keep track of the row/column swapping
        cols = range(A.ncols())
        rows = range(A.nrows())

        pivot_rows = []
        pivot_cols = []

        for p in range(A.nrows()):
            # find a pivot
            eligible_entries = [ (A[r,c].valuation(),r,c) for r in range(p,A.nrows()) for c in range(p,A.ncols()-protected_columns) ]
            if not eligible_entries:
                break
            eligible_entries.sort()
            _,r,c = eligible_entries[0]
            if A[r,c].is_zero():
                break

            # swap the pivot to p,p
            cols[c], cols[p] = cols[p], cols[c]
            rows[r], rows[p] = rows[p], rows[r]
            A.swap_rows(r,p)
            A.swap_columns(c,p)

            # remember the pivot
            pivot_rows.append(rows[p])
            pivot_cols.append(cols[p])

            # eliminate entries in column p
            A.rescale_row(p, ~A[p,p])
            for r in range(A.nrows()):
                if r == p: continue
                A.add_multiple_of_row(r,p,-A[r,p])

        return A, pivot_rows, pivot_cols

    def pivots(self):
        return self._complete_pivoting_echelonization()[1]

    def pivot_rows(self):
        return self._complete_pivoting_echelonization()[2]
