/*
    Copyright (C) 2020 William Hart
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "nmod_mat.h"
#include "fmpz.h"
#include "fmpz_mat.h"


/*
    Solve A X = B, A being m x n of rank r (determined mod p). We pick the
    r pivot rows and pivot columns of A, solve the nonsingular r x r
    system for the pivot-column entries of X (other entries zero) using
    the given square solver (Dixon lifting or multimodular), and then verify the remaining m - r
    equations exactly. If the verification fails, either the rank mod p
    was too small or the system is inconsistent; we try more primes until
    the product of primes exceeds the determinant bound (which certifies
    that the rank mod p was correct at least once).
*/
typedef int (*square_solver_fn)(fmpz_mat_t, fmpz_t, const fmpz_mat_t, const fmpz_mat_t);

static int
_fmpz_mat_can_solve_generic(fmpz_mat_t X, fmpz_t den,
                        const fmpz_mat_t A, const fmpz_mat_t B, square_solver_fn solve)
{
    ulong p;
    fmpz_t tested;
    nmod_mat_t LU;
    int result = 0, success;
    slong * perm, * pivots, * nonpiv;
    slong i, j, k, col, rank, m, n, cols;
    fmpz_mat_t Arank, Brank, Zrank, Anp, Bnp, T;
    fmpz_t det_bound;

    m = A->r;
    n = A->c;
    cols = B->c;

    if (m == 0 || cols == 0)
    {
        fmpz_mat_zero(X);
        fmpz_one(den);
        return 1;
    }

    if (n == 0)
    {
        fmpz_mat_zero(X);
        fmpz_one(den);
        return fmpz_mat_is_zero(B);
    }

    p = FMPZ_MAT_STARTING_PRIME();

    fmpz_init(det_bound);
    fmpz_init(tested);
    fmpz_one(tested);

    nmod_mat_init(LU, m, n, p);

    perm = flint_malloc(sizeof(slong) * m);
    pivots = flint_malloc(sizeof(slong) * n);
    nonpiv = flint_malloc(sizeof(slong) * m);

    /* Bound on all r x r minors of A: a prime is bad only if it divides
       some r x r minor, so once the product of the tested primes exceeds
       this bound at least one tested prime gave the true rank. */
    fmpz_mat_det_bound_submatrix(det_bound, A);

    while (1)
    {
        p = n_nextprime(p, 0);

        nmod_mat_set_mod(LU, p);
        fmpz_mat_get_nmod_mat(LU, A);

        for (i = 0; i < m; i++)
            perm[i] = i;

        rank = nmod_mat_lu(perm, LU, 0);

        col = 0;
        for (i = 0; i < rank; i++)
        {
            while (nmod_mat_entry(LU, i, col) == 0)
                col++;
            pivots[i] = col;
            col++;
        }

        fmpz_mat_init(Arank, rank, rank);
        fmpz_mat_init(Brank, rank, cols);
        fmpz_mat_init(Zrank, rank, cols);
        fmpz_mat_init(Anp, m - rank, rank);
        fmpz_mat_init(Bnp, m - rank, cols);
        fmpz_mat_init(T, m - rank, cols);

        for (i = 0; i < rank; i++)
        {
            for (k = 0; k < rank; k++)
                fmpz_set(fmpz_mat_entry(Arank, i, k), fmpz_mat_entry(A, perm[i], pivots[k]));
            for (j = 0; j < cols; j++)
                fmpz_set(fmpz_mat_entry(Brank, i, j), fmpz_mat_entry(B, perm[i], j));
        }
        for (i = rank; i < m; i++)
        {
            for (k = 0; k < rank; k++)
                fmpz_set(fmpz_mat_entry(Anp, i - rank, k), fmpz_mat_entry(A, perm[i], pivots[k]));
            for (j = 0; j < cols; j++)
                fmpz_set(fmpz_mat_entry(Bnp, i - rank, j), fmpz_mat_entry(B, perm[i], j));
        }

        success = solve(Zrank, den, Arank, Brank);

        if (success)
        {
            /* the r pivot equations hold by construction; check the rest */
            fmpz_mat_mul(T, Anp, Zrank);
            fmpz_mat_scalar_mul_fmpz(Bnp, Bnp, den);
            result = fmpz_mat_equal(T, Bnp);

            if (result)
            {
                fmpz_mat_zero(X);
                for (i = 0; i < rank; i++)
                    for (j = 0; j < cols; j++)
                        fmpz_swap(fmpz_mat_entry(X, pivots[i], j), fmpz_mat_entry(Zrank, i, j));
            }
        }

        fmpz_mat_clear(Arank);
        fmpz_mat_clear(Brank);
        fmpz_mat_clear(Zrank);
        fmpz_mat_clear(Anp);
        fmpz_mat_clear(Bnp);
        fmpz_mat_clear(T);

        if (result)
            break;
bad_prime:

        /* A prime p is "bad" only if the rank mod p is less than the
           true rank, which requires p to divide some r x r minor; once
           the product of tested primes exceeds the bound on the minors,
           at least one prime gave the true rank, so the system really is
           inconsistent. */
        fmpz_mul_ui(tested, tested, p);
        if (fmpz_cmp(tested, det_bound) > 0)
            break;
    }

    if (!result)
        fmpz_one(den);

    fmpz_clear(det_bound);
    nmod_mat_clear(LU);
    fmpz_clear(tested);
    flint_free(perm);
    flint_free(pivots);
    flint_free(nonpiv);

    return result;
}

int
fmpz_mat_can_solve_dixon_den(fmpz_mat_t X, fmpz_t den,
                        const fmpz_mat_t A, const fmpz_mat_t B)
{
    return _fmpz_mat_can_solve_generic(X, den, A, B, fmpz_mat_solve_dixon_den);
}

int
fmpz_mat_can_solve_multi_mod_den(fmpz_mat_t X, fmpz_t den,
                        const fmpz_mat_t A, const fmpz_mat_t B)
{
    return _fmpz_mat_can_solve_generic(X, den, A, B, fmpz_mat_solve_multi_mod_den);
}

/* the square subsystem is solved with fmpz_mat_solve, which picks the
   solver by shape */
int
_fmpz_mat_can_solve_auto_den(fmpz_mat_t X, fmpz_t den,
                        const fmpz_mat_t A, const fmpz_mat_t B)
{
    return _fmpz_mat_can_solve_generic(X, den, A, B, fmpz_mat_solve);
}
