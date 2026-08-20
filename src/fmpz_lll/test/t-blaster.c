/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Opus 4.8

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpz_lll.h"

/* Generate a random q-ary basis of dimension n:

       [ q          ]
       [    q       ]
       [       ...  ]
       [ a1 a2 ... 1]

   These are the lattices BLASter primarily targets. */
static void
_rand_qary_basis(fmpz_mat_t B, slong n, const fmpz_t q, flint_rand_t state)
{
    slong i, j;

    fmpz_mat_zero(B);

    for (i = 0; i < n - 1; i++)
        fmpz_set(fmpz_mat_entry(B, i, i), q);

    for (j = 0; j < n - 1; j++)
        fmpz_randm(fmpz_mat_entry(B, n - 1, j), state, q);

    fmpz_one(fmpz_mat_entry(B, n - 1, n - 1));
}

/* Check that the Gram determinant det(B B^T) is unchanged, i.e. that the
   applied transformation really was unimodular. */
static int
_gram_det_equal(const fmpz_mat_t A, const fmpz_mat_t C)
{
    slong n = fmpz_mat_nrows(A);
    slong m = fmpz_mat_ncols(A);
    fmpz_mat_t T, GA, GC;
    fmpz_t dA, dC;
    int eq;

    fmpz_mat_init(T, m, n);
    fmpz_mat_init(GA, n, n);
    fmpz_mat_init(GC, n, n);
    fmpz_init(dA);
    fmpz_init(dC);

    fmpz_mat_transpose(T, A);
    fmpz_mat_mul(GA, A, T);

    fmpz_mat_transpose(T, C);
    fmpz_mat_mul(GC, C, T);

    fmpz_mat_det(dA, GA);
    fmpz_mat_det(dC, GC);

    eq = fmpz_equal(dA, dC);

    fmpz_mat_clear(T);
    fmpz_mat_clear(GA);
    fmpz_mat_clear(GC);
    fmpz_clear(dA);
    fmpz_clear(dC);

    return eq;
}

TEST_FUNCTION_START(fmpz_lll_blaster, state)
{
    slong iter;
    fmpz_lll_t fl;

    fmpz_lll_context_init_default(fl);

    /* ---------------------------------------------------------------
       Trivial dimensions: n = 0 and n = 1 must succeed immediately.
       --------------------------------------------------------------- */
    {
        fmpz_mat_t B;

        fmpz_mat_init(B, 0, 0);
        if (fmpz_lll_blaster(B, NULL, fl) != 1)
        {
            flint_printf("FAIL: n = 0 not handled\n");
            fflush(stdout);
            flint_abort();
        }
        fmpz_mat_clear(B);

        fmpz_mat_init(B, 1, 3);
        fmpz_set_si(fmpz_mat_entry(B, 0, 0), 7);
        fmpz_set_si(fmpz_mat_entry(B, 0, 1), -3);
        fmpz_set_si(fmpz_mat_entry(B, 0, 2), 1);
        if (fmpz_lll_blaster(B, NULL, fl) != 1)
        {
            flint_printf("FAIL: n = 1 not handled\n");
            fflush(stdout);
            flint_abort();
        }
        fmpz_mat_clear(B);
    }

    /* ---------------------------------------------------------------
       Small fixed example: check reducedness and the U transform.
       --------------------------------------------------------------- */
    {
        fmpz_mat_t B, B0, U, UB;

        fmpz_mat_init(B, 2, 3);
        fmpz_mat_init(B0, 2, 3);
        fmpz_mat_init(U, 2, 2);
        fmpz_mat_init(UB, 2, 3);

        fmpz_set_si(fmpz_mat_entry(B, 0, 0), 1);
        fmpz_set_si(fmpz_mat_entry(B, 0, 1), 1);
        fmpz_set_si(fmpz_mat_entry(B, 0, 2), 1);
        fmpz_set_si(fmpz_mat_entry(B, 1, 0), -1);
        fmpz_set_si(fmpz_mat_entry(B, 1, 1), 0);
        fmpz_set_si(fmpz_mat_entry(B, 1, 2), 2);

        fmpz_mat_set(B0, B);
        fmpz_mat_one(U);

        if (fmpz_lll_blaster(B, U, fl) != 1)
        {
            flint_printf("FAIL: small example did not converge\n");
            fflush(stdout);
            flint_abort();
        }

        if (!fmpz_mat_is_reduced(B, fl->delta, fl->eta))
        {
            flint_printf("FAIL: small example not reduced\n");
            fmpz_mat_print_pretty(B);
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_mul(UB, U, B0);
        if (!fmpz_mat_equal(UB, B))
        {
            flint_printf("FAIL: small example U * B0 != B\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(B);
        fmpz_mat_clear(B0);
        fmpz_mat_clear(U);
        fmpz_mat_clear(UB);
    }

    /* ---------------------------------------------------------------
       Random q-ary lattices: reducedness, Gram determinant, and the
       unimodular transform B_out = U * B_in.
       --------------------------------------------------------------- */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        fmpz_mat_t B, B0, U, UB;
        fmpz_t q;
        slong n;
        flint_bitcnt_t bits;
        int res;

        n = 2 + n_randint(state, 24);
        bits = 4 + n_randint(state, 40);

        fmpz_init(q);
        fmpz_randbits(q, state, bits);
        fmpz_abs(q, q);
        if (fmpz_is_zero(q))
            fmpz_one(q);

        fmpz_mat_init(B, n, n);
        fmpz_mat_init(B0, n, n);
        fmpz_mat_init(U, n, n);
        fmpz_mat_init(UB, n, n);

        _rand_qary_basis(B, n, q, state);
        fmpz_mat_set(B0, B);
        fmpz_mat_one(U);

        res = fmpz_lll_blaster(B, U, fl);

        if (res == -1)
        {
            flint_printf("FAIL: numerical failure, n = %wd, bits = %wu\n",
                         n, bits);
            fflush(stdout);
            flint_abort();
        }

        if (!fmpz_mat_is_reduced(B, fl->delta, fl->eta))
        {
            flint_printf("FAIL: output not reduced, n = %wd, bits = %wu\n",
                         n, bits);
            fflush(stdout);
            flint_abort();
        }

        if (!_gram_det_equal(B, B0))
        {
            flint_printf("FAIL: Gram determinant changed, n = %wd, bits = %wu\n",
                         n, bits);
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_mul(UB, U, B0);
        if (!fmpz_mat_equal(UB, B))
        {
            flint_printf("FAIL: U * B0 != B, n = %wd, bits = %wu\n", n, bits);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(q);
        fmpz_mat_clear(B);
        fmpz_mat_clear(B0);
        fmpz_mat_clear(U);
        fmpz_mat_clear(UB);
    }

    /* ---------------------------------------------------------------
       Dimensions spanning the segment boundary (seg_len = 64), so that
       the inter-segment block reduction is actually exercised.
       --------------------------------------------------------------- */
    for (iter = 0; iter < 4 * flint_test_multiplier(); iter++)
    {
        fmpz_mat_t B, B0;
        fmpz_t q;
        slong n;

        n = 60 + n_randint(state, 80);   /* 60 .. 139: 1 to 3 segments */

        fmpz_init(q);
        fmpz_randbits(q, state, 20);
        fmpz_abs(q, q);
        if (fmpz_is_zero(q))
            fmpz_one(q);

        fmpz_mat_init(B, n, n);
        fmpz_mat_init(B0, n, n);

        _rand_qary_basis(B, n, q, state);
        fmpz_mat_set(B0, B);

        if (fmpz_lll_blaster(B, NULL, fl) == -1)
        {
            flint_printf("FAIL: numerical failure, multi-segment n = %wd\n", n);
            fflush(stdout);
            flint_abort();
        }

        if (!fmpz_mat_is_reduced(B, fl->delta, fl->eta))
        {
            flint_printf("FAIL: multi-segment output not reduced, n = %wd\n", n);
            fflush(stdout);
            flint_abort();
        }

        if (!_gram_det_equal(B, B0))
        {
            flint_printf("FAIL: multi-segment determinant changed, n = %wd\n", n);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(q);
        fmpz_mat_clear(B);
        fmpz_mat_clear(B0);
    }

    /* ---------------------------------------------------------------
       Rectangular bases (n rows in Z^m with m > n).
       --------------------------------------------------------------- */
    for (iter = 0; iter < 40 * flint_test_multiplier(); iter++)
    {
        fmpz_mat_t B, B0;
        slong n, m;
        flint_bitcnt_t bits;

        n = 2 + n_randint(state, 12);
        m = n + n_randint(state, 8);
        bits = 4 + n_randint(state, 20);

        fmpz_mat_init(B, n, m);
        fmpz_mat_init(B0, n, m);

        fmpz_mat_randbits(B, state, bits);

        /* Skip rank-deficient samples; BLASter assumes a full-rank basis. */
        {
            fmpz_mat_t T, G;
            fmpz_t d;
            int singular;

            fmpz_mat_init(T, m, n);
            fmpz_mat_init(G, n, n);
            fmpz_init(d);
            fmpz_mat_transpose(T, B);
            fmpz_mat_mul(G, B, T);
            fmpz_mat_det(d, G);
            singular = fmpz_is_zero(d);
            fmpz_mat_clear(T);
            fmpz_mat_clear(G);
            fmpz_clear(d);

            if (singular)
            {
                fmpz_mat_clear(B);
                fmpz_mat_clear(B0);
                continue;
            }
        }

        fmpz_mat_set(B0, B);

        if (fmpz_lll_blaster(B, NULL, fl) == -1)
        {
            flint_printf("FAIL: numerical failure, rectangular %wd x %wd\n", n, m);
            fflush(stdout);
            flint_abort();
        }

        if (!_gram_det_equal(B, B0))
        {
            flint_printf("FAIL: rectangular determinant changed, %wd x %wd\n",
                         n, m);
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(B);
        fmpz_mat_clear(B0);
    }

    TEST_FUNCTION_END(state);
}
