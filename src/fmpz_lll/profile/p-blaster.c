/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Opus 4.8

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Compares fmpz_lll_blaster against fmpz_lll_wrapper and fmpz_lll on
    random q-ary lattices, which are the lattices BLASter targets.

    algorithm 0 : fmpz_lll_blaster
    algorithm 1 : fmpz_lll_wrapper
    algorithm 2 : fmpz_lll
*/

#include "profiler.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpz_lll.h"
#include "ulong_extras.h"

typedef struct
{
    slong dim;
    flint_bitcnt_t bits;
    int algorithm;
} mat_lll_t;

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

void
sample(void *arg, ulong count)
{
    mat_lll_t *params = (mat_lll_t *) arg;
    slong dim = params->dim;
    flint_bitcnt_t bits = params->bits;
    int algorithm = params->algorithm;
    ulong i;
    fmpz_lll_t fl;
    fmpz_t q;
    fmpz_mat_t A, B;
    FLINT_TEST_INIT(state);

    fmpz_lll_context_init_default(fl);

    fmpz_init(q);
    fmpz_randbits(q, state, bits);
    fmpz_abs(q, q);
    if (fmpz_is_zero(q))
        fmpz_one(q);

    fmpz_mat_init(A, dim, dim);
    _rand_qary_basis(A, dim, q, state);
    fmpz_mat_init(B, dim, dim);

    for (i = 0; i < count; i++)
    {
        /* Reduction is destructive, so restore the input each round and
           only time the reduction itself. */
        fmpz_mat_set(B, A);

        prof_start();

        if (algorithm == 0)
            fmpz_lll_blaster(B, NULL, fl);
        else if (algorithm == 1)
            fmpz_lll_wrapper(B, NULL, fl);
        else
            fmpz_lll(B, NULL, fl);

        prof_stop();
    }

    fmpz_clear(q);
    fmpz_mat_clear(A);
    fmpz_mat_clear(B);
    FLINT_TEST_CLEAR(state);
}

int
main(void)
{
    double min, max;
    mat_lll_t params;
    slong dim;
    int alg;
    const char *names[3] = { "blaster", "wrapper", "fmpz_lll" };

    flint_printf("fmpz_lll_blaster profile\n");
    flint_printf("random q-ary lattices, q of 30 bits\n\n");
    flint_printf("%-6s  %-10s  %12s  %12s\n",
                 "dim", "algorithm", "min (us)", "max (us)");

    params.bits = 30;

    for (dim = 16; dim <= 512; dim *= 2)
    {
        params.dim = dim;

        for (alg = 0; alg < 3; alg++)
        {
            params.algorithm = alg;
            prof_repeat(&min, &max, sample, (void *) &params);

            flint_printf("%-6wd  %-10s  %12.1f  %12.1f\n",
                         dim, names[alg], min, max);
            fflush(stdout);
        }

        flint_printf("\n");
    }

    return 0;
}
