/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Compare fixed_log1p_bitwise_rs for several values of the reduction
   parameter r against arb_log at equal precision. */

#include "profiler.h"
#include "flint.h"
#include "mpn_extras.h"
#include "arb.h"
#include "fixed.h"

#define NUM_N 15
#define NUM_R 5

/* nanoseconds per call, by repetition scaling until >= 10 ms */
#define BENCH(ns, CODE) \
    do { \
        slong reps = 1; \
        timeit_t tm; \
        while (1) \
        { \
            slong ix; \
            timeit_start(tm); \
            for (ix = 0; ix < reps; ix++) { CODE; } \
            timeit_stop(tm); \
            if (tm->wall >= 10) \
                break; \
            reps *= 4; \
        } \
        (ns) = 1e6 * tm->wall / reps; \
    } while (0)

int main(void)
{
    flint_rand_t state;
    slong ns[NUM_N] = {1, 2, 3, 4, 5, 8, 12, 16, 24, 32, 48, 64, 96,
        128, 157};
    int rs[NUM_R] = {32, 64, 128, 192, 256};
    slong a, ri, pass;
    double tb[NUM_N][NUM_R], ta[NUM_N];

    flint_rand_init(state);

    for (a = 0; a < NUM_N; a++)
    {
        ta[a] = 1e30;
        for (ri = 0; ri < NUM_R; ri++)
            tb[a][ri] = 1e30;
    }

    /* warm the log table cache to the largest size */
    {
        ulong x[160], y[160];
        flint_mpn_rrandom(x, state, 157);
        fixed_log1p_bitwise_rs(y, x, 157, 256);
    }

    for (pass = 0; pass < 3; pass++)
    for (a = 0; a < NUM_N; a++)
    {
        slong n = ns[a], j;
        ulong x[160], y[160];
        arb_t xa, e;
        fmpz_t f;
        double t;

        for (j = 0; j < n; j++)
            x[j] = n_randlimb(state);

        for (ri = 0; ri < NUM_R; ri++)
        {
            BENCH(t, fixed_log1p_bitwise_rs(y, x, n, rs[ri]));
            tb[a][ri] = FLINT_MIN(tb[a][ri], t);
        }

        arb_init(xa);
        arb_init(e);
        fmpz_init(f);
        fmpz_set_ui_array(f, x, n);
        arb_set_fmpz(xa, f);
        arb_mul_2exp_si(xa, xa, -FLINT_BITS * n);
        arb_add_ui(xa, xa, 1, FLINT_BITS * n + 64);
        BENCH(t, arb_log(e, xa, FLINT_BITS * n));
        ta[a] = FLINT_MIN(ta[a], t);
        arb_clear(xa);
        arb_clear(e);
        fmpz_clear(f);
    }

    flint_printf(" bits   n |     r=32     r=64    r=128    r=192"
        "    r=256 |  arb_log  speedup\n");
    for (a = 0; a < NUM_N; a++)
    {
        double best = 1e30;
        slong bi = 0;

        for (ri = 0; ri < NUM_R; ri++)
            if (tb[a][ri] < best)
            {
                best = tb[a][ri];
                bi = ri;
            }

        flint_printf("%5wd %3wd |", FLINT_BITS * ns[a], ns[a]);
        for (ri = 0; ri < NUM_R; ri++)
            flint_printf(" %7.0f%s", tb[a][ri], ri == bi ? "*" : " ");
        flint_printf("| %8.0f   %.2fx\n", ta[a], ta[a] / best);
    }
    flint_printf("(ns per call; * marks the best r; r is clamped to"
        " FLINT_BITS n - 16 internally)\n");

    flint_rand_clear(state);
    return 0;
}
