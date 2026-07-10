/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Tune the default reduction parameters of fixed_exp_bitwise_rs and
   fixed_log1p_bitwise_rs (used when they are called with r = 0).

   For each function, the candidate reduction parameters form a
   ladder (16 for log only, then 32, 64, 128, 192, ...: every
   multiple of FLINT_BITS up to rmax), and for each consecutive pair
   the program finds the smallest n at which the larger parameter
   stops losing, by binary search under the assumption that the
   optimal r is nondecreasing in n (which holds to within noise: the
   margins near the crossovers are flat).  A candidate only counts
   once it is genuinely available, i.e. not clamped by
   r <= FLINT_BITS n - 16.

   Usage: tune-bitwise-r [nmax] [rmax]     (defaults 160, 320;
   nmax at most 32767 so that the crossovers fit the short tables)

   The output is a pair of tables per function in a form that can be
   pasted directly into src/fixed/exp_bitwise_rs.c and
   src/fixed/log1p_bitwise_rs.c. */

#include <stdlib.h>
#include "profiler.h"
#include "flint.h"
#include "mpn_extras.h"
#include "fixed.h"

/* nanoseconds per call, minimum of npasses repetition-scaled runs */
static double
time_call(void (*func)(nn_ptr, nn_srcptr, slong, int), nn_ptr y,
          nn_srcptr x, slong n, int r, int npasses)
{
    double best = 1e300;
    int pass;

    for (pass = 0; pass < npasses; pass++)
    {
        slong reps = 1;
        timeit_t tm;

        while (1)
        {
            slong ix;

            timeit_start(tm);
            for (ix = 0; ix < reps; ix++)
                func(y, x, n, r);
            timeit_stop(tm);
            if (tm->wall >= 5)
                break;
            reps *= 4;
        }

        best = FLINT_MIN(best, 1e6 * tm->wall / reps);
    }

    return best;
}

/* the smallest reduction parameter that is effective at size n
   (smaller requests are clamped up internally): r < 32 takes effect
   only in the n <= 4 code of fixed_log1p_bitwise_rs */
static int
r_floor(int r_first, slong n)
{
    if (r_first < 32 && n <= (FLINT_BITS == 64 ? 4 : 1))
        return r_first;
    return 32;
}

/* does r_new run at least about as fast as r_old at size n?
   Deterministic at the degenerate ends: when r_old is clamped up
   (into r_new or beyond) the two calls behave identically and r_new
   wins the tie; when r_new exceeds FLINT_BITS n - 16 it is not
   genuinely available and loses. */
static int
wins(void (*func)(nn_ptr, nn_srcptr, slong, int), nn_ptr y,
     nn_srcptr x, slong n, int r_new, int r_old, int r_first)
{
    double t_new, t_old;

    if ((slong) r_new > FLINT_BITS * n - 16)
        return 0;
    if (r_old < r_floor(r_first, n))
        return 1;

    t_new = time_call(func, y, x, n, r_new, 3);
    t_old = time_call(func, y, x, n, r_old, 3);

    return t_new <= 1.01 * t_old;
}

static void
tune_func(const char * name,
          void (*func)(nn_ptr, nn_srcptr, slong, int),
          int r_first, slong nmax, slong rmax)
{
    slong maxtab = rmax / FLINT_BITS + 3;
    int * r_tab;
    slong * n_tab;
    slong num, j, k;
    nn_ptr x, y;
    flint_rand_t state;

    r_tab = flint_malloc(maxtab * sizeof(int));
    n_tab = flint_malloc(maxtab * sizeof(slong));
    x = flint_malloc(nmax * sizeof(ulong));
    y = flint_malloc((nmax + 1) * sizeof(ulong));

    flint_rand_init(state);
    /* uniform random limbs: structured random (bit runs) would
       degenerate the greedy accept pattern and skew the reduction
       cost */
    for (j = 0; j < nmax; j++)
        x[j] = n_randlimb(state);

    r_tab[0] = r_first;
    n_tab[0] = 1;
    num = 1;

    while (num < maxtab)
    {
        int r_old = r_tab[num - 1];
        int r_new = (r_old < 32) ? 32 :
                    (r_old < FLINT_BITS) ? FLINT_BITS :
                    r_old + FLINT_BITS;
        slong lo, hi;

        if ((slong) r_new > rmax)
            break;

        /* smallest n in (n_tab[num-1], nmax] where r_new wins;
           the optimum is assumed nondecreasing in n */
        if (!wins(func, y, x, nmax, r_new, r_old, r_first))
            break;

        lo = n_tab[num - 1];        /* r_new loses (or ties clamped) */
        hi = nmax;                  /* r_new wins */
        while (hi - lo > 1)
        {
            slong mid = lo + (hi - lo) / 2;

            if (wins(func, y, x, mid, r_new, r_old, r_first))
                hi = mid;
            else
                lo = mid;
        }

        flint_printf("/* %s: r = %d from n = %wd */\n",
            name, r_new, hi);
        fflush(stdout);

        r_tab[num] = r_new;
        n_tab[num] = hi;
        num++;
    }

    flint_printf("\n/* smallest n at which each reduction parameter"
        " becomes optimal\n   (generated by"
        " src/fixed/tune/tune-bitwise-r.c) */\n");
    flint_printf("static const int _%s_r_tab[] = {", name);
    for (j = 0; j < num; j++)
        flint_printf("%s%d", j ? ", " : "", r_tab[j]);
    flint_printf("};\n");
    flint_printf("static const short _%s_n_tab[] = {", name);
    for (j = 0; j < num; j++)
        flint_printf("%s%wd", j ? ", " : "", n_tab[j]);
    flint_printf("};\n\n");

    /* sanity sweep: print the selected r over a few sizes */
    flint_printf("/* selection:");
    for (k = 1; k <= nmax; k = (k < 8) ? k + 1 : k * 3 / 2)
    {
        for (j = 0; j + 1 < num && k >= n_tab[j + 1]; j++)
            ;
        flint_printf(" n=%wd:r=%d", k, r_tab[j]);
    }
    flint_printf(" */\n\n");

    flint_rand_clear(state);
    flint_free(r_tab);
    flint_free(n_tab);
    flint_free(x);
    flint_free(y);
}

int
main(int argc, char * argv[])
{
    slong nmax = (argc > 1) ? atol(argv[1]) : 160;
    slong rmax = (argc > 2) ? atol(argv[2]) : 320;

    if (nmax < 2 || nmax > 32767 || rmax < 32)
    {
        flint_printf("usage: tune-bitwise-r [nmax (2..32767)]"
            " [rmax (>= 32)]\n");
        return 1;
    }

    flint_printf("/* tuning fixed_exp_bitwise_rs and"
        " fixed_log1p_bitwise_rs default r\n   (nmax = %wd,"
        " rmax = %wd) */\n\n", nmax, rmax);

    tune_func("fixed_exp_bitwise_rs", fixed_exp_bitwise_rs, 32,
        nmax, rmax);
    tune_func("fixed_log1p_bitwise_rs", fixed_log1p_bitwise_rs, 16,
        nmax, rmax);

    return 0;
}
