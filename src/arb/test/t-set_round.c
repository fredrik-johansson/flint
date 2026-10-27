/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "arb.h"

TEST_FUNCTION_START(arb_set_round, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x, y, z;
        arf_t lo, hi, r;
        slong prec, bits;
        int ok;

        arb_init(x);
        arb_init(y);
        arb_init(z);
        arf_init(lo);
        arf_init(hi);
        arf_init(r);

        bits = (n_randint(state, 10) == 0) ? 20000 : 1000;
        prec = 2 + n_randint(state, bits);
        arb_randtest_special(x, state, bits, 12);
        if (n_randint(state, 3) == 0)
            mag_zero(arb_radref(x));

        if (n_randint(state, 2))
        {
            arb_set_round(y, x, prec);
        }
        else
        {
            arb_set(y, x);
            arb_set_round(y, y, prec);
        }

        /* the reference: round with an ulp of error */
        arb_set(z, x);
        if (arf_set_round(arb_midref(z), arb_midref(x), prec, ARF_RND_DOWN))
            arf_mag_add_ulp(arb_radref(z), arb_radref(x), arb_midref(z), prec);

        ok = arb_overlaps(y, z);

        /* an infinite radius (or NaN midpoint) stays */
        if (mag_is_inf(arb_radref(x)))
            ok = ok && mag_is_inf(arb_radref(y));
        if (arf_is_nan(arb_midref(x)))
            ok = ok && arf_is_nan(arb_midref(y));

        if (arb_is_finite(x))
        {
            /* the endpoints of x lie in y */
            arf_set_mag(r, arb_radref(x));
            arf_add(hi, arb_midref(x), r, ARF_PREC_EXACT, ARF_RND_DOWN);
            arf_sub(lo, arb_midref(x), r, ARF_PREC_EXACT, ARF_RND_DOWN);
            ok = ok && arb_contains_arf(y, lo) && arb_contains_arf(y, hi);

            /* and y is not wider than the reference (up to mag rounding) */
            if (ok)
            {
                mag_t m;
                mag_init(m);
                mag_mul_ui(m, arb_radref(z), 260);
                mag_mul_2exp_si(m, m, -8);
                ok = mag_cmp(arb_radref(y), m) <= 0;
                mag_clear(m);
            }
        }

        /* for an exact LAGOM x: the midpoint is x rounded to nearest and
           the radius bounds |y - x| within a factor 1 + 2^-25 (plus the
           mag rounding) */
        if (ok && arb_is_exact(x) && arf_is_finite(arb_midref(x)) && !arf_is_zero(arb_midref(x)))
        {
            arf_t t;
            mag_t m;
            arf_init(t);
            mag_init(m);
            arf_set_round(t, arb_midref(x), prec, ARF_RND_NEAR);
            ok = arf_equal(t, arb_midref(y));
            arf_sub(t, arb_midref(y), arb_midref(x), ARF_PREC_EXACT, ARF_RND_DOWN);
            arf_get_mag(m, t);
            ok = ok && arf_cmpabs_mag(t, arb_radref(y)) <= 0;
            if (ok && !arf_is_zero(t))
            {
                /* not much larger than the error itself, or than 2^-30
                   ulp when tiny */
                mag_t u;
                mag_init(u);
                mag_mul_ui(u, m, 33);
                mag_mul_2exp_si(u, u, -5);
                arf_mag_set_ulp(m, arb_midref(y), prec + 29);
                mag_add(u, u, m);
                ok = mag_cmp(arb_radref(y), u) <= 0;
                mag_clear(u);
            }
            arf_clear(t);
            mag_clear(m);
        }

        /* arb_neg_round is the negation of arb_set_round */
        if (ok)
        {
            arb_t w;
            arb_init(w);
            if (n_randint(state, 2))
                arb_neg_round(w, x, prec);
            else
            {
                arb_set(w, x);
                arb_neg_round(w, w, prec);
            }
            arb_neg(w, w);
            ok = arb_equal(w, y) || (arf_is_nan(arb_midref(w)) && arf_is_nan(arb_midref(y)));
            arb_clear(w);
        }

        if (!ok)
        {
            flint_printf("FAIL!\n");
            flint_printf("prec = %wd\n\n", prec);
            flint_printf("x = "); arb_print(x); flint_printf("\n\n");
            flint_printf("y = "); arb_print(y); flint_printf("\n\n");
            flint_printf("z = "); arb_print(z); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(x);
        arb_clear(y);
        arb_clear(z);
        arf_clear(lo);
        arf_clear(hi);
        arf_clear(r);
    }

    /* arb_set_round_fmpz and arb_set_round_fmpz_2exp against the exact
       value */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t y;
        arf_t t, u;
        fmpz_t x, e;
        slong prec;
        int ok;

        arb_init(y);
        arf_init(t);
        arf_init(u);
        fmpz_init(x);
        fmpz_init(e);

        prec = 2 + n_randint(state, 300);
        fmpz_randtest(x, state, 1 + n_randint(state, 600));
        fmpz_randtest(e, state, n_randint(state, 3) == 0 ? 80 : 10);

        if (n_randint(state, 2))
        {
            arb_set_round_fmpz(y, x, prec);
            arf_set_fmpz(t, x);
        }
        else
        {
            arb_set_round_fmpz_2exp(y, x, e, prec);
            arf_set_fmpz_2exp(t, x, e);
        }

        ok = arb_contains_arf(y, t);
        arf_set_round(u, t, prec, ARF_RND_NEAR);
        /* rounded to nearest on the fast path (small exponents) */
        if (fmpz_bits(e) < FLINT_BITS - 8)
            ok = ok && arf_equal(u, arb_midref(y));

        if (!ok)
        {
            flint_printf("FAIL (fmpz)\n");
            flint_printf("prec = %wd\nx = ", prec); fmpz_print(x);
            flint_printf("\ne = "); fmpz_print(e);
            flint_printf("\ny = "); arb_printd(y, 30); flint_printf("\n");
            flint_abort();
        }

        arb_clear(y);
        arf_clear(t);
        arf_clear(u);
        fmpz_clear(x);
        fmpz_clear(e);
    }

    TEST_FUNCTION_END(state);
}
