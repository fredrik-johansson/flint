/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"

TEST_FUNCTION_START(arb_fmma, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, c, d, z, v;
        slong prec, bits;

        arb_init(a);
        arb_init(b);
        arb_init(c);
        arb_init(d);
        arb_init(z);
        arb_init(v);

        bits = (n_randint(state, 20) == 0) ? 20000 : 1000;
        prec = 2 + n_randint(state, bits);
        arb_randtest_special(a, state, bits, 100);
        arb_randtest_special(b, state, bits, 100);
        arb_randtest_special(c, state, bits, 100);
        arb_randtest_special(d, state, bits, 100);

        if (n_randint(state, 4) == 0)
            arb_neg(c, a), arb_set(d, b);

        /* reference: sum of the products with more precision */
        arb_mul(v, a, b, prec + 100);
        arb_addmul(v, c, d, prec + 100);

        switch (n_randint(state, 3))
        {
            case 0:
                arb_fmma(z, a, b, c, d, prec);
                break;
            case 1:
                arb_set(z, a);
                arb_fmma(z, z, b, c, d, prec);
                break;
            default:
                arb_set(z, d);
                arb_fmma(z, a, b, c, z, prec);
                break;
        }

        if (!arb_overlaps(z, v) || (arb_is_finite(v) && !arb_contains(z, v) &&
            arb_is_exact(a) && arb_is_exact(b) && arb_is_exact(c) && arb_is_exact(d) &&
            arb_rel_accuracy_bits(z) < prec - 4 && arb_is_finite(z)))
        {
            flint_printf("FAIL!\n");
            flint_printf("prec = %wd\n\n", prec);
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_printf("d = "); arb_print(d); flint_printf("\n\n");
            flint_printf("z = "); arb_print(z); flint_printf("\n\n");
            flint_printf("v = "); arb_print(v); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        arb_clear(d);
        arb_clear(z);
        arb_clear(v);
    }

    /* exact midpoints with small exponents: the result contains the exact
       value, its radius is within a factor 1 + 2^-20 of half an ulp plus
       a tiny multiple of the larger product (the deficits of truncated
       products), and squares go through the same code */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, c, d, z;
        arf_t s, t;
        mag_t u, w;
        slong prec, bits, emax;
        int ok;

        arb_init(a);
        arb_init(b);
        arb_init(c);
        arb_init(d);
        arb_init(z);
        arf_init(s);
        arf_init(t);
        mag_init(u);
        mag_init(w);

        bits = 1 + n_randint(state, (n_randint(state, 10) == 0) ? 3000 : 300);
        prec = 2 + n_randint(state, 400);
        arb_randtest(a, state, bits, 10);
        arb_randtest(b, state, bits, 10);
        arb_randtest(c, state, 1 + n_randint(state, bits), 10);
        arb_randtest(d, state, 1 + n_randint(state, bits), 10);
        mag_zero(arb_radref(a));
        mag_zero(arb_radref(b));
        mag_zero(arb_radref(c));
        mag_zero(arb_radref(d));

        switch (n_randint(state, 5))
        {
            case 0: arb_neg(c, a); arb_set(d, b); arf_mul_2exp_si(arb_midref(d), arb_midref(d), n_randint(state, 3)); break;
            case 1: arb_set(b, a); arb_set(d, c); break;   /* squares */
            default: break;
        }

        if (n_randint(state, 2) && arb_equal(a, b) && arb_equal(c, d))
            arb_sosq(z, a, c, prec);
        else if (n_randint(state, 2))
            arb_fmma(z, a, b, c, d, prec);
        else
        {
            arb_set(z, c);
            arb_fmma(z, a, b, z, d, prec);
        }

        arf_mul(s, arb_midref(a), arb_midref(b), ARF_PREC_EXACT, ARF_RND_DOWN);
        arf_mul(t, arb_midref(c), arb_midref(d), ARF_PREC_EXACT, ARF_RND_DOWN);
        emax = FLINT_MAX(arf_is_zero(s) ? WORD_MIN / 2 : ARF_EXP(s), arf_is_zero(t) ? WORD_MIN / 2 : ARF_EXP(t));
        arf_add(s, s, t, ARF_PREC_EXACT, ARF_RND_DOWN);

        ok = arb_contains_arf(z, s);

        if (ok && !arf_is_zero(s))
        {
            /* half an ulp of the midpoint, plus the slack */
            arf_mag_set_ulp(u, arb_midref(z), prec);
            mag_mul_2exp_si(u, u, -1);
            mag_mul_ui(u, u, 1048577);
            mag_mul_2exp_si(u, u, -20);
            mag_one(w);
            mag_mul_2exp_si(w, w, emax - prec - 20);
            mag_add(u, u, w);
            ok = mag_cmp(arb_radref(z), u) <= 0;
        }

        if (!ok)
        {
            flint_printf("FAIL (exact)\n");
            flint_printf("prec = %wd\n\n", prec);
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_printf("d = "); arb_print(d); flint_printf("\n\n");
            flint_printf("z = "); arb_print(z); flint_printf("\n\n");
            flint_printf("s = "); arf_print(s); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        arb_clear(d);
        arb_clear(z);
        arf_clear(s);
        arf_clear(t);
        mag_clear(u);
        mag_clear(w);
    }

    TEST_FUNCTION_END(state);
}
