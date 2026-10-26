/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arf.h"

/* Exact result, and the reference rounding for a strict mode */
static int
arf_fmma_naive(arf_t res, const arf_t a, const arf_t b, const arf_t c, const arf_t d, slong prec, arf_rnd_t rnd)
{
    arf_t t, u;
    int inexact;
    arf_init(t);
    arf_init(u);
    arf_mul(t, a, b, ARF_PREC_EXACT, ARF_RND_DOWN);
    arf_mul(u, c, d, ARF_PREC_EXACT, ARF_RND_DOWN);
    arf_add(t, t, u, ARF_PREC_EXACT, ARF_RND_DOWN);
    inexact = arf_set_round(res, t, prec, rnd);
    arf_clear(t);
    arf_clear(u);
    return inexact;
}

/* |z - v| <= tol/100 ulp of v at prec bits, and z == v when inexact == 0 */
static int
arf_relaxed_ok(const arf_t z, const arf_t v, slong prec, ulong tol, int inexact)
{
    arf_t err, t;
    fmpz_t e;
    int ok;

    if (arf_is_special(v) || prec == ARF_PREC_EXACT)
        return arf_equal(z, v);

    arf_init(err);
    arf_init(t);
    fmpz_init(e);

    arf_sub(err, z, v, ARF_PREC_EXACT, ARF_RND_DOWN);
    arf_abs(err, err);
    arf_mul_ui(err, err, 100, ARF_PREC_EXACT, ARF_RND_DOWN);
    fmpz_sub_ui(e, ARF_EXPREF(v), prec);
    arf_set_ui(t, tol);
    arf_mul_2exp_fmpz(t, t, e);
    ok = (arf_cmp(err, t) <= 0) && (inexact || arf_equal(z, v));

    arf_clear(err);
    arf_clear(t);
    fmpz_clear(e);
    return ok;
}

TEST_FUNCTION_START(arf_fmma, state)
{
    slong iter;

    for (iter = 0; iter < 20000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arf_t a, b, c, d, z, v, t;
        slong prec, bits;
        int r1, r2, relaxed, ok;
        arf_rnd_t rnd;

        arf_init(a);
        arf_init(b);
        arf_init(c);
        arf_init(d);
        arf_init(z);
        arf_init(v);
        arf_init(t);

        bits = (n_randint(state, 20) == 0) ? 30000 : 2000;
        arf_randtest_special(a, state, bits, 12);
        arf_randtest_special(b, state, bits, 12);
        arf_randtest_special(c, state, bits, 12);
        arf_randtest_special(d, state, bits, 12);

        if (n_randint(state, 4) == 0)
            arf_set(b, a);
        if (n_randint(state, 4) == 0)
            arf_set(d, c);

        /* near cancellation: c d ~ -a b */
        if (n_randint(state, 4) == 0 && arf_is_finite(a) && arf_is_finite(b))
        {
            arf_mul(t, a, b, ARF_PREC_EXACT, ARF_RND_DOWN);
            arf_neg(c, t);
            arf_one(d);
            if (n_randint(state, 2))
                arf_add_ui(c, c, 1, ARF_PREC_EXACT, ARF_RND_DOWN);
            if (n_randint(state, 2) && !arf_is_special(c))
                fmpz_add_si(ARF_EXPREF(c), ARF_EXPREF(c), n_randint(state, 5) - 2);
        }

        /* far apart */
        if (n_randint(state, 4) == 0 && arf_is_finite(a) && !arf_is_special(c))
            fmpz_add_si(ARF_EXPREF(c), ARF_EXPREF(a), n_randint(state, 3 * bits) - bits);

        prec = 2 + n_randint(state, bits);

        switch (n_randint(state, 7))
        {
            case 0: rnd = ARF_RND_DOWN; break;
            case 1: rnd = ARF_RND_UP; break;
            case 2: rnd = ARF_RND_FLOOR; break;
            case 3: rnd = ARF_RND_CEIL; break;
            case 4: rnd = ARF_RND_NEAR; break;
            case 5: rnd = ARF_RND_FAST; break;
            default: rnd = ARF_RND_ACCURATE; break;
        }
        relaxed = (rnd >= ARF_RND_FAST);

        r2 = arf_fmma_naive(v, a, b, c, d, ARF_PREC_EXACT, ARF_RND_DOWN);
        if (!relaxed)
            r2 = arf_set_round(v, v, prec, rnd);

        switch (n_randint(state, 4))
        {
            case 0:
                r1 = arf_fmma(z, a, b, c, d, prec, rnd);
                break;
            case 1:
                arf_set(z, a);
                r1 = arf_fmma(z, z, b, c, d, prec, rnd);
                break;
            case 2:
                arf_set(z, d);
                r1 = arf_fmma(z, a, b, c, z, prec, rnd);
                break;
            default:
                arf_set(z, b);
                r1 = arf_fmma(z, a, z, c, d, prec, rnd);
                break;
        }

        if (relaxed)
            ok = arf_relaxed_ok(z, v, prec, (rnd == ARF_RND_ACCURATE) ? 51 : 100, r1);
        else
            ok = arf_equal(z, v) && (r1 == r2);

        if (!ok)
        {
            flint_printf("FAIL!\n");
            flint_printf("prec = %wd, rnd = %d, r1 = %d, r2 = %d\n\n", prec, rnd, r1, r2);
            flint_printf("a = "); arf_print(a); flint_printf("\n\n");
            flint_printf("b = "); arf_print(b); flint_printf("\n\n");
            flint_printf("c = "); arf_print(c); flint_printf("\n\n");
            flint_printf("d = "); arf_print(d); flint_printf("\n\n");
            flint_printf("z = "); arf_print(z); flint_printf("\n\n");
            flint_printf("v = "); arf_print(v); flint_printf("\n\n");
            flint_abort();
        }

        arf_clear(a);
        arf_clear(b);
        arf_clear(c);
        arf_clear(d);
        arf_clear(z);
        arf_clear(v);
        arf_clear(t);
    }

    TEST_FUNCTION_END(state);
}
