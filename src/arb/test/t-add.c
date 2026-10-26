/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "arb.h"

/* the endpoints of a ball as arf values */
static void
_corners_add(arf_t lo, arf_t hi, const arb_t x)
{
    arf_t r;
    arf_init(r);
    arf_set_mag(r, arb_radref(x));
    arf_sub(lo, arb_midref(x), r, ARF_PREC_EXACT, ARF_RND_DOWN);
    arf_add(hi, arb_midref(x), r, ARF_PREC_EXACT, ARF_RND_DOWN);
    arf_clear(r);
}

/* random ball, usually of few limbs, with the exponent near e (if e != NULL) */
static void
_randball_add(arb_t x, flint_rand_t state, const fmpz * e)
{
    slong bits = n_randint(state, 8) == 0 ? 1 + n_randint(state, 3000) : 1 + n_randint(state, 200);
    arb_randtest(x, state, bits, 12);
    if (n_randint(state, 3) == 0)
        mag_zero(arb_radref(x));
    if (e != NULL && !arf_is_special(arb_midref(x)) && n_randint(state, 2))
        fmpz_add_si(ARF_EXPREF(arb_midref(x)), e, (slong) n_randint(state, 600) - 300);
}

TEST_FUNCTION_START(arb_add, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, c;
        fmpq_t x, y, z;

        arb_init(a);
        arb_init(b);
        arb_init(c);

        fmpq_init(x);
        fmpq_init(y);
        fmpq_init(z);

        arb_randtest(a, state, 1 + n_randint(state, 200), 10);
        arb_randtest(b, state, 1 + n_randint(state, 200), 10);
        arb_randtest(c, state, 1 + n_randint(state, 200), 10);

        arb_get_rand_fmpq(x, state, a, 1 + n_randint(state, 200));
        arb_get_rand_fmpq(y, state, b, 1 + n_randint(state, 200));

        arb_add(c, a, b, 2 + n_randint(state, 200));
        fmpq_add(z, x, y);

        if (!arb_contains_fmpq(c, z))
        {
            flint_printf("FAIL: containment\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("x = "); fmpq_print(x); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("y = "); fmpq_print(y); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_printf("z = "); fmpq_print(z); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);

        fmpq_clear(x);
        fmpq_clear(y);
        fmpq_clear(z);
    }

    /* aliasing of c and a */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b;
        fmpq_t x, y, z;

        arb_init(a);
        arb_init(b);

        fmpq_init(x);
        fmpq_init(y);
        fmpq_init(z);

        arb_randtest(a, state, 1 + n_randint(state, 200), 10);
        arb_randtest(b, state, 1 + n_randint(state, 200), 10);

        arb_get_rand_fmpq(x, state, a, 1 + n_randint(state, 200));
        arb_get_rand_fmpq(y, state, b, 1 + n_randint(state, 200));

        arb_add(a, a, b, 2 + n_randint(state, 200));
        fmpq_add(z, x, y);

        if (!arb_contains_fmpq(a, z))
        {
            flint_printf("FAIL: aliasing (c, a)\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("x = "); fmpq_print(x); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("y = "); fmpq_print(y); flint_printf("\n\n");
            flint_printf("z = "); fmpq_print(z); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);

        fmpq_clear(x);
        fmpq_clear(y);
        fmpq_clear(z);
    }

    /* aliasing of c and b */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b;
        fmpq_t x, y, z;

        arb_init(a);
        arb_init(b);

        fmpq_init(x);
        fmpq_init(y);
        fmpq_init(z);

        arb_randtest(a, state, 1 + n_randint(state, 200), 10);
        arb_randtest(b, state, 1 + n_randint(state, 200), 10);

        arb_get_rand_fmpq(x, state, a, 1 + n_randint(state, 200));
        arb_get_rand_fmpq(y, state, b, 1 + n_randint(state, 200));

        arb_add(b, a, b, 2 + n_randint(state, 200));
        fmpq_add(z, x, y);

        if (!arb_contains_fmpq(b, z))
        {
            flint_printf("FAIL: aliasing (c, b)\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("x = "); fmpq_print(x); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("y = "); fmpq_print(y); flint_printf("\n\n");
            flint_printf("z = "); fmpq_print(z); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);

        fmpq_clear(x);
        fmpq_clear(y);
        fmpq_clear(z);
    }

    /* one-pass rounding: the exact sums of the endpoints lie in the result,
       which is not much wider than the ball from a separately rounded
       midpoint */
    for (iter = 0; iter < 20000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x, y, r, w;
        arf_t xe[2], ye[2], p;
        slong prec, a, b;
        int ok, sub, inexact;

        arb_init(x); arb_init(y); arb_init(r); arb_init(w);
        for (a = 0; a < 2; a++) { arf_init(xe[a]); arf_init(ye[a]); }
        arf_init(p);

        prec = (n_randint(state, 8) == 0) ? 2 + n_randint(state, 3000) : 2 + n_randint(state, 300);
        if (n_randint(state, 50) == 0)
            prec = ARF_PREC_EXACT;
        sub = n_randint(state, 2);
        _randball_add(x, state, NULL);
        _randball_add(y, state, ARF_EXPREF(arb_midref(x)));
        if (n_randint(state, 4) == 0)
        {
            arb_set(y, x);
            if (n_randint(state, 2))
                arb_neg(y, y);
        }

        if (n_randint(state, 2))
        {
            arb_set(r, x);
            if (sub) arb_sub(r, r, y, prec); else arb_add(r, r, y, prec);
        }
        else
        {
            if (sub) arb_sub(r, x, y, prec); else arb_add(r, x, y, prec);
        }

        /* the reference: midpoint rounded by arf, radius bumped by an ulp */
        if (sub)
            inexact = arf_sub(arb_midref(w), arb_midref(x), arb_midref(y), prec, ARF_RND_DOWN);
        else
            inexact = arf_add(arb_midref(w), arb_midref(x), arb_midref(y), prec, ARF_RND_DOWN);
        mag_add(arb_radref(w), arb_radref(x), arb_radref(y));
        if (inexact)
            arf_mag_add_ulp(arb_radref(w), arb_radref(w), arb_midref(w), prec);

        _corners_add(xe[0], xe[1], x);
        _corners_add(ye[0], ye[1], y);
        ok = 1;
        for (a = 0; a < 2; a++)
            for (b = 0; b < 2; b++)
            {
                if (sub)
                    arf_sub(p, xe[a], ye[b], ARF_PREC_EXACT, ARF_RND_DOWN);
                else
                    arf_add(p, xe[a], ye[b], ARF_PREC_EXACT, ARF_RND_DOWN);
                ok = ok && arb_contains_arf(r, p);
            }

        if (ok)
        {
            mag_t m;
            mag_init(m);
            mag_mul_2exp_si(m, arb_radref(w), 2);
            ok = (mag_cmp(arb_radref(r), m) <= 0) || (prec != ARF_PREC_EXACT && arb_rel_accuracy_bits(r) >= prec - 2);
            mag_clear(m);
        }

        if (ok && prec == ARF_PREC_EXACT)
            ok = arb_equal(r, w);

        if (!ok)
        {
            flint_printf("FAIL: one-pass rounding\n");
            flint_printf("prec = %wd, sub = %d\n\n", prec, sub);
            flint_printf("x = "); arb_print(x); flint_printf("\n\n");
            flint_printf("y = "); arb_print(y); flint_printf("\n\n");
            flint_printf("r = "); arb_print(r); flint_printf("\n\n");
            flint_printf("w = "); arb_print(w); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(x); arb_clear(y); arb_clear(r); arb_clear(w);
        for (a = 0; a < 2; a++) { arf_clear(xe[a]); arf_clear(ye[a]); }
        arf_clear(p);
    }

    TEST_FUNCTION_END(state);
}
