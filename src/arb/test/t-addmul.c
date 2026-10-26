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

#ifndef mag_close
#define mag_close mag_close
int
mag_close(const mag_t am, const mag_t bm)
{
    arf_t t, a, b;
    int res1, res2;

    arf_init(t);
    arf_init(a);
    arf_init(b);

    arf_set_mag(a, am);
    arf_set_mag(b, bm);

    arf_mul_ui(t, b, 257, MAG_BITS, ARF_RND_UP);
    arf_mul_2exp_si(t, t, -8);
    res1 = arf_cmp(a, t) <= 0;

    arf_mul_ui(t, a, 257, MAG_BITS, ARF_RND_UP);
    arf_mul_2exp_si(t, t, -8);
    res2 = arf_cmp(b, t) <= 0;

    arf_clear(t);
    arf_clear(a);
    arf_clear(b);

    return res1 && res2;
}
#endif

#ifndef arb_equal_mid_close_mag
#define arb_equal_mid_close_mag arb_equal_mid_close_mag
int
arb_equal_mid_close_mag(const arb_t a, const arb_t b)
{
    return arf_equal(arb_midref(a), arb_midref(b)) &&
        (mag_close(arb_radref(a), arb_radref(b)) ||
            /* If a's and b's centres are infinite but their radii are finite, the radii don't need
               to be close: they represent signed infinity regardless. If their centres are NaN,
               then we should ignore their radii. */
            (arf_is_inf(arb_midref(a)) && arf_is_inf(arb_midref(b)) &&
                mag_is_finite(arb_radref(a)) && mag_is_finite(arb_radref(b))) ||
            (arf_is_nan(arb_midref(a)) && arf_is_nan(arb_midref(b))));
}
#endif

void
arb_addmul_naive(arb_t z, const arb_t x, const arb_t y, slong prec)
{
    arb_t t;
    arb_init(t);
    arb_mul(t, x, y, ARF_PREC_EXACT);
    arb_add(z, z, t, prec);
    arb_clear(t);
}

/* the same value enclosed with a comparable radius (the fused operation
   and the exact product followed by an addition round the midpoint at
   slightly different points; a truncated product adds its deficit) */
static int
_close_balls_addmul(const arb_t a, const arb_t b, slong prec)
{
    mag_t t, u;
    int ok;

    if (!arb_is_finite(a) || !arb_is_finite(b))
        return arb_equal_mid_close_mag(a, b);

    if (!arb_overlaps(a, b))
        return 0;

    mag_init(t);
    mag_init(u);
    arf_get_mag(u, arb_midref(a));
    mag_mul_2exp_si(u, u, 2 - prec);
    mag_mul_2exp_si(t, arb_radref(b), 2);
    mag_add(t, t, u);
    ok = (mag_cmp(arb_radref(a), t) <= 0);
    mag_mul_2exp_si(t, arb_radref(a), 2);
    mag_add(t, t, u);
    ok = ok && (mag_cmp(arb_radref(b), t) <= 0);
    mag_clear(t);
    mag_clear(u);
    return ok;
}

/* the endpoints of a ball as arf values */
static void
_corners_addmul(arf_t lo, arf_t hi, const arb_t x)
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
_randball_addmul(arb_t x, flint_rand_t state, const fmpz * e)
{
    slong bits = n_randint(state, 8) == 0 ? 1 + n_randint(state, 3000) : 1 + n_randint(state, 200);
    arb_randtest(x, state, bits, 12);
    if (n_randint(state, 3) == 0)
        mag_zero(arb_radref(x));
    if (e != NULL && !arf_is_special(arb_midref(x)) && n_randint(state, 2))
        fmpz_add_si(ARF_EXPREF(arb_midref(x)), e, (slong) n_randint(state, 600) - 300);
}

TEST_FUNCTION_START(arb_addmul, state)
{
    slong iter, iter2;

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
        arb_get_rand_fmpq(z, state, c, 1 + n_randint(state, 200));

        arb_addmul(c, a, b, 2 + n_randint(state, 200));
        fmpq_addmul(z, x, y);

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
        fmpq_set(z, x);

        arb_addmul(a, a, b, 2 + n_randint(state, 200));
        fmpq_addmul(z, x, y);

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
        fmpq_set(z, y);

        arb_addmul(b, a, b, 2 + n_randint(state, 200));
        fmpq_addmul(z, x, y);

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

    /* main test */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x, y, z, v;
        slong prec;

        arb_init(x);
        arb_init(y);
        arb_init(z);
        arb_init(v);

        for (iter2 = 0; iter2 < 100; iter2++)
        {
            arb_randtest_special(x, state, n_randint(state,2) ? 2000 : 200, 200);
            arb_randtest_special(y, state, n_randint(state,2) ? 2000 : 200, 200);
            arb_randtest_special(z, state, n_randint(state,2) ? 2000 : 200, 200);

            prec = 2 + n_randint(state, 2000);

            switch (n_randint(state, 5))
            {
            case 0:
                arb_set(v, z);
                arb_addmul(z, x, y, prec);
                arb_addmul_naive(v, x, y, prec);

                if (!_close_balls_addmul(z, v, prec))
                {
                    flint_printf("FAIL!\n");
                    flint_printf("x = "); arb_print(x); flint_printf("\n\n");
                    flint_printf("y = "); arb_print(y); flint_printf("\n\n");
                    flint_printf("z = "); arb_print(z); flint_printf("\n\n");
                    flint_printf("v = "); arb_print(v); flint_printf("\n\n");
                    flint_abort();
                }
                break;

            case 1:
                arb_set(y, x);
                arb_set(z, v);
                arb_addmul(z, x, y, prec);
                arb_addmul(v, x, x, prec);

                if (!arb_equal_mid_close_mag(z, v))
                {
                    flint_printf("FAIL (aliasing 1)!\n");
                    flint_printf("x = "); arb_print(x); flint_printf("\n\n");
                    flint_printf("z = "); arb_print(z); flint_printf("\n\n");
                    flint_printf("v = "); arb_print(v); flint_printf("\n\n");
                    flint_abort();
                }
                break;

            case 2:
                arb_set(v, x);
                arb_addmul(v, x, x, prec);
                arb_addmul(x, x, x, prec);

                if (!arb_equal_mid_close_mag(x, v))
                {
                    flint_printf("FAIL (aliasing 2)!\n");
                    flint_printf("x = "); arb_print(x); flint_printf("\n\n");
                    flint_printf("v = "); arb_print(v); flint_printf("\n\n");
                    flint_abort();
                }
                break;

            case 3:
                arb_set(v, x);
                arb_addmul(v, x, y, prec);
                arb_addmul(x, x, y, prec);

                if (!arb_equal_mid_close_mag(x, v))
                {
                    flint_printf("FAIL (aliasing 3)!\n");
                    flint_printf("x = "); arb_print(x); flint_printf("\n\n");
                    flint_printf("y = "); arb_print(y); flint_printf("\n\n");
                    flint_printf("v = "); arb_print(v); flint_printf("\n\n");
                    flint_abort();
                }
                break;

            default:
                arb_set(v, x);
                arb_addmul(v, x, y, prec);
                arb_addmul(x, y, x, prec);

                if (!arb_equal_mid_close_mag(x, v))
                {
                    flint_printf("FAIL (aliasing 4)!\n");
                    flint_printf("x = "); arb_print(x); flint_printf("\n\n");
                    flint_printf("y = "); arb_print(y); flint_printf("\n\n");
                    flint_printf("v = "); arb_print(v); flint_printf("\n\n");
                    flint_abort();
                }
                break;
            }
        }

        arb_clear(x);
        arb_clear(y);
        arb_clear(z);
        arb_clear(v);
    }

    /* one-pass rounding: the exact values at the endpoints lie in the
       result, which is not much wider than the ball from a separately
       rounded midpoint */
    for (iter = 0; iter < 20000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x, y, z, r, w;
        arf_t xe[2], ye[2], ze[2], p, q;
        mag_t xm, ym, zr;
        slong prec, a, b, c;
        int ok, sub, inexact;

        arb_init(x); arb_init(y); arb_init(z); arb_init(r); arb_init(w);
        for (a = 0; a < 2; a++) { arf_init(xe[a]); arf_init(ye[a]); arf_init(ze[a]); }
        arf_init(p); arf_init(q);
        mag_init(xm); mag_init(ym); mag_init(zr);

        prec = (n_randint(state, 8) == 0) ? 2 + n_randint(state, 3000) : 2 + n_randint(state, 300);
        if (n_randint(state, 50) == 0)
            prec = ARF_PREC_EXACT;
        sub = n_randint(state, 2);
        _randball_addmul(x, state, NULL);
        _randball_addmul(y, state, NULL);
        if (n_randint(state, 4) == 0)
            arb_set(y, x);
        if (n_randint(state, 4) == 0)
        {
            /* cancellation */
            arb_mul(z, x, y, prec == ARF_PREC_EXACT ? 300 : prec);
            if (!sub)
                arb_neg(z, z);
        }
        else
        {
            fmpz_t e;
            fmpz_init(e);
            if (!arf_is_special(arb_midref(x)) && !arf_is_special(arb_midref(y)))
                fmpz_add(e, ARF_EXPREF(arb_midref(x)), ARF_EXPREF(arb_midref(y)));
            _randball_addmul(z, state, e);
            fmpz_clear(e);
        }

        arb_set(r, z);
        if (sub) arb_submul(r, x, y, prec); else arb_addmul(r, x, y, prec);

        /* the reference: midpoint rounded by arf, radius bumped by an ulp */
        arb_set(w, z);
        arf_get_mag(xm, arb_midref(x));
        arf_get_mag(ym, arb_midref(y));
        mag_addmul(arb_radref(w), xm, arb_radref(y));
        mag_addmul(arb_radref(w), ym, arb_radref(x));
        mag_addmul(arb_radref(w), arb_radref(x), arb_radref(y));
        if (sub)
            inexact = arf_submul(arb_midref(w), arb_midref(x), arb_midref(y), prec, ARF_RND_DOWN);
        else
            inexact = arf_addmul(arb_midref(w), arb_midref(x), arb_midref(y), prec, ARF_RND_DOWN);
        if (inexact)
            arf_mag_add_ulp(arb_radref(w), arb_radref(w), arb_midref(w), prec);

        _corners_addmul(xe[0], xe[1], x);
        _corners_addmul(ye[0], ye[1], y);
        _corners_addmul(ze[0], ze[1], z);
        ok = 1;
        for (a = 0; a < 2; a++)
            for (b = 0; b < 2; b++)
                for (c = 0; c < 2; c++)
                {
                    arf_mul(p, xe[a], ye[b], ARF_PREC_EXACT, ARF_RND_DOWN);
                    if (sub)
                        arf_sub(q, ze[c], p, ARF_PREC_EXACT, ARF_RND_DOWN);
                    else
                        arf_add(q, ze[c], p, ARF_PREC_EXACT, ARF_RND_DOWN);
                    ok = ok && arb_contains_arf(r, q);
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
            ok = arf_equal(arb_midref(r), arb_midref(w));

        if (!ok)
        {
            flint_printf("FAIL: one-pass rounding\n");
            flint_printf("prec = %wd, sub = %d\n\n", prec, sub);
            flint_printf("x = "); arb_print(x); flint_printf("\n\n");
            flint_printf("y = "); arb_print(y); flint_printf("\n\n");
            flint_printf("z = "); arb_print(z); flint_printf("\n\n");
            flint_printf("r = "); arb_print(r); flint_printf("\n\n");
            flint_printf("w = "); arb_print(w); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(x); arb_clear(y); arb_clear(z); arb_clear(r); arb_clear(w);
        for (a = 0; a < 2; a++) { arf_clear(xe[a]); arf_clear(ye[a]); arf_clear(ze[a]); }
        arf_clear(p); arf_clear(q);
        mag_clear(xm); mag_clear(ym); mag_clear(zr);
    }

    TEST_FUNCTION_END(state);
}
