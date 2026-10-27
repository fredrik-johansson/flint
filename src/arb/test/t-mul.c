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
arb_mul_naive(arb_t z, const arb_t x, const arb_t y, slong prec)
{
    arf_t zm_exact, zm_rounded, zr, t, u;

    arf_init(zm_exact);
    arf_init(zm_rounded);
    arf_init(zr);
    arf_init(t);
    arf_init(u);

    arf_mul(zm_exact, arb_midref(x), arb_midref(y), ARF_PREC_EXACT, ARF_RND_DOWN);
    arf_set_round(zm_rounded, zm_exact, prec, ARB_RND);

    /* rounding error */
    if (arf_equal(zm_exact, zm_rounded))
    {
        arf_zero(zr);
    }
    else
    {
        fmpz_t e;
        fmpz_init(e);

        /* more accurate, but not what we are testing
        arf_sub(zr, zm_exact, zm_rounded, MAG_BITS, ARF_RND_UP);
        arf_abs(zr, zr); */

        fmpz_sub_ui(e, ARF_EXPREF(zm_rounded), prec);
        arf_one(zr);
        arf_mul_2exp_fmpz(zr, zr, e);
        fmpz_clear(e);
    }

    /* propagated error - note that (signed infinity) * nonzero should be a signed
       infinity, meaning the error should *not* propagate */
    if (!arb_is_exact(x) && !(
            arb_is_nonzero(x) && arf_is_inf(arb_midref(y)) && mag_is_finite(arb_radref(y))))
    {
        arf_set_mag(t, arb_radref(x));
        arf_abs(u, arb_midref(y));
        arf_addmul(zr, t, u, MAG_BITS, ARF_RND_UP);
    }

    if (!arb_is_exact(y) && !(
            arb_is_nonzero(y) && arf_is_inf(arb_midref(x)) && mag_is_finite(arb_radref(x))))
    {
        arf_set_mag(t, arb_radref(y));
        arf_abs(u, arb_midref(x));
        arf_addmul(zr, t, u, MAG_BITS, ARF_RND_UP);
    }

    if (!arb_is_exact(x) && !arb_is_exact(y))
    {
        arf_set_mag(t, arb_radref(x));
        arf_set_mag(u, arb_radref(y));
        arf_addmul(zr, t, u, MAG_BITS, ARF_RND_UP);
    }

    arf_set(arb_midref(z), zm_rounded);
    arf_get_mag(arb_radref(z), zr);

    arf_clear(zm_exact);
    arf_clear(zm_rounded);
    arf_clear(zr);
    arf_clear(t);
    arf_clear(u);
}

/* arb_mul rounds the midpoint to nearest at a precision that may be
   reduced when the propagated error dominates, and sets the radius from
   the actual rounding error, whereas the naive version rounds down and
   adds one ulp: we check that the balls overlap and that the radii agree
   up to two ulp and a small factor; rigour is covered by the containment
   tests above. */
int
arb_mul_close(const arb_t z, const arb_t v, slong prec)
{
    arf_t t, u;
    mag_t m;
    fmpz_t e;
    int ok;

    if (arf_is_special(arb_midref(v)) || arf_is_special(arb_midref(z)))
        return arb_equal_mid_close_mag(z, v);

    arf_init(t);
    arf_init(u);
    mag_init(m);
    fmpz_init(e);

    /* the balls overlap (the midpoint may be rounded to a lower
       precision when the propagated error dominates) */
    ok = arb_overlaps(z, v);

    fmpz_sub_ui(e, ARF_EXPREF(arb_midref(v)), prec - 1);
    arf_set_ui(u, 257);
    arf_mul_2exp_fmpz(u, u, e);
    arf_mul_2exp_si(u, u, -8);

    /* radii within two ulp and a factor 257/256 of each other */
    arf_set_mag(t, arb_radref(v));
    arf_add(t, t, u, MAG_BITS, ARF_RND_UP);
    arf_mul_ui(t, t, 257, MAG_BITS, ARF_RND_UP);
    arf_mul_2exp_si(t, t, -8);
    arf_set_mag(u, arb_radref(z));
    ok = ok && (arf_cmp(u, t) <= 0);

    arf_set_mag(t, arb_radref(z));
    fmpz_sub_ui(e, ARF_EXPREF(arb_midref(v)), prec - 1);
    arf_set_ui(u, 257);
    arf_mul_2exp_fmpz(u, u, e);
    arf_mul_2exp_si(u, u, -8);
    arf_add(t, t, u, MAG_BITS, ARF_RND_UP);
    arf_mul_ui(t, t, 257, MAG_BITS, ARF_RND_UP);
    arf_mul_2exp_si(t, t, -8);
    arf_set_mag(u, arb_radref(v));
    ok = ok && (arf_cmp(u, t) <= 0);

    arf_clear(t);
    arf_clear(u);
    mag_clear(m);
    fmpz_clear(e);

    return ok;
}

/* the endpoints of a ball as arf values */
static void
_corners_mul(arf_t lo, arf_t hi, const arb_t x)
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
_randball_mul(arb_t x, flint_rand_t state, const fmpz * e)
{
    slong bits = n_randint(state, 8) == 0 ? 1 + n_randint(state, 3000) : 1 + n_randint(state, 200);
    arb_randtest(x, state, bits, 12);
    if (n_randint(state, 3) == 0)
        mag_zero(arb_radref(x));
    if (e != NULL && !arf_is_special(arb_midref(x)) && n_randint(state, 2))
        fmpz_add_si(ARF_EXPREF(arb_midref(x)), e, (slong) n_randint(state, 600) - 300);
}

/* Squaring and a product by an equal copy use different truncated products
   (sqrhigh, mulhigh) with different deficits, so the measured rounding
   errors, and next to a tie the midpoints, may differ at the level of
   2^-MAG_BITS ulp: the balls overlap, and the radii agree up to an ulp. */
#ifndef arb_close_rounding
#define arb_close_rounding arb_close_rounding
int
arb_close_rounding(const arb_t z, const arb_t v, slong prec)
{
    arf_t a, u;
    int ok;

    if (arb_equal_mid_close_mag(z, v))
        return 1;

    if (arf_is_special(arb_midref(z)) || arf_is_special(arb_midref(v)) ||
        !mag_is_finite(arb_radref(z)) || !mag_is_finite(arb_radref(v)) ||
        !arb_overlaps(z, v))
        return 0;

    arf_init(a);
    arf_init(u);

    arf_set_mag(a, arb_radref(z));
    arf_set_mag(u, arb_radref(v));
    arf_sub(a, a, u, ARF_PREC_EXACT, ARF_RND_DOWN);

    arf_one(u);
    arf_mul_2exp_fmpz(u, u, ARF_EXPREF(arb_midref(v)));
    arf_mul_2exp_si(u, u, -prec);

    ok = (arf_cmpabs(a, u) <= 0);

    arf_clear(a);
    arf_clear(u);

    return ok;
}
#endif

TEST_FUNCTION_START(arb_mul, state)
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

        arb_mul(c, a, b, 2 + n_randint(state, 200));
        fmpq_mul(z, x, y);

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

        arb_mul(a, a, b, 2 + n_randint(state, 200));
        fmpq_mul(z, x, y);

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

        arb_mul(b, a, b, 2 + n_randint(state, 200));
        fmpq_mul(z, x, y);

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

            prec = 2 + n_randint(state, 2000);

            switch (n_randint(state, 5))
            {
            case 0:
                arb_mul(z, x, y, prec);
                arb_mul_naive(v, x, y, prec);

                if (!arb_mul_close(z, v, prec))
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
                arb_mul(z, x, y, prec);
                arb_mul(v, x, x, prec);

                if (!arb_close_rounding(z, v, prec))
                {
                    flint_printf("FAIL (aliasing 1)!\n");
                    flint_printf("x = "); arb_print(x); flint_printf("\n\n");
                    flint_printf("z = "); arb_print(z); flint_printf("\n\n");
                    flint_printf("v = "); arb_print(v); flint_printf("\n\n");
                    flint_abort();
                }
                break;

            case 2:
                arb_mul(v, x, x, prec);
                arb_mul(x, x, x, prec);

                if (!arb_equal(v, x))
                {
                    flint_printf("FAIL (aliasing 2)!\n");
                    flint_printf("x = "); arb_print(x); flint_printf("\n\n");
                    flint_printf("z = "); arb_print(z); flint_printf("\n\n");
                    flint_printf("v = "); arb_print(v); flint_printf("\n\n");
                    flint_abort();
                }
                break;

            case 3:
                arb_mul(v, x, y, prec);
                arb_mul(x, x, y, prec);

                if (!arb_equal(v, x))
                {
                    flint_printf("FAIL (aliasing 3)!\n");
                    flint_printf("x = "); arb_print(x); flint_printf("\n\n");
                    flint_printf("y = "); arb_print(y); flint_printf("\n\n");
                    flint_printf("v = "); arb_print(v); flint_printf("\n\n");
                    flint_abort();
                }
                break;

            default:
                arb_mul(v, x, y, prec);
                arb_mul(x, y, x, prec);

                if (!arf_equal(arb_midref(x), arb_midref(v))
                    || !mag_close(arb_radref(x), arb_radref(v)))
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

    /* one-pass rounding: the exact products of the endpoints lie in the
       result, which is not much wider than the ball from a separately
       rounded midpoint */
    for (iter = 0; iter < 20000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x, y, r, w;
        arf_t xe[2], ye[2], p;
        mag_t xm, ym, zr;
        slong prec, a, b;
        int ok, inexact;

        arb_init(x); arb_init(y); arb_init(r); arb_init(w);
        for (a = 0; a < 2; a++) { arf_init(xe[a]); arf_init(ye[a]); }
        arf_init(p);
        mag_init(xm); mag_init(ym); mag_init(zr);

        prec = (n_randint(state, 8) == 0) ? 2 + n_randint(state, 3000) : 2 + n_randint(state, 300);
        if (n_randint(state, 50) == 0)
            prec = ARF_PREC_EXACT;
        _randball_mul(x, state, NULL);
        _randball_mul(y, state, NULL);
        if (n_randint(state, 4) == 0)
            arb_set(y, x);

        if (n_randint(state, 2))
        {
            arb_set(r, x);
            arb_mul(r, r, y, prec);
        }
        else
            arb_mul(r, x, y, prec);

        /* the reference: midpoint rounded by arf, radius bumped by an ulp */
        arf_get_mag(xm, arb_midref(x));
        arf_get_mag(ym, arb_midref(y));
        mag_mul(zr, xm, arb_radref(y));
        mag_addmul(zr, ym, arb_radref(x));
        mag_addmul(zr, arb_radref(x), arb_radref(y));
        inexact = arf_mul(arb_midref(w), arb_midref(x), arb_midref(y), prec, ARF_RND_DOWN);
        if (inexact)
            arf_mag_add_ulp(arb_radref(w), zr, arb_midref(w), prec);
        else
            mag_set(arb_radref(w), zr);

        _corners_mul(xe[0], xe[1], x);
        _corners_mul(ye[0], ye[1], y);
        ok = 1;
        for (a = 0; a < 2; a++)
            for (b = 0; b < 2; b++)
            {
                arf_mul(p, xe[a], ye[b], ARF_PREC_EXACT, ARF_RND_DOWN);
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
            ok = arf_equal(arb_midref(r), arb_midref(w));

        if (!ok)
        {
            flint_printf("FAIL: one-pass rounding\n");
            flint_printf("prec = %wd\n\n", prec);
            flint_printf("x = "); arb_print(x); flint_printf("\n\n");
            flint_printf("y = "); arb_print(y); flint_printf("\n\n");
            flint_printf("r = "); arb_print(r); flint_printf("\n\n");
            flint_printf("w = "); arb_print(w); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(x); arb_clear(y); arb_clear(r); arb_clear(w);
        for (a = 0; a < 2; a++) { arf_clear(xe[a]); arf_clear(ye[a]); }
        arf_clear(p);
        mag_clear(xm); mag_clear(ym); mag_clear(zr);
    }

    TEST_FUNCTION_END(state);
}
