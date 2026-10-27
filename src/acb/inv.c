/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "arb/impl.h"
#include "mpn_extras/fp_small.h"

/*
    Fixed-point bounds on LAGOM mags for the error bounds below: the same
    MAG_BITS-bit mantissas as the mag_fast functions, with the lower bounds
    (truncated) that those lack, and a division rounded up.
*/

/* lower bound of |x| for a nonzero arf */
FLINT_FORCE_INLINE void
_mag_lower_arf(mag_t m, arf_srcptr x)
{
    ulong top;
    ARF_GET_TOP_LIMB(top, x);
    MAG_MAN(m) = top >> (FLINT_BITS - MAG_BITS);
    MAG_EXP(m) = ARF_EXP(x);
}

/* lower bound of |mid| - rad for a nonzero LAGOM midpoint and a finite
   LAGOM radius (zero if the difference is not positive): the top
   FLINT_BITS - 2 bits of the midpoint against the radius rounded up, the
   general routine when the difference nearly cancels */
static void
_mag_lower_arf_sub(mag_t r, arf_srcptr mid, mag_srcptr rad)
{
    slong shift, bc;
    ulong xm, rm, d;

    if (MAG_MAN(rad) == 0)
    {
        _mag_lower_arf(r, mid);
        return;
    }

    shift = ARF_EXP(mid) - MAG_EXP(rad);
    if (shift < 0)
    {
        MAG_MAN(r) = 0;
        MAG_EXP(r) = 0;
        return;
    }

    ARF_GET_TOP_LIMB(xm, mid);
    xm >>= 2;    /* |mid| >= xm 2^(exp - (FLINT_BITS - 2)) */

    /* rad in the same units, rounded up */
    if (shift <= FLINT_BITS - 2 - MAG_BITS)
        rm = MAG_MAN(rad) << (FLINT_BITS - 2 - MAG_BITS - shift);
    else if (shift < FLINT_BITS - 2)
        rm = (MAG_MAN(rad) >> (shift - (FLINT_BITS - 2 - MAG_BITS))) + 1;
    else
        rm = 1;

    /* the difference must keep MAG_BITS significant bits with a margin
       (always the arf route on 32-bit limbs, which have no room for it) */
    if (FLINT_BITS == 32 || xm <= rm || xm - rm < (UWORD(1) << (FLINT_BITS - 2 - 20)))
    {
        /* nearly cancelling (or not positive at this precision): the
           difference computed exactly enough with arf arithmetic */
        arf_t t;
        arf_init(t);
        arf_set_mag(t, rad);
        if (arf_sgn(mid) > 0)
            arf_sub(t, mid, t, MAG_BITS, ARF_RND_DOWN);
        else
        {
            arf_add(t, mid, t, MAG_BITS, ARF_RND_DOWN);
            arf_neg(t, t);
        }
        MAG_EXP(r) = 0;    /* r is uninitialised storage: a small exponent for the fmpz write */
        if (arf_sgn(t) <= 0)
            MAG_MAN(r) = 0;
        else
            arf_get_mag_lower(r, t);
        arf_clear(t);
        return;
    }

    d = xm - rm;
    bc = FLINT_BIT_COUNT(d);
    MAG_MAN(r) = d >> (bc - MAG_BITS);
    MAG_EXP(r) = ARF_EXP(mid) - (FLINT_BITS - 2) + bc;
}

/* lower bound of x y */
FLINT_FORCE_INLINE void
_mag_mul_lower(mag_t r, mag_srcptr x, mag_srcptr y)
{
    uint64_t p = (uint64_t) MAG_MAN(x) * MAG_MAN(y);

    if (p == 0)
    {
        MAG_MAN(r) = 0;
        MAG_EXP(r) = 0;
        return;
    }

    /* p >= 2^(2 MAG_BITS - 2): one or two leading zero bits */
    if (p >> (2 * MAG_BITS - 1))
    {
        MAG_MAN(r) = (ulong) (p >> MAG_BITS);
        MAG_EXP(r) = MAG_EXP(x) + MAG_EXP(y);
    }
    else
    {
        MAG_MAN(r) = (ulong) (p >> (MAG_BITS - 1));
        MAG_EXP(r) = MAG_EXP(x) + MAG_EXP(y) - 1;
    }
}

/* lower bound of x + y */
FLINT_FORCE_INLINE void
_mag_add_lower(mag_t r, mag_srcptr x, mag_srcptr y)
{
    slong shift;
    ulong s;

    if (MAG_MAN(x) == 0)
    {
        *r = *y;
        return;
    }
    if (MAG_MAN(y) == 0)
    {
        *r = *x;
        return;
    }
    if (MAG_EXP(x) < MAG_EXP(y))
    {
        mag_srcptr t = x;
        x = y;
        y = t;
    }

    shift = MAG_EXP(x) - MAG_EXP(y);
    s = MAG_MAN(x) + ((shift < MAG_BITS) ? (MAG_MAN(y) >> shift) : 0);
    MAG_EXP(r) = MAG_EXP(x);
    if (s >> MAG_BITS)
    {
        s >>= 1;
        MAG_EXP(r)++;
    }
    MAG_MAN(r) = s;
}

/* upper bound of x / y, y nonzero */
FLINT_FORCE_INLINE void
_mag_div_upper(mag_t r, mag_srcptr x, mag_srcptr y)
{
    ulong q;

    if (MAG_MAN(x) == 0)
    {
        MAG_MAN(r) = 0;
        MAG_EXP(r) = 0;
        return;
    }

    /* x / y = (q / 2^MAG_BITS) 2^(ex - ey) with q in (2^(MAG_BITS-1), 2^(MAG_BITS+1)) */
    q = (ulong) ((((uint64_t) MAG_MAN(x)) << MAG_BITS) / MAG_MAN(y)) + 1;
    MAG_EXP(r) = MAG_EXP(x) - MAG_EXP(y);
    if (q >> MAG_BITS)
    {
        q = (q + 1) >> 1;
        MAG_EXP(r)++;
    }
    MAG_MAN(r) = q;
    MAG_FAST_ADJUST_ONE_TOO_LARGE(r);
}

/*
    1 / (a + b i) for LAGOM parts with nonzero midpoints a, b and finite
    radii x, y: the midpoints a / s, -b / s with s a lower approximation of
    a^2 + b^2 at hprec = prec + 3 bits, each quotient rounded once with the
    measured error; the radius gets the error of s, |a| / s 2^(1 - hprec),
    and the propagated error of the inputs, bounded as in the general code
    below: the real part is exactly

        (a + x')/((a + x')^2 + (b + y')^2) - a/(a^2 + b^2) = P / Q with
        P = (b^2 - a^2) x' - a (x'^2 + y'^2 + 2 y' b),
        Q = (a^2 + b^2) ((a + x')^2 + (b + y')^2)

    for |x'| <= x, |y'| <= y, so |P| <= |a^2 - b^2| x + |a| (x^2 + y^2 +
    2 |b| y) and Q >= (a^2 + b^2) ((|a| - x)^2 + (|b| - y)^2), and
    likewise for the imaginary part. Returns 0 (without touching res) for
    other inputs.
*/
static int
_acb_inv_lagom(acb_t res, const acb_t z, slong prec)
{
    arf_srcptr a = arb_midref(acb_realref(z)), b = arb_midref(acb_imagref(z));
    mag_srcptr x = arb_radref(acb_realref(z)), y = arb_radref(acb_imagref(z));
    mag_t A, B, S, er, ei, t, u, v, w, n;
    arf_t s;
    slong hprec, qexp, acc, pr, pi, sexp;
    int small = 0, sinexact = 1;
#if FLINT_BITS == 64
    ulong a1 = 0, a0 = 0, b1 = 0, b0 = 0, s2 = 0, s1 = 0, s0 = 0;
#endif

    {
        const ulong M = UWORD(1) << (FLINT_BITS - 6);
        ulong range = (ulong) (ARF_EXP(a) + M) | (ulong) (ARF_EXP(b) + M)
                    | (ulong) (MAG_EXP(x) + M) | (ulong) (MAG_EXP(y) + M);
        if (range >= 2 * M || ARF_IS_SPECIAL(a) || ARF_IS_SPECIAL(b) || prec < 2 || prec == ARF_PREC_EXACT)
            return 0;
    }

#if FLINT_BITS == 64
    if (prec <= 2 * FLINT_BITS && ARF_SIZE(a) <= 2 && ARF_SIZE(b) <= 2)
    {
        /* The register path: s = a^2 + b^2 to three limbs (a lower
           approximation; the deficit is below one unit of the lowest
           limb) and the quotients by the three-limb division kernel. */
        ulong p3, p2, p1, p0, q3, q2, q1, q0, c, sticky, lead;
        slong ea = ARF_EXP(a), eb = ARF_EXP(b), d;

        a1 = ARF_NOPTR_D(a)[ARF_SIZE(a) - 1];
        a0 = (ARF_SIZE(a) == 2) ? ARF_NOPTR_D(a)[0] : 0;
        b1 = ARF_NOPTR_D(b)[ARF_SIZE(b) - 1];
        b0 = (ARF_SIZE(b) == 2) ? ARF_NOPTR_D(b)[0] : 0;

        /* the squares, aligned to the larger exponent: a^2 = P 2^(2 ea - 256)
           with P in [2^254, 2^256), and likewise b^2 */
        FLINT_MPN_SQR_2X2(p3, p2, p1, p0, a1, a0);
        FLINT_MPN_SQR_2X2(q3, q2, q1, q0, b1, b0);
        d = 2 * (ea - eb);
        if (d < 0)
        {
            ulong tmp;
            tmp = p3; p3 = q3; q3 = tmp;
            tmp = p2; p2 = q2; q2 = tmp;
            tmp = p1; p1 = q1; q1 = tmp;
            tmp = p0; p0 = q0; q0 = tmp;
            d = -d;
            sexp = 2 * eb;
        }
        else
            sexp = 2 * ea;

        /* (q3, q2, q1, q0) >> d, the bits shifted out going to sticky */
        if (d >= 4 * FLINT_BITS)
        {
            sticky = q3 | q2 | q1 | q0;
            q3 = q2 = q1 = q0 = 0;
        }
        else
        {
            ulong r = d % FLINT_BITS;
            sticky = 0;
            while (d >= FLINT_BITS)
            {
                sticky |= q0;
                q0 = q1; q1 = q2; q2 = q3; q3 = 0;
                d -= FLINT_BITS;
            }
            if (r != 0)
            {
                sticky |= q0 << (FLINT_BITS - r);
                q0 = (q0 >> r) | (q1 << (FLINT_BITS - r));
                q1 = (q1 >> r) | (q2 << (FLINT_BITS - r));
                q2 = (q2 >> r) | (q3 << (FLINT_BITS - r));
                q3 = q3 >> r;
            }
        }

        /* the sum, with carry c, in [2^254, 2^257) */
        p0 += q0; c = p0 < q0;
        p1 += c; c = p1 < c; p1 += q1; c += p1 < q1;
        p2 += c; c = p2 < c; p2 += q2; c += p2 < q2;
        p3 += c; c = p3 < c; p3 += q3; c += p3 < q3;

        /* normalise to three limbs: the value is (c, p3, p2, p1, p0)
           2^(sexp - 256) */
        if (c != 0)
        {
            /* c = 1: shift right by one, the top bit set */
            s2 = (p3 >> 1) | (UWORD(1) << (FLINT_BITS - 1));
            s1 = (p2 >> 1) | (p3 << (FLINT_BITS - 1));
            s0 = (p1 >> 1) | (p2 << (FLINT_BITS - 1));
            sticky |= (p1 << (FLINT_BITS - 1)) | p0;
            sexp = sexp + 1;
        }
        else
        {
            lead = flint_clz(p3);     /* 0 or 1 */
            if (lead == 0)
            {
                s2 = p3; s1 = p2; s0 = p1;
            }
            else
            {
                s2 = (p3 << 1) | (p2 >> (FLINT_BITS - 1));
                s1 = (p2 << 1) | (p1 >> (FLINT_BITS - 1));
                s0 = (p1 << 1) | (p0 >> (FLINT_BITS - 1));
                p0 = p0 << 1;
            }
            sticky |= p0;
            sexp = sexp - lead;
        }
        sinexact = (sticky != 0);

        /* s = (s2, s1, s0) 2^(sexp - 192) <= a^2 + b^2 < s (1 + 2^(2 - 192)) */
        MAG_MAN(S) = s2 >> (FLINT_BITS - MAG_BITS);
        MAG_EXP(S) = sexp;
        hprec = 3 * FLINT_BITS - 1;
        small = 1;
    }
    else
#endif
    {
        hprec = FLINT_MAX(prec + 3, MAG_BITS);
        arf_init(s);
        sinexact = arf_sosq(s, a, b, hprec, ARF_RND_DOWN);   /* s <= a^2 + b^2 < s (1 + 2^(1 - hprec)) */

        if (!ARF_IS_LAGOM(s) || ARF_IS_SPECIAL(s))
        {
            arf_clear(s);
            return 0;
        }

        _mag_lower_arf(S, s);               /* S <= s <= a^2 + b^2 */
        sexp = ARF_EXP(s);
    }

    mag_fast_init_set_arf(A, a);        /* |a| <= A */
    mag_fast_init_set_arf(B, b);

    /* the error of the denominator: |a / s| 2^(1 - hprec), |b / s| 2^(1 - hprec)
       (none when s is exact) */
    if (sinexact)
    {
        _mag_div_upper(er, A, S);
        MAG_EXP(er) += 1 - hprec;
        _mag_div_upper(ei, B, S);
        MAG_EXP(ei) += 1 - hprec;
    }
    else
    {
        MAG_MAN(er) = MAG_MAN(ei) = 0;
        MAG_EXP(er) = MAG_EXP(ei) = 0;
    }

    if (MAG_MAN(x) != 0 || MAG_MAN(y) != 0)
    {
        /* t = ((|a| - x)^2 + (|b| - y)^2) s, a lower bound of Q */
        _mag_lower_arf_sub(u, a, x);
        _mag_lower_arf_sub(v, b, y);
        if (MAG_MAN(u) == 0 && MAG_MAN(v) == 0)
        {
            /* the ball contains zero */
            if (!small)
                arf_clear(s);
            acb_indeterminate(res);
            return 1;
        }
        _mag_mul_lower(u, u, u);
        _mag_mul_lower(v, v, v);
        _mag_add_lower(t, u, v);
        _mag_mul_lower(t, t, S);

        /* w = max(A, B)^2 >= |a^2 - b^2|; v = x^2 + y^2 */
        if (MAG_EXP(A) > MAG_EXP(B) || (MAG_EXP(A) == MAG_EXP(B) && MAG_MAN(A) >= MAG_MAN(B)))
            mag_fast_mul(w, A, A);
        else
            mag_fast_mul(w, B, B);
        mag_fast_mul(v, x, x);
        mag_fast_addmul(v, y, y);

        /* n = A (v + 2 B y) + w x, the real numerator */
        mag_fast_mul(n, B, y);
        mag_fast_mul_2exp_si(n, n, 1);
        mag_fast_add(n, n, v);
        mag_fast_mul(n, n, A);
        mag_fast_addmul(n, w, x);
        _mag_div_upper(u, n, t);
        mag_fast_add(er, er, u);

        /* n = B (v + 2 A x) + w y, the imaginary numerator */
        mag_fast_mul(n, A, x);
        mag_fast_mul_2exp_si(n, n, 1);
        mag_fast_add(n, n, v);
        mag_fast_mul(n, n, B);
        mag_fast_addmul(n, w, y);
        _mag_div_upper(u, n, t);
        mag_fast_add(ei, ei, u);
    }

    /* the precision of each part lowered to what its error leaves (the
       quotients are below 2^qexp) */
    qexp = FLINT_MAX(ARF_EXP(a), ARF_EXP(b)) - sexp + 1;
    acc = qexp - MAG_EXP(er);
    pr = (MAG_MAN(er) != 0 && acc < prec - MAG_BITS) ? FLINT_MAX(2, acc + MAG_BITS) : prec;
    acc = qexp - MAG_EXP(ei);
    pi = (MAG_MAN(ei) != 0 && acc < prec - MAG_BITS) ? FLINT_MAX(2, acc + MAG_BITS) : prec;

#if FLINT_BITS == 64
    if (small)
    {
        /* each quotient exactly truncated to three limbs (within one unit
           if inexact: extra = inexact), rounded with its error; the
           operands are in registers, so aliasing of res and z is harmless */
        ulong q[3], sh;
        slong n = (prec > FLINT_BITS) ? 2 : 1, exp;
        int inexact, sa = ARF_SGNBIT(a), sb = ARF_SGNBIT(b);
        mag_srcptr erp = (MAG_MAN(er) != 0) ? er : NULL;
        mag_srcptr eip = (MAG_MAN(ei) != 0) ? ei : NULL;

        /* (no-ops: pr, pi <= prec <= n FLINT_BITS; they tell the compiler
           that the roundings below stay within q) */
        pr = FLINT_MIN(pr, n * FLINT_BITS);
        pi = FLINT_MIN(pi, n * FLINT_BITS);
#if FLINT_PREINVERT_LIMB_USE_NATIVE
#define _ACB_INV_FDIV(q, sh, a1, a0) _flint_mpn_fdiv_3_kernel((q) + 2, (q) + 1, (q), &(sh), a1, a0, s2, s1, s0)
#else
        ulong dinv = flint_mpn_preinv1(s2, s1);
#define _ACB_INV_FDIV(q, sh, a1, a0) _flint_mpn_fdiv_3_kernel_preinv((q) + 2, (q) + 1, (q), &(sh), a1, a0, s2, s1, s0, dinv)
#endif

        inexact = _ACB_INV_FDIV(q, sh, b1, b0);
        exp = ARF_EXP(b) - sexp + (slong) sh - 3 * FLINT_BITS;
        if (n == 1)
        {
            inexact |= (q[0] != 0);
            _arb_set_round_guarded(acb_imagref(res), q + 1, 1, sb ^ 1, exp + FLINT_BITS, pi, inexact, eip, NULL);
        }
        else
            _arb_set_round_guarded(acb_imagref(res), q, 2, sb ^ 1, exp, pi, inexact, eip, NULL);

        inexact = _ACB_INV_FDIV(q, sh, a1, a0);
#undef _ACB_INV_FDIV
        exp = ARF_EXP(a) - sexp + (slong) sh - 3 * FLINT_BITS;
        if (n == 1)
        {
            inexact |= (q[0] != 0);
            _arb_set_round_guarded(acb_realref(res), q + 1, 1, sa, exp + FLINT_BITS, pr, inexact, erp, NULL);
        }
        else
            _arb_set_round_guarded(acb_realref(res), q, 2, sa, exp, pr, inexact, erp, NULL);

        return 1;
    }
#endif

    /* the imaginary part first: it reads b, which the real part may
       overwrite when res = z */
    _arb_div_mid(acb_imagref(res), b, s, pi, (MAG_MAN(ei) != 0) ? ei : NULL);
    arf_neg(arb_midref(acb_imagref(res)), arb_midref(acb_imagref(res)));
    _arb_div_mid(acb_realref(res), a, s, pr, (MAG_MAN(er) != 0) ? er : NULL);

    arf_clear(s);
    return 1;
}

static void
_arb_arf_div_rounded_den(arb_t res, const arf_t x, const arf_t y, int yinexact, slong prec)
{
    int inexact = arf_div(arb_midref(res), x, y, prec, ARB_RND);

    if (yinexact && !arf_is_special(arb_midref(res)))
        arf_mag_set_ulp(arb_radref(res), arb_midref(res), prec - 1);
    else if (inexact)
        arf_mag_set_ulp(arb_radref(res), arb_midref(res), prec);
    else
        mag_zero(arb_radref(res));
}

static void
_arb_arf_div_rounded_den_add_err(arb_t res, const arf_t x, const arf_t y, int yinexact, slong prec)
{
    int inexact = arf_div(arb_midref(res), x, y, prec, ARB_RND);

    if (yinexact && !arf_is_special(arb_midref(res)))
        arf_mag_add_ulp(arb_radref(res), arb_radref(res), arb_midref(res), prec - 1);
    else if (inexact)
        arf_mag_add_ulp(arb_radref(res), arb_radref(res), arb_midref(res), prec);
}

void
acb_inv(acb_t res, const acb_t z, slong prec)
{
    mag_t am, bm;
    slong hprec;

#define a arb_midref(acb_realref(z))
#define b arb_midref(acb_imagref(z))
#define x arb_radref(acb_realref(z))
#define y arb_radref(acb_imagref(z))

    /* choose precision for the floating-point approximation of a^2+b^2 so
       that the double rounding result in less than
       2 ulp error; also use at least MAG_BITS bits since the
       value will be recycled for error bounds */
    hprec = FLINT_MAX(prec + 3, MAG_BITS);

    if (arb_is_zero(acb_imagref(z)))
    {
        arb_inv(acb_realref(res), acb_realref(z), prec);
        arb_zero(acb_imagref(res));
        return;
    }

    if (arb_is_zero(acb_realref(z)))
    {
        arb_inv(acb_imagref(res), acb_imagref(z), prec);
        arb_neg(acb_imagref(res), acb_imagref(res));
        arb_zero(acb_realref(res));
        return;
    }

    if (!acb_is_finite(z))
    {
        acb_indeterminate(res);
        return;
    }

    if (_acb_inv_lagom(res, z, prec))
        return;

    if (mag_is_zero(x) && mag_is_zero(y))
    {
        int inexact;

        arf_t a2b2;
        arf_init(a2b2);

        inexact = arf_sosq(a2b2, a, b, hprec, ARF_RND_DOWN);

        if (arf_is_special(a2b2))
        {
            acb_indeterminate(res);
        }
        else
        {
            _arb_arf_div_rounded_den(acb_realref(res), a, a2b2, inexact, prec);
            _arb_arf_div_rounded_den(acb_imagref(res), b, a2b2, inexact, prec);
            arf_neg(arb_midref(acb_imagref(res)), arb_midref(acb_imagref(res)));
        }

        arf_clear(a2b2);
        return;
    }

    mag_init(am);
    mag_init(bm);

    /* first bound |a|-x, |b|-y */
    arb_get_mag_lower(am, acb_realref(z));
    arb_get_mag_lower(bm, acb_imagref(z));

    if ((mag_is_zero(am) && mag_is_zero(bm)))
    {
        acb_indeterminate(res);
    }
    else
    {
        /*
        The propagated error in the real part is given exactly by

             (a+x')/((a+x')^2+(b+y'))^2 - a/(a^2+b^2) = P / Q,

             P = [(b^2-a^2) x' - a (x'^2+y'^2 + 2y'b)]
             Q = [(a^2+b^2)((a+x')^2+(b+y')^2)]

        where |x'| <= x and |y'| <= y, and analogously for the imaginary part.
        */
        mag_t t, u, v, w;
        arf_t a2b2;
        int inexact;

        mag_init(t);
        mag_init(u);
        mag_init(v);
        mag_init(w);

        arf_init(a2b2);

        inexact = arf_sosq(a2b2, a, b, hprec, ARF_RND_DOWN);

        /* compute denominator */
        /* t = (|a|-x)^2 + (|b|-x)^2 (lower bound) */
        mag_mul_lower(t, am, am);
        mag_mul_lower(u, bm, bm);
        mag_add_lower(t, t, u);
        /* u = a^2 + b^2 (lower bound) */
        arf_get_mag_lower(u, a2b2);
        /* t = ((|a|-x)^2 + (|b|-x)^2)(a^2 + b^2) (lower bound) */
        mag_mul_lower(t, t, u);

        /* compute numerator */
        /* real: |a^2-b^2| x  + |a| ((x^2 + y^2) + 2 |b| y)) */
        /* imag: |a^2-b^2| y  + |b| ((x^2 + y^2) + 2 |a| x)) */
        /* am, bm = upper bounds for a, b */
        arf_get_mag(am, a);
        arf_get_mag(bm, b);

        /* v = x^2 + y^2 */
        mag_mul(v, x, x);
        mag_addmul(v, y, y);

        /* u = |a| ((x^2 + y^2) + 2 |b| y) */
        mag_mul_2exp_si(u, bm, 1);
        mag_mul(u, u, y);
        mag_add(u, u, v);
        mag_mul(u, u, am);

        /* v = |b| ((x^2 + y^2) + 2 |a| x) */
        mag_mul_2exp_si(w, am, 1);
        mag_addmul(v, w, x);
        mag_mul(v, v, bm);

        /* w = |b^2 - a^2| (upper bound) */
        if (arf_cmpabs(a, b) >= 0)
            mag_mul(w, am, am);
        else
            mag_mul(w, bm, bm);

        mag_addmul(u, w, x);
        mag_addmul(v, w, y);

        mag_div(arb_radref(acb_realref(res)), u, t);
        mag_div(arb_radref(acb_imagref(res)), v, t);

        _arb_arf_div_rounded_den_add_err(acb_realref(res), a, a2b2, inexact, prec);
        _arb_arf_div_rounded_den_add_err(acb_imagref(res), b, a2b2, inexact, prec);
        arf_neg(arb_midref(acb_imagref(res)), arb_midref(acb_imagref(res)));

        mag_clear(t);
        mag_clear(u);
        mag_clear(v);
        mag_clear(w);

        arf_clear(a2b2);
    }

    mag_clear(am);
    mag_clear(bm);
#undef a
#undef b
#undef x
#undef y
}
