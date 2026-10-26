/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ARB_TERMS_H
#define ARB_TERMS_H

/*
    Sums of a few exact or truncated products, rounded once with the
    measured error: the machinery of the one-pass complex multiplication
    (acb/mul_onepass.c), shared with arb_fmma / arb_sosq (arb/fmma.c).
    Internal; include after arb/impl.h.

    A part is a sum of terms (-1)^sgn (p / B^n) 2^exp: exact products,
    or high products truncated to the working precision plus a guard limb
    whose deficit is charged to the radius. The terms are added exactly
    (in registers for short operands and precision, otherwise in a limb
    frame) and rounded to nearest once, the rounding error measured into
    the radius as in arb_add and arb_mul. Propagated errors are formed by
    _arb_mag_sum6 as one aligned sum of mag products.
*/

#include <stdint.h>
#include "mpn_extras.h"
#include "arb.h"

/* a term (-1)^sgn (p / B^n) 2^exp, p with a nonzero top limb having at
   most one leading zero bit; a lower approximation of the exact value with
   a deficit below 2^e units of p[0] (e < 0: exact) */
typedef struct
{
    nn_srcptr p;
    slong n;
    slong exp;
    slong e;
    int sgn;
}
_arb_term_struct;

/* exact product of {xp, xn} and {yp, yn} (top limbs normalised, so the
   product has at most one leading zero bit) into tmp (xn + yn limbs) */
FLINT_FORCE_INLINE void
_arb_term_exact(_arb_term_struct * T, nn_ptr tmp, nn_srcptr xp, slong xn, nn_srcptr yp, slong yn, slong exp, int sgn)
{
    slong n = xn + yn;

    if (xn == yn)
    {
        if (xn == 1)
            umul_ppmm(tmp[1], tmp[0], xp[0], yp[0]);
        else if (xp == yp)
            flint_mpn_sqr(tmp, xp, xn);
        else if (FLINT_MUL_USE_FUNC_TAB && FLINT_HAVE_MUL_N_FUNC(xn))
            FLINT_MPN_MUL_N_HARD(tmp, xp, yp, xn);
        else
            flint_mpn_mul_n(tmp, xp, yp, xn);
    }
    else if (xn == 2 && yn == 1)
        FLINT_MPN_MUL_2X1(tmp[2], tmp[1], tmp[0], xp[1], xp[0], yp[0]);
    else if (xn == 1 && yn == 2)
        FLINT_MPN_MUL_2X1(tmp[2], tmp[1], tmp[0], yp[1], yp[0], xp[0]);
    else if (xn > yn)
        flint_mpn_mul(tmp, xp, xn, yp, yn);
    else
        flint_mpn_mul(tmp, yp, yn, xp, xn);

    /* the top limb is nonzero (at least 2^(FLINT_BITS - 2)) */
    T->p = tmp;
    T->n = n;
    T->exp = exp;
    T->e = -1;
    T->sgn = sgn;
}

/* the top outn limbs of the product with a guard limb, normalised, into
   tmp (outn + 2 limbs) */
FLINT_FORCE_INLINE void
_arb_term_high(_arb_term_struct * T, nn_ptr tmp, nn_srcptr xp, slong xn, nn_srcptr yp, slong yn, slong outn, slong exp, int sgn)
{
    mp_limb_pair_t hi;

    if (xp == yp && xn == yn)
        hi = flint_mpn_sqrhigh_normalised(tmp + 1, xp, xn, outn);
    else
        hi = flint_mpn_mulhigh_normalised(tmp + 1, xp, xn, yp, yn, outn);

    tmp[0] = hi.m1;

    T->p = tmp;
    T->n = outn + 1;
    T->exp = exp - (slong) hi.m2;
    T->e = FLINT_BIT_COUNT(FLINT_MPN_MULHIGH_NORMALISED_GUARD_ERR_BOUND(outn));
    T->sgn = sgn;
}

/* whether the product of xn and yn limbs, of which the top outn are
   wanted, is to be computed in full (cheaper, or exact for free) */
FLINT_FORCE_INLINE int
_arb_product_full(slong xn, slong yn, slong outn)
{
    return xn + yn <= ARF_MUL_SMALL_LIMBS || outn >= xn + yn || _flint_mpn_mulhigh_prefer_full(xn, yn, outn);
}

/* whether the products x1 y1 and x2 y2 (truncated, with exponents within
   two bits of each other) agree to about 30 bits at the top, so that
   their difference would be dominated by the deficits */
static inline int
_arb_products_cancel(nn_srcptr x1, nn_srcptr y1, slong exp1, nn_srcptr x2, nn_srcptr y2, slong exp2)
{
    slong shift = exp1 - exp2;
    ulong h1, l1, h2, l2;

    if (shift < -2 || shift > 2)
        return 0;

    umul_ppmm(h1, l1, x1[0], y1[0]);
    umul_ppmm(h2, l2, x2[0], y2[0]);

    if (shift > 0) h2 >>= shift; else h1 >>= -shift;
    h1 = (h1 > h2) ? h1 - h2 : h2 - h1;
    return h1 < (UWORD(1) << (FLINT_BITS - 30));
}

/* the deficit of a truncated term as a mag */
FLINT_FORCE_INLINE void
_arb_term_deficit(mag_t d, const _arb_term_struct * T)
{
    MAG_MAN(d) = MAG_ONE_HALF;
    MAG_EXP(d) = T->exp - T->n * FLINT_BITS + T->e + 1;
}

/* res = (-1)^sgn {p, n} 2^(exp - n FLINT_BITS) rounded (n = 0: zero)
   with err (NULL for none) added to the radius */
FLINT_FORCE_INLINE void
_arb_round_term(arb_t res, nn_srcptr p, slong n, int sgn, slong exp, slong prec, mag_srcptr err)
{
    if (n == 0)
    {
        arf_zero(arb_midref(res));
        if (err == NULL)
            mag_zero(arb_radref(res));
        else
            mag_set(arb_radref(res), err);
    }
    else
        _arb_set_round_mpn_2exp_si_long(res, p, n, sgn, exp - n * FLINT_BITS, prec, err, NULL);
}

/* res = the single term T rounded, with its deficit (if any) and err */
static inline void
_arb_round_term1(arb_t res, const _arb_term_struct * T, slong prec, mag_srcptr err)
{
    mag_t tot;

    if (T->e >= 0)
    {
        _arb_term_deficit(tot, T);
        if (err != NULL)
            mag_fast_add(tot, tot, err);
        err = tot;
    }

    _arb_round_term(res, T->p, T->n, T->sgn, T->exp, prec, err);
}

/* frame limbs for the exact sums, on the stack up to this many */
#define ARB_SUM_FRAME_STACK 48

/*
    res = X + Y rounded to prec bits, radius = err (may be NULL) + the
    deficits + the rounding error, with prec lowered to the accuracy that
    err leaves.
*/
static inline void
_arb_round_sum2(arb_t res, const _arb_term_struct * X, const _arb_term_struct * Y, slong prec, mag_srcptr err)
{
    mag_t tot;
    slong shift, zn, zn_original, need;
    int sgn;
    ulong frame_stack[ARB_SUM_FRAME_STACK];
    nn_ptr frame;

    if (X->exp < Y->exp)
    {
        const _arb_term_struct * t = X;
        X = Y;
        Y = t;
    }

    if (err != NULL)
    {
        /* |X + Y| < 2^(X->exp + 1) */
        slong acc = X->exp + 1 - MAG_EXP(err);
        if (acc < prec - MAG_BITS)
            prec = FLINT_MAX(2, acc + MAG_BITS);
    }

    if (X->e >= 0)
    {
        _arb_term_deficit(tot, X);
        if (err != NULL)
            mag_fast_add(tot, tot, err);
        err = tot;
    }

    shift = X->exp - Y->exp;

    if (shift > prec + FLINT_BITS + 1 && shift > (X->n + 1) * FLINT_BITS)
    {
        /* X +/- eps: Y (including its deficit) lies more than FLINT_BITS
           bits below the precision */
        _arb_set_round_mpn_2exp_si_eps(res, X->p, X->n, X->sgn, X->exp - X->n * FLINT_BITS, X->sgn != Y->sgn, prec, err);
        return;
    }

    if (Y->e >= 0)
    {
        mag_t d;
        _arb_term_deficit(d, Y);
        if (err != NULL)
            mag_fast_add(tot, err, d);
        else
            *tot = *d;
        err = tot;
    }

    /* the frame on the stack, in the thread-local addition scratch, or
       allocated */
    need = _ARF_ADD_MPN_SIZE(X->n, Y->n, shift);
    frame = (need <= ARB_SUM_FRAME_STACK) ? frame_stack
          : (need <= ARF_ADD_TLS_ALLOC) ? _arf_add_tls_tmp(need) : flint_malloc(sizeof(ulong) * need);

    zn = _arf_add_mpn_raw(frame, &zn_original, &sgn, X->p, X->n, X->sgn, Y->p, Y->n, Y->sgn, shift);
    _arb_round_term(res, frame, zn, sgn, X->exp + (zn - zn_original) * FLINT_BITS, prec, err);

    if (need > ARF_ADD_TLS_ALLOC)
        flint_free(frame);
}

/*
    res = T[0] + ... + T[nterms - 1] (nterms >= 3) rounded likewise: the
    terms are added exactly in order of decreasing exponent, a term more
    than prec + 2 FLINT_BITS bits below the running sum being charged to
    the radius instead.
*/
static inline void
_arb_round_sum(arb_t res, const _arb_term_struct * T, slong nterms, slong prec, mag_srcptr err)
{
    mag_t def, tot;
    ulong frame_stack[2][ARB_SUM_FRAME_STACK];
    nn_ptr frame[2] = { frame_stack[0], frame_stack[1] };
    slong alloc[2] = { ARB_SUM_FRAME_STACK, ARB_SUM_FRAME_STACK };
    _arb_term_struct S;
    slong i, lead;
    int cur = 0, have_def = 0;

    /* leading term, and the deficits of the truncated products */
    lead = 0;
    for (i = 0; i < nterms; i++)
    {
        if (T[i].exp > T[lead].exp)
            lead = i;
        if (T[i].e >= 0)
        {
            mag_t d;
            _arb_term_deficit(d, T + i);
            if (have_def)
                mag_fast_add(def, def, d);
            else
            {
                *def = *d;
                have_def = 1;
            }
        }
    }

    if (err != NULL)
    {
        /* the sum is below 2^(exp + 2) */
        slong acc = T[lead].exp + 2 - MAG_EXP(err);
        if (acc < prec - MAG_BITS)
            prec = FLINT_MAX(2, acc + MAG_BITS);
    }

    S = T[lead];

    for (i = 0; i < nterms; i++)
    {
        slong shift, zn, zn_original, need;
        int sgn;
        nn_ptr dst;

        if (i == lead)
            continue;

        if (S.n == 0)
        {
            /* the previous terms cancelled exactly */
            S = T[i];
            continue;
        }

        shift = S.exp - T[i].exp;

        if (shift > prec + 2 * FLINT_BITS)
        {
            /* |T[i]| < 2^T[i].exp: charged to the radius */
            mag_t d;
            MAG_MAN(d) = MAG_ONE_HALF;
            MAG_EXP(d) = T[i].exp + 1;
            if (have_def)
                mag_fast_add(def, def, d);
            else
            {
                *def = *d;
                have_def = 1;
            }
            continue;
        }

        if (shift >= 0)
            need = _ARF_ADD_MPN_SIZE(S.n, T[i].n, shift);
        else
            need = _ARF_ADD_MPN_SIZE(T[i].n, S.n, -shift);

        if (need > alloc[cur])
        {
            if (frame[cur] != frame_stack[cur])
                flint_free(frame[cur]);
            frame[cur] = flint_malloc(sizeof(ulong) * need);
            alloc[cur] = need;
        }
        dst = frame[cur];

        if (shift >= 0)
        {
            zn = _arf_add_mpn_raw(dst, &zn_original, &sgn, S.p, S.n, S.sgn, T[i].p, T[i].n, T[i].sgn, shift);
            S.exp = S.exp + (zn - zn_original) * FLINT_BITS;
        }
        else
        {
            zn = _arf_add_mpn_raw(dst, &zn_original, &sgn, T[i].p, T[i].n, T[i].sgn, S.p, S.n, S.sgn, -shift);
            S.exp = T[i].exp + (zn - zn_original) * FLINT_BITS;
        }

        S.p = dst;
        S.n = zn;
        S.sgn = sgn;
        cur ^= 1;
    }

    if (have_def)
    {
        if (err != NULL)
        {
            mag_fast_add(tot, err, def);
            err = tot;
        }
        else
            err = def;
    }

    _arb_round_term(res, S.p, S.n, S.sgn, S.exp, prec, err);

    if (frame[0] != frame_stack[0])
        flint_free(frame[0]);
    if (frame[1] != frame_stack[1])
        flint_free(frame[1]);
}

/*
    res = (-1)^s0 (h0, l0) 2^(e0 - 2 FLINT_BITS) + (-1)^s1 (h1, l1) 2^(e1 -
    2 FLINT_BITS) for products of one-limb operands (h0, h1 >= 2^(FLINT_BITS
    - 2)), rounded in registers to prec <= 2 FLINT_BITS bits: the terms are
    aligned in three limbs, the leading one shifted down a bit so that the
    sum cannot carry out, the bits of the other below the three limbs being
    sticky (they only occur when it lies more than a limb below the top,
    where at most three bits can cancel).
*/
FLINT_FORCE_INLINE void
_arb_round_sum_1limb(arb_t res, ulong h0, ulong l0, slong e0, int s0,
    ulong h1, ulong l1, slong e1, int s1, slong prec, mag_srcptr err)
{
    ulong x2, x1, x0, y2, y1, y0, t2, t1, t0, sticky = 0;
    slong sh, exp3;
    int sgn;

    if (e0 < e1)
    {
        ulong t; slong e; int s;
        t = h0; h0 = h1; h1 = t;
        t = l0; l0 = l1; l1 = t;
        e = e0; e0 = e1; e1 = e;
        s = s0; s0 = s1; s1 = s;
    }

    /* the frame (x2, x1, x0) 2^exp3 holds X shifted down by one bit */
    exp3 = e0 + 1 - 3 * FLINT_BITS;
    x2 = h0 >> 1;
    x1 = (h0 << (FLINT_BITS - 1)) | (l0 >> 1);
    x0 = l0 << (FLINT_BITS - 1);

    sh = e0 - e1 + 1;
    if (sh < FLINT_BITS)
    {
        y2 = h1 >> sh;
        y1 = (h1 << (FLINT_BITS - sh)) | (l1 >> sh);
        y0 = l1 << (FLINT_BITS - sh);
    }
    else if (sh < 2 * FLINT_BITS)
    {
        sh -= FLINT_BITS;
        y2 = 0;
        if (sh == 0)
        {
            y1 = h1;
            y0 = l1;
        }
        else
        {
            y1 = h1 >> sh;
            y0 = (h1 << (FLINT_BITS - sh)) | (l1 >> sh);
            sticky = l1 << (FLINT_BITS - sh);
        }
    }
    else if (sh < 3 * FLINT_BITS)
    {
        sh -= 2 * FLINT_BITS;
        y2 = y1 = 0;
        if (sh == 0)
        {
            y0 = h1;
            sticky = l1;
        }
        else
        {
            y0 = h1 >> sh;
            sticky = (h1 << (FLINT_BITS - sh)) | (l1 != 0);
        }
    }
    else
    {
        y2 = y1 = y0 = 0;
        sticky = 1;
    }

    if (s0 == s1)
    {
        add_sssaaaaaa(t2, t1, t0, x2, x1, x0, y2, y1, y0);
        sgn = s0;
    }
    else if (x2 > y2 || (x2 == y2 && (x1 > y1 || (x1 == y1 && x0 >= y0))))
    {
        /* the sticky part of Y borrows a unit at the bottom, leaving a
           positive sticky remainder (Y then lies far below X) */
        sub_dddmmmsss(t2, t1, t0, x2, x1, x0, y2, y1, y0);
        if (sticky != 0)
            sub_dddmmmsss(t2, t1, t0, t2, t1, t0, 0, 0, 1);
        sgn = s0;

        if ((t2 | t1 | t0) == 0)
        {
            arf_zero(arb_midref(res));
            if (err == NULL)
                mag_zero(arb_radref(res));
            else
                mag_set(arb_radref(res), err);
            return;
        }
    }
    else
    {
        sub_dddmmmsss(t2, t1, t0, y2, y1, y0, x2, x1, x0);
        sgn = s1;
    }

    if (t2 == 0)
    {
        t2 = t1; t1 = t0; t0 = 0;
        exp3 -= FLINT_BITS;
        if (t2 == 0)
        {
            t2 = t1; t1 = 0;
            exp3 -= FLINT_BITS;
        }
    }

    _arb_set_round_3limbs(res, t2, t1, t0, sticky, sgn, exp3, prec, err, NULL);
}

/*
    res = T[0] + ... + T[nterms - 1] (nterms <= 3, terms of at most four
    limbs) rounded to prec <= 2 FLINT_BITS bits: the fixed-point
    accumulation of _arb_round_sum with a frame of four limbs in registers.
*/
FLINT_FORCE_INLINE void
_arb_round_acc4(arb_t res, const _arb_term_struct * T, slong nterms, slong prec, mag_srcptr err, int short_terms)
{
    mag_t tot;
    ulong s3 = 0, s2 = 0, s1 = 0, s0 = 0, serr = 0;
    ulong t[4];
    slong i, sum_exp, max_exp, zn;
    int negative;

    max_exp = T[0].exp;
    for (i = 1; i < nterms; i++)
        max_exp = FLINT_MAX(max_exp, T[i].exp);

    /* the frame (s3, s2, s1, s0) 2^(sum_exp - 4 FLINT_BITS) */
    sum_exp = max_exp + 3;

    for (i = 0; i < nterms; i++)
    {
        slong sh = sum_exp - T[i].exp, sl, sb, n = T[i].n;
        nn_srcptr p = T[i].p;
        ulong p3, p2, p1, p0, y3, y2, y1, y0, rest;

        if (sh >= 4 * FLINT_BITS)
        {
            serr++;
            continue;
        }

        /* the term top-aligned in four limbs (short_terms: at most two) */
        p3 = p[n - 1];
        p2 = (n >= 2) ? p[n - 2] : 0;
        p1 = (!short_terms && n >= 3) ? p[n - 3] : 0;
        p0 = (!short_terms && n >= 4) ? p[0] : 0;

        /* shifted right by sh = sl FLINT_BITS + sb bits (sh >= 3), the
           bits shifted out being sticky */
        sl = (ulong) sh / FLINT_BITS;
        sb = (ulong) sh % FLINT_BITS;

        if (sb == 0)
            rest = 0;
        else
        {
            rest = p0 << (FLINT_BITS - sb);
            p0 = (p1 << (FLINT_BITS - sb)) | (p0 >> sb);
            p1 = (p2 << (FLINT_BITS - sb)) | (p1 >> sb);
            p2 = (p3 << (FLINT_BITS - sb)) | (p2 >> sb);
            p3 = p3 >> sb;
        }

        switch (sl)
        {
            case 0: y3 = p3; y2 = p2; y1 = p1; y0 = p0; break;
            case 1: y3 = 0; y2 = p3; y1 = p2; y0 = p1; rest |= p0; break;
            case 2: y3 = y2 = 0; y1 = p3; y0 = p2; rest |= p1 | p0; break;
            default: y3 = y2 = y1 = 0; y0 = p3; rest |= p2 | p1 | p0; break;
        }

        serr += (rest != 0);

        if (T[i].sgn)
            sub_ddddmmmmssss(s3, s2, s1, s0, s3, s2, s1, s0, y3, y2, y1, y0);
        else
            add_ssssaaaaaaaa(s3, s2, s1, s0, s3, s2, s1, s0, y3, y2, y1, y0);
    }

    if (serr != 0)
    {
        /* serr < 2^2 units of 2^(sum_exp - 4 FLINT_BITS) */
        mag_t d;
        MAG_MAN(d) = serr << (MAG_BITS - 2);
        MAG_EXP(d) = sum_exp - 4 * FLINT_BITS + 2;
        if (err == NULL)
            *tot = *d;
        else
            mag_fast_add(tot, err, d);
        err = tot;
    }

    negative = (s3 >> (FLINT_BITS - 1));
    if (negative)
        sub_ddddmmmmssss(s3, s2, s1, s0, 0, 0, 0, 0, s3, s2, s1, s0);

    t[0] = s0; t[1] = s1; t[2] = s2; t[3] = s3;
    zn = 4;
    while (zn > 0 && t[zn - 1] == 0)
        zn--;

    if (zn == 0)
    {
        arf_zero(arb_midref(res));
        if (err == NULL)
            mag_zero(arb_radref(res));
        else
            mag_set(arb_radref(res), err);
    }
    else
        _arb_set_round_frame(res, t, zn, negative, sum_exp - 4 * FLINT_BITS, prec, err, 0);
}

/* one product of the error sum: p = man(x) man(y) < 2^(2 MAG_BITS) with
   x y <= p 2^(e - 2 MAG_BITS); a zero factor gives p = 0 and is sent far
   down in exponent so that it does not set the alignment */
#define ARB_MAG_PROD(p, e, x, y) \
    do { \
        (p) = (uint64_t) MAG_MAN(x) * MAG_MAN(y); \
        (e) = ((p) != 0) ? MAG_EXP(x) + MAG_EXP(y) : (WORD_MIN / 2); \
    } while (0)

/* an upper bound of x1 y1 + ... + x6 y6 for LAGOM mags, not all products
   zero: the products aligned to the largest exponent and summed in one
   word, each rounded up by a unit, the sum rounded up once */
FLINT_FORCE_INLINE void
_arb_mag_sum6(mag_t res, mag_srcptr x1, mag_srcptr y1, mag_srcptr x2, mag_srcptr y2,
    mag_srcptr x3, mag_srcptr y3, mag_srcptr x4, mag_srcptr y4,
    mag_srcptr x5, mag_srcptr y5, mag_srcptr x6, mag_srcptr y6)
{
    /* products of two 30-bit mantissas: 64-bit arithmetic on 32-bit limbs */
    uint64_t p1, p2, p3, p4, p5, p6, sum;
    slong e1, e2, e3, e4, e5, e6, E, bc;

    ARB_MAG_PROD(p1, e1, x1, y1);
    ARB_MAG_PROD(p2, e2, x2, y2);
    ARB_MAG_PROD(p3, e3, x3, y3);
    ARB_MAG_PROD(p4, e4, x4, y4);
    ARB_MAG_PROD(p5, e5, x5, y5);
    ARB_MAG_PROD(p6, e6, x6, y6);

    E = FLINT_MAX(FLINT_MAX(FLINT_MAX(e1, e2), FLINT_MAX(e3, e4)), FLINT_MAX(e5, e6));
    FLINT_ASSERT(E > WORD_MIN / 2);

    /* a term more than 2 MAG_BITS below contributes less than one unit
       (a zero one contributes a spurious unit, far below the result) */
    e1 = FLINT_MIN(E - e1, 63);
    e2 = FLINT_MIN(E - e2, 63);
    e3 = FLINT_MIN(E - e3, 63);
    e4 = FLINT_MIN(E - e4, 63);
    e5 = FLINT_MIN(E - e5, 63);
    e6 = FLINT_MIN(E - e6, 63);
    sum = (p1 >> e1) + (p2 >> e2) + (p3 >> e3) + (p4 >> e4) + (p5 >> e5) + (p6 >> e6) + 6;

    /* sum < 2^(2 MAG_BITS + 3); value <= sum 2^(E - 2 MAG_BITS) */
    bc = (sum >> 32) ? 32 + FLINT_BIT_COUNT((ulong) (sum >> 32)) : FLINT_BIT_COUNT((ulong) sum);
    if (bc > MAG_BITS)
    {
        sum = (sum >> (bc - MAG_BITS)) + 1;
        if (sum >> MAG_BITS)
        {
            sum = MAG_ONE_HALF;
            bc++;
        }
    }
    else
        sum <<= (MAG_BITS - bc);

    MAG_MAN(res) = (ulong) sum;
    MAG_EXP(res) = E - 2 * MAG_BITS + bc;
}

#endif
