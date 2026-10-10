/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "elem.h"

/*
    Arctangents.

    atan(x) for |x| < 1 is the kernel on |x| at enough limbs for the
    relative accuracy (the series take over for small arguments); for
    |x| > 1, pi/2 - atan(1/|x|) with one approximate reciprocal.

    atan2(y, x) evaluates atan(t) for the ratio t = min(|x|, |y|) / max(|x|, |y|)
    of the mantissas (one division, with relative error bounds) and adds
    the octant: theta = atan(t), pi/2 - atan(t), pi/2 + atan(t) or pi - atan(t),
    negated for y < 0. Only theta = atan(t) can be small, when x > 0 and
    |y| << x; then the kernel runs at as many more bits as t has leading
    zeros.
*/

/* res = (-1)^sgnbit (c pi/4 + sigma atan(v)) for v at w fraction limbs
   (within ev ulps, v < 1), c in {0, 1, 2, 3, 4}, sigma = +-1; the result
   must not cancel (c pi/4 >= 2 atan(v) when sigma = -1) */
static int
_nfloat_atan_finish(nfloat_ptr res, nn_srcptr v, slong w, ulong ev, int c, int sigma, int sgnbit, gr_ctx_t ctx)
{
    ulong y[NFLOAT_ELEM_MAX_LIMBS + 3];
    ulong T[NFLOAT_ELEM_MAX_LIMBS + 3];
    ulong err;

    /* atan' <= 1 */
    _mp_real_atan_kernel(y, &err, v, w);
    err += ev;

    if (c == 0)
    {
        y[w] = 0;
        return _nfloat_set_mpn_err(res, y, w + 1, -FLINT_BITS * w, err, err, sgnbit, ctx);
    }

    /* T = c floor(pi/4 B^w) +- atan(v), w + 1 limbs; pi/4 within one ulp
       below, c times */
    T[w] = mpn_mul_1(T, _mp_real_const_ptr(MP_REAL_CONST_ID_PI4, w), w, c);
    if (sigma > 0)
        T[w] += mpn_add_n(T, T, y, w);
    else
        T[w] -= mpn_sub_n(T, T, y, w);
    err += c;

    return _nfloat_set_mpn_err(res, T, w + 1, -FLINT_BITS * w, err, err, sgnbit, ctx);
}

#if FLINT_BITS == 64
/*
    In-register fast paths at one and two limbs (64-bit only): for
    2^-32 <= |x| < 1 (resp. 2^-40) the table-driven kernel of mp_real at
    two resp. three limbs on |x| exactly; for 1 < |x| < 2^63, pi/2 -
    atan(v) with v = 1/|x| = 2^-e / M from the quotient Q = floor(B^(w+1) / D)
    (1/M at w = 2 resp. 3 fraction limbs by hardware divisions) shifted
    right by e bits, below the true value by under two ulps of B^-w, and
    pi/2 as twice the floor of pi/4 at w limbs (two ulps below). The
    fast paths decline (return NFLOAT_ATAN_FAST_DECLINE) for other
    arguments and for |x| = 1.
*/

#define NFLOAT_ATAN_FAST_DECLINE (-1)

static int
_nfloat_atan_fast_1(nfloat_ptr res, ulong d, slong e, int sgnbit, gr_ctx_t ctx)
{
    ulong v1, v0, y1, y0, q1, q0, qi, r, t2, t1, t0;
    nn_srcptr P;

    if (e <= 0)
    {
        if (e < -32)
            return NFLOAT_ATAN_FAST_DECLINE;
        v1 = (e == 0) ? d : (d >> (-e));
        v0 = (e == 0) ? 0 : (d << (FLINT_BITS + e));
        _mp_real_small_atan_2(&y1, &y0, v1, v0);
        return _nfloat_small_set_1(res, 0, y1, y0, 0, MP_REAL_SMALL_ATAN_ERR_2, sgnbit, ctx);
    }

    if (e >= FLINT_BITS)
        return NFLOAT_ATAN_FAST_DECLINE;

    /* (qi, q1, q0) = floor(B^3 / d), i.e. 1/M at two fraction limbs */
    if (d == (UWORD(1) << (FLINT_BITS - 1)))
    {
        if (e == 1)
            return NFLOAT_ATAN_FAST_DECLINE;    /* |x| = 1 */
        qi = 2; q1 = q0 = 0;
    }
    else
    {
        qi = 1;
        r = -d;
        udiv_qrnnd(q1, r, r, UWORD(0), d);
        udiv_qrnnd(q0, r, r, UWORD(0), d);
        (void) r;
    }

    /* v = Q 2^-e < 1 (1 <= e < FLINT_BITS) */
    v1 = (q1 >> e) | (qi << (FLINT_BITS - e));
    v0 = (q0 >> e) | (q1 << (FLINT_BITS - e));

    _mp_real_small_atan_2(&y1, &y0, v1, v0);

    /* 2 P - atan(v) */
    P = _mp_real_const_ptr(MP_REAL_CONST_ID_PI4, 2);
    t2 = P[1] >> (FLINT_BITS - 1);
    t1 = (P[1] << 1) | (P[0] >> (FLINT_BITS - 1));
    t0 = P[0] << 1;
    sub_dddmmmsss(t2, t1, t0, t2, t1, t0, UWORD(0), y1, y0);

    return _nfloat_small_set_1(res, t2, t1, t0, 0, MP_REAL_SMALL_ATAN_ERR_2 + 4, sgnbit, ctx);
}

static int
_nfloat_atan_fast_2(nfloat_ptr res, ulong d1, ulong d0, slong e, int sgnbit, gr_ctx_t ctx)
{
    ulong v2, v1, v0, y2, y1, y0, q2, q1, q0, qi, r1, r0, t3, t2, t1, t0, dinv;
    nn_srcptr P;

    if (e <= 0)
    {
        if (e < -40)
            return NFLOAT_ATAN_FAST_DECLINE;
        if (e == 0)
        {
            v2 = d1; v1 = d0; v0 = 0;
        }
        else
        {
            unsigned int sh = -e;
            v2 = d1 >> sh;
            v1 = (d1 << (FLINT_BITS - sh)) | (d0 >> sh);
            v0 = d0 << (FLINT_BITS - sh);
        }
        _mp_real_small_atan_3(&y2, &y1, &y0, v2, v1, v0);
        return _nfloat_small_set_2(res, 0, y2, y1, y0, 0, MP_REAL_SMALL_ATAN_ERR_3, sgnbit, ctx);
    }

    if (e >= FLINT_BITS)
        return NFLOAT_ATAN_FAST_DECLINE;

    /* (qi, q2, q1, q0) = floor(B^5 / D), 1/M at three fraction limbs */
    if (d1 == (UWORD(1) << (FLINT_BITS - 1)) && d0 == 0)
    {
        if (e == 1)
            return NFLOAT_ATAN_FAST_DECLINE;    /* |x| = 1 */
        qi = 2; q2 = q1 = q0 = 0;
    }
    else
    {
        /* (B^2 - D) / D < 1 */
        qi = 1;
        sub_ddmmss(r1, r0, UWORD(0), UWORD(0), d1, d0);
        dinv = _MP_REAL_PREINV1(d1, d0);
        (void) dinv;
        _MP_REAL_UDIV_QR_3BY2(q2, r1, r0, r1, r0, UWORD(0), d1, d0, dinv);
        _MP_REAL_UDIV_QR_3BY2(q1, r1, r0, r1, r0, UWORD(0), d1, d0, dinv);
        _MP_REAL_UDIV_QR_3BY2(q0, r1, r0, r1, r0, UWORD(0), d1, d0, dinv);
        (void) r1;
        (void) r0;
    }

    v2 = (q2 >> e) | (qi << (FLINT_BITS - e));
    v1 = (q1 >> e) | (q2 << (FLINT_BITS - e));
    v0 = (q0 >> e) | (q1 << (FLINT_BITS - e));

    _mp_real_small_atan_3(&y2, &y1, &y0, v2, v1, v0);

    P = _mp_real_const_ptr(MP_REAL_CONST_ID_PI4, 3);
    t3 = P[2] >> (FLINT_BITS - 1);
    t2 = (P[2] << 1) | (P[1] >> (FLINT_BITS - 1));
    t1 = (P[1] << 1) | (P[0] >> (FLINT_BITS - 1));
    t0 = P[0] << 1;
    sub_ddddmmmmssss(t3, t2, t1, t0, t3, t2, t1, t0, UWORD(0), y2, y1, y0);

    return _nfloat_small_set_2(res, t3, t2, t1, t0, 0, MP_REAL_SMALL_ATAN_ERR_3 + 4, sgnbit, ctx);
}
#endif

int
nfloat_atan(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    ulong v[NFLOAT_ELEM_MAX_LIMBS + 2];
    slong n, e, w;
    int sgnbit;

    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_zero(res, ctx);
        if (NFLOAT_IS_NAN(x))
            return nfloat_nan(res, ctx);
        /* +-pi/2 */
        {
            ulong t[NFLOAT_ELEM_MAX_LIMBS + 2];
            w = _nfloat_elem_limbs(ctx, 0);
            t[w] = mpn_lshift(t, _mp_real_const_ptr(MP_REAL_CONST_ID_PI4, w), w, 1);
            return _nfloat_set_mpn_err(res, t, w + 1, -FLINT_BITS * w, 2, 2,
                NFLOAT_IS_NEG_INF(x), ctx);
        }
    }

    n = NFLOAT_CTX_NLIMBS(ctx);
    e = NFLOAT_EXP(x);
    sgnbit = NFLOAT_SGNBIT(x);

    /* atan(x) = x (1 - delta), 0 < delta < x^2 / 3 */
    if (2 * e <= -FLINT_BITS * n)
        return _nfloat_set_tiny_rel_err(res, x, -1, ctx);

#if FLINT_BITS == 64
    if (n <= 2)
    {
        int status = (n == 1) ? _nfloat_atan_fast_1(res, NFLOAT_D(x)[0], e, sgnbit, ctx)
                              : _nfloat_atan_fast_2(res, NFLOAT_D(x)[1], NFLOAT_D(x)[0], e, sgnbit, ctx);
        if (status != NFLOAT_ATAN_FAST_DECLINE)
            return status;
    }
#endif

    if (e <= 0)
    {
        int trunc;
        w = _nfloat_elem_limbs(ctx, -e);
        trunc = _nfloat_get_fixed(v, w, w, NFLOAT_D(x), n, e);
        return _nfloat_atan_finish(res, v, w, trunc, 0, 1, sgnbit, ctx);
    }

    w = _nfloat_elem_limbs(ctx, 0);

    if (e == 1 && NFLOAT_D(x)[n - 1] == (UWORD(1) << (FLINT_BITS - 1)) &&
        flint_mpn_zero_p(NFLOAT_D(x), n - 1))
    {
        /* |x| = 1: pi/4 */
        _nfloat_copy_limbs(v, _mp_real_const_ptr(MP_REAL_CONST_ID_PI4, w), w);
        return _nfloat_set_mpn_err(res, v, w, -FLINT_BITS * w, 1, 1, sgnbit, ctx);
    }

    {
        /* |x| > 1: v = 1/|x| = B^w 2^-e / M < 1: R = floor(B^(w + n) / D)
           or one more, in (B^w, 2 B^w], shifted right by e bits (the kernel
           costs the same for all v < 1, so this beats a reduction to
           (|x| - 1) / (|x| + 1) on [1, 2)) */
        ulong R[NFLOAT_ELEM_MAX_LIMBS + 3];
        slong ls;
        unsigned int bs;
        ulong ev;

        if (e > FLINT_BITS * (w + 1))
        {
            /* v below one ulp: pi/2 - v with 0 < v < 2^(1-e) */
            _nfloat_zero_limbs(v, w);
            return _nfloat_atan_finish(res, v, w, 1, 2, -1, sgnbit, ctx);
        }

        flint_mpn_invapprox(R, NFLOAT_D(x), n, w + n);
        /* (R, w + 2), the top limb zero */
        ls = e / FLINT_BITS;
        bs = e % FLINT_BITS;
        ev = 2;
        if (ls >= w + 1)
            _nfloat_zero_limbs(v, w);
        else
        {
            if (bs != 0)
                mpn_rshift(R + ls, R + ls, w + 2 - ls, bs);
            _nfloat_copy_limbs(v, R + ls, w + 1 - ls);
            _nfloat_zero_limbs(v + w + 1 - ls, ls);
        }

        if (ls < w + 1 && v[w] != 0)
        {
            /* v = 1 by the rounding of R (e = 1): one ulp below */
            FLINT_ASSERT(e == 1);
            for (ls = 0; ls < w; ls++)
                v[ls] = UWORD_MAX;
            ev++;
        }

        return _nfloat_atan_finish(res, v, w, ev, 2, -1, sgnbit, ctx);
    }
}

/* the magnitude comparison |a| < |b| for normal nfloats */
static int
_nfloat_cmpabs_normal(nfloat_srcptr a, nfloat_srcptr b, slong n)
{
    if (NFLOAT_EXP(a) != NFLOAT_EXP(b))
        return (NFLOAT_EXP(a) < NFLOAT_EXP(b)) ? -1 : 1;
    return mpn_cmp(NFLOAT_D(a), NFLOAT_D(b), n);
}

int
nfloat_atan2(nfloat_ptr res, nfloat_srcptr y, nfloat_srcptr x, gr_ctx_t ctx)
{
    ulong v[NFLOAT_ELEM_MAX_LIMBS + 2];
    ulong a[NFLOAT_ELEM_MAX_LIMBS + 1], b[NFLOAT_ELEM_MAX_LIMBS + 1];
    ulong q[NFLOAT_ELEM_MAX_LIMBS + 3];
    nfloat_srcptr num, den;
    slong n, w, wn, z, d;
    int c, sigma, ysgn, cmp;
    ulong eq;

    if (NFLOAT_IS_SPECIAL(x) || NFLOAT_IS_SPECIAL(y))
    {
        if (NFLOAT_IS_NAN(x) || NFLOAT_IS_NAN(y))
            return nfloat_nan(res, ctx);

        if (NFLOAT_IS_ZERO(y))
        {
            /* 0 for x >= 0, pi for x < 0 */
            if (NFLOAT_IS_ZERO(x) || NFLOAT_IS_POS_INF(x) ||
                (!NFLOAT_IS_SPECIAL(x) && !NFLOAT_SGNBIT(x)))
                return nfloat_zero(res, ctx);
            w = _nfloat_elem_limbs(ctx, 0);
            v[w] = mpn_lshift(v, _mp_real_const_ptr(MP_REAL_CONST_ID_PI4, w), w, 2);
            return _nfloat_set_mpn_err(res, v, w + 1, -FLINT_BITS * w, 4, 4, 0, ctx);
        }

        if (NFLOAT_IS_ZERO(x) && !NFLOAT_IS_SPECIAL(y))
        {
            /* +-pi/2 */
            w = _nfloat_elem_limbs(ctx, 0);
            v[w] = mpn_lshift(v, _mp_real_const_ptr(MP_REAL_CONST_ID_PI4, w), w, 1);
            return _nfloat_set_mpn_err(res, v, w + 1, -FLINT_BITS * w, 2, 2, NFLOAT_SGNBIT(y), ctx);
        }

        /* infinities */
        return nfloat_nan(res, ctx);
    }

    n = NFLOAT_CTX_NLIMBS(ctx);
    ysgn = NFLOAT_SGNBIT(y);
    cmp = _nfloat_cmpabs_normal(y, x, n);

    if (cmp == 0)
    {
        /* |y| = |x|: pi/4 or 3 pi/4 */
        w = _nfloat_elem_limbs(ctx, 0);
        _nfloat_zero_limbs(v, w);
        return _nfloat_atan_finish(res, v, w, 0, NFLOAT_SGNBIT(x) ? 3 : 1, 1, ysgn, ctx);
    }

    /* theta = c pi/4 + sigma atan(t), t = num / den < 1 */
    if (cmp < 0)
    {
        num = y;
        den = x;
        c = NFLOAT_SGNBIT(x) ? 4 : 0;
        sigma = NFLOAT_SGNBIT(x) ? -1 : 1;
    }
    else
    {
        num = x;
        den = y;
        c = 2;
        sigma = NFLOAT_SGNBIT(x) ? 1 : -1;
    }

    /* t = (Mnum / Mden) 2^-d, d >= 0 */
    d = NFLOAT_EXP(den) - NFLOAT_EXP(num);

    /* the leading zero bits of t: d - 1 or d */
    z = (d >= 1) ? d - 1 : 0;

    if (c == 0)
    {
        /* tiny theta = t (1 - delta), delta < t^2 / 3 < 2^(-2z) */
        if (2 * z >= FLINT_BITS * (n + 1) + 4)
        {
            /* the ratio of the mantissas at n + 1 limbs */
            wn = n + 1;
            _nfloat_zero_limbs(a, 1);
            _nfloat_copy_limbs(a + 1, NFLOAT_D(num), n);
            _nfloat_zero_limbs(b, 1);
            _nfloat_copy_limbs(b + 1, NFLOAT_D(den), n);
            eq = _nfloat_elem_div(q, a, 0, b, 0, wn);
            /* Q = (Mnum / Mden) B^wn, t = Q 2^(-d - FLINT_BITS wn); delta t
               is below one unit */
            return _nfloat_set_mpn_err(res, q, wn + 1,
                -d - FLINT_BITS * wn, eq + 1, eq, ysgn, ctx);
        }

        w = _nfloat_elem_limbs(ctx, z);
    }
    else
    {
        w = _nfloat_elem_limbs(ctx, 0);

        /* t below one ulp: theta = c pi/4 + sigma t with t < 2^-z */
        if (z >= FLINT_BITS * w + 2)
        {
            _nfloat_zero_limbs(v, w);
            return _nfloat_atan_finish(res, v, w, 1, c, sigma, ysgn, ctx);
        }
    }

    /* the mantissa ratio at wn >= w limbs, then t at w fraction limbs
       (wn < w would scale up the error count of the ratio) */
    wn = FLINT_MAX(_nfloat_elem_limbs(ctx, 0), w);
    if (wn >= n)
    {
        _nfloat_zero_limbs(a, wn - n);
        _nfloat_copy_limbs(a + wn - n, NFLOAT_D(num), n);
        _nfloat_zero_limbs(b, wn - n);
        _nfloat_copy_limbs(b + wn - n, NFLOAT_D(den), n);
        eq = _nfloat_elem_div(q, a, 0, b, 0, wn);
    }
    else
    {
        /* truncated mantissas (a reduced function precision) */
        eq = _nfloat_elem_div(q, NFLOAT_D(num) + n - wn, 1, NFLOAT_D(den) + n - wn, 1, wn);
    }

    /* Q = (Mnum / Mden) B^wn in (B^wn / 2, 2 B^wn), t = Q B^-wn 2^-d < 1:
       as a fraction at w limbs, shifted right by FLINT_BITS (wn - w) + d bits
       (one more ulp for the shifted-out bits, the error scaled down) */
    {
        slong sh = FLINT_BITS * (wn - w) + d;
        ulong ev;

        slong i;

        /* (v, w + 1) = floor(Q 2^-sh) */
        _nfloat_zero_limbs(v, w + 1);

        if (sh >= 0)
        {
            slong ls = sh / FLINT_BITS;
            unsigned int bs = sh % FLINT_BITS;

            for (i = 0; i <= w && ls + i <= wn; i++)
            {
                v[i] = q[ls + i] >> bs;
                if (bs != 0 && ls + i + 1 <= wn)
                    v[i] |= q[ls + i + 1] << (FLINT_BITS - bs);
            }

            ev = (sh >= FLINT_BITS) ? 2 : (eq >> sh) + 2;
        }
        else
        {
            slong ls = (-sh) / FLINT_BITS;
            unsigned int bs = (-sh) % FLINT_BITS;

            for (i = 0; i <= wn && ls + i <= w; i++)
            {
                v[ls + i] |= q[i] << bs;
                if (bs != 0 && ls + i + 1 <= w)
                    v[ls + i + 1] |= q[i] >> (FLINT_BITS - bs);
            }

            if (-sh >= FLINT_BITS - 4 || eq > (UWORD_MAX >> (-sh + 2)))
                return GR_UNABLE;
            ev = (eq << (-sh)) + 1;
        }

        /* t < 1 exactly; the approximation may reach 1 */
        if (v[w] != 0)
        {
            for (i = 0; i < w; i++)
                v[i] = UWORD_MAX;
            ev++;
        }

        return _nfloat_atan_finish(res, v, w, ev, c, sigma, ysgn, ctx);
    }
}
