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
    Trigonometric functions.

    For |x| >= 1, |x| is reduced modulo pi/2 in fixed point by the
    reduction of mp_real_sin_cos_bits (_mp_real_trig_reduce, on the
    mantissa realigned to the limb radix): |x| = a pi/2 + s v, s = +-1,
    v in [0, pi/4] at w fraction limbs within 3 ulps, and

        sin |x| = s sin v, cos v, -s sin v, -cos v   for a = 0, 1, 2, 3,
        cos |x| = cos v, -s sin v, -cos v, s sin v,
        tan |x| = s tan v (a even), -s / tan v (a odd).

    When the needed output is sin v or tan v and v has z leading zero
    bits, the reduction is repeated with z more bits so that the result
    keeps its relative accuracy. Arguments below 1 go to the kernels
    unreduced.

    sin(pi x), cos(pi x), tan(pi x) reduce |x| exactly by the reduction of
    mp_real_sin_cos_pi_bits (_mp_real_trig_pi_reduce): |x| = a/2 + s t
    with t in [0, 1/4], then v = pi t at enough limbs for t's leading
    zero bits.
*/

/* native reduction up to this exponent; beyond, mp_real balls */
#define NFLOAT_TRIG_MAX_EXP 65536

/* sin, cos, tan (op = 0, 1, 2) and of pi x (op = 3, 4, 5) by mp_real
   balls, for huge arguments and for more cancellation than the
   temporaries allow */
static int
_nfloat_trig_ball(mp_real_t r, slong * shift, const mp_real_t x, const mp_real_t y, int op, slong prec)
{
    (void) y;
    *shift = 0;

    switch (op)
    {
        case 0: mp_real_sin_cos_bits(r, NULL, x, prec); return 1;
        case 1: mp_real_sin_cos_bits(NULL, r, x, prec); return 1;
        case 2: return mp_real_tan_bits(r, x, prec);
        case 3: mp_real_sin_cos_pi_bits(r, NULL, x, prec); return 1;
        case 4: mp_real_sin_cos_pi_bits(NULL, r, x, prec); return 1;
        default: return mp_real_tan_pi_bits(r, x, prec);
    }
}

/* the outputs that are not NULL by the above (pi: 0 or 3) */
static int
_nfloat_trig_slow(nfloat_ptr rs, nfloat_ptr rc, nfloat_ptr rt, nfloat_srcptr x, int pi, gr_ctx_t ctx)
{
    slong min_prec = pi ? 0 : _nfloat_trig_ball_min_prec(x);
    int status = GR_SUCCESS;

    if (rs != NULL)
        status |= _nfloat_eval_ball(rs, _nfloat_trig_ball, x, NULL, pi + 0, min_prec, ctx);
    if (rc != NULL)
        status |= _nfloat_eval_ball(rc, _nfloat_trig_ball, x, NULL, pi + 1, min_prec, ctx);
    if (rt != NULL)
        status |= _nfloat_eval_ball(rt, _nfloat_trig_ball, x, NULL, pi + 2, min_prec, ctx);

    return status;
}

/* fraction limbs for a reduced argument with z leading zero bits; with
   div (the cotangent of the reduced argument), enough that the tangent
   normalizes to _nfloat_elem_limbs(ctx, 0) limbs with a left shift of at
   most NFLOAT_ELEM_MAX_NORM_SHIFT bits */
static slong
_nfloat_trig_limbs(gr_ctx_t ctx, slong z, int div)
{
    slong w = _nfloat_elem_limbs(ctx, z);

    if (div && z > NFLOAT_ELEM_MAX_NORM_SHIFT - 1)
        w = FLINT_MAX(w, _nfloat_elem_limbs(ctx, 0) + (z - NFLOAT_ELEM_MAX_NORM_SHIFT + 1 + FLINT_BITS - 1) / FLINT_BITS);

    return w;
}

/* the trigonometric functions of v (sin and cos for which = 0, tan for
   which = 1) with |x| = a pi/2 + s v (or pi |x| = ...), assembled into
   rs, rc (either may be NULL) resp. rt; v at w limbs within eps ulps */
static int
_nfloat_trig_eval(nfloat_ptr rs, nfloat_ptr rc, nfloat_ptr rt, nn_srcptr v, slong w,
    ulong eps, int a, int s, int xsgn, gr_ctx_t ctx)
{
    ulong ys[NFLOAT_ELEM_MAX_LIMBS + 2], yc[NFLOAT_ELEM_MAX_LIMBS + 2];
    ulong err;
    int status = GR_SUCCESS;

    if (rt == NULL)
    {
        int sneg = (s < 0), negs, negc;
        nfloat_ptr os = (a & 1) ? rc : rs;
        nfloat_ptr oc = (a & 1) ? rs : rc;

        _mp_real_sin_cos_kernel(ys, yc, &err, v, w);
        err += eps;

        /* sin |x| = s sin v, cos v, -s sin v, -cos v;
           cos |x| = cos v, -s sin v, -cos v, s sin v   (a = 0, 1, 2, 3) */
        if (a & 1)
        {
            negs = (a == 3);
            negc = sneg ^ (a == 1);
        }
        else
        {
            negs = sneg ^ (a == 2);
            negc = (a == 2);
        }

        /* os gets sin v, oc gets cos v; the sine of x is odd */
        if (os != NULL)
            status |= _nfloat_set_mpn_err(os, ys, w + 1, -FLINT_BITS * w, err, err,
                ((a & 1) ? negc : negs) ^ ((os == rs) ? xsgn : 0), ctx);
        if (oc != NULL)
            status |= _nfloat_set_mpn_err(oc, yc, w + 1, -FLINT_BITS * w, err, err,
                ((a & 1) ? negs : negc) ^ ((oc == rs) ? xsgn : 0), ctx);
        return status;
    }
    else
    {
        int neg;

        _mp_real_tan_kernel(ys, &err, v, w);
        err += 2 * eps;     /* tan' <= 2 on [0, pi/4] */

        if (!(a & 1))
        {
            /* s tan v */
            neg = (s < 0) ^ xsgn;
            return _nfloat_set_mpn_err(rt, ys, w + 1, -FLINT_BITS * w, err, err, neg, ctx);
        }
        else
        {
            /* -s / tan v: tan v normalized to wn limbs, then one division */
            ulong t[NFLOAT_ELEM_MAX_LIMBS + 1], one[NFLOAT_ELEM_MAX_LIMBS + 1];
            ulong q[NFLOAT_ELEM_MAX_LIMBS + 2];
            slong te, wn = _nfloat_elem_limbs(ctx, 0);
            ulong et, eq;

            neg = (s > 0) ^ xsgn;

            if (flint_mpn_zero_p(ys, w + 1))
                return GR_UNABLE;

            et = _nfloat_elem_normalize(t, &te, ys, w + 1, -FLINT_BITS * w, err, wn);
            if (et == UWORD_MAX)
                return GR_UNABLE;
            _nfloat_zero_limbs(one, wn);
            one[wn - 1] = UWORD(1) << (FLINT_BITS - 1);
            eq = _nfloat_elem_div(q, one, 0, t, et, wn);

            /* 1 / (T 2^te) = Q 2^(1 - 2 FLINT_BITS wn - te) */
            return _nfloat_set_mpn_err(rt, q, wn + 1, 1 - 2 * FLINT_BITS * wn - te, eq, eq, neg, ctx);
        }
    }
}

#if FLINT_BITS == 64
/*
    In-register fast paths for sin and cos at one and two limbs (64-bit
    only) for 2^-32 <= |x| < 2^32 (resp. 2^-40): |x| < 1 goes to the
    table-driven kernel of mp_real at two resp. three limbs unreduced
    (exactly), larger |x| is reduced as |x| = a pi/2 + s v with q from a
    one-limb 2/pi (one low at most, corrected once) and 2q times the
    floor P of pi/4 at three resp. four fraction limbs (an error of
    2(q + 1) units of B^-3 resp. B^-4), folded to v in [0, pi/4]. The
    kernel input drops the last limb of v (one ulp; sin', cos' <= 1).
    The fast paths decline (return NFLOAT_TRIG_FAST_DECLINE) when v <
    2^-32 resp. 2^-40, where sin v would lose too much relative
    accuracy.
*/

#define NFLOAT_TRIG_FAST_DECLINE (-1)

/* the outputs from sin v, cos v as fixed-point numbers Y 2^(-FLINT_BITS w)
   by the assignment of _nfloat_trig_eval */
#define NFLOAT_TRIG_FAST_SIGNS(negs_out, negc_out, a, sneg, xsgn) \
    do { \
        int __negs, __negc; \
        if ((a) & 1) { __negs = ((a) == 3); __negc = (sneg) ^ ((a) == 1); } \
        else { __negs = (sneg) ^ ((a) == 2); __negc = ((a) == 2); } \
        /* sin |x| from sin v (a even) or cos v (a odd), and conversely */ \
        negs_out = ((a) & 1) ? __negc : __negs; \
        negc_out = ((a) & 1) ? __negs : __negc; \
    } while (0)

static int
_nfloat_sin_cos_fast_1(nfloat_ptr rs, nfloat_ptr rc, ulong d, slong e, int xsgn, gr_ctx_t ctx)
{
    ulong v1, v0, s1, s0, g1, g0, c2, c1, c0;
    int a = 0, sneg = 0, nsv, ncv, status = GR_SUCCESS;
    nfloat_ptr os, oc;

    if (e <= 0)
    {
        if (e < -32)
            return NFLOAT_TRIG_FAST_DECLINE;
        v1 = (e == 0) ? d : (d >> (-e));
        v0 = (e == 0) ? 0 : (d << (FLINT_BITS + e));
    }
    else
    {
        nn_srcptr P = _mp_real_const_ptr(MP_REAL_CONST_ID_PI4, 3);
        ulong K = _mp_real_const_ptr(MP_REAL_CONST_ID_2_DIV_PI, 1)[0];
        ulong xi, x1, q, q2, h, l, h2, l2, p3, p2, p1, p0, t3, t2, t1, t0, u3, u2, u1, u0;

        /* |x| = (xi, x1) with the binary point after xi; 1 <= e <= 32 */
        xi = d >> (FLINT_BITS - e);
        x1 = d << e;

        /* q = floor(|x| 2/pi), or one less */
        umul_ppmm(h, l, x1, K);
        umul_ppmm(h2, l2, xi, K);
        add_ssaaaa(q, l2, h2, l2, UWORD(0), h);
        (void) l;

        /* t = |x| - 2q P >= 0 */
        q2 = 2 * q;
        umul_ppmm(p1, p0, q2, P[0]);
        umul_ppmm(h, l, q2, P[1]);
        add_ssaaaa(p2, p1, h, l, UWORD(0), p1);
        umul_ppmm(h, l, q2, P[2]);
        add_ssaaaa(p3, p2, h, l, UWORD(0), p2);
        sub_ddddmmmmssss(t3, t2, t1, t0, xi, x1, UWORD(0), UWORD(0), p3, p2, p1, p0);

        /* u = 2P (pi/2) */
        u3 = P[2] >> (FLINT_BITS - 1);
        u2 = (P[2] << 1) | (P[1] >> (FLINT_BITS - 1));
        u1 = (P[1] << 1) | (P[0] >> (FLINT_BITS - 1));
        u0 = P[0] << 1;

        if (t3 > u3 || (t3 == u3 && (t2 > u2 || (t2 == u2 && (t1 > u1 || (t1 == u1 && t0 >= u0))))))
        {
            sub_ddddmmmmssss(t3, t2, t1, t0, t3, t2, t1, t0, u3, u2, u1, u0);
            q++;
        }

        /* fold t > pi/4 */
        if (t3 != 0 || t2 > P[2] || (t2 == P[2] && (t1 > P[1] || (t1 == P[1] && t0 > P[0]))))
        {
            sub_ddddmmmmssss(t3, t2, t1, t0, u3, u2, u1, u0, t3, t2, t1, t0);
            q++;
            sneg = 1;
        }

        if (t2 < (UWORD(1) << 32))
            return NFLOAT_TRIG_FAST_DECLINE;

        a = q & 3;
        v1 = t2;
        v0 = t1;
    }

    _mp_real_small_sin_cos_g_2(&s1, &s0, &g1, &g0, v1, v0);

    NFLOAT_TRIG_FAST_SIGNS(nsv, ncv, a, sneg, xsgn);
    os = (a & 1) ? rc : rs;
    oc = (a & 1) ? rs : rc;

    if (os != NULL)
        status |= _nfloat_small_set_1(os, 0, s1, s0, 0, MP_REAL_SMALL_SIN_COS_ERR_2 + 2,
            nsv ^ ((os == rs) ? xsgn : 0), ctx);
    if (oc != NULL)
    {
        sub_dddmmmsss(c2, c1, c0, UWORD(1), UWORD(0), UWORD(0), UWORD(0), g1, g0);
        status |= _nfloat_small_set_1(oc, c2, c1, c0, 0, MP_REAL_SMALL_SIN_COS_ERR_2 + 2,
            ncv ^ ((oc == rs) ? xsgn : 0), ctx);
    }

    return status;
}

static int
_nfloat_sin_cos_fast_2(nfloat_ptr rs, nfloat_ptr rc, ulong d1, ulong d0, slong e, int xsgn, gr_ctx_t ctx)
{
    ulong v2, v1, v0, s2, s1, s0, g2, g1, g0, c3, c2, c1, c0;
    int a = 0, sneg = 0, nsv, ncv, status = GR_SUCCESS;
    nfloat_ptr os, oc;

    if (e <= 0)
    {
        if (e < -40)
            return NFLOAT_TRIG_FAST_DECLINE;
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
    }
    else
    {
        nn_srcptr P = _mp_real_const_ptr(MP_REAL_CONST_ID_PI4, 4);
        ulong K = _mp_real_const_ptr(MP_REAL_CONST_ID_2_DIV_PI, 1)[0];
        ulong xi, x2, x1, q, q2, h, l, h2, l2, p4, p3, p2, p1, p0;
        ulong t4, t3, t2, t1, t0, u4, u3, u2, u1, u0;

        /* |x| = (xi, x2, x1) with the binary point after xi; 1 <= e <= 32 */
        xi = d1 >> (FLINT_BITS - e);
        x2 = (d1 << e) | (d0 >> (FLINT_BITS - e));
        x1 = d0 << e;

        umul_ppmm(h, l, x2, K);
        umul_ppmm(h2, l2, xi, K);
        add_ssaaaa(q, l2, h2, l2, UWORD(0), h);
        (void) l;

        q2 = 2 * q;
        umul_ppmm(p1, p0, q2, P[0]);
        umul_ppmm(h, l, q2, P[1]);
        add_ssaaaa(p2, p1, h, l, UWORD(0), p1);
        umul_ppmm(h, l, q2, P[2]);
        add_ssaaaa(p3, p2, h, l, UWORD(0), p2);
        umul_ppmm(h, l, q2, P[3]);
        add_ssaaaa(p4, p3, h, l, UWORD(0), p3);
        sub_dddddmmmmmsssss(t4, t3, t2, t1, t0, xi, x2, x1, UWORD(0), UWORD(0), p4, p3, p2, p1, p0);

        u4 = P[3] >> (FLINT_BITS - 1);
        u3 = (P[3] << 1) | (P[2] >> (FLINT_BITS - 1));
        u2 = (P[2] << 1) | (P[1] >> (FLINT_BITS - 1));
        u1 = (P[1] << 1) | (P[0] >> (FLINT_BITS - 1));
        u0 = P[0] << 1;

        if (t4 > u4 || (t4 == u4 && (t3 > u3 || (t3 == u3 && (t2 > u2 || (t2 == u2
                && (t1 > u1 || (t1 == u1 && t0 >= u0))))))))
        {
            sub_dddddmmmmmsssss(t4, t3, t2, t1, t0, t4, t3, t2, t1, t0, u4, u3, u2, u1, u0);
            q++;
        }

        if (t4 != 0 || t3 > P[3] || (t3 == P[3] && (t2 > P[2] || (t2 == P[2]
                && (t1 > P[1] || (t1 == P[1] && t0 > P[0]))))))
        {
            sub_dddddmmmmmsssss(t4, t3, t2, t1, t0, u4, u3, u2, u1, u0, t4, t3, t2, t1, t0);
            q++;
            sneg = 1;
        }

        if (t3 < (UWORD(1) << 24))
            return NFLOAT_TRIG_FAST_DECLINE;

        a = q & 3;
        v2 = t3;
        v1 = t2;
        v0 = t1;
    }

    _mp_real_small_sin_cos_g_3(&s2, &s1, &s0, &g2, &g1, &g0, v2, v1, v0);

    NFLOAT_TRIG_FAST_SIGNS(nsv, ncv, a, sneg, xsgn);
    os = (a & 1) ? rc : rs;
    oc = (a & 1) ? rs : rc;

    if (os != NULL)
        status |= _nfloat_small_set_2(os, 0, s2, s1, s0, 0, MP_REAL_SMALL_SIN_COS_ERR_3 + 2,
            nsv ^ ((os == rs) ? xsgn : 0), ctx);
    if (oc != NULL)
    {
        sub_ddddmmmmssss(c3, c2, c1, c0, UWORD(1), UWORD(0), UWORD(0), UWORD(0), UWORD(0), g2, g1, g0);
        status |= _nfloat_small_set_2(oc, c3, c2, c1, c0, 0, MP_REAL_SMALL_SIN_COS_ERR_3 + 2,
            ncv ^ ((oc == rs) ? xsgn : 0), ctx);
    }

    return status;
}
#endif

/* which = 0: sin and/or cos; which = 1: tan */
static int
_nfloat_trig(nfloat_ptr rs, nfloat_ptr rc, nfloat_ptr rt, nfloat_srcptr x, gr_ctx_t ctx)
{
    ulong v[NFLOAT_ELEM_MAX_LIMBS + 1];
    ulong xbuf[NFLOAT_MAX_LIMBS + 1];
    mp_real_t X;
    slong n, e, w, z;
    int xsgn, a, s;
    int status;

    n = NFLOAT_CTX_NLIMBS(ctx);
    e = NFLOAT_EXP(x);
    xsgn = NFLOAT_SGNBIT(x);

    /* tiny x: sin x = x (1 - x^2/6 + ...), cos x = 1 - x^2/2 + ...,
       tan x = x (1 + x^2/3 + ...), x^2 < 2^(2e) <= 2^(-FLINT_BITS n) */
    if (2 * e <= -FLINT_BITS * n)
    {
        status = GR_SUCCESS;
        if (rt != NULL)
            return _nfloat_set_tiny_rel_err(rt, x, 1, ctx);
        if (rc != NULL)
            status |= _nfloat_one_plus_tiny(rc, 0, 1, -FLINT_BITS * n, ctx);
        if (rs != NULL)
            status |= _nfloat_set_tiny_rel_err(rs, x, -1, ctx);
        return status;
    }

#if FLINT_BITS == 64
    if (rt == NULL && e <= 32 && n <= 2)
    {
        status = (n == 1) ? _nfloat_sin_cos_fast_1(rs, rc, NFLOAT_D(x)[0], e, xsgn, ctx)
                          : _nfloat_sin_cos_fast_2(rs, rc, NFLOAT_D(x)[1], NFLOAT_D(x)[0], e, xsgn, ctx);
        if (status != NFLOAT_TRIG_FAST_DECLINE)
            return status;
    }
#endif

    if (e <= 0)
    {
        /* |x| < 1 unreduced: v = |x| with z = -e leading zero bits */
        int trunc;
        z = -e;
        w = _nfloat_elem_limbs(ctx, (rs != NULL || rt != NULL) ? z : 0);
        trunc = _nfloat_get_fixed(v, w, w, NFLOAT_D(x), n, e);
        return _nfloat_trig_eval(rs, rc, rt, v, w, trunc, 0, 1, xsgn, ctx);
    }

    if (e > NFLOAT_TRIG_MAX_EXP)
        return _nfloat_trig_slow(rs, rc, rt, x, 0, ctx);

    w = _nfloat_elem_limbs(ctx, 0);
    _nfloat_mp_real_view(X, xbuf, x, ctx);

    for (;;)
    {
        int need_rel, r;

        /* the reduction of mp_real_sin_cos_bits */
        r = _mp_real_trig_reduce(v, X, w);
        a = r & 3;
        s = (r & 4) ? -1 : 1;

        /* the output sin v (or tan v) needs relative accuracy */
        need_rel = (rt != NULL) || (rs != NULL && !(a & 1)) || (rc != NULL && (a & 1));

        if (!need_rel)
            break;

        z = _mp_real_elem_lzb(v, w);

        /* v within 3 ulps: z is accurate unless v is that small */
        if (z >= FLINT_BITS * w - 4)
            z = FLINT_BITS * w;

        if (_nfloat_trig_limbs(ctx, z, rt != NULL && (a & 1)) <= w)
            break;

        w = _nfloat_trig_limbs(ctx, z, rt != NULL && (a & 1));

        /* more cancellation than the temporaries allow */
        if (w > NFLOAT_ELEM_MAX_LIMBS)
            return _nfloat_trig_slow(rs, rc, rt, x, 0, ctx);
    }

    return _nfloat_trig_eval(rs, rc, rt, v, w, 3, a, s, xsgn, ctx);
}

int
nfloat_sin_cos(nfloat_ptr res1, nfloat_ptr res2, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_zero(res1, ctx) | nfloat_one(res2, ctx);
        return nfloat_nan(res1, ctx) | nfloat_nan(res2, ctx);
    }

    if (res1 == x || res2 == x)
    {
        ulong t[NFLOAT_MAX_ALLOC];
        nfloat_set(t, x, ctx);
        return _nfloat_trig(res1, res2, NULL, t, ctx);
    }

    return _nfloat_trig(res1, res2, NULL, x, ctx);
}

int
nfloat_sin(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_zero(res, ctx);
        return nfloat_nan(res, ctx);
    }

    return _nfloat_trig(res, NULL, NULL, x, ctx);
}

int
nfloat_cos(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_one(res, ctx);
        return nfloat_nan(res, ctx);
    }

    return _nfloat_trig(NULL, res, NULL, x, ctx);
}

int
nfloat_tan(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_zero(res, ctx);
        return nfloat_nan(res, ctx);
    }

    return _nfloat_trig(NULL, NULL, res, x, ctx);
}

/* sin(pi x), cos(pi x), tan(pi x) *******************************************/

static int
_nfloat_trig_pi(nfloat_ptr rs, nfloat_ptr rc, nfloat_ptr rt, nfloat_srcptr x, gr_ctx_t ctx)
{
    ulong T[NFLOAT_MAX_LIMBS + 1], xbuf[NFLOAT_MAX_LIMBS + 1];
    ulong v[NFLOAT_ELEM_MAX_LIMBS + 2];
    mp_real_t X;
    slong n, e, w, z, tl, te;
    int xsgn, a, s, r;
    int status = GR_SUCCESS;

    n = NFLOAT_CTX_NLIMBS(ctx);
    e = NFLOAT_EXP(x);
    xsgn = NFLOAT_SGNBIT(x);

    /* an even integer (all bits of weight >= 2) */
    if (e > FLINT_BITS * n)
    {
        if (rs != NULL)
            status |= nfloat_zero(rs, ctx);
        if (rc != NULL)
            status |= nfloat_one(rc, ctx);
        if (rt != NULL)
            status |= nfloat_zero(rt, ctx);
        return status;
    }

    /* tiny x: sin(pi x) = pi x (1 - delta), tan(pi x) = pi x (1 + delta),
       cos(pi x) = 1 - delta', delta < (pi x)^2 < 2^(2e + 4) */
    if (2 * e + 4 <= -FLINT_BITS * n - 2 * FLINT_BITS)
    {
        if (rc != NULL)
            status |= _nfloat_one_plus_tiny(rc, 0, 1, -FLINT_BITS * n, ctx);

        if (rs != NULL || rt != NULL)
        {
            /* pi |x| from (d, n) times pi at n + 2 limbs: within 3 units of
               B^-(n + 1) relative to the mantissa, delta below one more */
            ulong t[2 * NFLOAT_MAX_LIMBS + 4];
            ulong p[NFLOAT_MAX_LIMBS + 3];
            slong pn = n + 2;

            p[pn] = mpn_lshift(p, _mp_real_const_ptr(MP_REAL_CONST_ID_PI4, pn), pn, 2);
            flint_mpn_mul(t, p, pn + 1, NFLOAT_D(x), n);
            /* (t, 2n + 3) = pi D B^pn (to 4 ulps of the pi limbs); keep the
               top n + 2 limbs: units of B^(n + 1) */
            if (rs != NULL)
                status |= _nfloat_set_mpn_err(rs, t + n + 1, n + 2,
                    e - FLINT_BITS * n + FLINT_BITS * (n + 1) - FLINT_BITS * pn, 6, 6, xsgn, ctx);
            if (rt != NULL)
                status |= _nfloat_set_mpn_err(rt, t + n + 1, n + 2,
                    e - FLINT_BITS * n + FLINT_BITS * (n + 1) - FLINT_BITS * pn, 6, 6, xsgn, ctx);
        }
        return status;
    }

    /* |x| = a/2 + s t exactly, t in [0, 1/4] (the reduction of
       mp_real_sin_cos_pi_bits) */
    _nfloat_mp_real_view(X, xbuf, x, ctx);
    r = _mp_real_trig_pi_reduce(T, &tl, &te, X);
    a = r & 3;
    s = (r & 4) ? -1 : 1;

    if (tl == 0)
    {
        /* x a multiple of 1/2: sin, cos in {0, 1, -1}, tan in {0, pole} */
        if (rt != NULL)
            return (a & 1) ? nfloat_nan(rt, ctx) : nfloat_zero(rt, ctx);

        if (a & 1)
        {
            if (rs != NULL)
            {
                status |= nfloat_one(rs, ctx);
                NFLOAT_SGNBIT(rs) = (a == 3) ^ xsgn;
            }
            if (rc != NULL)
                status |= nfloat_zero(rc, ctx);
        }
        else
        {
            if (rs != NULL)
                status |= nfloat_zero(rs, ctx);
            if (rc != NULL)
            {
                status |= nfloat_one(rc, ctx);
                NFLOAT_SGNBIT(rc) = (a == 2);
            }
        }
        return status;
    }

    /* 2^-(z + 1) <= t < 2^-z: relative accuracy for t with z leading
       zero bits */
    z = (slong) flint_clz(T[tl - 1]) - FLINT_BITS * te;
    w = _nfloat_trig_limbs(ctx, z, rt != NULL && (a & 1));
    if (w > NFLOAT_ELEM_MAX_LIMBS)
        return _nfloat_trig_slow(rs, rc, rt, x, 3, ctx);

    /* v = pi t within 3 ulps */
    _mp_real_trig_pi_v_fixed(v, T, tl, te, w);
    v[w] = 0;

    return _nfloat_trig_eval(rs, rc, rt, v, w, 3, a, s, xsgn, ctx);
}

int
nfloat_sin_pi(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_zero(res, ctx);
        return nfloat_nan(res, ctx);
    }

    return _nfloat_trig_pi(res, NULL, NULL, x, ctx);
}

int
nfloat_cos_pi(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_one(res, ctx);
        return nfloat_nan(res, ctx);
    }

    return _nfloat_trig_pi(NULL, res, NULL, x, ctx);
}

int
nfloat_sin_cos_pi(nfloat_ptr res1, nfloat_ptr res2, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_zero(res1, ctx) | nfloat_one(res2, ctx);
        return nfloat_nan(res1, ctx) | nfloat_nan(res2, ctx);
    }

    if (res1 == x || res2 == x)
    {
        ulong t[NFLOAT_MAX_ALLOC];
        nfloat_set(t, x, ctx);
        return _nfloat_trig_pi(res1, res2, NULL, t, ctx);
    }

    return _nfloat_trig_pi(res1, res2, NULL, x, ctx);
}

int
nfloat_tan_pi(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_zero(res, ctx);
        return nfloat_nan(res, ctx);
    }

    return _nfloat_trig_pi(NULL, NULL, res, x, ctx);
}
