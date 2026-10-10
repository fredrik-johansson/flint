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
    Exponentials.

    exp(x) for |x| < 2^(FLINT_BITS - 3) (larger arguments overflow or underflow the
    exponent range) is reduced as |x| = q log 2 + t, t in [0, log 2), in
    fixed point with one integral and nk + 1 fraction limbs, nk the
    kernel limbs: q is a lower approximation of |x| / log 2 from two
    limbs of 1/log 2 (exact or one low, then corrected; the reduction
    step is that of mp_real/exp.c), and log 2 is
    the floor L of log 2 B^(nk + 1), so that the reduction errs by
    (q + 1)(log 2 - L) < 2^(FLINT_BITS - 2) B^-(nk + 1), a quarter ulp of B^-nk. For
    x < 0, exp(x) = 2^-(q+1) exp(L - t'), t' = |x| - q log 2. Then
    exp(x) = 2^k exp(t) with the kernel on t at nk fraction limbs.
    Positive x < 1 go to the kernel directly, and small negative x to
    the alternating series of exp(-|x|).
*/

/* exp(x) for the normal nfloat x = (-1)^sgnbit (d, n) 2^(e - FLINT_BITS n),
   e <= NFLOAT_EXP_MAX_ARG_EXP, with nk kernel limbs: sets (y, nk + 1)
   (y needs nk + 2 limbs) and *err, and returns k such that
   exp(x) = Y 2^(k - FLINT_BITS nk) within *err units, Y = (y, nk + 1).
   Y is in [B^nk, 2 B^nk) except on the paths with k = 0 for x < 1. */
slong
_nfloat_exp_mpn(nn_ptr y, ulong * err, nn_srcptr d, slong n, slong e, int sgnbit, slong nk)
{
    ulong X[NFLOAT_ELEM_MAX_LIMBS + 3];
    nn_srcptr L;
    ulong q;
    slong N, k;
    int trunc;

    if (e <= 0)
    {
        slong z = -e;

        if (!sgnbit)
        {
            /* x in (0, 1): truncating x costs one ulp, times exp' < e */
            trunc = _nfloat_get_fixed(X, nk, nk, d, n, e);
            _mp_real_exp_kernel(y, err, X, nk);
            *err += 3 * trunc;
            return 0;
        }

        /* the series of exp(-|x|) where it is faster (exp' <= 1) */
        (void) z;
        trunc = _nfloat_get_fixed(X, nk, nk, d, n, e);
        if (_mp_real_exp_neg_series(y, err, X, nk))
        {
            *err += trunc;
            return 0;
        }
    }

    /* X = |x| with one integral and N fraction limbs (truncated) */
    N = nk + 1;
    _nfloat_get_fixed(X, N + 1, N, d, n, e);
    L = _mp_real_const_ptr(MP_REAL_CONST_ID_LOG2, N);

    /* q <= |x| / log 2, exact or one less, and t = X - q L resp.
       (q + 1) L - |x| from the two's complement (mp_real/impl.h) */
    q = _mp_real_exp_quotient(X[N], X[N - 1]);
    if (sgnbit)
        mpn_neg(X, X, N + 1);
    k = _mp_real_exp_reduce_step(X, L, N, q, sgnbit);

    FLINT_ASSERT(X[N] == 0);

    /* t within two ulps of B^-nk (the reduction, the truncation of x and
       of the low fraction limb), times exp' < 2 */
    _mp_real_exp_kernel(y, err, X + 1, nk);
    *err += 4;
    return k;
}

#if FLINT_BITS == 64
/*
    In-register fast paths at one and two limbs for |x| < 2^10: the
    reduction |x| = q log 2 + t (t in [0, log 2), q from
    _mp_real_exp_quotient, one more log 2 if q was one low) on the
    mantissa with one integral and two resp. three fraction limbs plus
    one more for the product q L, L the floor of log 2 (the reduction
    errs by q 2^-192 resp. q 2^-256), then the table-driven kernel of
    mp_real at two resp. three limbs on the top fraction limbs of t (or
    of L - t for x < 0, exp(x) = 2^-(q+1) exp(L - t)): one more ulp for
    the dropped limb, twice through exp' < 2.
*/
static int
_nfloat_exp_fast_1(nfloat_ptr res, ulong d, slong e, int sgnbit, gr_ctx_t ctx)
{
    nn_srcptr L = _mp_real_const_ptr(MP_REAL_CONST_ID_LOG2, 3);
    ulong xi, x1, x0, q, p3, p2, p1, p0, h, l, t3, t2, t1, t0, y2, y1, y0;
    slong k;

    /* |x| = (xi, x1, x0) with the binary point after xi; -64 <= e <= 10 */
    if (e > 0)
    {
        xi = d >> (FLINT_BITS - e);
        x1 = d << e;
        x0 = 0;
    }
    else if (e == 0)
    {
        xi = 0; x1 = d; x0 = 0;
    }
    else if (e > -FLINT_BITS)
    {
        xi = 0;
        x1 = d >> (-e);
        x0 = d << (FLINT_BITS + e);
    }
    else
    {
        xi = 0; x1 = 0; x0 = d;
    }

    q = _mp_real_exp_quotient(xi, x1);
    umul_ppmm(p1, p0, q, L[0]);
    umul_ppmm(h, l, q, L[1]);
    add_ssaaaa(p2, p1, h, l, UWORD(0), p1);
    umul_ppmm(h, l, q, L[2]);
    add_ssaaaa(p3, p2, h, l, UWORD(0), p2);
    sub_ddddmmmmssss(t3, t2, t1, t0, xi, x1, x0, UWORD(0), p3, p2, p1, p0);

    if (t3 != 0 || t2 > L[2] || (t2 == L[2] && (t1 > L[1] || (t1 == L[1] && t0 >= L[0]))))
    {
        sub_ddddmmmmssss(t3, t2, t1, t0, t3, t2, t1, t0, UWORD(0), L[2], L[1], L[0]);
        q++;
    }

    if (sgnbit)
    {
        sub_dddmmmsss(t2, t1, t0, L[2], L[1], L[0], t2, t1, t0);
        k = -(slong) (q + 1);
    }
    else
        k = (slong) q;

    _mp_real_small_exp_2(&y2, &y1, &y0, t2, t1);
    return _nfloat_small_set_1(res, y2, y1, y0, k, MP_REAL_SMALL_EXP_ERR_2 + 3, 0, ctx);
}

static int
_nfloat_exp_fast_2(nfloat_ptr res, ulong d1, ulong d0, slong e, int sgnbit, gr_ctx_t ctx)
{
    nn_srcptr L = _mp_real_const_ptr(MP_REAL_CONST_ID_LOG2, 4);
    ulong xi, x2, x1, x0, q, p4, p3, p2, p1, p0, h, l, t4, t3, t2, t1, t0;
    ulong y3, y2, y1, y0;
    slong k;
    unsigned int s;

    /* |x| = (xi, x2, x1, x0) with the binary point after xi (truncated
       below B^-3); -128 <= e <= 10 */
    xi = x2 = x1 = x0 = 0;
    if (e > 0)
    {
        xi = d1 >> (FLINT_BITS - e);
        x2 = (d1 << e) | (d0 >> (FLINT_BITS - e));
        x1 = d0 << e;
    }
    else if (e == 0)
    {
        x2 = d1; x1 = d0;
    }
    else if (e > -FLINT_BITS)
    {
        s = -e;
        x2 = d1 >> s;
        x1 = (d1 << (FLINT_BITS - s)) | (d0 >> s);
        x0 = d0 << (FLINT_BITS - s);
    }
    else if (e == -FLINT_BITS)
    {
        x1 = d1; x0 = d0;
    }
    else if (e > -2 * FLINT_BITS)
    {
        s = -e - FLINT_BITS;
        x1 = d1 >> s;
        x0 = (d1 << (FLINT_BITS - s)) | (d0 >> s);
    }
    else
        x0 = d1;

    q = _mp_real_exp_quotient(xi, x2);
    umul_ppmm(p1, p0, q, L[0]);
    umul_ppmm(h, l, q, L[1]);
    add_ssaaaa(p2, p1, h, l, UWORD(0), p1);
    umul_ppmm(h, l, q, L[2]);
    add_ssaaaa(p3, p2, h, l, UWORD(0), p2);
    umul_ppmm(h, l, q, L[3]);
    add_ssaaaa(p4, p3, h, l, UWORD(0), p3);
    sub_dddddmmmmmsssss(t4, t3, t2, t1, t0, xi, x2, x1, x0, UWORD(0), p4, p3, p2, p1, p0);

    if (t4 != 0 || t3 > L[3] || (t3 == L[3] && (t2 > L[2] || (t2 == L[2]
            && (t1 > L[1] || (t1 == L[1] && t0 >= L[0]))))))
    {
        sub_dddddmmmmmsssss(t4, t3, t2, t1, t0, t4, t3, t2, t1, t0, UWORD(0), L[3], L[2], L[1], L[0]);
        q++;
    }

    if (sgnbit)
    {
        sub_ddddmmmmssss(t3, t2, t1, t0, L[3], L[2], L[1], L[0], t3, t2, t1, t0);
        k = -(slong) (q + 1);
    }
    else
        k = (slong) q;

    _mp_real_small_exp_3(&y3, &y2, &y1, &y0, t3, t2, t1);
    return _nfloat_small_set_2(res, y3, y2, y1, y0, k, MP_REAL_SMALL_EXP_ERR_3 + 3, 0, ctx);
}
#endif

int
nfloat_exp(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    ulong y[NFLOAT_ELEM_MAX_LIMBS + 2];
    ulong err;
    slong n, e, nk, k;
    int sgnbit;

    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_one(res, ctx);
        if (NFLOAT_IS_POS_INF(x))
            return nfloat_pos_inf(res, ctx);
        if (NFLOAT_IS_NEG_INF(x))
            return nfloat_zero(res, ctx);
        return nfloat_nan(res, ctx);
    }

    n = NFLOAT_CTX_NLIMBS(ctx);
    e = NFLOAT_EXP(x);
    sgnbit = NFLOAT_SGNBIT(x);

    if (e > NFLOAT_EXP_MAX_ARG_EXP)
        return sgnbit ? _nfloat_underflow(res, 0, ctx) : _nfloat_overflow(res, 0, ctx);

    /* |x| < 2^e: exp(x) = 1 + delta, |delta| < 2^(e + 1) <= 2^(-FLINT_BITS n) */
    if (e < -FLINT_BITS * n
        || (e < -NFLOAT_CTX_FUNC_PREC(ctx) && !NFLOAT_CTX_HAS_DIRECTED_ROUNDING(ctx)))
        return _nfloat_one_plus_tiny(res, 0, sgnbit, -FLINT_BITS * n, ctx);

#if FLINT_BITS == 64
    if (n == 1 && e <= 10)
        return _nfloat_exp_fast_1(res, NFLOAT_D(x)[0], e, sgnbit, ctx);
    if (n == 2 && e <= 10)
        return _nfloat_exp_fast_2(res, NFLOAT_D(x)[1], NFLOAT_D(x)[0], e, sgnbit, ctx);
#endif

    nk = _nfloat_elem_limbs(ctx, 0);
    k = _nfloat_exp_mpn(y, &err, NFLOAT_D(x), n, e, sgnbit, nk);
    return _nfloat_set_mpn_err(res, y, nk + 1, k - FLINT_BITS * nk, err, err, 0, ctx);
}

int
nfloat_expm1(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    ulong y[NFLOAT_ELEM_MAX_LIMBS + 2];
    ulong err;
    slong n, e, nk, k;
    int sgnbit;

    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_zero(res, ctx);
        if (NFLOAT_IS_POS_INF(x))
            return nfloat_pos_inf(res, ctx);
        if (NFLOAT_IS_NEG_INF(x))
            return nfloat_neg_one(res, ctx);
        return nfloat_nan(res, ctx);
    }

    n = NFLOAT_CTX_NLIMBS(ctx);
    e = NFLOAT_EXP(x);
    sgnbit = NFLOAT_SGNBIT(x);

    if (e > NFLOAT_EXP_MAX_ARG_EXP)
    {
        if (!sgnbit)
            return _nfloat_overflow(res, 0, ctx);
        /* -1 + exp(x), exp(x) < 2^(-2^(FLINT_BITS - 3)) */
        return _nfloat_one_plus_tiny(res, 1, 1, -FLINT_BITS * n, ctx);
    }

    /* expm1(x) = x (1 + delta), delta = x/2 + x^2/6 + ..., |delta| < |x| */
    if (e < -FLINT_BITS * n
        || (e < -NFLOAT_CTX_FUNC_PREC(ctx) && !NFLOAT_CTX_HAS_DIRECTED_ROUNDING(ctx)))
        return _nfloat_set_tiny_rel_err(res, x, sgnbit ? -1 : 1, ctx);

    if (e <= 0)
    {
        /* |x| < 2^-z: an absolute accuracy of 2^-(p + z) */
        nk = _nfloat_elem_limbs(ctx, -e);
        k = _nfloat_exp_mpn(y, &err, NFLOAT_D(x), n, e, sgnbit, nk);

        if (!sgnbit)
        {
            /* exp(x) - 1 in [0, e - 1) */
            FLINT_ASSERT(k == 0);
            y[nk] -= 1;
            return _nfloat_set_mpn_err(res, y, nk + 1, -FLINT_BITS * nk, err, err, 0, ctx);
        }
        else
        {
            /* 1 - Y 2^k B^-nk, k in {0, -1, -2}: shift Y down by -k bits
               (one more unit for the dropped bits) */
            if (k != 0)
            {
                int dropped;
                FLINT_ASSERT(k >= -2);
                dropped = (y[0] & ((UWORD(1) << (-k)) - 1)) != 0;
                mpn_rshift(y, y, nk + 1, (unsigned int) (-k));
                err = (err >> (-k)) + 1 + dropped;
            }
            mpn_neg(y, y, nk + 1);
            y[nk] += 1;
            return _nfloat_set_mpn_err(res, y, nk + 1, -FLINT_BITS * nk, err, err, 1, ctx);
        }
    }
    else
    {
        nk = _nfloat_elem_limbs(ctx, 0);
        k = _nfloat_exp_mpn(y, &err, NFLOAT_D(x), n, e, sgnbit, nk);

        if (!sgnbit)
        {
            /* Y 2^(k - FLINT_BITS nk) - 1 with k >= 1: subtract the bit
               FLINT_BITS nk - k of Y, or a fraction of a unit */
            ulong elo = err;
            slong b = FLINT_BITS * nk - k;

            if (b >= 0)
                mpn_sub_1(y + b / FLINT_BITS, y + b / FLINT_BITS,
                    nk + 1 - b / FLINT_BITS, UWORD(1) << (b % FLINT_BITS));
            else
                elo++;

            return _nfloat_set_mpn_err(res, y, nk + 1, k - FLINT_BITS * nk, elo, err, 0, ctx);
        }
        else
        {
            /* 1 - exp(x), exp(x) < 1/2 (k <= -1): Y 2^k at nk fraction
               limbs (one more unit for the dropped bits) */
            slong s = -k;

            FLINT_ASSERT(k <= -1);

            if (s >= FLINT_BITS * (nk + 1))
            {
                _nfloat_zero_limbs(y, nk + 1);
                err = 1;
            }
            else
            {
                slong ls = s / FLINT_BITS, i;
                unsigned int bs = s % FLINT_BITS;
                int dropped = 0;

                for (i = 0; i < ls; i++)
                    dropped |= (y[i] != 0);
                if (bs != 0)
                {
                    dropped |= ((y[ls] << (FLINT_BITS - bs)) != 0);
                    mpn_rshift(y, y + ls, nk + 1 - ls, bs);
                }
                else
                    flint_mpn_copyi(y, y + ls, nk + 1 - ls);
                _nfloat_zero_limbs(y + nk + 1 - ls, ls);
                err = (s >= FLINT_BITS) ? 1 : (err >> s) + 1;
                err += dropped;
            }

            mpn_neg(y, y, nk + 1);
            y[nk] += 1;
            return _nfloat_set_mpn_err(res, y, nk + 1, -FLINT_BITS * nk, err, err, 1, ctx);
        }
    }
}

/*
    exp2(x) = 2^k exp(f log 2) with k = floor(x), f = x - k in [0, 1):
    the reduction is exact apart from the truncation of f (at nk + 1
    fraction limbs) and the product f L, L the floor of log 2 B^(nk + 1).
*/
int
nfloat_exp2(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    ulong X[NFLOAT_ELEM_MAX_LIMBS + 3];
    ulong T[NFLOAT_ELEM_MAX_LIMBS + 2];
    ulong y[NFLOAT_ELEM_MAX_LIMBS + 2];
    nn_srcptr L;
    ulong err;
    slong n, e, nk, N, k;
    int sgnbit;

    if (NFLOAT_IS_SPECIAL(x))
        return nfloat_exp(res, x, ctx);

    n = NFLOAT_CTX_NLIMBS(ctx);
    e = NFLOAT_EXP(x);
    sgnbit = NFLOAT_SGNBIT(x);

    if (e > NFLOAT_EXP_MAX_ARG_EXP)
        return sgnbit ? _nfloat_underflow(res, 0, ctx) : _nfloat_overflow(res, 0, ctx);

    /* 2^x = 1 + delta, |delta| < 2^e */
    if (e < -FLINT_BITS * n
        || (e < -NFLOAT_CTX_FUNC_PREC(ctx) && !NFLOAT_CTX_HAS_DIRECTED_ROUNDING(ctx)))
        return _nfloat_one_plus_tiny(res, 0, sgnbit, -FLINT_BITS * n, ctx);

    nk = _nfloat_elem_limbs(ctx, 0);
    N = nk + 1;

    /* X = x with one integral limb in two's complement and N fraction
       limbs, truncated towards -inf */
    {
        int trunc = _nfloat_get_fixed(X, N + 1, N, NFLOAT_D(x), n, e);

        if (sgnbit)
        {
            mpn_neg(X, X, N + 1);
            if (trunc)
                mpn_sub_1(X, X, N + 1, 1);
        }

        k = (slong) X[N];

        /* an integer: 2^k exactly */
        if (!trunc && flint_mpn_zero_p(X, N))
        {
            nfloat_one(res, ctx);
            NFLOAT_EXP(res) = k + 1;
            NFLOAT_HANDLE_UNDERFLOW_OVERFLOW(res, ctx);
            return GR_SUCCESS;
        }

        /* t = f L in [0, log 2): the product truncated (2 ulps of B^-N),
           f's truncation (1 ulp) times log 2 */
        L = _mp_real_const_ptr(MP_REAL_CONST_ID_LOG2, N);
        flint_mpn_mulhigh_n(T, X, L, N);
    }

    /* t within two ulps of B^-nk, times exp' < 2 */
    _mp_real_exp_kernel(y, &err, T + 1, nk);
    err += 4;
    return _nfloat_set_mpn_err(res, y, nk + 1, k - FLINT_BITS * nk, err, err, 0, ctx);
}
