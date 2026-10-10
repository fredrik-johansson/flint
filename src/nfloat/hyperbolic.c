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
    Hyperbolic functions from E = exp(|x|) and its reciprocal:
    cosh = (E + 1/E) / 2, sinh = (E - 1/E) / 2, tanh = (E - 1/E) / (E + 1/E),
    composed with relative error bounds. For |x| < 1, E - 1/E is computed
    in fixed point with absolute errors at a precision extended by the
    -e + 1 bits cancelling.
*/

/* E = exp(|x|) at w limbs */
static void
_nfloat_exp_abs_rf(_nfloat_rf_t E, nfloat_srcptr x, slong w, gr_ctx_t ctx)
{
    ulong y[NFLOAT_ELEM_MAX_LIMBS + 2];
    ulong err;
    slong k;

    k = _nfloat_exp_mpn(y, &err, NFLOAT_D(x), NFLOAT_CTX_NLIMBS(ctx), NFLOAT_EXP(x), 0, w);
    _nfloat_rf_set_mpn(E, y, w + 1, k - FLINT_BITS * w, err, w);
}

/* For 0 < |x| < 1, with E = exp(|x|): sets (T, wb + 1) = E - 1/E and
   (S, wb + 1) = E + 1/E in units of B^-wb, both within the returned
   number of units. The difference is computed in fixed point, so that
   wb limbs with -e extra bits give a relative accuracy for E - 1/E ~ 2x
   (converting E and 1/E to relative errors before subtracting them would
   blow up the error counts). */
static ulong
_nfloat_exp_pm_fixed(nn_ptr T, nn_ptr S, nfloat_srcptr x, slong wb, gr_ctx_t ctx)
{
    ulong y[NFLOAT_ELEM_MAX_LIMBS + 2], r[NFLOAT_ELEM_MAX_LIMBS + 2];
    _nfloat_rf_t E, R;
    ulong ey, er;
    slong k;

    k = _nfloat_exp_mpn(y, &ey, NFLOAT_D(x), NFLOAT_CTX_NLIMBS(ctx), NFLOAT_EXP(x), 0, wb);
    FLINT_ASSERT(k == 0);
    (void) k;

    /* R = 1/E in (1/e, 1), R->exp in {0, -1} */
    _nfloat_rf_set_mpn(E, y, wb + 1, -FLINT_BITS * wb, ey, wb);
    _nfloat_rf_inv(R, E, wb);
    FLINT_ASSERT(R->exp == 0 || R->exp == -1);

    if (R->exp == 0)
    {
        _nfloat_copy_limbs(r, R->d, wb);
        er = R->err;
    }
    else
    {
        mpn_rshift(r, R->d, wb, 1);
        er = R->err + 1;
    }
    r[wb] = 0;

    if (T != NULL)
        mpn_sub_n(T, y, r, wb + 1);
    if (S != NULL)
        mpn_add_n(S, y, r, wb + 1);

    return ey + er;
}

/* sinh (want & 1), cosh (want & 2) */
static int
_nfloat_sinh_cosh(nfloat_ptr rs, nfloat_ptr rc, nfloat_srcptr x, gr_ctx_t ctx)
{
    _nfloat_rf_t E, R, T;
    slong n, e, w;
    int sgnbit, status = GR_SUCCESS;

    n = NFLOAT_CTX_NLIMBS(ctx);
    e = NFLOAT_EXP(x);
    sgnbit = NFLOAT_SGNBIT(x);

    if (e > NFLOAT_EXP_MAX_ARG_EXP)
    {
        if (rs != NULL)
            status |= _nfloat_overflow(rs, sgnbit, ctx);
        if (rc != NULL)
            status |= _nfloat_overflow(rc, 0, ctx);
        return status;
    }

    /* sinh x = x (1 + delta), cosh x = 1 + delta', delta, delta' < x^2 */
    if (2 * e <= -FLINT_BITS * n)
    {
        if (rs != NULL)
            status |= _nfloat_set_tiny_rel_err(rs, x, 1, ctx);
        if (rc != NULL)
            status |= _nfloat_one_plus_tiny(rc, 0, 0, -FLINT_BITS * n, ctx);
        return status;
    }

    if (e <= 0)
    {
        ulong T[NFLOAT_ELEM_MAX_LIMBS + 2], S[NFLOAT_ELEM_MAX_LIMBS + 2];
        ulong err;

        w = _nfloat_elem_limbs(ctx, -e + 1);
        err = _nfloat_exp_pm_fixed((rs != NULL) ? T : NULL, (rc != NULL) ? S : NULL, x, w, ctx);

        if (rs != NULL)
            status |= _nfloat_set_mpn_err(rs, T, w + 1, -FLINT_BITS * w - 1, err, err, sgnbit, ctx);
        if (rc != NULL)
            status |= _nfloat_set_mpn_err(rc, S, w + 1, -FLINT_BITS * w - 1, err, err, 0, ctx);
        return status;
    }

    w = _nfloat_elem_limbs(ctx, 0);

    _nfloat_exp_abs_rf(E, x, w, ctx);
    _nfloat_rf_inv(R, E, w);

    if (rs != NULL)
    {
        _nfloat_rf_sub(T, E, R, w);
        T->exp--;
        status |= _nfloat_set_rf(rs, T, sgnbit, w, ctx);
    }

    if (rc != NULL)
    {
        _nfloat_rf_add(T, E, R, w);
        T->exp--;
        status |= _nfloat_set_rf(rc, T, 0, w, ctx);
    }

    return status;
}

int
nfloat_sinh(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
        return nfloat_set(res, x, ctx) | (NFLOAT_IS_NAN(x) ? nfloat_nan(res, ctx) : GR_SUCCESS);

    return _nfloat_sinh_cosh(res, NULL, x, ctx);
}

int
nfloat_cosh(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_one(res, ctx);
        if (NFLOAT_IS_INF(x))
            return nfloat_pos_inf(res, ctx);
        return nfloat_nan(res, ctx);
    }

    return _nfloat_sinh_cosh(NULL, res, x, ctx);
}

int
nfloat_sinh_cosh(nfloat_ptr res1, nfloat_ptr res2, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
        return nfloat_sinh(res1, x, ctx) | nfloat_cosh(res2, x, ctx);

    if (res1 == x || res2 == x)
    {
        ulong t[NFLOAT_MAX_ALLOC];
        nfloat_set(t, x, ctx);
        return _nfloat_sinh_cosh(res1, res2, t, ctx);
    }

    return _nfloat_sinh_cosh(res1, res2, x, ctx);
}

int
nfloat_tanh(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    _nfloat_rf_t E, R, N, D;
    slong n, e, w;
    int sgnbit;

    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_zero(res, ctx);
        if (NFLOAT_IS_POS_INF(x))
            return nfloat_one(res, ctx);
        if (NFLOAT_IS_NEG_INF(x))
            return nfloat_neg_one(res, ctx);
        return nfloat_nan(res, ctx);
    }

    n = NFLOAT_CTX_NLIMBS(ctx);
    e = NFLOAT_EXP(x);
    sgnbit = NFLOAT_SGNBIT(x);

    /* tanh x = x (1 - delta), delta < x^2 / 3 */
    if (2 * e <= -FLINT_BITS * n)
        return _nfloat_set_tiny_rel_err(res, x, -1, ctx);

    /* |tanh x| = 1 - 2 / (E^2 + 1) with 2 / (E^2 + 1) < 2 exp(-2 |x|) <
       2^(1 - 2.88 |x|) <= 2^(-FLINT_BITS n) for |x| >= 2^(e - 1) >= 23 n + 1 */
    if (e - 1 >= (slong) FLINT_BIT_COUNT(23 * n + 1))
        return _nfloat_one_plus_tiny(res, sgnbit, 1, -FLINT_BITS * n, ctx);

    w = _nfloat_elem_limbs(ctx, 0);

    if (e <= 0)
    {
        ulong T[NFLOAT_ELEM_MAX_LIMBS + 2], S[NFLOAT_ELEM_MAX_LIMBS + 2];
        ulong err;
        /* E - 1/E in [2^e, 2^(e+2)) B^wb, with -e + 1 more bits for the
           cancellation, and enough that its normalization to w limbs
           shifts left by at most NFLOAT_ELEM_MAX_NORM_SHIFT bits */
        slong wb = _nfloat_elem_limbs(ctx, -e + 1);
        if (-e > NFLOAT_ELEM_MAX_NORM_SHIFT)
            wb = FLINT_MAX(wb, w + (-e - NFLOAT_ELEM_MAX_NORM_SHIFT + FLINT_BITS - 1) / FLINT_BITS);

        /* (E - 1/E) / (E + 1/E), the difference in fixed point */
        err = _nfloat_exp_pm_fixed(T, S, x, wb, ctx);
        _nfloat_rf_set_mpn(N, T, wb + 1, -FLINT_BITS * wb, err, w);
        _nfloat_rf_set_mpn(D, S, wb + 1, -FLINT_BITS * wb, err, w);
    }
    else
    {
        _nfloat_exp_abs_rf(E, x, w, ctx);
        _nfloat_rf_inv(R, E, w);
        _nfloat_rf_sub(N, E, R, w);
        _nfloat_rf_add(D, E, R, w);
    }

    _nfloat_rf_div(N, N, D, w);
    return _nfloat_set_rf(res, N, sgnbit, w, ctx);
}
