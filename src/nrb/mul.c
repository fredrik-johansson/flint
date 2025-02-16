/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <mpfr.h>
#include <stdio.h>
#include "gr.h"
#include "fmpz.h"
#include "arf.h"
#include "arb.h"
#include "nfloat.h"
#include "gr_generic.h"
#include "gr_special.h"
#include "nrb.h"

int
_nrb_fix_range(nrb_ptr res, nrb_ctx_t ctx)
{
    if (NRB_ERR(res) > NRB_MAX_ERR)
    {
        int e;
        slong err_exp, trim_limbs, n, new_n;
        double err = NRB_ERR(res);
        double err2;

        if (err == D_INF)
            return nrb_zero_pm_inf(res, ctx);

        n = NRB_N(res);

        err2 = frexp(err, &e);
        err_exp = e;


        if (n == 0)
        {
            NRB_ERR(res) = err2;
            NRB_EXP(res) += err_exp;
        }
        else
        {
            trim_limbs = err_exp / FLINT_BITS - 1;

            // flint_printf("trim_limbs %wd\n", trim_limbs);

            if (trim_limbs >= n)
            {
                NRB_N(res) = 0;
                NRB_EXP(res) = NRB_EXP(res) - n * FLINT_BITS + err_exp;
                /* todo: check/tighten */
                NRB_ERR(res) = (err2 + ULP_N1) * NRB_CORRECTION_A;
            }
            else
            {
                NRB_N(res) = new_n = n - trim_limbs;
                flint_mpn_copyi(NRB_D(res), NRB_D(res) + trim_limbs, new_n);
                /* todo: check/tighten */
                NRB_ERR(res) = d_mul_2exp_inrange2(err * NRB_CORRECTION_A, -trim_limbs * FLINT_BITS);
            }
        }
    }

    /* todo: NRB_ERR(res) < NRB_MIN_ERR */

    if (NRB_EXP(res) > NFLOAT_MAX_EXP)
        return nrb_zero_pm_inf(res, ctx);

    if (NRB_EXP(res) < NFLOAT_MIN_EXP)
    {
        int e;
        slong err_exp;
        double err = NRB_ERR(res);

        frexp(err, &e);
        err_exp = e;

        NRB_EXP(res) = 1 + FLINT_MAX(NFLOAT_MIN_EXP, NRB_EXP(res) + err_exp);
        NRB_ERR(res) = 1.0;
        NRB_N(res) = 0;
    }

    NRB_ASSERT_VALID(res, ctx);

    return GR_SUCCESS;
}

#define NRB_RETURN_FIX_RANGE(res, ctx) \
    do { \
        if (FLINT_UNLIKELY(NRB_ERR(res) > NRB_MAX_ERR || \
                           NRB_EXP(res) > NFLOAT_MAX_EXP || \
                           NRB_EXP(res) < NFLOAT_MIN_EXP)) \
            return _nrb_fix_range(res, ctx); \
        else \
            return GR_SUCCESS; \
    } while (0)


int nrb_mul_2exp_si(nrb_ptr res, nrb_srcptr x, slong y, nrb_ctx_t ctx)
{
    if (NRB_IS_SPECIAL(x))
    {
        if (NRB_ERR(x) == 0.0)
            return nrb_zero(res, ctx);
        else if (NRB_ERR(x) == D_INF)
            return nrb_zero_pm_inf(res, ctx);
        else
            return GR_UNABLE;
    }

    if (y < NRB_MIN_EXP || y > NRB_MAX_EXP)
        return GR_UNABLE;

    nrb_set(res, x, ctx);
    NRB_EXP(res) += y;

    NRB_RETURN_FIX_RANGE(res, ctx);
}


int nrb_add(nrb_ptr res, nrb_srcptr x, nrb_srcptr y, nrb_ctx_t ctx)
{
    return GR_UNABLE;
}

int nrb_sub(nrb_ptr res, nrb_srcptr x, nrb_srcptr y, nrb_ctx_t ctx)
{
    return GR_UNABLE;
}

static int
_mulhigh_check_exact_product(nn_srcptr x, nn_srcptr y, slong n)
{
    slong val = 0;

    while (x[0] == 0)
    {
        val += FLINT_BITS;
        x++;
    }

    val += flint_ctz(x[0]);

    while (y[0] == 0)
    {
        val += FLINT_BITS;
        y++;
    }

    val += flint_ctz(y[0]);

    return val >= n * FLINT_BITS;
}

static mp_limb_pair_t
_flint_mpn_mulhigh_normalised3(nn_ptr rp, nn_srcptr xp, slong xn, nn_srcptr yp, slong yn, slong n)
{
    mp_limb_pair_t ret;

    FLINT_ASSERT(n >= 1);

    if (rp == xp || rp == yp)
    {
        nn_ptr t;
        TMP_INIT;
        TMP_START;
        t = TMP_ALLOC(sizeof(ulong) * n);
        ret = _flint_mpn_mulhigh_normalised3(t, xp, xn, yp, yn, n);
        flint_mpn_copyi(rp, t, n);
        TMP_END;
        return ret;
    }

    if (xp == yp)
        ret.m1 = flint_mpn_sqrhigh(rp, xp + xn - n, n);
    else
        ret.m1 = flint_mpn_mulhigh_n(rp, xp + xn - n, yp + yn - n, n);

    if (LIMB_MSB_IS_SET(rp[n - 1]))
    {
        ret.m2 = 0;
    }
    else
    {
        ret.m2 = 1;
        mpn_lshift(rp, rp, n, 1);
        rp[0] |= (ret.m1 >> (FLINT_BITS - 1));
        ret.m1 <<= 1;
    }

    return ret;
}

/* handles aliasing */
FLINT_FORCE_INLINE
mp_limb_pair_t flint_mpn_mulhigh_normalised2(nn_ptr rp, nn_srcptr xp, nn_srcptr yp, slong n)
{
    FLINT_ASSERT(n >= 1);

    if (FLINT_HAVE_MULHIGH_NORMALISED_FUNC(n))
        return flint_mpn_mulhigh_normalised_func_tab[n](rp, xp, yp);
    else
        return _flint_mpn_mulhigh_normalised3(rp, xp, n, yp, n, n);
}

/* handles aliasing, allowing xp and yp to be shifted */
FLINT_FORCE_INLINE
mp_limb_pair_t flint_mpn_mulhigh_normalised3(nn_ptr rp, nn_srcptr xp, slong xn, nn_srcptr yp, slong yn, slong n)
{
    FLINT_ASSERT(n >= 1);

    if (FLINT_HAVE_MULHIGH_NORMALISED_FUNC(n))
        return flint_mpn_mulhigh_normalised_func_tab[n](rp, xp + xn - n, yp + yn - n);
    else
        return _flint_mpn_mulhigh_normalised3(rp, xp, xn, yp, yn, n);
}


#if FLINT_BITS == 64
#define ULP_N1 0x1.0p-64
#define ULP_N2 0x1.0p-128
#define ULP_N3 0x1.0p-192
#define ULP_N4 0x1.0p-256
#else
#define ULP_N1 0x1.0p-32
#define ULP_N2 0x1.0p-64
#define ULP_N3 0x1.0p-96
#define ULP_N4 0x1.0p-128
#endif

#define NRB_CORRECTION_A  (1.0 + 0x1.0p-50)
#define NRB_CORRECTION_B  0x1.0000000000000p-1016

int
_nrb_mul_special(nrb_ptr res, nrb_srcptr x, nrb_srcptr y, gr_ctx_t ctx)
{
    mp_limb_pair_t mul_res;
    double xerr, yerr, err;
    slong n, xn, yn, correction;
    slong xexp, yexp;
    nn_srcptr xp, yp;

    xn = NRB_N(x);
    yn = NRB_N(y);

    xerr = NRB_ERR(x);
    yerr = NRB_ERR(y);

    xexp = NRB_EXP(x);
    yexp = NRB_EXP(y);

    xp = NRB_D(x);
    yp = NRB_D(y);

    if (FLINT_UNLIKELY(xn == 0 || yn == 0))
    {
        /* 0 * anything -> 0 */
        if ((xn == 0 && xerr == 0.0) || (yn == 0 && yerr == 0.0))
        {
            return nrb_zero(res, ctx);
        }

        /* (+/- inf) * nonzero -> (+/- inf) */
        if (xerr == D_INF || yerr == D_INF)
            return nrb_zero_pm_inf(res, ctx);

        /* (x +/- xerr) * (x +/- yerr) */
        if (yn != 0)
            yerr = yp[yn - 1] * ULP_N1 + d_mul_2exp(yerr, -yn * FLINT_BITS);
        if (xn != 0)
            xerr = xp[xn - 1] * ULP_N1 + d_mul_2exp(xerr, -xn * FLINT_BITS);

        NRB_ERR(res) = (xerr * yerr) * NRB_CORRECTION_A;
        NRB_EXP(res) = xexp + yexp;
        NRB_N(res) = 0;
    }
    else
    {
        n = FLINT_MIN(xn, yn);

        double xm, ym;
        double erra, errb, errc;
        slong expa, expb, expc, rexp;

        xm = (double) xp[xn - 1];
        ym = (double) yp[yn - 1];

        /* truncation errors */
        if (xn < n)
            xerr += (xp[xn - n - 1] + 1.0) * ULP_N1;
        if (yn < n)
            yerr += (yp[xn - n - 1] + 1.0) * ULP_N1;

        erra = xm * yerr;
        expa = xexp + yexp - FLINT_BITS - yn * FLINT_BITS;
        errb = ym * xerr;
        expb = xexp + yexp - FLINT_BITS - xn * FLINT_BITS;
        errc = xerr * yerr;
        expc = xexp + yexp - xn * FLINT_BITS - yn * FLINT_BITS;

        rexp = xexp + yexp - n * FLINT_BITS;

        /* todo: clamp shifts to avoid spurious subnormals? */
        err = 0.0;
        err += d_mul_2exp(erra, -(rexp - expa));
        err += d_mul_2exp(errb, -(rexp - expb));
        err += d_mul_2exp(errc, -(rexp - expc));

        mul_res = flint_mpn_mulhigh_normalised3(NRB_D(res), xp, xn, yp, yn, n);
        correction = mul_res.m2;
        err *= (1.0 + correction);

        if (mul_res.m1 != 0 || !_mulhigh_check_exact_product(xp + xn - n, yp + yn - n, n))
        {
            err += (mul_res.m1 + 2.0 * n) * ULP_N1;
        }

        err *= NRB_CORRECTION_A;
        NRB_SGNBIT(res) = NRB_SGNBIT(x) ^ NRB_SGNBIT(y);
        NRB_N(res) = n;
        NRB_EXP(res) = xexp + yexp - correction;
        NRB_ERR(res) = err;
    }

    NRB_RETURN_FIX_RANGE(res, ctx);
}

int
nrb_mul(nrb_ptr res, nrb_srcptr x, nrb_srcptr y, nrb_ctx_t ctx)
{
    mp_limb_pair_t mul_res;
    double xerr, yerr, err;
    slong xexp, yexp;
    slong n, correction;
    int nx, ny;
    nn_srcptr xp, yp;

    n = NRB_CTX_NLIMBS(ctx);
    nx = NRB_N(x);
    ny = NRB_N(y);

    // flint_printf("DO MUL %wd %wd %wd\n", n, nx, ny);

    if (FLINT_UNLIKELY(nx != n || ny != n))
        return _nrb_mul_special(res, x, y, ctx);

    xp = NRB_D(x);
    yp = NRB_D(y);

    xerr = NRB_ERR(x);
    yerr = NRB_ERR(y);

    xexp = NRB_EXP(x);
    yexp = NRB_EXP(y);

    if (n == 1)
    {
        ulong r1, r0, x0, y0;

        x0 = NRB_D(x)[0];
        y0 = NRB_D(y)[0];

        umul_ppmm(r1, r0, x0, y0);

        if (LIMB_MSB_IS_SET(r1))
        {
            NRB_D(res)[0] = r1;
            correction = 0;
            err = (x0 * yerr + y0 * xerr + xerr * yerr + r0) * ULP_N1;
        }
        else
        {
            NRB_D(res)[0] = (r1 << 1) | (r0 >> (FLINT_BITS - 1));
            correction = 1;
            err = (x0 * yerr + y0 * xerr + xerr * yerr) * (2.0 * ULP_N1);
            err += (r0 << 1) * ULP_N1;
        }
    }
    else if (n == 2)
    {
        ulong r3, r2, r1, r0, x1, x0, y1, y0;

        x0 = NRB_D(x)[0];
        x1 = NRB_D(x)[1];
        y0 = NRB_D(y)[0];
        y1 = NRB_D(y)[1];

        err = (x1 * yerr + y1 * xerr) * ULP_N1 + (xerr * yerr) * ULP_N2;

        FLINT_MPN_MUL_2X2(r3, r2, r1, r0, x1, x0, y1, y0);

        if (LIMB_MSB_IS_SET(r3))
        {
            NRB_D(res)[0] = r2;
            NRB_D(res)[1] = r3;
            correction = 0;
            err += r1 * ULP_N1 + r0 * ULP_N2;
        }
        else
        {
            NRB_D(res)[0] = (r2 << 1) | (r1 >> (FLINT_BITS - 1));
            NRB_D(res)[1] = (r3 << 1) | (r2 >> (FLINT_BITS - 1));
            correction = 1;
            err = (2.0 * err) + (((r1 << 1) + (r0 != 0)) * ULP_N1);
        }
    }
    else
    {
        double xm, ym;

        xm = (double) xp[n - 1];
        ym = (double) yp[n - 1];

        err = (xm * yerr + ym * xerr) * ULP_N1 + d_mul_2exp(xerr * yerr, -n * FLINT_BITS);

        mul_res = flint_mpn_mulhigh_normalised2(NRB_D(res), xp, yp, n);
        correction = mul_res.m2;
        err *= (1.0 + correction);

        if (mul_res.m1 != 0 || !_mulhigh_check_exact_product(xp, yp, n))
        {
            err += (mul_res.m1 + 2.0 * n) * ULP_N1;
        }
    }

    err *= NRB_CORRECTION_A;
    NRB_SGNBIT(res) = NRB_SGNBIT(x) ^ NRB_SGNBIT(y);
    NRB_N(res) = n;
    NRB_EXP(res) = xexp + yexp - correction;
    NRB_ERR(res) = err;

    NRB_RETURN_FIX_RANGE(res, ctx);
}
