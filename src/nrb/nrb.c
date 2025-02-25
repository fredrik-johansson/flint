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

int _nrb_is_valid(nrb_srcptr x, nrb_ctx_t ctx)
{
    int err_in_normal_range;
    int exp_in_normal_range;

    if (!(NRB_ERR(x) >= 0))
    {
        NRB_EXPLAIN("Error %g not in range\n", NRB_ERR(x));
        return 0;
    }

    if (NRB_N(x) < 0 || NRB_N(x) > NRB_CTX_NLIMBS(ctx))
    {
        NRB_EXPLAIN("n %wd not in range\n", (slong) NRB_N(x));
        return 0;
    }

    err_in_normal_range = (NRB_ERR(x) >= NRB_MIN_ERR && NRB_ERR(x) <= NRB_MAX_ERR);
    exp_in_normal_range = (NRB_EXP(x) >= NRB_MIN_EXP && NRB_EXP(x) <= NRB_MAX_EXP);

    if (NRB_N(x) == 0)
    {
        int ok = (NRB_ERR(x) == 0.0) || (NRB_ERR(x) == D_INF) || (err_in_normal_range && exp_in_normal_range);

        if (!ok)
        {
            NRB_EXPLAIN("Special value: either error %g ", NRB_ERR(x));
            NRB_EXPLAIN("or exponent %wd not in range\n", NRB_EXP(x));
        }

        return ok;
    }

    if (!(NRB_ERR(x) == 0.0 || err_in_normal_range))
    {
        NRB_EXPLAIN("Error %g not in range for normal value\n", NRB_ERR(x));
        return 0;
    }

    if (!exp_in_normal_range)
    {
        NRB_EXPLAIN("Exponent %wd not in range for normal value\n", NRB_EXP(x));
        return 0;
    }

    if (NRB_SGNBIT(x) != 0 && NRB_SGNBIT(x) != 1)
    {
        NRB_EXPLAIN("Sign bit %d\n", NRB_SGNBIT(x));
        return 0;
    }

    if (!LIMB_MSB_IS_SET(NRB_D(x)[NRB_N(x) - 1]))
    {
        NRB_EXPLAIN("Limb %wu not normalised\n", NRB_D(x)[NRB_N(x) - 1]);
        return 0;
    }

    return 1;
}


int
nrb_one(nrb_ptr res, nrb_ctx_t ctx)
{
    slong n = NRB_CTX_NLIMBS(ctx);

    NRB_ERR(res) = 0.0;
    NRB_EXP(res) = 1;
    NRB_SGNBIT(res) = 0;
    NRB_N(res) = n;

    flint_mpn_zero(NRB_D(res), n - 1);
    NRB_D(res)[n - 1] = UWORD(1) << (FLINT_BITS - 1);
    return GR_SUCCESS;
}

int
nrb_neg_one(nrb_ptr res, nrb_ctx_t ctx)
{
    slong n = NRB_CTX_NLIMBS(ctx);

    NRB_ERR(res) = 0.0;
    NRB_EXP(res) = 1;
    NRB_SGNBIT(res) = 1;
    NRB_N(res) = n;

    flint_mpn_zero(NRB_D(res), n - 1);
    NRB_D(res)[n - 1] = UWORD(1) << (FLINT_BITS - 1);
    return GR_SUCCESS;
}

int
nrb_set_ui(nrb_ptr res, ulong x, nrb_ctx_t ctx)
{
    NRB_ERR(res) = 0.0;

    if (x == 0)
    {
        NRB_N(res) = 0;
    }
    else
    {
        slong n = NRB_CTX_NLIMBS(ctx);
        unsigned int l = flint_clz(x);

        NRB_EXP(res) = FLINT_BITS - l;
        NRB_SGNBIT(res) = 0;
        NRB_N(res) = n;

        flint_mpn_zero(NRB_D(res), n - 1);
        NRB_D(res)[n - 1] = x << l;
    }

    return GR_SUCCESS;
}

int
nrb_set_si(nrb_ptr res, slong x, nrb_ctx_t ctx)
{
    NRB_ERR(res) = 0.0;

    if (x == 0)
    {
        NRB_N(res) = 0;
    }
    else
    {
        slong n = NRB_CTX_NLIMBS(ctx);
        ulong y = FLINT_UABS(x);
        unsigned int l = flint_clz(y);

        NRB_EXP(res) = FLINT_BITS - l;
        NRB_SGNBIT(res) = x < 0;
        NRB_N(res) = n;

        flint_mpn_zero(NRB_D(res), n - 1);
        NRB_D(res)[n - 1] = y << l;
    }

    return GR_SUCCESS;
}

int
nrb_zero_pm_inf(nrb_ptr res, nrb_ctx_t ctx)
{
    NRB_ERR(res) = D_INF;
    NRB_N(res) = 0;
    return GR_SUCCESS;
}

int
nrb_zero_pm_2exp_si(nrb_ptr res, slong e, nrb_ctx_t ctx)
{
    if (e >= NRB_MAX_EXP)
        return nrb_zero_pm_inf(res, ctx);

    NRB_ERR(res) = 1.0;
    NRB_N(res) = 0;
    NRB_EXP(res) = FLINT_MAX(e, NFLOAT_MIN_EXP);

    return GR_SUCCESS;
}

int
nrb_zero_pm_one(nrb_ptr res, nrb_ctx_t ctx)
{
    NRB_ERR(res) = 1.0;
    NRB_N(res) = 0;
    NRB_EXP(res) = 0;

    return GR_SUCCESS;
}


int
nrb_write_debug(gr_stream_t out, nrb_srcptr x, nrb_ctx_t ctx)
{
    slong i;

    gr_stream_write(out, "\n{");

    gr_stream_write(out, "error ");
    {
        char t[64];
        flint_sprintf(t, "%g", NRB_ERR(x));
        gr_stream_write(out, t);
    }

    gr_stream_write(out, ", exp ");
    gr_stream_write_si(out, NRB_EXP(x));
    gr_stream_write(out, ", sgnbit ");
    gr_stream_write_si(out, NRB_SGNBIT(x));
    gr_stream_write(out, ", n ");
    gr_stream_write_si(out, NRB_N(x));
    gr_stream_write(out, ", digits ");

    for (i = 0; i < NRB_N(x); i++)
    {
        gr_stream_write_ui(out, NRB_D(x)[i]);
        if (i + 1 < NRB_N(x))
            gr_stream_write(out, " ");
    }

    gr_stream_write(out, "}\n");

    return GR_SUCCESS;
}

int
nrb_print_debug(nrb_srcptr x, nrb_ctx_t ctx)
{
    gr_stream_t out;
    gr_stream_init_file(out, stdout);
    nrb_write_debug(out, x, ctx);
    fflush(stdout);

    return GR_SUCCESS;
}

int
nrb_write(gr_stream_t out, nrb_srcptr x, nrb_ctx_t ctx)
{
    gr_ctx_t arb_ctx;
    arb_t t;
    int status;

    gr_ctx_init_real_arb(arb_ctx, NRB_CTX_PREC(ctx));
    arb_init(t);
    nrb_get_arb(t, x, ctx);
    status = gr_write(out, t, arb_ctx);
    arb_clear(t);
    gr_ctx_clear(arb_ctx);
    return status;
}

int
nrb_randtest(nrb_ptr res, flint_rand_t state, nrb_ctx_t ctx)
{
    arb_t t;
    int status;

    arb_init(t);
    arb_randtest(t, state, NRB_CTX_PREC(ctx), n_randint(state, 2) ? 2 : 10);
    status = nrb_set_arb(res, t, ctx);
    arb_clear(t);

    return status;
}

int
nrb_randtest_ebits(nrb_ptr res, flint_rand_t state, slong ebits, nrb_ctx_t ctx)
{
    arb_t t;
    int status;

    arb_init(t);
    arb_randtest(t, state, NRB_CTX_PREC(ctx), ebits);
    status = nrb_set_arb(res, t, ctx);
    arb_clear(t);

    return status;
}

truth_t nrb_is_zero(nrb_srcptr x, nrb_ctx_t ctx)
{
    double err;
    int e;
    slong n = NRB_N(x);
    double d;

    err = NRB_ERR(x);

    if (err == 0.0)
        return (n == 0) ? T_TRUE : T_FALSE;

    if (err == D_INF || n == 0)
        return T_UNKNOWN;

    if (err < d_2exp_inrange(FLINT_BITS - 1))
        return T_FALSE;

    err = d_pos_normal_frexp(err, &e);

    if (e < FLINT_BITS * n)
        return T_FALSE;

    if (e == FLINT_BITS * n)
    {
        d = NRB_D(x)[n - 1] * ULP_N1;

#if FLINT_BITS == 32
            if (n >= 2)
                d += NRB_D(x)[n - 2] * ULP_N2;
#endif

        if (err < d * (1.0 - 0x1.0p-50))
            return T_FALSE;
    }

    return T_UNKNOWN;
}

truth_t nrb_is_one(nrb_srcptr x, nrb_ctx_t ctx)
{
    gr_ctx_t arb_ctx;
    arb_t t;
    truth_t ans;

    gr_ctx_init_real_arb(arb_ctx, NRB_CTX_PREC(ctx));
    arb_init(t);
    nrb_get_arb(t, x, ctx);
    ans = gr_is_one(t, arb_ctx);
    arb_clear(t);
    gr_ctx_clear(arb_ctx);
    return ans;
}

truth_t nrb_is_neg_one(nrb_srcptr x, nrb_ctx_t ctx)
{
    gr_ctx_t arb_ctx;
    arb_t t;
    truth_t ans;

    gr_ctx_init_real_arb(arb_ctx, NRB_CTX_PREC(ctx));
    arb_init(t);
    nrb_get_arb(t, x, ctx);
    ans = gr_is_neg_one(t, arb_ctx);
    arb_clear(t);
    gr_ctx_clear(arb_ctx);
    return ans;
}

truth_t nrb_equal(nrb_srcptr x, nrb_srcptr y, nrb_ctx_t ctx)
{
    gr_ctx_t arb_ctx;
    arb_t t, u;
    truth_t ans;

    gr_ctx_init_real_arb(arb_ctx, NRB_CTX_PREC(ctx));
    arb_init(t);
    arb_init(u);
    nrb_get_arb(t, x, ctx);
    nrb_get_arb(u, y, ctx);
    ans = gr_equal(t, u, arb_ctx);
    arb_clear(t);
    arb_clear(u);
    gr_ctx_clear(arb_ctx);
    return ans;
}

int nrb_set_fmpz(nrb_ptr res, const fmpz_t x, nrb_ctx_t ctx)
{
    if (!COEFF_IS_MPZ(*x))
    {
        return nrb_set_si(res, *x, ctx);
    }
    else
    {
        arb_t t;
        arb_init(t);
        arb_set_fmpz(t, x);
        int status = nrb_set_arb(res, t, ctx);
        arb_clear(t);
        return status;

        /* todo: optimize */
/*
        slong bits, width;

        bits = width = fmpz_bits(x);

        if (width > NRB_CTX_PREC(ctx))
            width -= fmpz_val2(x);

        if (width <= NRB_CTX_PREC(ctx))
        {
        }
        else
        {
        }
        NRB_ASSERT_VALID(res, ctx);

        return GR_SUCCESS;
*/
    }
}

void
nrb_swap(nrb_ptr x, nrb_ptr y, nrb_ctx_t ctx)
{
    slong i, n = NRB_CTX_DATA_NLIMBS(ctx);

    for (i = 0; i < n; i++)
        FLINT_SWAP(ulong, NRB_DATA(x)[i], NRB_DATA(y)[i]);
}

int
nrb_set(nrb_ptr res, nrb_srcptr x, nrb_ctx_t ctx)
{
    slong i, n = NRB_CTX_DATA_NLIMBS(ctx);

    for (i = 0; i < n; i++)
        NRB_DATA(res)[i] = NRB_DATA(x)[i];

    return GR_SUCCESS;
}

void
nrb_set_shallow(nrb_ptr res, nrb_srcptr x, nrb_ctx_t ctx)
{
    slong i, n = NRB_CTX_DATA_NLIMBS(ctx);

    for (i = 0; i < n; i++)
        NRB_DATA(res)[i] = NRB_DATA(x)[i];
}

int
nrb_neg(nrb_ptr res, nrb_srcptr x, nrb_ctx_t ctx)
{
    if (res != x)
    {
        slong i, n = NRB_CTX_DATA_NLIMBS(ctx);

        for (i = 0; i < n; i++)
            NRB_DATA(res)[i] = NRB_DATA(x)[i];
    }

    NRB_SGNBIT(res) = !NRB_SGNBIT(res);

    return GR_SUCCESS;
}

int
nrb_abs(nrb_ptr res, nrb_srcptr x, nrb_ctx_t ctx)
{
    if (res != x)
    {
        slong i, n = NRB_CTX_DATA_NLIMBS(ctx);

        for (i = 0; i < n; i++)
            NRB_DATA(res)[i] = NRB_DATA(x)[i];
    }

    NRB_SGNBIT(res) = 0;

    return GR_SUCCESS;
}

#define WRAP_ARB_FUNC_1(func) \
    arb_t t; \
    arb_init(t); \
    nrb_get_arb(t, x, ctx); \
    func(t, t, NRB_CTX_PREC(ctx)); \
    int status = nrb_set_arb(res, t, ctx); \
    arb_clear(t); \
    return status;

#define WRAP_ARB_FUNC_2(func) \
    arb_t t, u; \
    arb_init(t); \
    arb_init(u); \
    nrb_get_arb(t, x, ctx); \
    nrb_get_arb(u, y, ctx); \
    func(t, t, u, NRB_CTX_PREC(ctx)); \
    int status = nrb_set_arb(res, t, ctx); \
    arb_clear(t); \
    arb_clear(u); \
    return status;

int nrb_inv(nrb_ptr res, nrb_srcptr x, nrb_ctx_t ctx) { WRAP_ARB_FUNC_1(arb_inv) }
int nrb_div(nrb_ptr res, nrb_srcptr x, nrb_srcptr y, nrb_ctx_t ctx) { WRAP_ARB_FUNC_2(arb_div) }
int nrb_sqrt(nrb_ptr res, nrb_srcptr x, nrb_ctx_t ctx) { WRAP_ARB_FUNC_1(arb_sqrt) }
int nrb_rsqrt(nrb_ptr res, nrb_srcptr x, nrb_ctx_t ctx) { WRAP_ARB_FUNC_1(arb_rsqrt) }
int nrb_pow(nrb_ptr res, nrb_srcptr x, nrb_srcptr y, nrb_ctx_t ctx) { WRAP_ARB_FUNC_2(arb_pow) }
int nrb_exp(nrb_ptr res, nrb_srcptr x, nrb_ctx_t ctx) { WRAP_ARB_FUNC_1(arb_exp) }
int nrb_log(nrb_ptr res, nrb_srcptr x, nrb_ctx_t ctx) { WRAP_ARB_FUNC_1(arb_log) }
int nrb_sin(nrb_ptr res, nrb_srcptr x, nrb_ctx_t ctx) { WRAP_ARB_FUNC_1(arb_sin) }
int nrb_cos(nrb_ptr res, nrb_srcptr x, nrb_ctx_t ctx) { WRAP_ARB_FUNC_1(arb_cos) }

int
nrb_get_interval_arf(arf_t a, arf_t b, nrb_srcptr x, nrb_ctx_t ctx, slong prec)
{
    if (NRB_ERR(x) == D_INF)
    {
        arf_neg_inf(a);
        arf_pos_inf(b);
    }
    else
    {
        arf_t mid, rad;
        slong n = NRB_N(x);

        arf_init(mid);
        arf_init(rad);

        if (n == 0)
        {
            arf_zero(mid);
        }
        else
        {
            arf_set_mpn(mid, NRB_D(x), n, NRB_SGNBIT(x));
            arf_mul_2exp_si(mid, mid, NRB_EXP(x) - n * FLINT_BITS);
        }

        arf_set_d(rad, NRB_ERR(x));
        arf_mul_2exp_si(rad, rad, NRB_EXP(x) - n * FLINT_BITS);

        arf_sub(a, mid, rad, prec, ARF_RND_FLOOR);
        arf_add(b, mid, rad, prec, ARF_RND_CEIL);

        arf_clear(mid);
        arf_clear(rad);
    }

    return GR_SUCCESS;
}

int
nrb_get_arb(arb_t res, nrb_srcptr x, nrb_ctx_t ctx)
{
    if (NRB_ERR(x) == D_INF)
    {
        arb_zero_pm_inf(res);
    }
    else
    {
        slong n = NRB_N(x);

        if (n == 0)
        {
            arf_zero(arb_midref(res));
        }
        else
        {
            arf_set_mpn(arb_midref(res), NRB_D(x), n, NRB_SGNBIT(x));
            arf_mul_2exp_si(arb_midref(res), arb_midref(res), NRB_EXP(x) - n * FLINT_BITS);
        }

        if (NRB_ERR(x) == 0.0)
        {
            mag_zero(arb_radref(res));
        }
        else
        {
            mag_set_d(arb_radref(res), NRB_ERR(x));
            mag_mul_2exp_si(arb_radref(res), arb_radref(res), NRB_EXP(x) - n * FLINT_BITS);
        }
    }

    return GR_SUCCESS;
}

int nrb_get_mid_arf(arf_t res, nrb_srcptr x, nrb_ctx_t ctx)
{
    if (NRB_ERR(x) == D_INF)
    {
        arf_zero(res);
    }
    else
    {
        slong n = NRB_N(x);

        if (n == 0)
        {
            arf_zero(res);
        }
        else
        {
            arf_set_mpn(res, NRB_D(x), n, NRB_SGNBIT(x));
            arf_mul_2exp_si(res, res, NRB_EXP(x) - n * FLINT_BITS);
        }
    }

    return GR_SUCCESS;
}

int nrb_get_rad_arf(arf_t res, nrb_srcptr x, nrb_ctx_t ctx)
{
    if (NRB_ERR(x) == D_INF)
    {
        arf_pos_inf(res);
    }
    else if (NRB_ERR(x) == 0.0)
    {
        arf_zero(res);
    }
    else
    {
        arf_set_d(res, NRB_ERR(x));
        arf_mul_2exp_si(res, res, NRB_EXP(x) - NRB_N(x) * FLINT_BITS);
    }

    return GR_SUCCESS;
}

int
nrb_set_arb(nrb_ptr res, const arb_t x, nrb_ctx_t ctx)
{
    slong xexp, rexp;
    double err;
    slong n = NRB_CTX_NLIMBS(ctx);
    slong xn;
    nn_srcptr xp;

    if (!arf_is_finite(arb_midref(x)))
    {
        if (arf_is_nan(arb_midref(x)))
            return GR_UNABLE;
        else
            return GR_DOMAIN;
    }

    if (!mag_is_finite(arb_radref(x)))
        return nrb_zero_pm_inf(res, ctx);

    if (!ARB_IS_LAGOM(x))
        return GR_UNABLE;

//    arb_print(x); flint_printf("\n");

    if (arf_is_zero(arb_midref(x)))
    {
        /* 0 */
        if (mag_is_zero(arb_radref(x)))
            return nrb_zero(res, ctx);

        /* 0 +/- r */
        NRB_EXP(res) = MAG_EXP(arb_radref(x)) + MAG_BITS;
        NRB_N(res) = 0;
        NRB_ERR(res) = d_mul_2exp_inrange2(MAG_MAN(arb_radref(x)), -MAG_BITS);
    }
    else
    {
        ARF_GET_MPN_READONLY(xp, xn, arb_midref(x));
        xexp = ARF_EXP(arb_midref(x));

        if (mag_is_zero(arb_radref(x)))
        {
            err = 0.0;
        }
        else
        {
            /* scale err to n-limb ulp of x */

            rexp = MAG_EXP(arb_radref(x)) - xexp + n * FLINT_BITS;
            err = MAG_MAN(arb_radref(x));

            /* sanity checks to make sure we don't overflow the double
               exponent range; the actual range will be fixed later
               when we call _nrb_fix_range */
            rexp = FLINT_MAX(rexp, -128);
            if (rexp > 768)
            {
                slong trim_limbs;

                trim_limbs = rexp / FLINT_BITS - 1;

                if (trim_limbs >= n)
                {
                    /* XXX: checkme */
                    NRB_EXP(res) = MAG_EXP(arb_radref(x));
                    NRB_N(res) = 0;
                    NRB_SGNBIT(res) = 0;
                    NRB_ERR(res) = d_mul_2exp_inrange2(err * NRB_CORRECTION_A, -MAG_BITS);
                    return _nrb_fix_range(res, ctx);
                }
                else
                {
                    n -= trim_limbs;
                    rexp -= trim_limbs * FLINT_BITS;
                    err = d_mul_2exp_inrange2(err, rexp - MAG_BITS);
                }
            }
            else
            {
                err = d_mul_2exp_inrange2(err, rexp - MAG_BITS);
            }
        }

        if (xn <= n)
        {
            flint_mpn_zero(NRB_D(res), n - xn);
            flint_mpn_copyi(NRB_D(res) + n - xn, xp, xn);
        }
        else
        {
            flint_mpn_copyi(NRB_D(res), xp + xn - n, n);
            err += ((double) xp[xn - n - 1] + 1.0) * ULP_N1;
            err *= NRB_CORRECTION_A;
        }

        NRB_ERR(res) = err;
        NRB_N(res) = n;
        NRB_SGNBIT(res) = ARF_SGNBIT(arb_midref(x));
        NRB_EXP(res) = xexp;
    }

    return _nrb_fix_range(res, ctx);
}

int
nrb_inplace_add_error_d_2exp_si(nrb_ptr res, double err, slong err_exp, nrb_ctx_t ctx)
{
    slong xexp, rexp, n;
    double xerr;

    FLINT_ASSERT(err >= 0.0);

    if (err == 0.0 || NRB_ERR(res) == D_INF)
        return GR_SUCCESS;

    if (err == D_INF || err_exp > NRB_MAX_EXP)
        return nrb_zero_pm_inf(res, ctx);

    err_exp = FLINT_MAX(err_exp, NRB_MIN_EXP);

    if (err < NRB_MIN_ERR)
    {
        int e;
        err = frexp(err, &e);
        err_exp += e;
    }
    else
    {
        /* todo: safely avoid this in common cases */
        int e;
        err = d_pos_normal_frexp(err, &e);
        err_exp += e;
    }

    n = NRB_N(res);
    xerr = NRB_ERR(res);
    xexp = NRB_EXP(res);

    if (n == 0 && xerr == 0.0)
    {
        NRB_ERR(res) = err;
        NRB_EXP(res) = err_exp;
        return GR_SUCCESS;
    }

    /* scale err to n-limb ulp of x */

    rexp = err_exp - xexp + n * FLINT_BITS;

    /* sanity checks to make sure we don't overflow the double
        exponent range; the actual range will be fixed later
        when we call _nrb_fix_range */
    rexp = FLINT_MAX(rexp, -128);
    if (rexp > 768)
    {
        slong trim_limbs;

        trim_limbs = rexp / FLINT_BITS - 1;

        if (trim_limbs >= n)
        {
            FLINT_ASSERT(xerr < d_2exp_inrange(768 - 32));

            NRB_EXP(res) = err_exp;
            NRB_N(res) = 0;
            NRB_SGNBIT(res) = 0;
            NRB_ERR(res) = err * NRB_CORRECTION_A;
            return _nrb_fix_range(res, ctx);
        }
        else
        {
            n -= trim_limbs;
            rexp -= trim_limbs * FLINT_BITS;
            err = d_mul_2exp(err, rexp);
            /* todo: avoid denormal here */
            xerr = d_mul_2exp(xerr, -trim_limbs * FLINT_BITS);
            NRB_ERR(res) = (xerr + err) * NRB_CORRECTION_A;
            flint_mpn_copyi(NRB_D(res), NRB_D(res) + trim_limbs, n);
            NRB_N(res) = n;
        }
    }
    else
    {
        err = d_mul_2exp(err, rexp);
        NRB_ERR(res) = (xerr + err) * NRB_CORRECTION_A;
    }

    return _nrb_fix_range(res, ctx);
}
