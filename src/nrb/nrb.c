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
nrb_write_debug(gr_stream_t out, nrb_srcptr x, gr_ctx_t ctx)
{
    slong i;

    gr_stream_write(out, "\n{");

    gr_stream_write(out, "rad_exp ");
    if (NRB_RAD_EXP(x) == NFLOAT_EXP_ZERO)
        gr_stream_write(out, "ZERO");
    else if (NRB_RAD_EXP(x) == NFLOAT_EXP_POS_INF)
        gr_stream_write(out, "INF");
    else
        gr_stream_write_si(out, NRB_RAD_EXP(x));

    gr_stream_write(out, ", rad_man ");
    gr_stream_write_ui(out, NRB_RAD_D(x));
    gr_stream_write(out, ", exp ");
    if (NRB_EXP(x) == NFLOAT_EXP_ZERO)
        gr_stream_write(out, "ZERO");
    else
        gr_stream_write_si(out, NRB_EXP(x));
    gr_stream_write(out, ", sgnbit ");
    gr_stream_write_ui(out, NRB_SGNBIT(x));
    gr_stream_write(out, ", digits ");

    for (i = 0; i < NFLOAT_CTX_NLIMBS(ctx); i++)
    {
        gr_stream_write_ui(out, NRB_D(x)[i]);
        if (i + 1 < NFLOAT_CTX_NLIMBS(ctx))
            gr_stream_write(out, " ");
    }

    gr_stream_write(out, "}\n");

    return GR_SUCCESS;
}

int
nrb_print_debug(nrb_srcptr x, gr_ctx_t ctx)
{
    gr_stream_t out;
    gr_stream_init_file(out, stdout);
    nrb_write_debug(out, x, ctx);
    fflush(stdout);

    return GR_SUCCESS;
}

int
nrb_write(gr_stream_t out, nrb_srcptr x, gr_ctx_t ctx)
{
    gr_ctx_t arb_ctx;
    arb_t t;
    int status;

    gr_ctx_init_real_arb(arb_ctx, NFLOAT_CTX_PREC(ctx));
    arb_init(t);
    nrb_get_arb(t, x, ctx);
    status = gr_write(out, t, arb_ctx);
    arb_clear(t);
    gr_ctx_clear(arb_ctx);
    return status;
}

int
nrb_randtest(nrb_ptr res, flint_rand_t state, gr_ctx_t ctx)
{
    arb_t t;
    int status;

    arb_init(t);
    arb_randtest(t, state, NFLOAT_CTX_PREC(ctx), n_randint(state, 2) ? 2 : 10);
    status = nrb_set_arb(res, t, ctx);
    arb_clear(t);
    return status;
}

truth_t nrb_is_zero(nrb_srcptr x, gr_ctx_t ctx)
{
    gr_ctx_t arb_ctx;
    arb_t t;
    truth_t ans;

    gr_ctx_init_real_arb(arb_ctx, NFLOAT_CTX_PREC(ctx));
    arb_init(t);
    nrb_get_arb(t, x, ctx);
    ans = gr_is_zero(t, arb_ctx);
    arb_clear(t);
    gr_ctx_clear(arb_ctx);
    return ans;
}

truth_t nrb_is_one(nrb_srcptr x, gr_ctx_t ctx)
{
    gr_ctx_t arb_ctx;
    arb_t t;
    truth_t ans;

    gr_ctx_init_real_arb(arb_ctx, NFLOAT_CTX_PREC(ctx));
    arb_init(t);
    nrb_get_arb(t, x, ctx);
    ans = gr_is_one(t, arb_ctx);
    arb_clear(t);
    gr_ctx_clear(arb_ctx);
    return ans;
}

truth_t nrb_is_neg_one(nrb_srcptr x, gr_ctx_t ctx)
{
    gr_ctx_t arb_ctx;
    arb_t t;
    truth_t ans;

    gr_ctx_init_real_arb(arb_ctx, NFLOAT_CTX_PREC(ctx));
    arb_init(t);
    nrb_get_arb(t, x, ctx);
    ans = gr_is_neg_one(t, arb_ctx);
    arb_clear(t);
    gr_ctx_clear(arb_ctx);
    return ans;
}

truth_t nrb_equal(nrb_srcptr x, nrb_srcptr y, gr_ctx_t ctx)
{
    gr_ctx_t arb_ctx;
    arb_t t, u;
    truth_t ans;

    gr_ctx_init_real_arb(arb_ctx, NFLOAT_CTX_PREC(ctx));
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

int nrb_set_fmpz(nrb_ptr res, const fmpz_t x, gr_ctx_t ctx)
{
    if (!COEFF_IS_MPZ(*x))
    {
        return nrb_set_si(res, *x, ctx);
    }
    else
    {
        /* todo: optimize */
        slong bits, width;

        bits = width = fmpz_bits(x);

        if (width > NFLOAT_CTX_PREC(ctx))
            width -= fmpz_val2(x);

        nfloat_set_fmpz(NRB_MID(res), x, ctx);

        if (width <= NFLOAT_CTX_PREC(ctx))
        {
            NRB_ZERO_RAD(res);
        }
        else
        {
            NRB_RAD_EXP(res) = bits - NFLOAT_CTX_PREC(ctx) + 1;
            NRB_RAD_D(res) = MAG1_ONE_HALF;
        }

        NRB_ASSERT_VALID(res, ctx);

        return GR_SUCCESS;
    }
}

int
nrb_zero_pm_inf(nrb_ptr res, gr_ctx_t ctx)
{
    MAG1_INF(NRB_RAD(res));
    NRB_ZERO_MID(res);
    return GR_SUCCESS;
}

/* Enclose [-2^exp, 2^exp] = [0.5*2^(exp+1)] */
int
nrb_zero_pm_2exp_si(nrb_ptr res, slong exp, gr_ctx_t ctx)
{
    if (exp >= NFLOAT_MAX_EXP)
        return nrb_zero_pm_inf(res, ctx);

    NRB_RAD_EXP(res) = FLINT_MAX(exp + 1, NFLOAT_MIN_EXP);
    NRB_RAD_D(res) = MAG_ONE_HALF;
    NRB_ZERO_MID(res);

    NRB_ASSERT_VALID(res, ctx);

    return GR_SUCCESS;
}

/* Enclose (-1)^sgnbit * {[0, 2^exp] = [0.5*2^exp +/- 0.5*2^exp]} */
int
nrb_zero_to_2exp_si_sgnbit(nrb_ptr res, slong exp, int sgnbit, gr_ctx_t ctx)
{
    if (exp > NFLOAT_MAX_EXP)
        return nrb_zero_pm_inf(res, ctx);

    exp = FLINT_MAX(exp, NFLOAT_MIN_EXP);

    NRB_RAD_EXP(res) = exp;
    NRB_RAD_D(res) = MAG_ONE_HALF;

    NRB_EXP(res) = exp;
    NRB_SGNBIT(res) = sgnbit;
    flint_mpn_zero(NRB_D(res), NFLOAT_CTX_NLIMBS(ctx) - 1);
    NRB_D(res)[NFLOAT_CTX_NLIMBS(ctx) - 1] = UWORD(1) << (FLINT_BITS - 1);

    NRB_ASSERT_VALID(res, ctx);

    return GR_SUCCESS;
}

void
nrb_swap(nrb_ptr x, nrb_ptr y, gr_ctx_t ctx)
{
    slong i, n = NRB_CTX_DATA_NLIMBS(ctx);

    for (i = 0; i < n; i++)
        FLINT_SWAP(ulong, NRB_DATA(x)[i], NRB_DATA(y)[i]);
}

int
nrb_set(nrb_ptr res, nrb_srcptr x, gr_ctx_t ctx)
{
    slong i, n = NRB_CTX_DATA_NLIMBS(ctx);

    for (i = 0; i < n; i++)
        NRB_DATA(res)[i] = NRB_DATA(x)[i];

    return GR_SUCCESS;
}

void
nrb_set_shallow(nrb_ptr res, nrb_srcptr x, gr_ctx_t ctx)
{
    slong i, n = NRB_CTX_DATA_NLIMBS(ctx);

    for (i = 0; i < n; i++)
        NRB_DATA(res)[i] = NRB_DATA(x)[i];
}

int
nrb_neg(nrb_ptr res, nrb_srcptr x, gr_ctx_t ctx)
{
    if (res != x)
    {
        slong i, n = NRB_CTX_DATA_NLIMBS(ctx);

        for (i = 0; i < n; i++)
            NRB_DATA(res)[i] = NRB_DATA(x)[i];
    }

    if (!NFLOAT_IS_ZERO(NRB_MID(res)))
        NRB_SGNBIT(res) = !NRB_SGNBIT(res);

    return GR_SUCCESS;
}

/* todo: improve enclosure */
int
nrb_abs(nrb_ptr res, nrb_srcptr x, gr_ctx_t ctx)
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

int
nrb_mid_get_arf(arf_t res, nrb_srcptr x, gr_ctx_t ctx)
{
    return nfloat_get_arf(res, NRB_MID(x), ctx);
}

int
nrb_rad_get_mag(mag_t res, nrb_srcptr x, gr_ctx_t ctx)
{
    if (NRB_IS_UNBOUNDED(x))
    {
        mag_inf(res);
    }
    else if (NRB_RAD_IS_ZERO(x))
    {
        mag_zero(res);
    }
    else
    {
        MAG_MAN(res) = NRB_RAD_D(x);
        fmpz_set_si(MAG_EXPREF(res), NRB_RAD_EXP(x));
    }

    return GR_SUCCESS;
}

int
nrb_get_arb(arb_t res, nrb_srcptr x, gr_ctx_t ctx)
{
    nrb_mid_get_arf(arb_midref(res), x, ctx);
    nrb_rad_get_mag(arb_radref(res), x, ctx);
    return GR_SUCCESS;
}

int
nrb_rad_get_arf(arf_t res, nrb_srcptr x, gr_ctx_t ctx)
{
    mag_t t;
    mag_init(t);
    nrb_rad_get_mag(t, x, ctx);
    arf_set_mag(res, t);
    /* exponent is small; no need to clear */
    return GR_SUCCESS;
}

int
nrb_set_arf(nrb_ptr res, const arf_t x, gr_ctx_t ctx)
{
    if (ARF_IS_SPECIAL(x))
    {
        if (arf_is_zero(x))
            return nrb_zero(res, ctx);
        else
            return GR_UNABLE;
    }

    if (fmpz_cmp_si(ARF_EXPREF(x), NFLOAT_MAX_EXP) > 0)
        return nrb_zero_pm_inf(res, ctx);

    if (fmpz_cmp_si(ARF_EXPREF(x), NFLOAT_MIN_EXP) < 0)
        return nrb_zero_to_2exp_si_sgnbit(res, NFLOAT_MIN_EXP, ARF_SGNBIT(x), ctx);

    /* nfloat_set_arf simply truncates the mantissa, so the
       following will neither overflow nor underflow */
    GR_MUST_SUCCEED(nfloat_set_arf(NRB_MID(res), x, ctx));

    if (arf_bits(x) <= NFLOAT_CTX_PREC(ctx))
    {
        NRB_ZERO_RAD(res);
    }
    else
    {
        /* 1 ulp error; radius can underflow but not overflow */
        slong exp = ARF_EXP(x) - NFLOAT_CTX_PREC(ctx);

        FLINT_ASSERT(exp < NFLOAT_MAX_EXP);

        NRB_RAD_EXP(res) = FLINT_MAX(exp + 1, NFLOAT_MIN_EXP);
        NRB_RAD_D(res) = MAG_ONE_HALF;
    }

    NRB_ASSERT_VALID(res, ctx);

    return GR_SUCCESS;
}

int
nrb_rad_add_mag(nrb_ptr res, const mag_t x, gr_ctx_t ctx)
{
    if (mag_is_zero(x))
        return GR_SUCCESS;

    if (mag_is_inf(x) || fmpz_cmp_si(MAG_EXPREF(x), NFLOAT_MAX_EXP) > 0)
        return nrb_zero_pm_inf(res, ctx);  /* deliberately modifies midpoint too */

    mag1_t x1;

    if (fmpz_cmp_si(MAG_EXPREF(x), NFLOAT_MIN_EXP) < 0)
    {
        MAG1_EXP(x1) = NFLOAT_MIN_EXP + 1;
        MAG1_MAN(x1) = MAG_ONE_HALF;
    }
    else
    {
        MAG1_EXP(x1) = MAG_EXP(x);
        MAG1_MAN(x1) = MAG_MAN(x);
    }

    mag1_add(NRB_RAD(res), NRB_RAD(res), x1);

    NRB_ASSERT_VALID(res, ctx);

    return GR_SUCCESS;
}

int
nrb_set_arb(nrb_ptr res, const arb_t x, gr_ctx_t ctx)
{
    int status;

    status = nrb_set_arf(res, arb_midref(x), ctx);
    status |= nrb_rad_add_mag(res, arb_radref(x), ctx);

    NRB_ASSERT_VALID(res, ctx);

    return status;
}

int
nrb_add(nrb_ptr res, nrb_srcptr x, nrb_srcptr y, gr_ctx_t ctx)
{
    mag1_t rad;

    mag1_add(rad, NRB_RAD(x), NRB_RAD(y));

    if (MAG1_IS_INF(rad))
        return nrb_zero_pm_inf(res, ctx);

    int subtract = (NRB_SGNBIT(x) != NRB_SGNBIT(y));

    if (nfloat_add(NRB_MID(res), NRB_MID(x), NRB_MID(y), ctx) != GR_SUCCESS)
    {
        if (!subtract)   /* overflow */
            return nrb_zero_pm_inf(res, ctx);
        else             /* underflow */
            mag1_add_2exp_si(rad, rad, NFLOAT_MIN_EXP + 1);
    }
    else if (!NFLOAT_IS_ZERO(NRB_MID(res)))
    {
        mag1_add_2exp_si(rad, rad, NRB_EXP(res) - NFLOAT_CTX_PREC(ctx) + 1);
    }

    *NRB_RAD(res) = *rad;
    NRB_ASSERT_VALID(res, ctx);

    return GR_SUCCESS;
}

int
nrb_sub(nrb_ptr res, nrb_srcptr x, nrb_srcptr y, gr_ctx_t ctx)
{
    mag1_t rad;

    mag1_add(rad, NRB_RAD(x), NRB_RAD(y));

    if (MAG1_IS_INF(rad))
        return nrb_zero_pm_inf(res, ctx);

    int subtract = (NRB_SGNBIT(x) == NRB_SGNBIT(y));

    if (nfloat_sub(NRB_MID(res), NRB_MID(x), NRB_MID(y), ctx) != GR_SUCCESS)
    {
        if (!subtract)   /* overflow */
            return nrb_zero_pm_inf(res, ctx);
        else             /* underflow */
            mag1_add_2exp_si(rad, rad, NFLOAT_MIN_EXP + 1);
    }
    else if (!NFLOAT_IS_ZERO(NRB_MID(res)))
    {
        mag1_add_2exp_si(rad, rad, NRB_EXP(res) - NFLOAT_CTX_PREC(ctx) + 1);
    }

    *NRB_RAD(res) = *rad;
    NRB_ASSERT_VALID(res, ctx);

    return GR_SUCCESS;
}

#define _MAG1_ADJUST_ONE_TOO_LARGE_1(xexp, xman) \
    do { \
        ulong __t = (xman) >> MAG1_BITS; \
        (xman) = ((xman) >> __t) + (__t); \
        (xexp) += __t; \
    } while (0)


/* does not check for overflow */
static void
_nrb_unsafe_mid_get_mag1(mag1_t res, nrb_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_ZERO(NRB_MID(x)))
    {
        MAG1_ZERO(res);
    }
    else
    {
        MAG1_EXP(res) = NRB_EXP(x);
        MAG1_MAN(res) = (NRB_D(x)[NFLOAT_CTX_NLIMBS(ctx) - 1] >> (FLINT_BITS - MAG1_BITS)) + UWORD(1);
        _MAG1_ADJUST_ONE_TOO_LARGE_1(MAG1_EXP(res), MAG1_MAN(res));
    }
}

static int
nrb_fix_rad_underflow_overflow(nrb_ptr res, gr_ctx_t ctx)
{
    if (!MAG1_IS_SPECIAL(NRB_RAD(res)))
    {
        if (MAG1_EXP(NRB_RAD(res)) < NFLOAT_MIN_EXP)
        {
            MAG1_EXP(NRB_RAD(res)) = NFLOAT_MIN_EXP;
            MAG1_MAN(NRB_RAD(res)) = MAG1_ONE_HALF;
            return GR_SUCCESS;
        }
        else if (MAG1_EXP(NRB_RAD(res)) > NFLOAT_MAX_EXP)
        {
            return nrb_zero_pm_inf(res, ctx);
        }
    }

    return GR_SUCCESS;
}

/* am * br + ar * bm + ar * br = am * br + ar (bm + br) */

int
nrb_mul(nrb_ptr res, nrb_srcptr x, nrb_srcptr y, gr_ctx_t ctx)
{
    mag1_t rad;
    mag1_t xm;
    mag1_t ym;
    slong xyexp;

    if (NRB_IS_UNBOUNDED(x) || NRB_IS_UNBOUNDED(y))
        return nrb_zero_pm_inf(res, ctx);

    _nrb_unsafe_mid_get_mag1(xm, x, ctx);
    _nrb_unsafe_mid_get_mag1(ym, y, ctx);

    mag1_unsafe_mul(rad, NRB_RAD(x), NRB_RAD(y));
    mag1_unsafe_addmul(rad, xm, NRB_RAD(y));
    mag1_unsafe_addmul(rad, NRB_RAD(x), ym);

    xyexp = NRB_EXP(x) + NRB_EXP(y);

    if (nfloat_mul(NRB_MID(res), NRB_MID(x), NRB_MID(y), ctx) != GR_SUCCESS)
    {
        if (xyexp > 0)   /* overflow */
            return nrb_zero_pm_inf(res, ctx);
        else             /* underflow */
            mag1_unsafe_add_2exp_si(rad, rad, xyexp);
    }
    else if (!NFLOAT_IS_ZERO(NRB_MID(res)))
    {
        mag1_unsafe_add_2exp_si(rad, rad, NRB_EXP(res) - NFLOAT_CTX_PREC(ctx) + 1);
    }

    *NRB_RAD(res) = *rad;
    int status = nrb_fix_rad_underflow_overflow(res, ctx);
    NRB_ASSERT_VALID(res, ctx);

    return status;
}

/* assumes x is finite */
void _mag1_mul_2exp_si(mag1_t res, const mag1_t x, slong yexp)
{
    if (MAG1_IS_ZERO(x))
    {
        *res = *x;
    }
    else
    {
        if (yexp >= NFLOAT_MAX_EXP)
        {
            MAG1_INF(res);
        }
        else
        {
            slong exp;

            exp = MAG1_EXP(x) + FLINT_MAX(yexp, NFLOAT_MIN_EXP);

            if (yexp < NFLOAT_MIN_EXP)
            {
                MAG1_EXP(res) = NFLOAT_MIN_EXP;
                MAG1_MAN(res) = MAG1_ONE_HALF;
            }
            else
            {
                MAG1_EXP(res) = exp;
                MAG1_MAN(res) = MAG_MAN(x);
            }
        }
    }
}

/* todo */
int nrb_mul_2exp_si(nrb_ptr res, nrb_srcptr x, slong yexp, gr_ctx_t ctx)
{
    if (NRB_IS_UNBOUNDED(x))
        return nrb_set(res, x, ctx);

    if (nfloat_mul_2exp_si(NRB_MID(res), NRB_MID(x), yexp, ctx) != GR_SUCCESS)
        return GR_UNABLE;

    _mag1_mul_2exp_si(NRB_RAD(res), NRB_RAD(x), yexp);
    return GR_SUCCESS;
}


int nrb_sin(nrb_ptr res, nrb_srcptr x, gr_ctx_t ctx)
{
    arb_t t;
    arb_init(t);
    nrb_get_arb(t, x, ctx);
    arb_sin(t, t, NFLOAT_CTX_PREC(ctx));
    int status = nrb_set_arb(res, t, ctx);
    arb_clear(t);
    return status;
}

int nrb_cos(nrb_ptr res, nrb_srcptr x, gr_ctx_t ctx)
{
    arb_t t;
    arb_init(t);
    nrb_get_arb(t, x, ctx);
    arb_cos(t, t, NFLOAT_CTX_PREC(ctx));
    int status = nrb_set_arb(res, t, ctx);
    arb_clear(t);
    return status;
}
