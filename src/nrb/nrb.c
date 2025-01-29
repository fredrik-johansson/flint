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

        NRB_ASSERT_VALID(res);

        return GR_SUCCESS;
    }
}

int
nrb_zero_pm_inf(nrb_ptr res, gr_ctx_t ctx)
{
    NRB_RAD_EXP(res) = NFLOAT_EXP_POS_INF;
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

    NRB_ASSERT_VALID(res);

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

    NRB_ASSERT_VALID(res);

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

    NRB_ASSERT_VALID(res);

    return GR_SUCCESS;
}

int
nrb_rad_add_mag(nrb_ptr res, const mag_t x, gr_ctx_t ctx)
{
    if (mag_is_zero(x))
        return GR_SUCCESS;

    if (mag_is_inf(x) || fmpz_cmp_si(MAG_EXPREF(x), NFLOAT_MAX_EXP) > 0)
        return nrb_zero_pm_inf(res, ctx);  /* deliberately modifies midpoint too */

    slong exp;
    ulong man;

    if (fmpz_cmp_si(MAG_EXPREF(x), NFLOAT_MIN_EXP) < 0)
    {
        exp = NFLOAT_MIN_EXP + 1;
        man = MAG_ONE_HALF;
    }
    else
    {
        exp = MAG_EXP(x);
        man = MAG_MAN(x);
        _mag1_add(&NRB_RAD_EXP(res), &NRB_RAD_D(res),
                NRB_RAD_EXP(res), NRB_RAD_D(res),
                exp, man);
    }

    NRB_ASSERT_VALID(res);

    return GR_SUCCESS;
}

int
nrb_set_arb(nrb_ptr res, const arb_t x, gr_ctx_t ctx)
{
    int status;

    status = nrb_set_arf(res, arb_midref(x), ctx);
    status |= nrb_rad_add_mag(res, arb_radref(x), ctx);

    NRB_ASSERT_VALID(res);

    return status;
}

int
nrb_add(nrb_ptr res, nrb_srcptr x, nrb_srcptr y, gr_ctx_t ctx)
{
    slong rad_exp;
    ulong rad_man;

    _mag1_add(&rad_exp, &rad_man, NRB_RAD_EXP(x), NRB_RAD_D(x), NRB_RAD_EXP(y), NRB_RAD_D(y));

    if (rad_exp == NFLOAT_EXP_POS_INF)
        return nrb_zero_pm_inf(res, ctx);

    int subtract = (NRB_SGNBIT(x) != NRB_SGNBIT(y));

    /* add midpoints */
    if (nfloat_add(NRB_MID(res), NRB_MID(x), NRB_MID(y), ctx) != GR_SUCCESS)
    {
        if (!subtract)   /* overflow */
            return nrb_zero_pm_inf(res, ctx);
        else             /* underflow */
            _mag1_add_2exp(&rad_exp, &rad_man, rad_exp, rad_man, NFLOAT_MIN_EXP + 1);
    }
    else if (!NFLOAT_IS_ZERO(NRB_MID(res)))
    {
        /* add possible rounding error */
        slong dexp;

        /* xxx? */
        dexp = NRB_EXP(res) - FLINT_BITS * NFLOAT_CTX_NLIMBS(ctx) + 1;
        dexp = FLINT_MAX(dexp, NFLOAT_MIN_EXP);

        _mag1_add_2exp(&rad_exp, &rad_man, rad_exp, rad_man, dexp);
    }

    NRB_RAD_EXP(res) = rad_exp;
    NRB_RAD_D(res) = rad_man;

    NRB_ASSERT_VALID(res);

    return GR_SUCCESS;
}

int
nrb_sub(nrb_ptr res, nrb_srcptr x, nrb_srcptr y, gr_ctx_t ctx)
{
    slong rad_exp;
    ulong rad_man;

    _mag1_add(&rad_exp, &rad_man, NRB_RAD_EXP(x), NRB_RAD_D(x), NRB_RAD_EXP(y), NRB_RAD_D(y));

    if (rad_exp == NFLOAT_EXP_POS_INF)
        return nrb_zero_pm_inf(res, ctx);

    int subtract = (NRB_SGNBIT(x) == NRB_SGNBIT(y));

    /* add midpoints */
    if (nfloat_sub(NRB_MID(res), NRB_MID(x), NRB_MID(y), ctx) != GR_SUCCESS)
    {
        if (!subtract)   /* overflow */
            return nrb_zero_pm_inf(res, ctx);
        else             /* underflow */
            _mag1_add_2exp(&rad_exp, &rad_man, rad_exp, rad_man, NFLOAT_MIN_EXP + 1);
    }
    else if (!NFLOAT_IS_ZERO(NRB_MID(res)))
    {
        /* add possible rounding error */
        slong dexp;

        /* xxx? */
        dexp = NRB_EXP(res) - FLINT_BITS * NFLOAT_CTX_NLIMBS(ctx) + 1;
        dexp = FLINT_MAX(dexp, NFLOAT_MIN_EXP);

        _mag1_add_2exp(&rad_exp, &rad_man, rad_exp, rad_man, dexp);
    }

    NRB_RAD_EXP(res) = rad_exp;
    NRB_RAD_D(res) = rad_man;

    NRB_ASSERT_VALID(res);

    return GR_SUCCESS;
}

int
nrb_mul(nrb_ptr res, nrb_srcptr x, nrb_srcptr y, gr_ctx_t ctx)
{
    return GR_UNABLE;
}

