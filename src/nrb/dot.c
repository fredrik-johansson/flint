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
#include "nrb.h"


void
_nrb_vec_dot_addmul(ulong * s, slong sexp, double * serr, slong sulp_exp, nrb_srcptr xi, nrb_srcptr yi, int subtract, slong snlimbs)
{
    slong xexp, delta;
    slong delta_limbs, delta_bits;
    slong xn, yn, tn, tn2;
    int xsgnbit;
    ulong t[NFLOAT_MAX_LIMBS + 1];
    double xerr, yerr, terr;
    double xm, ym;
    nn_srcptr xd, yd;

    xn = NRB_N(xi);
    yn = NRB_N(yi);

    xexp = NRB_EXP(xi) + NRB_EXP(yi);
    xsgnbit = NRB_SGNBIT(xi) ^ NRB_SGNBIT(yi) ^ subtract;
    delta = sexp - xexp;

    xerr = NRB_ERR(xi);
    yerr = NRB_ERR(yi);
    terr = 0.0;

    xd = NRB_D(xi);
    yd = NRB_D(yi);

    xm = (double) xd[xn - 1];
    ym = (double) yd[yn - 1];

    /* Todo: when terr == 0, check for possible exact product */

    if (delta < FLINT_BITS)
    {
        tn = FLINT_MIN(xn, yn);

        /*
            Various possibilities:

            |ssss|ssss|ssss|ssss|ssss|
                |tttt|tttt|tttt|tttt|tttt|cccc|      cccc = correction limb in high product

            |ssss|ssss|ssss|ssss|ssss|
                |tttt|tttt|tttt|tttt|cccc|

            |ssss|ssss|ssss|ssss|ssss|
                |tttt|tttt|tttt|cccc|

            |ssss|ssss|ssss|ssss|ssss|
                |tttt|tttt|cccc|
        */

        if (xn != tn)
        {
            xerr = _nrb_trim_error(xerr, xd, xn, tn);
            xd += (xn - tn);
        }

        if (yn != tn)
        {
            yerr = _nrb_trim_error(yerr, yd, yn, tn);
            yd += (yn - tn);
        }

        /* Error in ulps of the product excluding correction limb */
        terr = (xm * yerr + ym * xerr) * ULP_N1 + d_mul_2exp_inrange(xerr * yerr, FLINT_MAX(-tn * FLINT_BITS, -NRB_MAX_ERROR_RIGHT_SHIFT2));

        t[0] = flint_mpn_mulhigh_n(t + 1, xd, yd, tn);
        tn2 = tn + 1;

        /* Error in ulps including correction limb */
        terr *= ULP_N1;
        terr += (2 * tn);


        /* Error rescaled to s. */
        terr = d_mul_2exp_inrange(terr, -delta + (tn - sn));
        *serr += terr;




        if (tn2 > snlimbs)
            flint_abort();

        mpn_rshift(t, t, tn2, delta);

        if (xsgnbit)
            mpn_sub_n(s + snlimbs - tn2, s + snlimbs - tn2, t, tn2);
        else
            mpn_add_n(s + snlimbs - tn2, s + snlimbs - tn2, t, tn2);


#if 0
        if (nlimbs == 3)
        {
            if (xsgnbit)
                sub_ddddmmmmssss(s[3], s[2], s[1], s[0], s[3], s[2], s[1], s[0],
                    t[3] >> delta,
                    (t[2] >> delta) | (t[3] << (FLINT_BITS - delta)),
                    (t[1] >> delta) | (t[2] << (FLINT_BITS - delta)),
                    (t[0] >> delta) | (t[1] << (FLINT_BITS - delta)));
            else
                add_ssssaaaaaaaa(s[3], s[2], s[1], s[0], s[3], s[2], s[1], s[0],
                    t[3] >> delta,
                    (t[2] >> delta) | (t[3] << (FLINT_BITS - delta)),
                    (t[1] >> delta) | (t[2] << (FLINT_BITS - delta)),
                    (t[0] >> delta) | (t[1] << (FLINT_BITS - delta)));
        }
        else if (nlimbs == 4)
        {
            if (xsgnbit)
                sub_dddddmmmmmsssss(s[4], s[3], s[2], s[1], s[0], s[4], s[3], s[2], s[1], s[0],
                    t[4] >> delta,
                    (t[3] >> delta) | (t[4] << (FLINT_BITS - delta)),
                    (t[2] >> delta) | (t[3] << (FLINT_BITS - delta)),
                    (t[1] >> delta) | (t[2] << (FLINT_BITS - delta)),
                    (t[0] >> delta) | (t[1] << (FLINT_BITS - delta)));
            else
                add_sssssaaaaaaaaaa(s[4], s[3], s[2], s[1], s[0], s[4], s[3], s[2], s[1], s[0],
                    t[4] >> delta,
                    (t[3] >> delta) | (t[4] << (FLINT_BITS - delta)),
                    (t[2] >> delta) | (t[3] << (FLINT_BITS - delta)),
                    (t[1] >> delta) | (t[2] << (FLINT_BITS - delta)),
                    (t[0] >> delta) | (t[1] << (FLINT_BITS - delta)));
        }
#endif

    }
    else
    {
        tn = FLINT_MIN(xn, yn);

        delta_limbs = delta / FLINT_BITS;
        delta_bits = delta % FLINT_BITS;

        /* alternative criterion: if (delta < (FLINT_BITS * nlimbs) + 2 * pad_bits) */
        if (delta_limbs < snlimbs)
        {
            /* todo: squaring case */
            flint_mpn_mulhigh_n(t, xd + delta_limbs - 1, yd + delta_limbs - 1, snlimbs - delta_limbs);

            if (delta_bits != 0)
                mpn_rshift(t, t, snlimbs - delta_limbs, delta_bits);

            if (xsgnbit)
                mpn_sub(s, s, snlimbs, t, snlimbs - delta_limbs);
            else
                mpn_add(s, s, snlimbs, t, snlimbs - delta_limbs);
        }
        else
        {
            /* Skip term. */
        }
    }
}


int
__nrb_vec_dot(nrb_ptr res, nrb_srcptr initial, int subtract, nrb_srcptr x, slong sizeof_xstep, nrb_srcptr y, slong sizeof_ystep, slong len, nrb_ctx_t ctx)
{
    slong i, xexp, sexp, sulp_exp, ulp_exp, sn, xn, yn, pad_bits;
    int xsgnbit;
    ulong s[NFLOAT_MAX_LIMBS + 1];
    ulong t[NFLOAT_MAX_LIMBS + 1];
    nrb_srcptr xi, yi;
    slong snlimbs;
    int status;

    sexp = WORD_MIN;
    sulp_exp = WORD_MIN;

    if (initial != NULL)
    {
        xn = NRB_N(initial);

        if (FLINT_UNLIKELY(xn == 0))
        {
            if (NRB_ERR(initial) == D_INF)
            {
                return nrb_zero_pm_inf(res, ctx);
            }
            else if (NRB_ERR(initial) != 0.0)
            {
                sexp = NFLOAT_EXP(initial);
                sulp_exp = xexp - xn * FLINT_BITS;
            }
        }
        else
        {
            sexp = NFLOAT_EXP(initial);
            sulp_exp = xexp - xn * FLINT_BITS;
        }
    }

    /* todo: flag whether we need to handle zeros subsequently */

    for (i = 0, xi = x, yi = y; i < len; i++, xi = (char *) xi + sizeof_xstep, yi = (char *) yi + sizeof_ystep)
    {
        xn = NRB_N(xi);
        yn = NRB_N(yi);

        if (FLINT_UNLIKELY(xn == 0 || yn == 0))
        {
            if (xn == 0 && NRB_ERR(xi) == 0.0)
                continue;
            if (yn == 0 && NRB_ERR(yi) == 0.0)
                continue;
            if (NRB_ERR(xi) == D_INF || NRB_ERR(yi) == D_INF)
                return nrb_zero_pm_inf(res, ctx);
        }

        xexp = NFLOAT_EXP(xi) + NFLOAT_EXP(yi);
        ulp_exp = xexp - FLINT_MIN(xn, yn) * FLINT_BITS;

        sexp = FLINT_MAX(sexp, xexp);
        sulp_exp = FLINT_MAX(sulp_exp, ulp_exp);
    }

    if (sexp == WORD_MIN)
        return nfloat_zero(res, ctx);

    pad_bits = FLINT_BIT_COUNT(len + 1) + 1;
    sexp += pad_bits;

    snlimbs = (sexp - sulp_exp + FLINT_BITS - 1) / FLINT_BITS;
    snlimbs = FLINT_MAX(snlimbs, 2);
    snlimbs = FLINT_MIN(snlimbs, NRB_CTX_NLIMBS(ctx) + 1);

    sulp_exp = sexp - snlimbs * FLINT_BITS;

    flint_mpn_zero(s, snlimbs);

    if (initial != NULL)
        return GR_UNABLE;
/*
    if (initial != NULL && !NFLOAT_IS_ZERO(initial))
        _nrb_vec_dot_set_initial(s, sexp, initial, subtract, nlimbs);
*/

    for (i = 0, xi = x, yi = y; i < len; i++, xi = (char *) xi + sizeof_xstep, yi = (char *) yi + sizeof_ystep)
    {
        _nrb_vec_dot_addmul(s, sexp, xi, yi, 0, nlimbs);
    }

    if (LIMB_MSB_IS_SET(s[snlimbs - 1]))
    {
        mpn_neg(s, s, snlimbs);
        xsgnbit = !subtract;
    }
    else
    {
        xsgnbit = subtract;
    }

    status = nrb_set_mpn_2exp(res, s, snlimbs, sexp, xsgnbit, ctx);
    status |= nrb_inplace_add_error_d_2exp_si(res, err, sulp_exp, ctx);
    return status;
}
