/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef NRB_H
#define NRB_H

#ifdef NRB_INLINES_C
#define NRB_INLINE
#else
#define NRB_INLINE static inline
#endif

#include "arb_types.h"
#include "nfloat.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef void * nrb_ptr;
typedef const void * nrb_srcptr;

/* 2 radius limbs + 2 nfloat header limbs */
#define NRB_HEADER_LIMBS 4
#define NRB_CTX_DATA_NLIMBS(ctx) (NFLOAT_CTX_NLIMBS(ctx) + NRB_HEADER_LIMBS)

/* This must be the same as MAG_BITS. */
#define NRB_RAD_BITS 30

#define NRB_DATA(x) ((nn_ptr) (x))
#define NRB_RAD(x) NRB_DATA(x)
#define NRB_MID(x) (NRB_DATA(x) + 2)

#define NRB_RAD_EXP(x) (((slong *) (x))[0])
#define NRB_RAD_D(x) (((nn_ptr) (x))[1])

#define NRB_EXP(x) NFLOAT_EXP(NRB_MID(x))
#define NRB_SGNBIT(x) NFLOAT_SGNBIT(NRB_MID(x))
#define NRB_D(x) NFLOAT_D(NRB_MID(x))

#define NRB_ZERO_RAD(res) do { NRB_RAD_EXP(res) = NFLOAT_EXP_ZERO; } while (0)
#define NRB_ZERO_MID(res) do { NRB_EXP(res) = NFLOAT_EXP_ZERO; NRB_SGNBIT(res) = 0; } while (0)

#define NRB_IS_UNBOUNDED(res) (NRB_RAD_EXP(x) == NFLOAT_EXP_POS_INF)
#define NRB_RAD_IS_ZERO(res) (NRB_RAD_EXP(x) == NFLOAT_EXP_ZERO)

int nrb_ctx_init(gr_ctx_t ctx, slong prec);
int nrb_ctx_write(gr_stream_t out, gr_ctx_t ctx);

int nrb_write(gr_stream_t out, nrb_srcptr x, gr_ctx_t ctx);
int nrb_randtest(nrb_ptr res, flint_rand_t state, gr_ctx_t ctx);

truth_t nrb_is_zero(nrb_srcptr x, gr_ctx_t ctx);
truth_t nrb_is_one(nrb_srcptr x, gr_ctx_t ctx);
truth_t nrb_is_neg_one(nrb_srcptr x, gr_ctx_t ctx);
truth_t nrb_equal(nrb_srcptr x, nrb_srcptr y, gr_ctx_t ctx);

NRB_INLINE void
nrb_init(nrb_ptr res, gr_ctx_t ctx)
{
    NRB_ZERO_RAD(res);
    NRB_ZERO_MID(res);
}

#define NRB_DEBUG 1

int nrb_write_debug(gr_stream_t out, nrb_srcptr x, gr_ctx_t ctx);
int nrb_print_debug(nrb_srcptr x, gr_ctx_t ctx);

#define NRB_ASSERT_VALID(res) \
    do { \
        if (NRB_DEBUG) \
        { \
            int mid_ok = (NRB_EXP(res) == NFLOAT_EXP_ZERO || \
                     (NFLOAT_MIN_EXP <= NRB_EXP(res) && NRB_EXP(res) <= NFLOAT_MAX_EXP && \
                        (LIMB_MSB_IS_SET(NRB_D(res)[NFLOAT_CTX_NLIMBS(ctx) - 1])))); \
            int rad_ok = (NRB_RAD_EXP(res) == NFLOAT_EXP_ZERO || \
                     NRB_RAD_EXP(res) == NFLOAT_EXP_POS_INF || \
                     (NFLOAT_MIN_EXP <= NRB_RAD_EXP(res) && NRB_RAD_EXP(res) <= NFLOAT_MAX_EXP && \
                      (UWORD(1) << (NRB_RAD_BITS - 1)) <= NRB_RAD_D(res) && NRB_RAD_D(res) <= (UWORD(1) << NRB_RAD_BITS) - 1)); \
            if (!mid_ok || !rad_ok) \
                nrb_print_debug(res, ctx); \
            FLINT_ASSERT(mid_ok); \
            FLINT_ASSERT(rad_ok); \
        } \
    } while (0)


NRB_INLINE void
nrb_clear(nrb_ptr res, gr_ctx_t ctx)
{
    NRB_ASSERT_VALID(res);
}

NRB_INLINE int
nrb_zero(nrb_ptr res, gr_ctx_t ctx)
{
    NRB_ZERO_RAD(res);
    NRB_ZERO_MID(res);
    NRB_ASSERT_VALID(res);
    return GR_SUCCESS;
}

NRB_INLINE int
nrb_one(nrb_ptr res, gr_ctx_t ctx)
{
    NRB_ZERO_RAD(res);
    nfloat_one(NRB_MID(res), ctx);
    NRB_ASSERT_VALID(res);
    return GR_SUCCESS;
}

NRB_INLINE int
nrb_neg_one(nrb_ptr res, gr_ctx_t ctx)
{
    NRB_ZERO_RAD(res);
    nfloat_neg_one(NRB_MID(res), ctx);
    NRB_ASSERT_VALID(res);
    return GR_SUCCESS;
}

NRB_INLINE int
nrb_set_ui(nrb_ptr res, ulong x, gr_ctx_t ctx)
{
    NRB_ZERO_RAD(res);
    nfloat_set_ui(NRB_MID(res), x, ctx);
    NRB_ASSERT_VALID(res);
    return GR_SUCCESS;
}

NRB_INLINE int
nrb_set_si(nrb_ptr res, ulong x, gr_ctx_t ctx)
{
    NRB_ZERO_RAD(res);
    nfloat_set_si(NRB_MID(res), x, ctx);
    NRB_ASSERT_VALID(res);
    return GR_SUCCESS;
}

int nrb_set_fmpz(nrb_ptr res, const fmpz_t x, gr_ctx_t ctx);

int nrb_zero_pm_inf(nrb_ptr res, gr_ctx_t ctx);
int nrb_zero_pm_2exp_si(nrb_ptr res, slong exp, gr_ctx_t ctx);
int nrb_zero_to_2exp_si_sgnbit(nrb_ptr res, slong exp, int sgnbit, gr_ctx_t ctx);

void nrb_swap(nrb_ptr x, nrb_ptr y, gr_ctx_t ctx);
int nrb_set(nrb_ptr res, nrb_srcptr x, gr_ctx_t ctx);
void nrb_set_shallow(nrb_ptr res, nrb_srcptr x, gr_ctx_t ctx);

int nrb_mid_get_arf(arf_t res, nrb_srcptr x, gr_ctx_t ctx);
int nrb_rad_get_mag(mag_t res, nrb_srcptr x, gr_ctx_t ctx);
int nrb_get_arb(arb_t res, nrb_srcptr x, gr_ctx_t ctx);
int nrb_rad_get_arf(arf_t res, nrb_srcptr x, gr_ctx_t ctx);
int nrb_set_arf(nrb_ptr res, const arf_t x, gr_ctx_t ctx);
int nrb_rad_add_mag(nrb_ptr res, const mag_t x, gr_ctx_t ctx);
int nrb_set_arb(nrb_ptr res, const arb_t x, gr_ctx_t ctx);

int nrb_abs(nrb_ptr res, nrb_srcptr x, gr_ctx_t ctx);
int nrb_neg(nrb_ptr res, nrb_srcptr x, gr_ctx_t ctx);
int nrb_add(nrb_ptr res, nrb_srcptr x, nrb_srcptr y, gr_ctx_t ctx);
int nrb_sub(nrb_ptr res, nrb_srcptr x, nrb_srcptr y, gr_ctx_t ctx);
int nrb_mul(nrb_ptr res, nrb_srcptr x, nrb_srcptr y, gr_ctx_t ctx);


/* Magnitude arithmetic */

#define MAG1_ONE_HALF (UWORD(1) << (NRB_RAD_BITS - 1))

#define _MAG1_ADJUST_ONE_TOO_LARGE(xexp, xman) \
    do { \
        ulong __t = (xman) >> NRB_RAD_BITS; \
        (xman) = ((xman) >> __t) + (__t & (xman)); \
        (xexp) += __t; \
    } while (0)

/* todo: separate for finite */
NRB_INLINE void
_mag1_add(slong * rexp, ulong * rman, slong xexp, ulong xman, slong yexp, ulong yman)
{
    if (xexp == NFLOAT_EXP_ZERO || yexp == NFLOAT_EXP_POS_INF)
    {
        *rexp = yexp;
        *rman = yman;
    }
    else if (yexp == NFLOAT_EXP_ZERO || xexp == NFLOAT_EXP_POS_INF)
    {
        *rexp = xexp;
        *rman = xman;
    }
    else
    {
        slong shift;
        slong exp;
        ulong man;

        FLINT_ASSERT(xexp >= NFLOAT_MIN_EXP && xexp <= NFLOAT_MAX_EXP);
        FLINT_ASSERT(yexp >= NFLOAT_MIN_EXP && yexp <= NFLOAT_MAX_EXP);

        shift = xexp - yexp;

        if (shift < 0)
        {
            FLINT_SWAP(slong, xexp, yexp);
            FLINT_SWAP(ulong, xman, yman);
        }

        exp = xexp;

        if (shift == 0)
        {
            man = xman + yman;
            _MAG1_ADJUST_ONE_TOO_LARGE(exp, man);
        }
        else if (shift < NRB_RAD_BITS)
        {
            man = xman + (yman >> shift) + UWORD(1);
        }
        else
        {
            man = xman + UWORD(1);
        }

        /* todo: combine the adjustments */
        _MAG1_ADJUST_ONE_TOO_LARGE(exp, man);

        if (exp > NFLOAT_MAX_EXP)
        {
            exp = NFLOAT_EXP_POS_INF;
            man = 0;
        }

        *rexp = exp;
        *rman = man;
    }
}

NRB_INLINE void
_mag1_add_2exp(slong * rexp, ulong * rman, slong xexp, ulong xman, slong yexp)
{
    if (xexp == NFLOAT_EXP_ZERO)
    {
        /* XXX: overflow (cannot happen as currently called, but ...) */
        *rexp = yexp + 1;
        *rman = MAG1_ONE_HALF;
    }
    else if (xexp == NFLOAT_EXP_POS_INF)
    {
        *rexp = xexp;
        *rman = xman;
    }
    else
    {
        slong shift;
        slong exp;
        ulong man;

        FLINT_ASSERT(xexp >= NFLOAT_MIN_EXP && xexp <= NFLOAT_MAX_EXP);
        FLINT_ASSERT(yexp >= NFLOAT_MIN_EXP && yexp <= NFLOAT_MAX_EXP);

        shift = xexp - yexp;

        if (shift > 0)
        {
            exp = xexp;

            if (shift >= NRB_RAD_BITS)
                man = xman + UWORD(1);
            else
                man = xman + (UWORD(1) << (NRB_RAD_BITS - shift));
        }
        else
        {
            shift = -shift;

            exp = yexp + 1;

            if (shift >= NRB_RAD_BITS)
                man = MAG1_ONE_HALF + 1;
            else
                man = MAG1_ONE_HALF + (xman >> (shift + 1)) + UWORD(1);
        }

        /* todo: combine the adjustments */
        _MAG1_ADJUST_ONE_TOO_LARGE(exp, man);

        if (exp > NFLOAT_MAX_EXP)
        {
            exp = NFLOAT_EXP_POS_INF;
            man = 0;
        }

        *rexp = exp;
        *rman = man;
    }
}


#ifdef __cplusplus
}
#endif

#endif
