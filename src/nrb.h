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
#define MAG1_BITS 30
#define MAG1_ONE_HALF (UWORD(1) << (MAG1_BITS - 1))
#define MAG1_ONE_MINUS_EPS ((UWORD(1) << MAG1_BITS) - 1)

#if FLINT_BITS == 64
# define MAG1_FIXMUL(x, y) (((x) * (y)) >> MAG1_BITS)
#else
# define MAG1_FIXMUL(x, y) ((ulong) (((unsigned long long int) (x) * (unsigned long long int) (y)) >> MAG1_BITS))
#endif

typedef struct
{
    ulong man;
    slong exp;
}
mag1_struct;

typedef mag1_struct mag1_t[1];

#define NRB_DATA(x) ((nn_ptr) (x))
#define NRB_MID(x) (NRB_DATA(x) + 2)
#define NRB_RAD(x) ((mag1_struct *) NRB_DATA(x))

#define NRB_RAD_EXP(x) (NRB_RAD(x)->exp)
#define NRB_RAD_D(x) (NRB_RAD(x)->man)

#define NRB_EXP(x) NFLOAT_EXP(NRB_MID(x))
#define NRB_SGNBIT(x) NFLOAT_SGNBIT(NRB_MID(x))
#define NRB_D(x) NFLOAT_D(NRB_MID(x))

#define MAG1_MAN(x) ((x)->man)
#define MAG1_EXP(x) ((x)->exp)
#define MAG1_IS_SPECIAL(x) ((x)->man == UWORD(0))
#define MAG1_IS_INF(x) ((x)->exp == NFLOAT_EXP_POS_INF)
#define MAG1_IS_ZERO(x) ((x)->exp == NFLOAT_EXP_ZERO)
#define MAG1_ZERO(x) do { MAG1_MAN(x) = 0; MAG1_EXP(x) = NFLOAT_EXP_ZERO; } while (0)
#define MAG1_INF(x) do { MAG1_MAN(x) = 0; MAG1_EXP(x) = NFLOAT_EXP_POS_INF; } while (0)
#define _MAG1_EXP_INRANGE(exp) (NFLOAT_MIN_EXP <= (exp) && (exp) <= NFLOAT_MAX_EXP)
#define _MAG1_MAN_INRANGE(man) (MAG1_ONE_HALF <= (man) && (man) <= MAG1_ONE_MINUS_EPS)

#define NRB_IS_UNBOUNDED(x) MAG1_IS_INF(NRB_RAD(x))
#define NRB_IS_BOUNDED(x) (!MAG1_IS_INF(NRB_RAD(x)))
#define NRB_RAD_IS_ZERO(x) MAG1_IS_ZERO(NRB_RAD(x))
#define NRB_RAD_IS_SPECIAL(x) MAG1_IS_SPECAL(NRB_RAD(x))

#define NRB_ZERO_RAD(x) MAG1_ZERO(NRB_RAD(x))
#define NRB_ZERO_MID(x) do { NRB_EXP(x) = NFLOAT_EXP_ZERO; NRB_SGNBIT(x) = 0; } while (0)

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

int nrb_print_debug(nrb_srcptr x, gr_ctx_t ctx);

NRB_INLINE int
_mag1_is_valid(const mag1_t x)
{
    if (MAG1_EXP(x) == NFLOAT_EXP_ZERO && MAG1_MAN(x) == 0)
        return 1;
    if (MAG1_EXP(x) == NFLOAT_EXP_POS_INF && MAG1_MAN(x) == 0)
        return 1;
    if (_MAG1_EXP_INRANGE(MAG1_EXP(x)) && _MAG1_MAN_INRANGE(MAG1_MAN(x)))
        return 1;
    return 0;
}

#if FLINT_WANT_ASSERT
#define NRB_ASSERT_VALID(res, ctx) \
    do { \
        if () \
        { \
            int mid_ok = (NRB_EXP(res) == NFLOAT_EXP_ZERO || \
                     (NFLOAT_MIN_EXP <= NRB_EXP(res) && NRB_EXP(res) <= NFLOAT_MAX_EXP && \
                        (LIMB_MSB_IS_SET(NRB_D(res)[NFLOAT_CTX_NLIMBS(ctx) - 1])))); \
            int rad_ok = _mag1_is_valid(NRB_RAD(res)); \
            if (!mid_ok || !rad_ok) \
                nrb_print_debug(res, ctx); \
            FLINT_ASSERT(mid_ok); \
            FLINT_ASSERT(rad_ok); \
        } \
    } while (0)
#else
#define NRB_ASSERT_VALID(res, ctx)
#endif


NRB_INLINE void
nrb_clear(nrb_ptr res, gr_ctx_t ctx)
{
    NRB_ASSERT_VALID(res, ctx);
}

NRB_INLINE int
nrb_zero(nrb_ptr res, gr_ctx_t ctx)
{
    NRB_ZERO_RAD(res);
    NRB_ZERO_MID(res);
    NRB_ASSERT_VALID(res, ctx);
    return GR_SUCCESS;
}

NRB_INLINE int
nrb_one(nrb_ptr res, gr_ctx_t ctx)
{
    NRB_ZERO_RAD(res);
    nfloat_one(NRB_MID(res), ctx);
    NRB_ASSERT_VALID(res, ctx);
    return GR_SUCCESS;
}

NRB_INLINE int
nrb_neg_one(nrb_ptr res, gr_ctx_t ctx)
{
    NRB_ZERO_RAD(res);
    nfloat_neg_one(NRB_MID(res), ctx);
    NRB_ASSERT_VALID(res, ctx);
    return GR_SUCCESS;
}

NRB_INLINE int
nrb_set_ui(nrb_ptr res, ulong x, gr_ctx_t ctx)
{
    NRB_ZERO_RAD(res);
    nfloat_set_ui(NRB_MID(res), x, ctx);
    NRB_ASSERT_VALID(res, ctx);
    return GR_SUCCESS;
}

NRB_INLINE int
nrb_set_si(nrb_ptr res, ulong x, gr_ctx_t ctx)
{
    NRB_ZERO_RAD(res);
    nfloat_set_si(NRB_MID(res), x, ctx);
    NRB_ASSERT_VALID(res, ctx);
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

int nrb_mul_2exp_si(nrb_ptr res, nrb_srcptr x, slong yexp, gr_ctx_t ctx);

int nrb_sin(nrb_ptr res, nrb_srcptr x, gr_ctx_t ctx);
int nrb_cos(nrb_ptr res, nrb_srcptr x, gr_ctx_t ctx);

/* Magnitude arithmetic */

#define _MAG1_ADJUST_ONE_TOO_LARGE(xexp, xman) \
    do { \
        ulong __t = (xman) >> MAG1_BITS; \
        (xman) = ((xman) >> __t) + (__t & (xman)); \
        (xexp) += __t; \
    } while (0)

#define _MAG1_ADJUST_ONE_TOO_SMALL(xexp, xman) \
    do { \
        ulong __t = !(xman >> (MAG1_BITS - 1)); \
        (xman) = ((xman) << __t); \
        (xexp) -= __t; \
    } while (0)

NRB_INLINE void
mag1_add(mag1_t res, const mag1_t x, const mag1_t y)
{
    slong exp, xexp, yexp;
    ulong man, xman, yman;
    slong shift;

    xexp = x->exp;
    xman = x->man;
    yexp = y->exp;
    yman = y->man;

    if (MAG1_IS_SPECIAL(x) || MAG1_IS_SPECIAL(y))
    {
        if (xexp == NFLOAT_EXP_ZERO || yexp == NFLOAT_EXP_POS_INF)
            *res = *y;
        else
            *res = *x;
        return;
    }

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
        man = (man >> 1) + (man & 1);
        exp++;
    }
    else if (shift < MAG1_BITS)
    {
        man = xman + (yman >> shift) + UWORD(1);
    }
    else
    {
        man = xman + UWORD(1);
    }

    _MAG1_ADJUST_ONE_TOO_LARGE(exp, man);

    if (exp > NFLOAT_MAX_EXP)
    {
        exp = NFLOAT_EXP_POS_INF;
        man = 0;
    }

    MAG1_EXP(res) = exp;
    MAG1_MAN(res) = man;
}

NRB_INLINE void
mag1_add_2exp_si(mag1_t res, const mag1_t x, slong yexp)
{
    slong exp, xexp;
    ulong man, xman;
    slong shift;

    xexp = x->exp;
    xman = x->man;

    if (yexp >= NFLOAT_MAX_EXP || xexp == NFLOAT_EXP_POS_INF)
    {
        MAG1_INF(res);
        return;
    }

    yexp = FLINT_MAX(yexp, NFLOAT_MIN_EXP);

    if (xexp == NFLOAT_EXP_ZERO)
    {
        exp = yexp + 1;
        man = MAG1_ONE_HALF;
    }
    else
    {
        FLINT_ASSERT(xexp >= NFLOAT_MIN_EXP && xexp <= NFLOAT_MAX_EXP);

        shift = xexp - yexp;

        if (shift > 0)
        {
            exp = xexp;

            if (shift >= MAG1_BITS)
                man = xman + UWORD(1);
            else
                man = xman + (UWORD(1) << (MAG1_BITS - shift));
        }
        else
        {
            shift = -shift;

            exp = yexp + 1;

            if (shift >= MAG1_BITS)
                man = MAG1_ONE_HALF + 1;
            else
                man = MAG1_ONE_HALF + (xman >> (shift + 1)) + UWORD(1);
        }

        _MAG1_ADJUST_ONE_TOO_LARGE(exp, man);

        if (exp > NFLOAT_MAX_EXP)
        {
            exp = NFLOAT_EXP_POS_INF;
            man = 0;
        }
    }

    MAG1_EXP(res) = exp;
    MAG1_MAN(res) = man;
}

NRB_INLINE void
mag1_unsafe_add_2exp_si(mag1_t res, const mag1_t x, slong yexp)
{
    slong exp, xexp;
    ulong man, xman;
    slong shift;

    xexp = x->exp;
    xman = x->man;

    if (xexp == NFLOAT_EXP_ZERO)
    {
        exp = yexp + 1;
        man = MAG1_ONE_HALF;
    }
    else
    {
        shift = xexp - yexp;

        if (shift > 0)
        {
            exp = xexp;

            if (shift >= MAG1_BITS)
                man = xman + UWORD(1);
            else
                man = xman + (UWORD(1) << (MAG1_BITS - shift));
        }
        else
        {
            shift = -shift;

            exp = yexp + 1;

            if (shift >= MAG1_BITS)
                man = MAG1_ONE_HALF + 1;
            else
                man = MAG1_ONE_HALF + (xman >> (shift + 1)) + UWORD(1);
        }

        _MAG1_ADJUST_ONE_TOO_LARGE(exp, man);
    }

    MAG1_EXP(res) = exp;
    MAG1_MAN(res) = man;
}

/* assumes x and y both finite; does not check result for overflow or underflow */
NRB_INLINE void
mag1_unsafe_mul(mag1_t res, const mag1_t x, const mag1_t y)
{
    if (MAG1_MAN(x) == 0 || MAG1_MAN(y) == 0)
    {
        MAG1_ZERO(res);
    }
    else
    {
        MAG1_MAN(res) = MAG1_FIXMUL(MAG1_MAN(x), MAG1_MAN(y)) + UWORD(1);
        MAG1_EXP(res) = MAG1_EXP(x) + MAG1_EXP(y);
        _MAG1_ADJUST_ONE_TOO_SMALL(MAG1_EXP(res), MAG1_MAN(res));
    }
}

NRB_INLINE void
mag1_unsafe_nonzero_mul(mag1_t res, const mag1_t x, const mag1_t y)
{
    MAG1_MAN(res) = MAG1_FIXMUL(MAG1_MAN(x), MAG1_MAN(y)) + UWORD(1);
    MAG1_EXP(res) = MAG1_EXP(x) + MAG1_EXP(y);
    _MAG1_ADJUST_ONE_TOO_SMALL(MAG1_EXP(res), MAG1_MAN(res));
}

NRB_INLINE void
mag1_unsafe_mul_2exp_si(mag1_t res, const mag1_t x, slong yexp)
{
    MAG1_MAN(res) = MAG1_MAN(x);
    MAG1_EXP(res) = MAG1_EXP(x) + (MAG1_MAN(x) == 0) ? 0 : yexp;
}

NRB_INLINE void
mag1_unsafe_addmul(mag1_t res, const mag1_t x, const mag1_t y)
{
    if (MAG1_MAN(res) == 0)
    {
        mag1_unsafe_mul(res, x, y);
    }
    else if (MAG1_MAN(x) == 0 || MAG1_MAN(y) == 0)
    {
        return;
    }
    else
    {
        slong shift, exp;

        /* x*y < 2^e */
        exp = MAG1_EXP(x) + MAG1_EXP(y);
        shift = MAG1_EXP(res) - exp;

        if (shift >= 0)
        {
            if (shift >= MAG1_BITS)
                MAG1_MAN(res)++;
            else
                MAG1_MAN(res) = MAG1_MAN(res) + (MAG1_FIXMUL(MAG1_MAN(x), MAG1_MAN(y)) >> shift) + 1;
        }
        else
        {
            shift = -shift;
            MAG1_EXP(res) = exp;

            if (shift >= MAG1_BITS)
                MAG1_MAN(res) = MAG1_FIXMUL(MAG1_MAN(x), MAG1_MAN(y)) + 2;
            else
                MAG1_MAN(res) = MAG1_FIXMUL(MAG1_MAN(x), MAG1_MAN(y)) + (MAG1_MAN(res) >> shift) + 2;

            _MAG1_ADJUST_ONE_TOO_SMALL(MAG1_EXP(res), MAG1_MAN(res));
        }

        _MAG1_ADJUST_ONE_TOO_LARGE(MAG1_EXP(res), MAG1_MAN(res));
    }
}

#ifdef __cplusplus
}
#endif

#endif
