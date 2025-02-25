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

#include "double_extras.h"
#include "arb_types.h"
#include "nfloat.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    slong nlimbs;
    int flags;
}
_nrb_ctx_struct;

#define NRB_MIN_EXP NFLOAT_MIN_EXP
#define NRB_MAX_EXP NFLOAT_MAX_EXP

/* An error (measured in ulps) must either be 0.0, +inf, or a finite
   number in the following range. If the error becomes too large, we
   must trim the number of limbs. If the error becomes too small, we
   can increase the number of limbs (up to the maximum limit)
   or perturb the error.

   The range is chosen to allow perform a wide range of
   manipulations without risking either underflow or overflow.

   The max error is currently equivalent to several limbs. Setting
   this to a smaller value would result in more aggressive trimming.
   However, at lower precision, trimming is not actually desirable
   since this can create more overhead than just powering through
   the arithmetic.

   The main reason why we want to avoid underflow is not actually
   correctness (though that does require some attention); rather,
   it is very important to avoid creating denormal numbers since
   this destroys performance on many CPUs.
*/

#define NRB_MAX_ERR 0x1.0000000000000p+448
#define NRB_MIN_ERR 0x1.0000000000000p-128

#define NRB_MAX_ERR_EXP 449
#define NRB_MIN_ERR_EXP -127

#define NRB_MIN_ERR2 0x1.0000000000000p-96
#define NRB_MIN_ERR2_EXP -95


/* We can divide a normalized error by 2^NRB_MAX_ERROR_RIGHT_SHIFT
   without risking denormals. */
#define NRB_MAX_ERROR_RIGHT_SHIFT 768


typedef struct
{
    /* Scaled to ulp: absolute error is err * 2^(exp - n*FLINT_BITS) */
    double err;
    /* For normal values, satisfies NFLOAT_MIN_EXP <= exp <= NFLOAT_MIN_EXP. */
    slong exp;
    /* 0 or 1 */
    int sgnbit;
    /* Used limbs. At low precision, this will always be either 0
       (for a special value) or equal to N = NRB_CTX_NLIMBS(ctx). At higher
       precision, we allow a sliding precision with 0 <= n <= N
       in order to work around the limited exponent range of a double
       and to improve efficiency for wide balls. */
    int n;
}
nrb_head_struct;

#if FLINT_BITS == 64
#define NRB_HEADER_LIMBS 3
#else
#define NRB_HEADER_LIMBS 5
#endif

#define NRB_CTX(ctx) ((_nrb_ctx_struct *)(ctx))
#define NRB_CTX_FLAGS(ctx) (NRB_CTX(ctx)->flags)
#define NRB_CTX_NLIMBS(ctx) (NRB_CTX(ctx)->nlimbs)
#define NRB_CTX_PREC(ctx) ((NRB_CTX(ctx)->nlimbs) * FLINT_BITS)

#if FLINT_BITS == 64
#define NRB_CTX_DATA_NLIMBS(ctx) (NRB_HEADER_LIMBS + NRB_CTX_NLIMBS(ctx))
#else
/* Must pad to an even number of limbs on 32-bit so that the double is 64-bit aligned. */
#define NRB_CTX_DATA_NLIMBS(ctx) (NRB_HEADER_LIMBS + NRB_CTX_NLIMBS(ctx) + !(NRB_CTX_NLIMBS(ctx) & 1))
#endif

typedef gr_ctx_struct nrb_ctx_struct;
typedef nrb_ctx_struct nrb_ctx_t[1];

#define NRB_ERR(x) (((nrb_head_struct *) (x))->err)
#define NRB_EXP(x) (((nrb_head_struct *) (x))->exp)
#define NRB_SGNBIT(x) (((nrb_head_struct *) (x))->sgnbit)
#define NRB_N(x) (((nrb_head_struct *) (x))->n)
#define NRB_DATA(x) (((nn_ptr) (x)))
#define NRB_D(x) (NRB_DATA(x) + NRB_HEADER_LIMBS)

#define NRB_IS_UNBOUNDED(x) (NRB_ERR(x) == D_INF)
#define NRB_IS_EXACT(x) (NRB_ERR(x) == 0.0)
#define NRB_IS_SPECIAL(x) (NRB_N(x) == 0)
#define NRB_IS_ZERO(x) (NRB_IS_SPECIAL(x) && NRB_IS_EXACT(x))

typedef void * nrb_ptr;
typedef const void * nrb_srcptr;

typedef struct { ulong head[NRB_HEADER_LIMBS]; ulong d[64 / FLINT_BITS]; } nrb64_struct;
typedef struct { ulong head[NRB_HEADER_LIMBS]; ulong d[128 / FLINT_BITS]; } nrb128_struct;
typedef struct { ulong head[NRB_HEADER_LIMBS]; ulong d[192 / FLINT_BITS]; } nrb192_struct;
typedef struct { ulong head[NRB_HEADER_LIMBS]; ulong d[256 / FLINT_BITS]; } nrb256_struct;
typedef struct { ulong head[NRB_HEADER_LIMBS]; ulong d[384 / FLINT_BITS]; } nrb384_struct;
typedef struct { ulong head[NRB_HEADER_LIMBS]; ulong d[512 / FLINT_BITS]; } nrb512_struct;
typedef struct { ulong head[NRB_HEADER_LIMBS]; ulong d[1024 / FLINT_BITS]; } nrb1024_struct;
typedef struct { ulong head[NRB_HEADER_LIMBS]; ulong d[2048 / FLINT_BITS]; } nrb2048_struct;
typedef struct { ulong head[NRB_HEADER_LIMBS]; ulong d[4096 / FLINT_BITS]; } nrb4096_struct;

typedef nrb64_struct nrb64_t[1];
typedef nrb128_struct nrb128_t[1];
typedef nrb192_struct nrb192_t[1];
typedef nrb256_struct nrb256_t[1];
typedef nrb384_struct nrb384_t[1];
typedef nrb512_struct nrb512_t[1];
typedef nrb1024_struct nrb1024_t[1];
typedef nrb2048_struct nrb2048_t[1];
typedef nrb4096_struct nrb4096_t[1];



int nrb_ctx_init(gr_ctx_t ctx, slong prec, int flags);
int nrb_ctx_write(gr_stream_t out, nrb_ctx_t ctx);

NRB_INLINE
int nrb_ctx_set_real_prec(nrb_ctx_t ctx, slong prec)
{
    return GR_UNABLE;
}

NRB_INLINE
int nrb_ctx_get_real_prec(slong * res, nrb_ctx_t ctx)
{
    *res = NRB_CTX_PREC(ctx);
    return GR_SUCCESS;
}

int nrb_write(gr_stream_t out, nrb_srcptr x, nrb_ctx_t ctx);
int nrb_randtest(nrb_ptr res, flint_rand_t state, nrb_ctx_t ctx);
int nrb_randtest_ebits(nrb_ptr res, flint_rand_t state, slong ebits, nrb_ctx_t ctx);

truth_t nrb_is_zero(nrb_srcptr x, nrb_ctx_t ctx);
truth_t nrb_is_one(nrb_srcptr x, nrb_ctx_t ctx);
truth_t nrb_is_neg_one(nrb_srcptr x, nrb_ctx_t ctx);
truth_t nrb_equal(nrb_srcptr x, nrb_srcptr y, nrb_ctx_t ctx);

NRB_INLINE void nrb_init(nrb_ptr res, nrb_ctx_t ctx)
{
    NRB_ERR(res) = 0.0;
    NRB_SGNBIT(res) = 0;
    NRB_N(res) = 0;
}

int nrb_print_debug(nrb_srcptr x, nrb_ctx_t ctx);

#define NRB_DEBUG 1

#if NRB_DEBUG
#define NRB_EXPLAIN(s, f) flint_printf(s, f)
#else
#define NRB_EXPLAIN(s, f)
#endif

int _nrb_fix_range(nrb_ptr res, nrb_ctx_t ctx);

int _nrb_is_valid(nrb_srcptr x, nrb_ctx_t ctx);

#if NRB_DEBUG
#include <assert.h>

#define NRB_ASSERT_VALID(x, ctx) \
    do { \
        assert (_nrb_is_valid(x, ctx)); \
    } while (0)
#else
#define NRB_ASSERT_VALID(res, ctx)
#endif


NRB_INLINE void
nrb_clear(nrb_ptr res, nrb_ctx_t ctx)
{
    NRB_ASSERT_VALID(res, ctx);
}

NRB_INLINE int
nrb_zero(nrb_ptr res, nrb_ctx_t ctx)
{
    NRB_ERR(res) = 0.0;
    NRB_N(res) = 0;
    NRB_ASSERT_VALID(res, ctx);
    return GR_SUCCESS;
}

int nrb_one(nrb_ptr res, nrb_ctx_t ctx);
int nrb_neg_one(nrb_ptr res, nrb_ctx_t ctx);
int nrb_set_ui(nrb_ptr res, ulong x, nrb_ctx_t ctx);
int nrb_set_si(nrb_ptr res, slong x, nrb_ctx_t ctx);
int nrb_set_fmpz(nrb_ptr res, const fmpz_t x, nrb_ctx_t ctx);

int nrb_zero_pm_inf(nrb_ptr res, nrb_ctx_t ctx);
int nrb_zero_pm_2exp_si(nrb_ptr res, slong exp, nrb_ctx_t ctx);

void nrb_swap(nrb_ptr x, nrb_ptr y, nrb_ctx_t ctx);
int nrb_set(nrb_ptr res, nrb_srcptr x, nrb_ctx_t ctx);
void nrb_set_shallow(nrb_ptr res, nrb_srcptr x, nrb_ctx_t ctx);

/*
int nrb_get_rad_mag(mag_t res, nrb_srcptr x, nrb_ctx_t ctx);
int nrb_set_arf(nrb_ptr res, const arf_t x, nrb_ctx_t ctx);
*/
int nrb_get_arb(arb_t res, nrb_srcptr x, nrb_ctx_t ctx);
int nrb_set_arb(nrb_ptr res, const arb_t x, nrb_ctx_t ctx);

int nrb_get_interval_arf(arf_t a, arf_t b, nrb_srcptr x, nrb_ctx_t ctx, slong prec);
int nrb_get_mid_arf(arf_t res, nrb_srcptr x, nrb_ctx_t ctx);
int nrb_get_rad_arf(arf_t res, nrb_srcptr x, nrb_ctx_t ctx);

int nrb_inplace_add_error_d_2exp_si(nrb_ptr res, double err, slong err_exp, nrb_ctx_t ctx);

int nrb_abs(nrb_ptr res, nrb_srcptr x, nrb_ctx_t ctx);
int nrb_neg(nrb_ptr res, nrb_srcptr x, nrb_ctx_t ctx);
int nrb_add(nrb_ptr res, nrb_srcptr x, nrb_srcptr y, nrb_ctx_t ctx);
int nrb_sub(nrb_ptr res, nrb_srcptr x, nrb_srcptr y, nrb_ctx_t ctx);
int nrb_mul(nrb_ptr res, nrb_srcptr x, nrb_srcptr y, nrb_ctx_t ctx);

int nrb_mul_2exp_si(nrb_ptr res, nrb_srcptr x, slong yexp, nrb_ctx_t ctx);

int nrb_inv(nrb_ptr res, nrb_srcptr x, nrb_ctx_t ctx);
int nrb_div(nrb_ptr res, nrb_srcptr x, nrb_srcptr y, nrb_ctx_t ctx);
int nrb_sqrt(nrb_ptr res, nrb_srcptr x, nrb_ctx_t ctx);
int nrb_rsqrt(nrb_ptr res, nrb_srcptr x, nrb_ctx_t ctx);
int nrb_pow(nrb_ptr res, nrb_srcptr x, nrb_srcptr y, nrb_ctx_t ctx);
int nrb_exp(nrb_ptr res, nrb_srcptr x, nrb_ctx_t ctx);
int nrb_log(nrb_ptr res, nrb_srcptr x, nrb_ctx_t ctx);
int nrb_sin(nrb_ptr res, nrb_srcptr x, nrb_ctx_t ctx);
int nrb_cos(nrb_ptr res, nrb_srcptr x, nrb_ctx_t ctx);

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

/* Correction factor to ensure that error bound calculations in
   e.g. add and mul round up. On 64-bit we mostly need to compensate
   for possible downward rounding in the 53-bit double arithmetic.
   On 32-bit we also need to account for the fact that approximating
   |x| by its leading limb gives a lower bound with only ~32 bit
   accuracy. (Todo: we could take two limbs on 32-bit, and thereby
   improve accuracy.) */
#if FLINT_BITS == 64
#define NRB_CORRECTION_A  (1.0 + 0x1.0p-50)
#else
#define NRB_CORRECTION_A  (1.0 + 0x1.0p-28)
#endif

#define NRB_CORRECTION_B  0x1.0p-64

#ifdef __cplusplus
}
#endif

#endif
