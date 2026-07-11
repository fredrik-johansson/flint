/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FIXED_H
#define FIXED_H

#include "flint.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Efficient low-level fixed-point real arithmetic.

   A fixed-point number (x, n) is an unsigned n-limb fraction
   x[0], ..., x[n-1] representing sum x[i] 2^(FLINT_BITS (i - n)),
   i.e. 0 <= x < 1 with unit in the last place (ulp)
   2^(-FLINT_BITS n).  Outputs of size n + 1 additionally carry an
   integer (units) limb at index n.

   The series evaluation functions below require an argument reduced
   below 2^-32 (checked with a FLINT_ASSERT) and dispatch internally
   on the top limb of x: nonzero selects the 32-bit reduction range,
   zero the 64-bit-or-higher one.  Error bounds are in ulp; for exp
   and the non-alternating (hyperbolic) functions they are one-sided
   (the result never exceeds the true value), for the alternating
   functions two-sided. */

#define FIXED_EXP_RS_MAX_ERR(n) 10
/* for r = 0 (tuned default) the bound holds with the selected r,
   which never exceeds 512 */
#define FIXED_EXP_BITWISE_RS_MAX_ERR(n, r) \
    (9 * (slong) ((r) ? (r) : 512) + 100)
#define FIXED_LOG1P_BITWISE_RS_MAX_ERR(n, r) \
    (3 * (slong) ((r) ? (r) : 512) + 64)
#define FIXED_SIN_COS_BITWISE_RS_MAX_ERR(n, r) \
    (6 * (slong) ((r) ? (r) : 512) + 128)
#define FIXED_ATAN_BITWISE_RS_MAX_ERR(n, r) \
    (4 * (slong) ((r) ? (r) : 512) + 64)
#define FIXED_SIN_RS_MAX_ERR(n) 15
#define FIXED_COS_RS_MAX_ERR(n) 15
#define FIXED_SIN_COS_RS_MAX_ERR(n) 15
#define FIXED_SINH_RS_MAX_ERR(n) 15
#define FIXED_COSH_RS_MAX_ERR(n) 15
#define FIXED_SINH_COSH_RS_MAX_ERR(n) 15
#define FIXED_ATAN_RS_MAX_ERR(n) 15
#define FIXED_ATANH_RS_MAX_ERR(n) 15

/* exp((x, n)) -> (res, n + 1) */
void fixed_exp_rs(nn_ptr res, nn_srcptr x, slong n);

/* exp((x, n)) -> (res, n + 1) for any 0 <= x < 1 (n >= 2 on 32-bit
   limbs; r = 0 selects a tuned default), using bitwise
   argument reduction with a runtime-cached table of log(1 + 2^-i)
   followed by Taylor evaluation below 2^-r (r >= 32 a tuning
   parameter) and shift-and-add reconstruction.  All work happens at
   the output precision, so the error bound grows linearly with r;
   callers wanting sub-ulp accuracy should pad the precision by one
   limb themselves. */
void fixed_exp_bitwise_rs(nn_ptr res, nn_srcptr x, slong n, int r);

/* log(1 + (x, n)) -> (res, n) for any 0 <= x < 1, by the dual
   reduction: greedily multiply P by factors 1 + 2^-i (each one
   shift-and-add) while P (1 + 2^-i) <= 1 + x, then
   log((1+x)/P) = 2 atanh(((1+x) - P)/((1+x) + P)) with a single
   division, then add the tabulated logarithms.  The error bound
   grows linearly with r as for fixed_exp_bitwise_rs.  Requires
   r = 0 (which selects a tuned default) or r >= 16; values below 32 shorten the reduction further and are
   effective in the specialized code for n <= 4. */
void fixed_log1p_bitwise_rs(nn_ptr res, nn_srcptr x, slong n, int r);

/* sin and cos of (x, n) in [0, 1) -> (ysin, n + 1), (ycos, n + 1)
   (either may be NULL), by greedy reduction with the doubled angles
   2 atan(2^-i) and unimodular reconstruction; requires n >= 2 on
   32-bit limbs and r = 0 (tuned default) or r >= 16. */
void fixed_sin_cos_bitwise_rs(nn_ptr ysin, nn_ptr ycos, nn_srcptr x,
    slong n, int r);

/* atan((x, n)) -> (res, n) for x in [0, 1), by greedy vectoring;
   r = 0 selects a tuned default; otherwise r >= 16. */
void fixed_atan_bitwise_rs(nn_ptr res, nn_srcptr x, slong n, int r);

/* sin, cos, sinh, cosh of (x, n) -> (res, n + 1); the combined
   versions allow either output to be NULL */
void fixed_sin_rs(nn_ptr res, nn_srcptr x, slong n);
void fixed_cos_rs(nn_ptr res, nn_srcptr x, slong n);
void fixed_sin_cos_rs(nn_ptr ysin, nn_ptr ycos, nn_srcptr x, slong n);
void fixed_sinh_rs(nn_ptr res, nn_srcptr x, slong n);
void fixed_cosh_rs(nn_ptr res, nn_srcptr x, slong n);
void fixed_sinh_cosh_rs(nn_ptr ysinh, nn_ptr ycosh, nn_srcptr x, slong n);

/* atan, atanh of (x, n) -> (res, n) */
void fixed_atan_rs(nn_ptr res, nn_srcptr x, slong n);
void fixed_atanh_rs(nn_ptr res, nn_srcptr x, slong n);

/* Fallbacks used on all architectures and precisions: rectangular
   splitting at constant full precision with coefficients generated on
   the fly, requiring only x < 2^-32 (the number of terms is chosen
   from the actual leading zero bits of x).  Exposed for testing. */
void _fixed_exp_rs_fallback(nn_ptr res, nn_srcptr x, slong n);

/* Internal: exp for the wider range x < 2^-16 (n <= 5), used by the
   small reduction parameters of fixed_exp_bitwise_rs.  The hardcoded
   series family exists only for 64-bit limbs. */
#if FLINT_BITS == 64
void _fixed_exp_rs16(nn_ptr res, nn_srcptr x, slong n);

/* Internal: the fully specialized exp series, one per n <= 5, each
   built for the reduction parameter hardcoded alongside it in
   fixed_exp_bitwise_rs and using the smallest number of terms N with
   N r + log2(N!) >= 64 n. */
void _fixed_exp_rs_opt_1(nn_ptr res, nn_srcptr x);
void _fixed_exp_rs_opt_2(nn_ptr res, nn_srcptr x);
void _fixed_exp_rs_opt_3(nn_ptr res, nn_srcptr x);
void _fixed_exp_rs_opt_5(nn_ptr res, nn_srcptr x);
#endif

/* Internal: thread-local cached table of L_i = log(1 + 2^-i)
   truncated to _fixed_exp_logs_n fraction limbs, entries
   i = 0.._fixed_exp_logs_r, shared by fixed_exp_bitwise_rs and
   fixed_log1p_bitwise_rs. */
extern FLINT_TLS_PREFIX nn_ptr _fixed_exp_logs;
extern FLINT_TLS_PREFIX slong _fixed_exp_logs_n;
extern FLINT_TLS_PREFIX slong _fixed_exp_logs_r;
void _fixed_exp_logs_ensure(slong nc, slong rc);

/* Internal: number of slots the used array of _fixed_bitwise_reduce
   must provide: each index i = istart..r is used at most once, plus
   the window-boundary and final steps, which may repeat an index a
   bounded number of times to absorb the truncation creep of the
   table (see exp_bitwise_rs.c). */
#define FIXED_BITWISE_REDUCE_USED_ALLOC(r) \
    ((r) + 2 * ((r) / FLINT_BITS) + 12)

/* Internal: shared greedy table-subtraction reduction (see
   exp_bitwise_rs.c); returns the number of indices recorded in
   used, which must have room for
   FIXED_BITWISE_REDUCE_USED_ALLOC(r) entries. */
slong _fixed_bitwise_reduce(nn_ptr t, slong wn, int r, slong istart,
    nn_srcptr tab, slong tabn, slong * used);

/* Internal: thread-local cached table of the doubled angles
   alpha_i = 2 atan(2^-i) (entry 0 is unused and zeroed: pi/2 does
   not fit the fraction format), shared by fixed_sin_cos_bitwise_rs
   and fixed_atan_bitwise_rs. */
extern FLINT_TLS_PREFIX nn_ptr _fixed_atans;
extern FLINT_TLS_PREFIX slong _fixed_atans_n;
extern FLINT_TLS_PREFIX slong _fixed_atans_r;
void _fixed_atans_ensure(slong nc, slong rc);
void _fixed_sin_cos_rs_fallback(nn_ptr ysin, nn_ptr ycos, nn_srcptr x,
    slong n, int alternating);
void _fixed_atan_rs_fallback(nn_ptr res, nn_srcptr x, slong n,
    int alternating);

/* Internal: atanh for the wider range x < 2^-16, used by the small
   reduction parameters of fixed_log1p_bitwise_rs. */
void _fixed_atanh_rs16(nn_ptr res, nn_srcptr x, slong n);

/* Internal: atan and sin/cos for the wider range x < 2^-16, used by
   the small reduction parameters of the bitwise trigonometric
   functions.  The sin/cos routine requires both outputs. */
void _fixed_atan_rs16(nn_ptr res, nn_srcptr x, slong n);
void _fixed_sin_cos_rs16(nn_ptr ysin, nn_ptr ycos, nn_srcptr x, slong n);

#ifdef __cplusplus
}
#endif

#endif
