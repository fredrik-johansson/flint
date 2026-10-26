/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "longlong.h"
#include "mpn_extras.h"
#include "gr.h"
#include "nfloat.h"
#include "impl.h"
#include "../mpn_extras/fp_small.h"

/*
    Division.

    With normalized mantissas A = a / B^n and Y = b / B^n in [1/2, 1), the
    quotient A / Y lies in (1/2, 2). Setting s = [a >= b], the integer
    floor(a B^n / (2^s b)) has exactly n limbs with the top bit set, and
    the exponent of the result is exp(x) - exp(y) + s. All kernels below
    produce this truncated mantissa together with an inexactness flag
    (the quotient lies strictly between the truncation and the
    truncation plus one ulp), which is all that directed rounding needs.

    * 1 and 2 limbs: the dividend is shifted right by s bits (branch-free)
      before a 2/1 resp. two 3/2 divisions, so that the quotient comes out
      normalized; the remainder gives the inexactness. Where hardware
      division is fast (FLINT_PREINVERT_LIMB_USE_NATIVE) it is used
      directly, otherwise a precomputed inverse.

    * Divisors with a single nonzero limb (e.g. integers): a chain of 2/1
      divisions, O(n).

    * Divisors of two limbs: the exact register-based division
      (flint_mpn_divhigh_exact), with remainder.

    * Longer: without directed rounding, the schoolbook approximate
      division _flint_mpn_divapprox_basecase_preinv1 at n limbs, which
      gives floor(q) or floor(q) + 1 of the exact (scaled) quotient q,
      i.e. an error below 1 ulp. With directed rounding, the shared
      flint_mpn_divhigh computes one guard limb; its truncation is
      certified (and inexact) when the guard limb is more than the error
      bound away from 0 and B - 1, and otherwise (probability about
      2^(5 - FLINT_BITS) for random input, but always for exact quotients)
      flint_mpn_divhigh_exact gives the exact truncation and remainder.
*/

/* Special values; x or y is special. Division by zero gives NaN. */
static int
_nfloat_div_special(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
{
    if (NFLOAT_IS_NAN(x) || NFLOAT_IS_NAN(y) || NFLOAT_IS_ZERO(y))
        return nfloat_nan(res, ctx);

    if (NFLOAT_IS_ZERO(x))
        return nfloat_zero(res, ctx);

    if (NFLOAT_IS_INF(x))
    {
        if (NFLOAT_IS_INF(y))
            return nfloat_nan(res, ctx);

        if (NFLOAT_IS_NEG_INF(x) ^ NFLOAT_SGNBIT(y))
            return nfloat_neg_inf(res, ctx);
        else
            return nfloat_pos_inf(res, ctx);
    }

    /* finite / inf */
    return nfloat_zero(res, ctx);
}

/* Writes exponent and sign of a nonzero result, rounding up the (already
   written) mantissa if needed. */
FLINT_FORCE_INLINE int
_nfloat_div_finish(nfloat_ptr res, slong exp, int sgnbit, int inexact, slong n, gr_ctx_t ctx)
{
    if (inexact && nfloat_should_round_up(sgnbit, ctx))
        NFLOAT_MANT_INCREMENT(NFLOAT_D(res), n, exp);

    NFLOAT_CHECK_EXP_RANGE(res, exp, sgnbit, ctx);
    NFLOAT_EXP(res) = exp;
    NFLOAT_SGNBIT(res) = sgnbit;
    return GR_SUCCESS;
}

/******************************************************************************/
/*  One limb                                                                  */
/******************************************************************************/

FLINT_FORCE_INLINE int
_nfloat_div_1_finish(nfloat_ptr res, ulong q, ulong r, slong exp, int sgnbit, gr_ctx_t ctx)
{
    if (r != 0 && nfloat_should_round_up(sgnbit, ctx))
    {
        q++;
        if (q == 0)
        {
            q = UWORD(1) << (FLINT_BITS - 1);
            exp++;
        }
    }

    NFLOAT_CHECK_EXP_RANGE(res, exp, sgnbit, ctx);
    NFLOAT_EXP(res) = exp;
    NFLOAT_SGNBIT(res) = sgnbit;
    NFLOAT_D(res)[0] = q;
    return GR_SUCCESS;
}

/******************************************************************************/
/*  Two limbs                                                                 */
/******************************************************************************/

#define _nfloat_div_2_kernel _flint_mpn_fdiv_2_kernel
#define _nfloat_div_2_kernel_preinv _flint_mpn_fdiv_2_kernel_preinv

FLINT_FORCE_INLINE int
_nfloat_div_2_finish(nfloat_ptr res, ulong q1, ulong q0, int inexact, slong exp, int sgnbit, gr_ctx_t ctx)
{
    if (inexact && nfloat_should_round_up(sgnbit, ctx))
    {
        add_ssaaaa(q1, q0, q1, q0, 0, 1);
        if (q1 == 0)
        {
            q1 = UWORD(1) << (FLINT_BITS - 1);
            exp++;
        }
    }

    NFLOAT_CHECK_EXP_RANGE(res, exp, sgnbit, ctx);
    NFLOAT_EXP(res) = exp;
    NFLOAT_SGNBIT(res) = sgnbit;
    NFLOAT_D(res)[0] = q0;
    NFLOAT_D(res)[1] = q1;
    return GR_SUCCESS;
}

/******************************************************************************/
/*  General length                                                            */
/******************************************************************************/

/* (q, n + 1) = floor(a B / d) for a normalized limb d; returns the
   remainder. GMP's assembly code is faster for long operands. */
static ulong
_nfloat_divrem_1(nn_ptr q, nn_srcptr a, slong n, ulong d, ulong dinv)
{
    if (n >= 8)
    {
        return mpn_divrem_1(q, 1, a, n, d);
    }
    else
    {
        ulong t[8];
        t[0] = 0;
        _nfloat_copy_limbs(t + 1, a, n);
        return flint_mpn_divrem_1_preinv(q, t, n + 1, d, dinv, 0);
    }
}

/* Normalizes the quotient (q, n + 1) (top limb 0 or 1) into the mantissa of
   res and finishes. q may be clobbered. */
FLINT_FORCE_INLINE int
_nfloat_div_mpn_finish(nfloat_ptr res, nn_ptr q, int inexact, slong exp, int sgnbit, slong n, gr_ctx_t ctx)
{
    if (q[n] != 0)
    {
        inexact |= (q[0] & 1);
        mpn_rshift(q, q, n + 1, 1);
        exp++;
    }

    _nfloat_copy_limbs(NFLOAT_D(res), q, n);
    return _nfloat_div_finish(res, exp, sgnbit, inexact, n, ctx);
}

/* Division of mantissas with a single-limb divisor (normalized limb d). */
static int
_nfloat_div_mpn_1(nfloat_ptr res, nn_srcptr a, ulong d, ulong dinv, slong exp, int sgnbit, slong n, gr_ctx_t ctx)
{
    ulong q[NFLOAT_MAX_LIMBS + 1];
    int inexact;

    inexact = (_nfloat_divrem_1(q, a, n, d, dinv) != 0);
    return _nfloat_div_mpn_finish(res, q, inexact, exp, sgnbit, n, ctx);
}

/* Sets (t, off + n) = floor(a B^off / 2^s) for s in {0, 1}, off >= 1,
   branch-free. */
FLINT_FORCE_INLINE void
_nfloat_shifted_numerator(nn_ptr t, nn_srcptr a, slong n, slong off, ulong s)
{
    ulong mask = -s;
    slong i;

    _nfloat_zero_limbs(t, off - 1);
    t[off - 1] = (a[0] << (FLINT_BITS - 1)) & mask;
    for (i = 0; i < n - 1; i++)
        t[off + i] = (a[i] >> s) | ((a[i + 1] << (FLINT_BITS - 1)) & mask);
    t[off + n - 1] = a[n - 1] >> s;
}

/* s = [a >= b], b of length bn <= n with b[bn - 1] != 0 compared as the
   top limbs of an n-limb number */
FLINT_FORCE_INLINE ulong
_nfloat_div_cmp(nn_srcptr a, nn_srcptr b, slong bn, slong n)
{
    ulong s = (a[n - 1] > b[bn - 1]);
    if (FLINT_UNLIKELY(a[n - 1] == b[bn - 1]))
        s = (mpn_cmp(a + n - bn, b, bn) >= 0);
    return s;
}

/* Divisor of bn >= 3 limbs: schoolbook approximate division with the 3/2
   inverse (for n <= NFLOAT_MAX_LIMBS this beats the register-based,
   divide-and-conquer and short divisions), giving the quotient or one
   more. */
static int
_nfloat_div_mpn_basecase(nfloat_ptr res, nn_srcptr a, nn_srcptr b, slong bn, slong exp, int sgnbit, slong n, gr_ctx_t ctx)
{
    ulong t[2 * NFLOAT_MAX_LIMBS + 1];
    ulong q[NFLOAT_MAX_LIMBS + 1];
    ulong s, dinv, FLINT_SET_BUT_UNUSED(qh);
    int inexact;

    FLINT_ASSERT(bn >= 3);

    s = _nfloat_div_cmp(a, b, bn, n);
    exp += s;

    if (!NFLOAT_CTX_HAS_DIRECTED_ROUNDING(ctx))
    {
        /* Error < 1 ulp; inexactness is irrelevant. */
        dinv = flint_mpn_preinv1(b[bn - 1], b[bn - 2]);
        _nfloat_shifted_numerator(t, a, n, bn, s);

        /* the output must not overlap the divisor */
        if (FLINT_UNLIKELY(b >= NFLOAT_D(res) && b < NFLOAT_D(res) + n))
        {
            qh = _flint_mpn_divapprox_basecase_preinv1(q, t, n + bn, b, bn, dinv);
            _nfloat_copy_limbs(NFLOAT_D(res), q, n);
        }
        else
        {
            qh = _flint_mpn_divapprox_basecase_preinv1(NFLOAT_D(res), t, n + bn, b, bn, dinv);
        }

        FLINT_ASSERT(qh == 0);
        return _nfloat_div_finish(res, exp, sgnbit, 0, n, ctx);
    }
    else
    {
        /* One guard limb: (q, n + 1), the truncation certified unless the
           guard limb is within the error bound of 0 or B - 1. */
        ulong se;

        flint_mpn_divhigh(q, a, n, b, bn, n);

        if (FLINT_LIKELY(q[0] > FLINT_MPN_DIVHIGH_GUARD_ERR_BOUND &&
                         q[0] < UWORD_MAX - FLINT_MPN_DIVHIGH_GUARD_ERR_BOUND))
        {
            _nfloat_copy_limbs(NFLOAT_D(res), q + 1, n);
            return _nfloat_div_finish(res, exp, sgnbit, 1, n, ctx);
        }

        inexact = flint_mpn_divhigh_exact(NFLOAT_D(res), &se, a, n, b, bn, n - 1);
        FLINT_ASSERT(se == s);
        return _nfloat_div_finish(res, exp, sgnbit, inexact, n, ctx);
    }
}

/* Divisor of two limbs (n >= 3): exact register-based division */
FLINT_STATIC_NOINLINE int
_nfloat_div_mpn_bn2(nfloat_ptr res, nn_srcptr a, nn_srcptr b, slong exp, int sgnbit, slong n, gr_ctx_t ctx)
{
    ulong s;
    int inexact;

    inexact = flint_mpn_divhigh_exact(NFLOAT_D(res), &s, a, n, b, 2, n - 1);
    return _nfloat_div_finish(res, exp + s, sgnbit, inexact, n, ctx);
}

/* Divisor with a single nonzero limb */
FLINT_STATIC_NOINLINE int
_nfloat_div_mpn_1_noinv(nfloat_ptr res, nn_srcptr a, ulong d, slong exp, int sgnbit, slong n, gr_ctx_t ctx)
{
    ulong dinv;
    NFLOAT_INVERT_LIMB(dinv, d);
    return _nfloat_div_mpn_1(res, a, d, dinv, exp, sgnbit, n, ctx);
}

/* Division of normalized n-limb mantissas, n >= 2. Also used by inv
   (a = 1/2). exp = exp(x) - exp(y). */
static int
_nfloat_div_mpn(nfloat_ptr res, nn_srcptr a, nn_srcptr b, slong exp, int sgnbit, slong n, gr_ctx_t ctx)
{
    slong bn = n;

    /* Strip zero limbs of the divisor. */
    if (b[0] == 0)
    {
        do {
            b++;
            bn--;
        } while (b[0] == 0);

        if (bn == 1)
            return _nfloat_div_mpn_1_noinv(res, a, b[0], exp, sgnbit, n, ctx);
    }

    if (FLINT_LIKELY(bn >= 3))
        return _nfloat_div_mpn_basecase(res, a, b, bn, exp, sgnbit, n, ctx);
    else
        return _nfloat_div_mpn_bn2(res, a, b, exp, sgnbit, n, ctx);
}

/******************************************************************************/
/*  Scalar functions                                                          */
/******************************************************************************/

FLINT_FORCE_INLINE int
_nfloat_div_1(nfloat_ptr res, ulong a, ulong b, slong exp, int sgnbit, gr_ctx_t ctx)
{
    ulong q, r, s;
    q = _nfloat_div_1_kernel(&r, &s, a, b);
    return _nfloat_div_1_finish(res, q, r, exp + s, sgnbit, ctx);
}

FLINT_FORCE_INLINE int
_nfloat_div_2(nfloat_ptr res, ulong a1, ulong a0, ulong b1, ulong b0, slong exp, int sgnbit, gr_ctx_t ctx)
{
    ulong q1, q0, s;
    int inexact;
    inexact = _nfloat_div_2_kernel(&q1, &q0, &s, a1, a0, b1, b0);
    return _nfloat_div_2_finish(res, q1, q0, inexact, exp + s, sgnbit, ctx);
}

int
nfloat_div(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
{
    slong n, exp;
    int sgnbit;

    if (NFLOAT_IS_SPECIAL(x) || NFLOAT_IS_SPECIAL(y))
        return _nfloat_div_special(res, x, y, ctx);

    n = NFLOAT_CTX_NLIMBS(ctx);
    sgnbit = NFLOAT_SGNBIT(x) ^ NFLOAT_SGNBIT(y);
    exp = NFLOAT_EXP(x) - NFLOAT_EXP(y);

    if (n == 1)
        return _nfloat_div_1(res, NFLOAT_D(x)[0], NFLOAT_D(y)[0], exp, sgnbit, ctx);
    else if (n == 2)
        return _nfloat_div_2(res, NFLOAT_D(x)[1], NFLOAT_D(x)[0],
                    NFLOAT_D(y)[1], NFLOAT_D(y)[0], exp, sgnbit, ctx);
    else
        return _nfloat_div_mpn(res, NFLOAT_D(x), NFLOAT_D(y), exp, sgnbit, n, ctx);
}

int
nfloat_inv(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    slong n, exp;
    int sgnbit;
    ulong half = UWORD(1) << (FLINT_BITS - 1);

    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_INF(x))
            return nfloat_zero(res, ctx);
        return nfloat_nan(res, ctx);
    }

    /* 1 = (1/2) 2^1 */
    n = NFLOAT_CTX_NLIMBS(ctx);
    sgnbit = NFLOAT_SGNBIT(x);
    exp = 1 - NFLOAT_EXP(x);

    if (n == 1)
    {
        return _nfloat_div_1(res, half, NFLOAT_D(x)[0], exp, sgnbit, ctx);
    }
    else if (n == 2)
    {
        return _nfloat_div_2(res, half, 0, NFLOAT_D(x)[1], NFLOAT_D(x)[0], exp, sgnbit, ctx);
    }
    else
    {
        ulong a[NFLOAT_MAX_LIMBS];
        _nfloat_zero_limbs(a, n - 1);
        a[n - 1] = half;
        return _nfloat_div_mpn(res, a, NFLOAT_D(x), exp, sgnbit, n, ctx);
    }
}

/* Division by c = d 2^(-norm) with d normalized, i.e. by the nfloat with
   mantissa d / B and exponent FLINT_BITS - norm. */
FLINT_FORCE_INLINE int
_nfloat_div_limb(nfloat_ptr res, nfloat_srcptr x, ulong d, ulong dinv, slong norm, int sgnbit, gr_ctx_t ctx)
{
    slong n = NFLOAT_CTX_NLIMBS(ctx);
    slong exp = NFLOAT_EXP(x) - (FLINT_BITS - norm);

    if (n == 1)
    {
        ulong q, r, s;
        q = _nfloat_div_1_kernel_preinv(&r, &s, NFLOAT_D(x)[0], d, dinv);
        return _nfloat_div_1_finish(res, q, r, exp + s, sgnbit, ctx);
    }

    return _nfloat_div_mpn_1(res, NFLOAT_D(x), d, dinv, exp, sgnbit, n, ctx);
}

static int
_nfloat_div_ui_special(nfloat_ptr res, nfloat_srcptr x, ulong c, int csgnbit, gr_ctx_t ctx)
{
    if (c == 0 || NFLOAT_IS_NAN(x))
        return nfloat_nan(res, ctx);
    if (NFLOAT_IS_ZERO(x))
        return nfloat_zero(res, ctx);
    if (NFLOAT_IS_NEG_INF(x) ^ csgnbit)
        return nfloat_neg_inf(res, ctx);
    else
        return nfloat_pos_inf(res, ctx);
}

int
nfloat_div_ui(nfloat_ptr res, nfloat_srcptr x, ulong c, gr_ctx_t ctx)
{
    ulong d, dinv;
    slong norm;

    if (NFLOAT_IS_SPECIAL(x) || c == 0)
        return _nfloat_div_ui_special(res, x, c, 0, ctx);

    norm = flint_clz(c);
    d = c << norm;
    NFLOAT_INVERT_LIMB(dinv, d);
    return _nfloat_div_limb(res, x, d, dinv, norm, NFLOAT_SGNBIT(x), ctx);
}

int
nfloat_div_si(nfloat_ptr res, nfloat_srcptr x, slong c, gr_ctx_t ctx)
{
    ulong d, dinv, uc;
    slong norm;
    int csgnbit = (c < 0);

    uc = FLINT_UABS(c);

    if (NFLOAT_IS_SPECIAL(x) || c == 0)
        return _nfloat_div_ui_special(res, x, uc, csgnbit, ctx);

    norm = flint_clz(uc);
    d = uc << norm;
    NFLOAT_INVERT_LIMB(dinv, d);
    return _nfloat_div_limb(res, x, d, dinv, norm, NFLOAT_SGNBIT(x) ^ csgnbit, ctx);
}

/******************************************************************************/
/*  Vector functions                                                          */
/******************************************************************************/

int
_nfloat_vec_div(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, slong len, gr_ctx_t ctx)
{
    slong i, n = NFLOAT_CTX_NLIMBS(ctx);
    slong sz = NFLOAT_CTX_DATA_NLIMBS(ctx);
    int status = GR_SUCCESS;
    nn_ptr R = (nn_ptr) res;
    nn_srcptr X = (nn_srcptr) x;
    nn_srcptr Y = (nn_srcptr) y;

    if (n == 1)
    {
        for (i = 0; i < len; i++)
        {
            nn_ptr r = R + i * sz;
            nn_srcptr u = X + i * sz;
            nn_srcptr v = Y + i * sz;

            if (FLINT_UNLIKELY(NFLOAT_IS_SPECIAL(u) || NFLOAT_IS_SPECIAL(v)))
                status |= _nfloat_div_special(r, u, v, ctx);
            else
                status |= _nfloat_div_1(r, NFLOAT_D(u)[0], NFLOAT_D(v)[0],
                    NFLOAT_EXP(u) - NFLOAT_EXP(v), NFLOAT_SGNBIT(u) ^ NFLOAT_SGNBIT(v), ctx);
        }
    }
    else if (n == 2)
    {
        for (i = 0; i < len; i++)
        {
            nn_ptr r = R + i * sz;
            nn_srcptr u = X + i * sz;
            nn_srcptr v = Y + i * sz;

            if (FLINT_UNLIKELY(NFLOAT_IS_SPECIAL(u) || NFLOAT_IS_SPECIAL(v)))
                status |= _nfloat_div_special(r, u, v, ctx);
            else
                status |= _nfloat_div_2(r, NFLOAT_D(u)[1], NFLOAT_D(u)[0],
                    NFLOAT_D(v)[1], NFLOAT_D(v)[0],
                    NFLOAT_EXP(u) - NFLOAT_EXP(v), NFLOAT_SGNBIT(u) ^ NFLOAT_SGNBIT(v), ctx);
        }
    }
    else
    {
        for (i = 0; i < len; i++)
        {
            nn_ptr r = R + i * sz;
            nn_srcptr u = X + i * sz;
            nn_srcptr v = Y + i * sz;

            if (FLINT_UNLIKELY(NFLOAT_IS_SPECIAL(u) || NFLOAT_IS_SPECIAL(v)))
                status |= _nfloat_div_special(r, u, v, ctx);
            else
                status |= _nfloat_div_mpn(r, NFLOAT_D(u), NFLOAT_D(v),
                    NFLOAT_EXP(u) - NFLOAT_EXP(v), NFLOAT_SGNBIT(u) ^ NFLOAT_SGNBIT(v), n, ctx);
        }
    }

    return status;
}

/* Division by a single-limb mantissa d (normalized) with exponent
   shift FLINT_BITS - norm and sign csgnbit. */
static int
_nfloat_vec_div_limb(nfloat_ptr res, nfloat_srcptr x, slong len, ulong d, slong norm, int csgnbit, gr_ctx_t ctx)
{
    slong i, sz = NFLOAT_CTX_DATA_NLIMBS(ctx);
    int status = GR_SUCCESS;
    ulong dinv;

    NFLOAT_INVERT_LIMB(dinv, d);

    for (i = 0; i < len; i++)
    {
        nn_ptr r = (nn_ptr) res + i * sz;
        nn_srcptr u = (nn_srcptr) x + i * sz;

        if (FLINT_UNLIKELY(NFLOAT_IS_SPECIAL(u)))
            status |= _nfloat_div_ui_special(r, u, 1, csgnbit, ctx);
        else
            status |= _nfloat_div_limb(r, u, d, dinv, norm, NFLOAT_SGNBIT(u) ^ csgnbit, ctx);
    }

    return status;
}

/*
    Division of (x, n) by a fixed (b, n) with n >= 3 and at least two
    nonzero limbs in b, via the approximate reciprocal R = floor(B^(2n+1) / b)
    or one more, which lies in (B^(n+1), 2 B^(n+1)), so that R = B^(n+1) + R1
    with R1 < B^(n+1). For a mantissa X, the quotient with n + 1 fraction
    limbs is q = X R_exact / B^n; we compute Q = X B + H, H the high n + 1
    limbs of (X B) R1 (a lower bound with error < 2 units). Then
    Q <= X R / B^n, Q >= X R / B^n - 2 and |X R - X R_exact| / B^n < 1, so
    q lies in (Q - 1, Q + 3). The truncation of Q to n limbs is certified
    (and inexact) when the guard bits D satisfy 2 <= D <= 2^k - 4; the
    error is otherwise below 1 ulp anyway, so only directed rounding
    needs the exact fallback.
*/
static int
_nfloat_vec_div_mpn_recip(nfloat_ptr res, nfloat_srcptr x, slong len, nfloat_srcptr c, gr_ctx_t ctx)
{
    slong i, n = NFLOAT_CTX_NLIMBS(ctx);
    slong sz = NFLOAT_CTX_DATA_NLIMBS(ctx);
    ulong one[NFLOAT_MAX_LIMBS];
    ulong R[NFLOAT_MAX_LIMBS + 3];
    ulong Xp[NFLOAT_MAX_LIMBS + 1];
    ulong Q[NFLOAT_MAX_LIMBS + 2];
    ulong H[NFLOAT_MAX_LIMBS + 1];
    nn_srcptr b = NFLOAT_D(c);
    slong cexp = NFLOAT_EXP(c);
    int csgnbit = NFLOAT_SGNBIT(c);
    int directed = NFLOAT_CTX_HAS_DIRECTED_ROUNDING(ctx);
    int status = GR_SUCCESS;

    /* (R, n + 3) = floor(B^(n-1) B^(n+2) / b) or one more */
    _nfloat_zero_limbs(one, n - 1);
    one[n - 1] = 1;
    flint_mpn_divapprox_fraction(R, one, n, b, n, n + 2);
    FLINT_ASSERT(R[n + 2] == 0 && R[n + 1] == 1);

    Xp[0] = 0;

    for (i = 0; i < len; i++)
    {
        nn_ptr r = (nn_ptr) res + i * sz;
        nn_srcptr u = (nn_srcptr) x + i * sz;
        slong exp;
        int sgnbit, certified;

        if (FLINT_UNLIKELY(NFLOAT_IS_SPECIAL(u)))
        {
            status |= _nfloat_div_special(r, u, c, ctx);
            continue;
        }

        sgnbit = NFLOAT_SGNBIT(u) ^ csgnbit;
        exp = NFLOAT_EXP(u) - cexp;

        _nfloat_copy_limbs(Xp + 1, NFLOAT_D(u), n);
        flint_mpn_mulhigh_n(H, Xp, R, n + 1);
        Q[0] = H[0];
        Q[n + 1] = mpn_add_n(Q + 1, Xp + 1, H + 1, n);

        if (!directed)
        {
            status |= _nfloat_div_mpn_finish(r, Q + 1, 0, exp, sgnbit, n, ctx);
            continue;
        }

        if (Q[n + 1] != 0)
            certified = ((Q[1] & 1) || Q[0] >= 2) && (!(Q[1] & 1) || Q[0] <= UWORD_MAX - 3);
        else
            certified = (Q[0] >= 2) && (Q[0] <= UWORD_MAX - 3);

        if (FLINT_LIKELY(certified))
            status |= _nfloat_div_mpn_finish(r, Q + 1, 1, exp, sgnbit, n, ctx);
        else
            status |= _nfloat_div_mpn(r, NFLOAT_D(u), b, exp, sgnbit, n, ctx);
    }

    return status;
}

int
_nfloat_vec_div_scalar(nfloat_ptr res, nfloat_srcptr x, slong len, nfloat_srcptr c, gr_ctx_t ctx)
{
    slong i, n = NFLOAT_CTX_NLIMBS(ctx);
    slong sz = NFLOAT_CTX_DATA_NLIMBS(ctx);
    int status = GR_SUCCESS;
    nn_srcptr b;
    slong cexp;
    int csgnbit;

    if (len <= 1 || NFLOAT_IS_SPECIAL(c))
    {
        for (i = 0; i < len; i++)
            status |= nfloat_div((nn_ptr) res + i * sz, (nn_srcptr) x + i * sz, c, ctx);
        return status;
    }

    /* the scalar may alias an entry of res */
    {
        ulong t[NFLOAT_MAX_ALLOC];
        NFLOAT_EXP(t) = NFLOAT_EXP(c);
        NFLOAT_SGNBIT(t) = NFLOAT_SGNBIT(c);
        _nfloat_copy_limbs(NFLOAT_D(t), NFLOAT_D(c), n);
        c = t;

        b = NFLOAT_D(c);
        cexp = NFLOAT_EXP(c);
        csgnbit = NFLOAT_SGNBIT(c);

        if (n == 1)
        {
            return _nfloat_vec_div_limb(res, x, len, b[0], FLINT_BITS - cexp, csgnbit, ctx);
        }
        else if (n == 2)
        {
            ulong b1 = b[1], b0 = b[0], dinv;

            if (b0 == 0)
                return _nfloat_vec_div_limb(res, x, len, b1, FLINT_BITS - cexp, csgnbit, ctx);

            dinv = flint_mpn_preinv1(b1, b0);

            for (i = 0; i < len; i++)
            {
                nn_ptr r = (nn_ptr) res + i * sz;
                nn_srcptr u = (nn_srcptr) x + i * sz;
                ulong q1, q0, s;
                int inexact;

                if (FLINT_UNLIKELY(NFLOAT_IS_SPECIAL(u)))
                {
                    status |= _nfloat_div_special(r, u, c, ctx);
                    continue;
                }

                inexact = _nfloat_div_2_kernel_preinv(&q1, &q0, &s,
                    NFLOAT_D(u)[1], NFLOAT_D(u)[0], b1, b0, dinv);
                status |= _nfloat_div_2_finish(r, q1, q0, inexact,
                    NFLOAT_EXP(u) - cexp + s, NFLOAT_SGNBIT(u) ^ csgnbit, ctx);
            }

            return status;
        }
        else
        {
            if (flint_mpn_zero_p(b, n - 1))
                return _nfloat_vec_div_limb(res, x, len, b[n - 1], FLINT_BITS - cexp, csgnbit, ctx);

            return _nfloat_vec_div_mpn_recip(res, x, len, c, ctx);
        }
    }
}

int
_nfloat_vec_div_scalar_ui(nfloat_ptr res, nfloat_srcptr x, slong len, ulong c, gr_ctx_t ctx)
{
    slong i, sz = NFLOAT_CTX_DATA_NLIMBS(ctx);
    int status = GR_SUCCESS;

    if (c == 0)
    {
        for (i = 0; i < len; i++)
            status |= nfloat_div_ui((nn_ptr) res + i * sz, (nn_srcptr) x + i * sz, c, ctx);
        return status;
    }

    return _nfloat_vec_div_limb(res, x, len, c << flint_clz(c), flint_clz(c), 0, ctx);
}

int
_nfloat_vec_div_scalar_si(nfloat_ptr res, nfloat_srcptr x, slong len, slong c, gr_ctx_t ctx)
{
    slong i, sz = NFLOAT_CTX_DATA_NLIMBS(ctx);
    int status = GR_SUCCESS;
    ulong uc;

    if (c == 0)
    {
        for (i = 0; i < len; i++)
            status |= nfloat_div_si((nn_ptr) res + i * sz, (nn_srcptr) x + i * sz, c, ctx);
        return status;
    }

    uc = FLINT_UABS(c);
    return _nfloat_vec_div_limb(res, x, len, uc << flint_clz(uc), flint_clz(uc), c < 0, ctx);
}
