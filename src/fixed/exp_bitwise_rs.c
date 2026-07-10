/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "mpn_extras.h"
#include "arb.h"
#include "fixed.h"

/* exp via bitwise argument reduction: subtract in turn each
   L_i = log(1 + 2^-i), i = 0, 1, ..., r, for which L_i <= x, evaluate
   the Taylor series on the reduced argument below 2^-r, and multiply
   back the used (1 + 2^-i) factors, each of which is a single
   shift-and-add.

   The classical single-pass reduction is sound on [0, 1): since
   (1 + 2^-i)^2 > 1 + 2^(1-i) we have L_{i-1} < 2 L_i, so by induction
   the remainder is below L_i after step i (base case 1 < 2 log 2),
   and each L_i is used at most once.  The stored logarithms are
   truncations, which lets the remainder creep at most one guard ulp
   per step above the exact invariant; a single extra conditional
   subtraction of L_r after the loop therefore guarantees a remainder
   strictly below L_r < 2^-r.

   All work happens at the output precision of n limbs -- callers
   needing sub-ulp accuracy are expected to pad the precision
   themselves, so padding internally would double up.  The error
   therefore grows linearly with r: at most r + 2 table entries are
   subtracted, each a truncation short of the true logarithm by less
   than two ulp (multiplying the result by exp of the deficit), and
   each of the at most r + 2 reconstruction shifts truncates below
   one ulp, all amplified by the remaining product < e; with the
   series error this gives less than 8.3 (r + 2) + 80 ulp, bounded by
   FIXED_EXP_BITWISE_RS_MAX_ERR(n, r) = 9 r + 100.  (With one guard
   limb these terms all drop below one output ulp -- measured 1.00
   ulp flat for r up to 320 -- so a caller-side padded limb absorbs
   the bound with 50+ bits to spare.)  r is clamped to
   FLINT_BITS n - 16: this keeps the reduction slop (r + 2 ulp)
   below L_r, so the single extra subtraction suffices, and past
   that point the series would degenerate anyway.

   The log table is generated on demand with arb (to be replaced by a
   native computation later) and cached per thread at the largest
   precision and index range requested so far; smaller requests read
   the top limbs of the cached entries, which equal the truncations
   at the smaller precision exactly. */

FLINT_TLS_PREFIX nn_ptr _fixed_exp_logs = NULL;
FLINT_TLS_PREFIX slong _fixed_exp_logs_n = 0;   /* limbs per entry */
FLINT_TLS_PREFIX slong _fixed_exp_logs_r = 0;   /* entries 0..r */
static FLINT_TLS_PREFIX int _fixed_exp_logs_cleanup_registered = 0;

static void
_fixed_exp_logs_cleanup(void)
{
    flint_free(_fixed_exp_logs);
    _fixed_exp_logs = NULL;
    _fixed_exp_logs_n = 0;
    _fixed_exp_logs_r = 0;
    _fixed_exp_logs_cleanup_registered = 0;
}

#define LOGS_L(i) (_fixed_exp_logs + (i) * _fixed_exp_logs_n)

void
_fixed_exp_logs_ensure(slong nc, slong rc)
{
    arb_t v, t;
    arf_t lb;
    fmpz_t f;
    slong i, prec;

    if (nc <= _fixed_exp_logs_n && rc <= _fixed_exp_logs_r)
        return;

    nc = FLINT_MAX(nc, _fixed_exp_logs_n);
    rc = FLINT_MAX(rc, _fixed_exp_logs_r);

    flint_free(_fixed_exp_logs);
    _fixed_exp_logs = flint_malloc((rc + 1) * nc * sizeof(ulong));
    _fixed_exp_logs_n = nc;
    _fixed_exp_logs_r = rc;

    prec = FLINT_BITS * nc + 64;

    arb_init(v);
    arb_init(t);
    arf_init(lb);
    fmpz_init(f);

    for (i = 0; i <= rc; i++)
    {
        /* L_i = log(1 + 2^-i), truncated (floor) to nc fraction
           limbs */
        arb_one(v);
        arb_mul_2exp_si(v, v, -i);
        arb_log1p(t, v, prec);

        arb_get_lbound_arf(lb, t, prec);
        arf_mul_2exp_si(lb, lb, FLINT_BITS * nc);
        arf_get_fmpz(f, lb, ARF_RND_FLOOR);
        FLINT_ASSERT(fmpz_sgn(f) >= 0);
        fmpz_get_ui_array(LOGS_L(i), nc, f);
    }

    arb_clear(v);
    arb_clear(t);
    arf_clear(lb);
    fmpz_clear(f);

    if (!_fixed_exp_logs_cleanup_registered)
    {
        flint_register_cleanup_function(_fixed_exp_logs_cleanup);
        _fixed_exp_logs_cleanup_registered = 1;
    }
}


/* Use the sinh series (half the terms) plus a squaring and a square
   root once the direct exp series gets long enough.  Measured
   crossovers on x86-64: for r < 64 (series through the pre32 path)
   sinh wins from about 45 terms (n ~ 24 at r = 32); for r >= 64 the
   windowed pre64 series is more efficient per term and mpn_sqrtrem
   sets a higher floor, moving the crossover to about 128 terms
   (n ~ 2r).  The margins near the crossovers are within a few
   percent, so the exact placement is not critical. */
#define EXP_USE_SINH(wn, r) \
    (FLINT_BITS * (wn) >= (((r) >= 64) ? 128 : 45) * (slong) (r))

/* exp(t) of the reduced argument t < 2^-r into y (wn + 1 limbs:
   wn fraction limbs and a units limb); the series functions pick the
   32- or 64-bit internal range from the top limb of t.  When the
   direct series would need many terms, evaluate sinh(t) instead --
   the odd series has half the terms -- and reconstruct
   exp(t) = sinh(t) + sqrt(1 + sinh(t)^2), costing one squaring and
   one square root.  The sinh error (<= 15 ulp at the guarded
   precision), the truncated squaring and the floor of the integer
   square root all sit below one guard ulp of the final result. */
static void
_fixed_exp_reduced(nn_ptr y, nn_srcptr t, slong wn, int r, int use_sinh)
{
    if (!use_sinh)
    {
        fixed_exp_rs(y, t, wn);
    }
    else
    {
        nn_ptr s, u2, rt, rem;
        TMP_INIT;

        TMP_START;
        s = TMP_ALLOC(((wn + 1) + (2 * wn + 1) + (wn + 1) + (wn + 2))
            * sizeof(ulong));
        u2 = s + (wn + 1);
        rt = u2 + (2 * wn + 1);
        rem = rt + (wn + 1);

        fixed_sinh_rs(s, t, wn);

        /* u2 = (1 + sinh(t)^2) 2^(2 FLINT_BITS wn) */
        flint_mpn_zero(u2, wn);
        flint_mpn_sqrhigh(u2 + wn, s, wn);
        u2[2 * wn] = 1;

        /* cosh(t) = sqrt(1 + sinh(t)^2) */
        mpn_sqrtrem(rt, rem, u2, 2 * wn + 1);

        /* exp(t) = sinh(t) + cosh(t) */
        mpn_add_n(y, rt, s, wn + 1);

        TMP_END;
    }
}


/* multiply y (ylen limbs, fraction plus a units limb) by the factors
   (1 + 2^-used[j]) for j0 <= j < num: each is a shift and an add.
   The used indices are increasing; process them in chunks sharing the
   same limb offset.  The value stays below e, so no carry ever leaves
   the units limb.  sh is scratch of ylen limbs. */
static void
_fixed_exp_recon(nn_ptr y, nn_ptr sh, slong ylen, const slong * used,
    slong j, slong num)
{
    if (j == 0 && num > 0 && used[0] == 0)
    {
        mpn_add_n(y, y, y, ylen);               /* factor 2 */
        j = 1;
    }

    for (; j < num && used[j] < FLINT_BITS; j++)
    {
        mpn_rshift(sh, y, ylen, used[j]);
        mpn_add_n(y, y, sh, ylen);
    }

    while (j < num)
    {
        slong q = used[j] / FLINT_BITS;
        slong stop = (q + 1) * FLINT_BITS;
        slong len = ylen - q;

        for (; j < num && used[j] < stop; j++)
        {
            slong b = used[j] - q * FLINT_BITS;
            ulong cy;

            if (b != 0)
                mpn_rshift(sh, y + q, len, b);
            else
                flint_mpn_copyi(sh, y + q, len);

            cy = mpn_add_n(y, y, sh, len);
            mpn_add_1(y + len, y + len, q, cy);
        }
    }
}


/* Greedy bitwise reduction of t (wn limbs) by L_0, ..., L_r,
   recording the used indices; returns their number.  The subtractions
   are decided one FLINT_BITS window at a time on a single-limb model
   h of the remainder: h is kept a lower bound on the remainder's
   current window limb, with upper slack e growing by one per model
   subtraction (each table limb is a truncation), so h > lt proves
   L_i <= t and h <= lt - e - 1 proves t < L_i.  Decisions inside the
   ambiguity band [lt - e, lt] -- a vanishing fraction of random
   steps, plus the last few steps of each window where lt itself is
   small -- fall back to the exact compare-subtract, as does the first
   step of each window (the remainder may still carry one bit above
   it).  Decided subtractions are recorded and applied to the
   full-precision remainder in unconditional batches, so the bulk of
   the loop touches one limb per step and the full-width work is
   branch-free.  The decisions agree exactly with the plain
   full-precision greedy. */
static slong
_fixed_exp_reduce(nn_ptr t, slong wn, int r, slong * used)
{
    slong num = 0, bj = 0, nc = _fixed_exp_logs_n;
    slong i, c;

    for (c = 0; FLINT_BITS * c <= r; c++)
    {
        slong i0 = FLINT_BITS * c;
        slong i1 = FLINT_MIN((slong) r, FLINT_BITS * (c + 1) - 1);
        ulong h, e;
        nn_srcptr Lp;

        /* apply pending subtractions from the previous window */
        for (; bj < num; bj++)
            mpn_sub_n(t, t, LOGS_L(used[bj]) + (nc - wn), wn);

        /* exact step at the window boundary */
        Lp = LOGS_L(i0) + (nc - wn);
        if (mpn_cmp(t, Lp, wn) >= 0)
        {
            mpn_sub_n(t, t, Lp, wn);
            used[num++] = i0;
        }
        bj = num;

        h = t[wn - 1 - c];
        e = 0;

        for (i = i0 + 1; i <= i1; i++)
        {
            ulong lt = _fixed_exp_logs[i * nc + (nc - 1 - c)];

            if (h - lt + e <= 2 * e)
            {
                /* ambiguity band: apply pending and decide exactly */
                for (; bj < num; bj++)
                    mpn_sub_n(t, t, LOGS_L(used[bj]) + (nc - wn), wn);
                Lp = LOGS_L(i) + (nc - wn);
                if (mpn_cmp(t, Lp, wn) >= 0)
                {
                    mpn_sub_n(t, t, Lp, wn);
                    used[num++] = i;
                }
                bj = num;
                h = t[wn - 1 - c];
                e = 0;
            }
            else
            {
                ulong sub = (h > lt);

                h -= (lt + 1) & (-sub);
                used[num] = i;
                num += sub;
                e += sub;
            }
        }
    }

    for (; bj < num; bj++)
        mpn_sub_n(t, t, LOGS_L(used[bj]) + (nc - wn), wn);

    /* the truncated table lets the remainder creep marginally above
       L_r; one extra subtraction restores t < L_r < 2^-r */
    {
        nn_srcptr Lr = LOGS_L(r) + (nc - wn);

        if (mpn_cmp(t, Lr, wn) >= 0)
        {
            mpn_sub_n(t, t, Lr, wn);
            used[num++] = r;
        }
    }

    return num;
}

#if FLINT_BITS == 64
#include "exp_bitwise_rs_small.inc"
#endif

/* smallest n at which each reduction parameter becomes optimal
   (generated by src/fixed/tune/tune-bitwise-r.c; x86-64 defaults) */
static const int _fixed_exp_bitwise_rs_r_tab[] = {32, 64, 128, 192, 256, 320};
static const short _fixed_exp_bitwise_rs_n_tab[] = {1, 12, 39, 115, 137, 154};

void
fixed_exp_bitwise_rs(nn_ptr res, nn_srcptr x, slong n, int r)
{
    slong num, wn;
    nn_ptr t, y, sh;
    slong * used;
    TMP_INIT;

    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT(r == 0 || r >= 32);
    /* on 32-bit limbs, n = 1 would clamp r to 16 below, violating
       the 2^-32 contract of the series functions */
    FLINT_ASSERT(FLINT_BITS == 64 || n >= 2);

    if (r == 0)
    {
        /* tuned default; the largest tabulated value serves all
           larger n */
        slong j;

        for (j = 0; j + 1 < (slong) (sizeof(_fixed_exp_bitwise_rs_n_tab)
                / sizeof(short))
                && n >= _fixed_exp_bitwise_rs_n_tab[j + 1]; j++)
            ;
        r = _fixed_exp_bitwise_rs_r_tab[j];
    }

#if FLINT_BITS == 64
#ifndef FIXED_EXP_BITWISE_NO_SMALL
    if (n <= 7)
    {
        _fixed_exp_bitwise_rs_small(res, x, n, r);
        return;
    }
#endif
#endif

    r = FLINT_MIN((slong) r, FLINT_BITS * n - 16);

    wn = n;

    _fixed_exp_logs_ensure(wn, r);

    TMP_START;
    t = TMP_ALLOC((wn + 2 * (n + 1)) * sizeof(ulong));
    y = t + wn;
    sh = y + (n + 1);
    used = TMP_ALLOC((r + 2) * sizeof(slong));

    flint_mpn_copyi(t, x, n);

    num = _fixed_exp_reduce(t, wn, r, used);

    /* exp of the reduced argument */
    _fixed_exp_reduced(y, t, wn, r, EXP_USE_SINH(wn, r));

    _fixed_exp_recon(y, sh, n + 1, used, 0, num);

    flint_mpn_copyi(res, y, n + 1);

    TMP_END;
}
