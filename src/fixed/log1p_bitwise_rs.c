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
#include "fixed.h"

/* log(1 + x) by the dual of the bitwise exp reduction (the L-mode
   BKM recurrence): greedily multiply P (starting from 1) by factors
   1 + 2^-i, i = 1..r, while P (1 + 2^-i) <= X = 1 + x.  This greedy
   is exactly the exp reduction applied to log(X/P), so each factor
   is used at most once (L_{i-1} < 2 L_i), one extra conditional step
   at i = r absorbs the truncation creep, and afterwards
   X/P in [1, 1 + 2^-r).  The same cached table of
   L_i = log(1 + 2^-i) serves both directions.

   The residual is obtained by fusing the normalization with the
   atanh transformation: log(X/P) = 2 atanh((X - P)/(X + P)), where
   the numerator is the deficit D = X - P (maintained exactly through
   the reduction) and the single division by S = X + P in [2, 4)
   yields t = D/S < 2^-r, meeting the atanh series contract directly;
   the odd atanh series needs half the terms of log(1 + u).  Finally
   log X = sum of used L_i + 2 atanh(t).

   Decisions run one FLINT_BITS window at a time on a single-limb
   model, as in exp: h is the deficit's window limb and
   lt = MPN_RIGHT_SHIFT_LOW(P[n], P[n-1], i mod FLINT_BITS) is the
   window limb of P >> i, computable from P's top two limbs alone
   (independently of the window, since P < 2).  Because accepted
   factors update D immediately, the model resyncs exactly after
   every accept, so the ambiguity band degenerates to the single
   value h == lt; rejected steps (about half) cost a handful of
   register operations.

   P grows gradually: multiplying by 1 + 2^-i moves the exact least
   significant bit down by i, so the significant fraction length of P
   is the sum of the accepted i divided by FLINT_BITS, capped at n
   limbs (only past that cap does P get truncated at all).  All
   accept-side operations (the shift, the additions into P and the
   subtraction from D, and the table additions into the accumulator)
   run on the significant slices with a carry/borrow ripple across
   the remaining limbs.

   Error accounting (output ulp 2^(-FLINT_BITS n), all one-sided
   sources): at most r + 2 table entries are added, each short of the
   true logarithm by < 2 ulp; each accepted factor truncates P by
   < 1 ulp, entering the result as log(Pi/P) < (r + 2) ulp in total;
   the division floors once; the doubled atanh error contributes
   <= 30 ulp.  Total < 3 (r + 2) + 40, bounded by
   FIXED_LOG1P_BITWISE_RS_MAX_ERR(n, r) = 3 r + 64. */

#define LOGS_L(i) (_fixed_exp_logs + (i) * _fixed_exp_logs_n)

#if FLINT_BITS == 64

#ifndef FIXED_LOG1P_BITWISE_NO_SMALL
#include "log1p_bitwise_rs_small.inc"
#endif


/* n = 1: everything in registers.  With r <= 48 there is a single
   reduction window; each step is a funnel shift of P's top limbs, a
   branch-free compare-select, and masked updates of the deficit, the
   product and the (single-limb) logarithm accumulator.  The final
   division is a 3-by-2 limb schoolbook step, and at t < 2^-16 the
   atanh series collapses to t + t^3/3 (the t^5 term is below the
   64-bit ulp already for r = 16). */
static void
_fixed_log1p_bitwise_rs_1(nn_ptr res, nn_srcptr x, int r)
{
    ulong p0 = 0, d0 = x[0], a0 = 0;
    ulong t0, z, w, q, rr, h, l, sn1, sn0, n2, n1, s1, s0, cy;
    slong i, nc;
    int k;

    r = FLINT_MIN(r, FLINT_BITS * 1 - 16);

    _fixed_exp_logs_ensure(1, r);
    nc = _fixed_exp_logs_n;

    for (i = 1; i <= r + 1; i++)
    {
        slong ii = FLINT_MIN(i, (slong) r);
        ulong v = MPN_RIGHT_SHIFT_LOW(UWORD(1), p0, (int) ii);
        ulong lt = _fixed_exp_logs[ii * nc + (nc - 1)];
        ulong m = -(ulong) (v <= d0);

        d0 -= v & m;
        p0 += v & m;      /* no carry: P + v <= X < 2 */
        a0 += lt & m;     /* no carry: acc < log 2 */
    }

    /* t = d0 2^64 / S, S = X + P = (s1 : s0) with s1 in {2, 3}:
       one normalized 3-by-2 division step with two corrections */
    add_ssaaaa(s1, s0, UWORD(1), x[0], UWORD(1), p0);
    if (d0 == 0)
    {
        res[0] = a0;
        return;
    }
    k = flint_clz(s1);
    sn1 = (s1 << k) | (s0 >> (FLINT_BITS - k));
    sn0 = s0 << k;
    n2 = d0 >> (FLINT_BITS - k);
    n1 = d0 << k;
    udiv_qrnnd(q, rr, n2, n1, sn1);
    umul_ppmm(h, l, q, sn0);
    while (h > rr || (h == rr && l > 0))
    {
        q--;
        sub_ddmmss(h, l, h, l, UWORD(0), sn0);
        cy = rr + sn1;
        if (cy < rr)
            break;      /* remainder overflowed: q is now correct */
        rr = cy;
    }
    t0 = q;

    /* atanh(t) = t + t^3/3 + O(t^5), t < 2^-16 */
    umul_ppmm(z, l, t0, t0);
    umul_ppmm(w, l, t0, z);
    res[0] = a0 + 2 * (t0 + w / 3);
}

/* n = 2: two-limb register version.  Window 0 (i < 64) works on
   two-limb funnel shifts of P and branch-free two-limb
   compare-selects; window 1 (i >= 64) reduces to a single
   significant limb.  The division falls back to mpn_tdiv_qr and the
   atanh series runs four terms (enough down to r = 16). */
static void
_fixed_log1p_bitwise_rs_2(nn_ptr res, nn_srcptr x, int r)
{
    ulong p1 = 0, p0 = 0, d1 = x[1], d0 = x[0], a1 = 0, a0 = 0;
    ulong t2[2], w[2], S[3], nd[4], rem[4];
    slong i, nc;

    r = FLINT_MIN(r, FLINT_BITS * 2 - 16);

    _fixed_exp_logs_ensure(2, r);
    nc = _fixed_exp_logs_n;

#define STEP2(ii, v1, v0) \
    do { \
        nn_srcptr Lq = _fixed_exp_logs + (ii) * nc + (nc - 2); \
        ulong bw, e1, e0, m; \
        sub_dddmmmsss(bw, e1, e0, UWORD(0), d1, d0, \
            UWORD(0), v1, v0); \
        m = ~bw;                          /* accept iff no borrow */ \
        d1 = (d1 & bw) | (e1 & m); \
        d0 = (d0 & bw) | (e0 & m); \
        add_ssaaaa(p1, p0, p1, p0, (v1) & m, (v0) & m); \
        add_ssaaaa(a1, a0, a1, a0, Lq[1] & m, Lq[0] & m); \
    } while (0)

    for (i = 1; i <= FLINT_MIN((slong) r, FLINT_BITS - 1); i++)
    {
        ulong v1 = MPN_RIGHT_SHIFT_LOW(UWORD(1), p1, (int) i);
        ulong v0 = MPN_RIGHT_SHIFT_LOW(p1, p0, (int) i);
        STEP2(i, v1, v0);
    }

    if (r >= FLINT_BITS)
    {
        STEP2(FLINT_BITS, UWORD(1), p1);      /* i = 64: v = P >> 64 */

        for (i = FLINT_BITS + 1; i <= (slong) r; i++)
        {
            ulong v0 = MPN_RIGHT_SHIFT_LOW(UWORD(1), p1,
                (int) (i - FLINT_BITS));
            STEP2(i, UWORD(0), v0);
        }
    }

    /* one extra step at i = r absorbs the truncation creep */
    if (r > FLINT_BITS)
    {
        ulong v0 = MPN_RIGHT_SHIFT_LOW(UWORD(1), p1,
            (int) (r - FLINT_BITS));
        STEP2(r, UWORD(0), v0);
    }
    else if (r == FLINT_BITS)
        STEP2(r, UWORD(1), p1);           /* v = P >> 64 exactly */
    else
    {
        ulong v1 = MPN_RIGHT_SHIFT_LOW(UWORD(1), p1, (int) r);
        ulong v0 = MPN_RIGHT_SHIFT_LOW(p1, p0, (int) r);
        STEP2(r, v1, v0);
    }

#undef STEP2

    /* t = D 2^128 / (X + P) via mpn_tdiv_qr; S = X + P in [2, 4) */
    add_sssaaaaaa(S[2], S[1], S[0], UWORD(1), x[1], x[0],
        UWORD(1), p1, p0);
    nd[0] = 0;
    nd[1] = 0;
    nd[2] = d0;
    nd[3] = d1;
    t2[0] = t2[1] = 0;
    if (d1 != 0 || d0 != 0)
        mpn_tdiv_qr(t2, rem, 0, nd, 4, S, 3);

    /* generated atanh series for t < 2^-16 (single denominator
       division by truncated-inverse multiplication) */
    _fixed_atanh_rs16(w, t2, 2);

    /* res = acc + 2 atanh(t) */
    add_ssaaaa(res[1], res[0], a1, a0,
        (w[1] << 1) | (w[0] >> (FLINT_BITS - 1)), w[0] << 1);
}

#endif /* FLINT_BITS == 64 */

/* smallest n at which each reduction parameter becomes optimal
   (generated by src/fixed/tune/tune-bitwise-r.c; x86-64 defaults) */
static const int _fixed_log1p_bitwise_rs_r_tab[] = {16, 32, 64, 128};
static const short _fixed_log1p_bitwise_rs_n_tab[] = {1, 4, 18, 64};

void
fixed_log1p_bitwise_rs(nn_ptr res, nn_srcptr x, slong n, int r)
{
#if FLINT_BITS == 64
    int r0;
#endif
    slong i, c, wn, nc, pl, lsb_bits, qr, ds;
    nn_ptr P, D, acc, sh, t, nd;
    TMP_INIT;

    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT(r == 0 || r >= 16);

#if FLINT_BITS == 64
    r0 = r;
#endif

    if (r == 0)
    {
        /* tuned default; the largest tabulated value serves all
           larger n */
        slong j;

        for (j = 0; j + 1 < (slong) (sizeof(_fixed_log1p_bitwise_rs_n_tab)
                / sizeof(short))
                && n >= _fixed_log1p_bitwise_rs_n_tab[j + 1]; j++)
            ;
        r = _fixed_log1p_bitwise_rs_r_tab[j];
    }

    /* r < 32 (a shorter reduction with a slightly longer series) is
       available in the register implementations below, whose inline
       atanh evaluations only need t < 2^-16; the general path uses
       fixed_atanh_rs and hence at least r = 32 */
#if FLINT_BITS == 64
#ifndef FIXED_LOG1P_BITWISE_NO_SMALL
    /* r = 0: the specialized sizes run with a compile-time constant r
       and the hand-written atanh series built for it.  n = 1 and n = 2
       keep their own register code, whose atanh evaluation is already
       inline. */
    if (n <= 4 && r0 == 0 && _fixed_log1p_bitwise_rs_opt(res, x, n))
        return;
#endif

    if (n == 1)
    {
        _fixed_log1p_bitwise_rs_1(res, x, (int) FLINT_MAX(r, 16));
        return;
    }
    if (n == 2)
    {
        _fixed_log1p_bitwise_rs_2(res, x, (int) FLINT_MAX(r, 16));
        return;
    }
#ifndef FIXED_LOG1P_BITWISE_NO_SMALL
    if (n <= 7)
    {
        switch (n)
        {
            case 3: _fixed_log1p_bitwise_rs_3(res, x, (int) r); return;
            case 4: _fixed_log1p_bitwise_rs_4(res, x, (int) r); return;
            case 5: _fixed_log1p_bitwise_rs_5(res, x, (int) r); return;
            case 6: _fixed_log1p_bitwise_rs_6(res, x, (int) r); return;
            default: _fixed_log1p_bitwise_rs_7(res, x, (int) r); return;
        }
    }
#endif
#endif /* FLINT_BITS == 64 */

    r = FLINT_MAX(r, 32);
    r = FLINT_MIN((slong) r, FLINT_BITS * n - 16);

    wn = n + 1;

    _fixed_exp_logs_ensure(n, r);
    nc = _fixed_exp_logs_n;

    TMP_START;
    P = TMP_ALLOC((wn + wn + n + n + n + 2 * n + 2) * sizeof(ulong));
    sh = P + wn;
    D = sh + wn;
    acc = D + n;
    t = acc + n;
    nd = t + n;

    /* P = 1, deficit D = X - P = x, accumulator = 0 */
    flint_mpn_zero(P, n);
    P[n] = 1;
    flint_mpn_copyi(D, x, n);
    flint_mpn_zero(acc, n);
    pl = 0;             /* significant fraction limbs of P */
    lsb_bits = 0;       /* depth of P's exact lowest fraction bit */

/* materialize v = trunc(P >> ii) into sh as (P >> b) (limb offset q
   applied at use); z_ = index below which (P >> b) is exactly zero */
#define MAKE_SH(ii, q_, b_, z_) \
    do { \
        (q_) = (ii) / FLINT_BITS; \
        (b_) = (ii) - (q_) * FLINT_BITS; \
        (z_) = FLINT_MAX(0, n - pl - 1); \
        if ((b_) != 0) \
            mpn_rshift(sh + (z_), P + (z_), wn - (z_), b_); \
        else \
            flint_mpn_copyi(sh + (z_), P + (z_), wn - (z_)); \
    } while (0)

/* accept factor 1 + 2^-ii, with sh already holding P >> b:
   D -= v, P += v, acc += L_ii, on the significant slices */
#define APPLY(ii, q_, b_, z_) \
    do { \
        slong len_ = FLINT_MIN(n, wn - (q_)); \
        slong j0_ = FLINT_MAX(0, (z_) - (q_)); \
        ulong cy_; \
        if (j0_ < len_) \
        { \
            cy_ = mpn_sub_n(D + j0_, D + j0_, sh + (q_) + j0_, \
                len_ - j0_); \
            if (cy_ != 0) \
                mpn_sub_1(D, D, j0_, cy_); \
            cy_ = mpn_add_n(P + j0_, P + j0_, sh + (q_) + j0_, \
                len_ - j0_); \
            if (cy_ != 0) \
                mpn_add_1(P + len_, P + len_, wn - len_, cy_); \
        } \
        cy_ = mpn_add(acc, acc, n, LOGS_L(ii) + (nc - n), n - (q_)); \
        (void) cy_;                 /* acc < log 2 < 1 */ \
        lsb_bits += (ii); \
        pl = FLINT_MIN(n, (lsb_bits + FLINT_BITS - 1) / FLINT_BITS); \
    } while (0)

/* exact compare-accept of factor ii */
#define EXACT_STEP(ii) \
    do { \
        slong q_, b_, z_, len_, j0_, k_; \
        int ge_; \
        MAKE_SH(ii, q_, b_, z_); \
        len_ = FLINT_MIN(n, wn - q_); \
        j0_ = FLINT_MAX(0, z_ - q_); \
        /* v's limbs are sh[q_ + j] for j0_ <= j < len_, zero \
           elsewhere on the n-limb grid; note sh[n] = P[n] >> b_ \
           is nonzero only for b_ = 0, in which case v's top limb \
           n - q_ is included in len_ */ \
        ge_ = 0; \
        for (k_ = n - 1; k_ >= len_; k_--) \
            if (D[k_] != 0) \
            { \
                ge_ = 1;            /* D has a higher limb set */ \
                break; \
            } \
        if (!ge_ && j0_ < len_) \
        { \
            int c_ = mpn_cmp(D + j0_, sh + q_ + j0_, len_ - j0_); \
            /* on tie, D's lower limbs are >= v's (which are zero) */ \
            ge_ = (c_ >= 0); \
        } \
        else if (!ge_) \
            ge_ = 1;                /* v truncates to zero: cannot \
                                       happen for i <= r (clamped) */ \
        if (ge_) \
            APPLY(ii, q_, b_, z_); \
    } while (0)

    for (c = 0; FLINT_BITS * c <= r; c++)
    {
        slong i0 = (c == 0) ? 1 : FLINT_BITS * c;
        slong i1 = FLINT_MIN((slong) r, FLINT_BITS * (c + 1) - 1);
        ulong h;

        if (i0 > (slong) r)
            break;

        /* exact step at the window boundary */
        EXACT_STEP(i0);

        h = D[n - 1 - c];

        for (i = i0 + 1; i <= i1; i++)
        {
            slong b = i - FLINT_BITS * c;
            ulong lt;

            /* unlike exp (where t < L_i < 2^-i), the deficit obeys
               only D < P 2^-i < 2^(1-i): a bit above the window can
               persist, but then v = P >> i < 2^(1-i) <= D is a
               certain accept */
            if (c > 0 && D[n - c] != 0)
            {
                slong q_, b_, z_;
                MAKE_SH(i, q_, b_, z_);
                APPLY(i, q_, b_, z_);
                h = D[n - 1 - c];
                continue;
            }

            lt = MPN_RIGHT_SHIFT_LOW(P[n], P[n - 1], (int) b);

            if (h > lt)
            {
                slong q_, b_, z_;
                MAKE_SH(i, q_, b_, z_);
                APPLY(i, q_, b_, z_);
                h = D[n - 1 - c];
            }
            else if (h == lt)
            {
                /* the lower limbs decide */
                EXACT_STEP(i);
                h = D[n - 1 - c];
            }
            /* else certain reject: D < (h+1) u <= lt u <= P >> i */
        }
    }

    /* the truncated table and P let the residual creep marginally
       above 2^-r; one extra step restores X/P < 1 + 2^-r */
    EXACT_STEP(r);

#undef MAKE_SH
#undef APPLY
#undef EXACT_STEP

    /* t = D / (X + P) < 2^-r: the single division fuses the
       normalization by P with the atanh transformation */
    {
        flint_mpn_copyi(sh, P, wn);
        sh[n] += 1 + mpn_add_n(sh, sh, x, n);   /* S = X + P in [2, 4) */

        qr = 0;
        while (qr < n && D[n - 1 - qr] == 0)
            qr++;
        ds = 2 * n - qr;

        flint_mpn_zero(t, n);
        if (ds >= wn && qr < n)
        {
            flint_mpn_zero(nd, n);
            flint_mpn_copyi(nd + n, D, n - qr);
            /* quotient has ds - wn + 1 = n - qr limbs */
            mpn_tdiv_qr(t, nd, 0, nd, ds, sh, wn);
        }
    }

    /* log X = acc + 2 atanh(t); on 32-bit limbs n = 1 clamps r to
       16 above, where the wider-range atanh applies */
    if (r < 32)
        _fixed_atanh_rs16(sh, t, n);
    else
        fixed_atanh_rs(sh, t, n);
    mpn_lshift(sh, sh, n, 1);
    mpn_add_n(res, acc, sh, n);

    TMP_END;
}
