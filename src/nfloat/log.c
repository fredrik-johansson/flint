/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "elem.h"

/*
    Logarithms.

    For u = M 2^E, M in [1/2, 1), log u = (E - 1) log 2 + log1p(2M - 1),
    with the log1p kernel on v = 2M - 1 in [0, 1) at nk fraction limbs
    and log 2 the floor L of log 2 B^(nk + 1), whose error times
    |E - 1| < 2^(FLINT_BITS - 2) stays below a quarter ulp of B^-nk. Near 1 (E = 0 or 1)
    with |u - 1| < 2^-z the sum cancels up to z + 1 bits, which the
    kernel absorbs by running at z + 1 more bits; from enough leading zero
    bits (by kernel limbs, the crossover of mp_real/log.c),
    log(1 +- s) = +-2 atanh(s / (2 +- s)) by the series of mp_real
    instead. Close enough to 1, log(1 + D) = D (1 + delta) with a tiny
    relative error.

    log2 u = (E - 1) + log1p(2M - 1) / log 2, the division a
    multiplication by 1/log 2 (the mp_real constant), exact in its
    integral part.
*/

/* T = log2(e) * (Y, w) = Y + Y R at the same scale, R the fraction of
   1/log 2: (T, w + 1) with err ulps updated (the value 1.443.. times
   the error, plus the truncations) */
static void
_nfloat_mul_ilog2(nn_ptr T, nn_srcptr Y, slong w, ulong * err)
{
    ulong t[NFLOAT_ELEM_MAX_LIMBS + 1];

    flint_mpn_mulhigh_n(t, Y, _mp_real_const_ptr(MP_REAL_CONST_ID_INV_LOG2_FRAC, w), w);
    T[w] = mpn_add_n(T, Y, t, w);
    *err = *err + (*err >> 1) + 3;
}

/*
    log(1 + s) (neg = 0) resp. log(1 - s) (neg = 1) for the fraction
    (s, w), s < 1 resp. s <= 1/2, s known within serr ulps (as a lower
    bound for neg = 0, an upper bound for neg = 1), w at least the kernel
    limbs for z + 1 extra bits where s < 2^-z. Sets (y, w + 1), returns the
    exponent of its unit and sets *err.
*/
static slong
_nfloat_log1pm_fixed(nn_ptr y, ulong * err, nn_srcptr s, ulong serr, slong w, int neg)
{
    slong z = _mp_real_elem_lzb(s, w);

    FLINT_ASSERT(z != WORD_MAX);

    if (z >= _mp_real_log_series_min_z(w, !neg) && z >= 8)
    {
        /* +-2 atanh(t), t = s / (2 +- s) within serr + 1 ulps, atanh' < 1 +
           2^-15 */
        ulong t[NFLOAT_ELEM_MAX_LIMBS];
        ulong ea;

        _mp_real_elem_half_ratio(t, s, w, neg);
        if (z >= 32)
            _mp_real_atanh_rs(y, &ea, t, w);
        else
        {
            _mp_real_series_atan(y, t, w, (flint_bitcnt_t) z, MP_REAL_SERIES_ATANH);
            ea = 6;
        }
        y[w] = 0;
        *err = ea + serr + 2;
        return 1 - FLINT_BITS * w;
    }
    else if (!neg)
    {
        /* the kernel (log1p' <= 1) */
        _mp_real_log1p_kernel(y, err, s, w);
        y[w] = 0;
        *err += serr;
        return -FLINT_BITS * w;
    }
    else
    {
        /* -log(1 - s) = log 2 - log1p(1 - 2s), 1 - 2s within 2 serr ulps,
           log1p' <= 1; F at w fraction limbs, L at w + 1 */
        ulong v[NFLOAT_ELEM_MAX_LIMBS + 1];
        ulong T[NFLOAT_ELEM_MAX_LIMBS + 2];
        ulong ek;

        v[w] = mpn_lshift(v, s, w, 1);
        mpn_neg(v, v, w);

        T[0] = 0;
        _mp_real_log1p_kernel(T + 1, &ek, v, w);
        mpn_sub_n(T, _mp_real_const_ptr(MP_REAL_CONST_ID_LOG2, w + 1), T, w + 1);
        _nfloat_copy_limbs(y, T + 1, w);
        y[w] = 0;
        /* L within one ulp of B^-(w+1), the dropped limb one ulp */
        *err = ek + 2 * serr + 2;
        return -FLINT_BITS * w;
    }
}

/* assembles the result into res, or stores it in fix */
static int
_nfloat_log_out(nfloat_ptr res, _nfloat_fix_t fix, nn_ptr y, slong len, slong e, ulong elo, ulong ehi, int sgnbit, gr_ctx_t ctx)
{
    if (fix == NULL)
        return _nfloat_set_mpn_err(res, y, len, e, elo, ehi, sgnbit, ctx);

    _nfloat_copy_limbs(fix->d, y, len);
    fix->len = len;
    fix->e = e;
    fix->err = FLINT_MAX(elo, ehi);
    fix->sgnbit = sgnbit;
    return GR_SUCCESS;
}

/* res (or fix) = log(u) resp. log2(u) for u = (d, dn) 2^(E - FLINT_BITS dn)
   normalized, u != 1, at extra more bits of precision; ehi_extra units are
   added to the upper error (for log1p of huge arguments) */
static int
_nfloat_log_mpn(nfloat_ptr res, _nfloat_fix_t fix, nn_srcptr d, slong dn, slong E, int base2, slong extra, int ehi_extra, gr_ctx_t ctx)
{
    ulong y[NFLOAT_ELEM_MAX_LIMBS + 3];
    ulong s[NFLOAT_ELEM_MAX_LIMBS + 2];
    ulong err;
    slong n = NFLOAT_CTX_NLIMBS(ctx);
    slong w, z, ue;
    int trunc, neg;

    if (E == 0 || E == 1)
    {
        /* s = |u - 1| with u - 1 = (2M - 1) for E = 1, -(1 - M) for E = 0 */
        neg = (E == 0);

        if (!neg)
        {
            /* the fraction bits of 2M, exactly */
            ulong t[NFLOAT_ELEM_MAX_LIMBS + 2];
            slong j;
            _nfloat_get_fixed(t, dn + 1, dn, d, dn, 1);
            for (j = dn - 1; j >= 0 && t[j] == 0; j--)
                ;
            if (j < 0)
            {
                if (fix != NULL)
                {
                    fix->len = 0;
                    return GR_SUCCESS;
                }
                return nfloat_zero(res, ctx);
            }
            z = FLINT_BITS * (dn - 1 - j) + flint_clz(t[j]);
        }
        else
        {
            /* 1 - M: the leading ones of M */
            slong j;
            z = 0;
            for (j = dn - 1; j >= 0 && d[j] == UWORD_MAX; j--)
                z += FLINT_BITS;
            if (j >= 0)
                z += flint_clz(~d[j]);
            z--;
        }

        /* log(1 + D) = D (1 + delta), |delta| < |D| < 2^-z (in the
           default rounding mode, also when D is below the function
           precision) */
        if (!base2 && fix == NULL && (z >= FLINT_BITS * n || (z > NFLOAT_CTX_FUNC_PREC(ctx)
                && !NFLOAT_CTX_HAS_DIRECTED_ROUNDING(ctx))))
        {
            ulong D[NFLOAT_MAX_ALLOC];
            int status;

            /* D exactly as an nfloat (at most FLINT_BITS n bits) */
            FLINT_ASSERT(dn == n);
            if (!neg)
            {
                _nfloat_get_fixed(s, dn + 1, dn, d, dn, 1);
                status = nfloat_set_mpn_2exp(D, s, dn, 0, 0, ctx);
            }
            else
            {
                _nfloat_copy_limbs(s, d, dn);
                mpn_neg(s, s, dn);
                status = nfloat_set_mpn_2exp(D, s, dn, 0, 1, ctx);
            }

            FLINT_ASSERT(status == GR_SUCCESS);
            (void) status;
            /* log(1 + D) < D for D > 0; |log(1 + D)| > |D| for D < 0 */
            return _nfloat_set_tiny_rel_err(res, D, neg ? 1 : -1, ctx);
        }

        /* s at w fraction limbs for z + 1 extra bits */
        w = _nfloat_elem_limbs(ctx, FLINT_MIN(z, FLINT_BITS * n + 64) + 1 + 2 * base2 + extra);

        if (!neg)
        {
            trunc = _nfloat_get_fixed(s, w + 1, w, d, dn, 1);
            s[w] = 0;
        }
        else
        {
            /* 1 - M rounded up: the negated truncation of M */
            trunc = _nfloat_get_fixed(s, w, w, d, dn, 0);
            mpn_neg(s, s, w);
        }

        ue = _nfloat_log1pm_fixed(y, &err, s, trunc, w, neg);

        if (base2)
        {
            ulong T[NFLOAT_ELEM_MAX_LIMBS + 3];
            _nfloat_mul_ilog2(T, y, w + 1, &err);
            return _nfloat_log_out(res, fix, T, w + 2, ue, err, err + ehi_extra, neg, ctx);
        }

        return _nfloat_log_out(res, fix, y, w + 1, ue, err, err + ehi_extra, neg, ctx);
    }
    else
    {
        /* (E - 1) log 2 + F, F = log1p(2M - 1), in fixed point with one
           integral limb and w + 1 fraction limbs: T = |c| L +- F */
        ulong T[NFLOAT_ELEM_MAX_LIMBS + 4];
        ulong cc, cy;
        slong c = E - 1;

        w = _nfloat_elem_limbs(ctx, extra);

        trunc = _nfloat_get_fixed(s, w + 1, w, d, dn, 1);   /* 2M at w limbs */
        T[0] = 0;
        _mp_real_log1p_kernel(T + 1, &err, s, w);
        T[w + 1] = 0;
        err += 2 * trunc;
        cc = (c < 0) ? -(ulong) c : (ulong) c;

        if (base2)
        {
            /* (E - 1) + F / log 2: the integer part is exact */
            ulong G[NFLOAT_ELEM_MAX_LIMBS + 3];

            _nfloat_mul_ilog2(G, T + 1, w, &err);
            /* G = F / log 2 < 1 (up to the error) at w fraction limbs */
            if (c > 0)
            {
                G[w] += cc;
                return _nfloat_log_out(res, fix, G, w + 1, -FLINT_BITS * w, err, err + ehi_extra, 0, ctx);
            }
            else
            {
                /* |c| - G with |c| >= 2 */
                mpn_neg(G, G, w + 1);
                G[w] += cc;
                return _nfloat_log_out(res, fix, G, w + 1, -FLINT_BITS * w, err, err + ehi_extra, 1, ctx);
            }
        }
        else
        {
            nn_srcptr L = _mp_real_const_ptr(MP_REAL_CONST_ID_LOG2, w + 1);

            if (c > 0)
            {
                MP_REAL_ADDMUL_1(cy, T, L, w + 1, cc);
                T[w + 1] += cy;
                neg = 0;
            }
            else
            {
                /* F - |c| L < 0, negated */
                MP_REAL_SUBMUL_1(cy, T, L, w + 1, cc);
                T[w + 1] -= cy;
                mpn_neg(T, T, w + 2);
                neg = 1;
            }

            /* |c| (log 2 - L) < 2^(FLINT_BITS - 2) B^-(w+1) and the dropped limb */
            err += 2;
            return _nfloat_log_out(res, fix, T + 1, w + 1, -FLINT_BITS * w, err, err + ehi_extra, neg, ctx);
        }
    }
}

#if FLINT_BITS == 64
/*
    In-register fast paths at one and two limbs (64-bit only): log u =
    c L +- F with c = E - 1, F = log1p(2M - 1) by the table-driven kernel
    of mp_real at two resp. three limbs (v = 2M - 1 exact) and L the
    floor of log 2 at three resp. four fraction limbs, whose error times
    |c| < 2^62 stays below a quarter ulp of the kernel. The sum, with one
    integral limb, drops its lowest limb (one ulp). Near 1 (E = 0 or 1)
    the cancellation costs as many bits as |log u| has leading zeros,
    so the fast paths require |u - 1| >= 2^-32 resp. 2^-40, where the
    relative error stays below 2^-90 resp. 2^-146; the general code
    handles the rest.
*/
static int
_nfloat_log_fast_1(nfloat_ptr res, ulong d, slong E, gr_ctx_t ctx)
{
    nn_srcptr L = _mp_real_const_ptr(MP_REAL_CONST_ID_LOG2, 3);
    ulong f1, f0, cc, p3, p2, p1, p0, h, l;
    slong c = E - 1;

    _mp_real_small_log1p_2(&f1, &f0, d << 1, 0);

    cc = (c < 0) ? -(ulong) c : (ulong) c;
    umul_ppmm(p1, p0, cc, L[0]);
    umul_ppmm(h, l, cc, L[1]);
    add_ssaaaa(p2, p1, h, l, UWORD(0), p1);
    umul_ppmm(h, l, cc, L[2]);
    add_ssaaaa(p3, p2, h, l, UWORD(0), p2);
    (void) p0;

    if (c >= 0)
        add_sssaaaaaa(p3, p2, p1, p3, p2, p1, UWORD(0), f1, f0);
    else
        sub_dddmmmsss(p3, p2, p1, p3, p2, p1, UWORD(0), f1, f0);

    return _nfloat_small_set_1(res, p3, p2, p1, 0, MP_REAL_SMALL_LOG_ERR_2 + 2, c < 0, ctx);
}

static int
_nfloat_log_fast_2(nfloat_ptr res, ulong d1, ulong d0, slong E, gr_ctx_t ctx)
{
    nn_srcptr L = _mp_real_const_ptr(MP_REAL_CONST_ID_LOG2, 4);
    ulong f2, f1, f0, cc, p4, p3, p2, p1, p0, h, l;
    slong c = E - 1;

    _mp_real_small_log1p_3(&f2, &f1, &f0, (d1 << 1) | (d0 >> (FLINT_BITS - 1)), d0 << 1, 0);

    cc = (c < 0) ? -(ulong) c : (ulong) c;
    umul_ppmm(p1, p0, cc, L[0]);
    umul_ppmm(h, l, cc, L[1]);
    add_ssaaaa(p2, p1, h, l, UWORD(0), p1);
    umul_ppmm(h, l, cc, L[2]);
    add_ssaaaa(p3, p2, h, l, UWORD(0), p2);
    umul_ppmm(h, l, cc, L[3]);
    add_ssaaaa(p4, p3, h, l, UWORD(0), p3);
    (void) p0;

    if (c >= 0)
        add_ssssaaaaaaaa(p4, p3, p2, p1, p4, p3, p2, p1, UWORD(0), f2, f1, f0);
    else
        sub_ddddmmmmssss(p4, p3, p2, p1, p4, p3, p2, p1, UWORD(0), f2, f1, f0);

    return _nfloat_small_set_2(res, p4, p3, p2, p1, 0, MP_REAL_SMALL_LOG_ERR_3 + 2, c < 0, ctx);
}
#endif

static int
_nfloat_log(nfloat_ptr res, nfloat_srcptr x, int base2, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x) || NFLOAT_SGNBIT(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_neg_inf(res, ctx);
        if (NFLOAT_IS_POS_INF(x))
            return nfloat_pos_inf(res, ctx);
        return nfloat_nan(res, ctx);
    }

    if (base2 && NFLOAT_D(x)[NFLOAT_CTX_NLIMBS(ctx) - 1] == (UWORD(1) << (FLINT_BITS - 1))
        && flint_mpn_zero_p(NFLOAT_D(x), NFLOAT_CTX_NLIMBS(ctx) - 1))
        return nfloat_set_si(res, NFLOAT_EXP(x) - 1, ctx);

#if FLINT_BITS == 64
    if (!base2)
    {
        slong E = NFLOAT_EXP(x);

        if (NFLOAT_CTX_NLIMBS(ctx) == 1)
        {
            ulong d = NFLOAT_D(x)[0];
            if ((E != 0 && E != 1) || (E == 1 && (d << 1) >= (UWORD(1) << 32))
                    || (E == 0 && d <= UWORD_MAX - (UWORD(1) << 32)))
                return _nfloat_log_fast_1(res, d, E, ctx);
        }
        else if (NFLOAT_CTX_NLIMBS(ctx) == 2)
        {
            ulong d1 = NFLOAT_D(x)[1];
            if ((E != 0 && E != 1) || (E == 1 && (d1 << 1) >= (UWORD(1) << 25))
                    || (E == 0 && d1 <= UWORD_MAX - (UWORD(1) << 24)))
                return _nfloat_log_fast_2(res, d1, NFLOAT_D(x)[0], E, ctx);
        }
    }
#endif

    return _nfloat_log_mpn(res, NULL, NFLOAT_D(x), NFLOAT_CTX_NLIMBS(ctx), NFLOAT_EXP(x), base2, 0, 0, ctx);
}

int
nfloat_log(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    return _nfloat_log(res, x, 0, ctx);
}

int
nfloat_log2(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    return _nfloat_log(res, x, 1, ctx);
}

int
nfloat_log1p(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    slong n, e;
    int sgnbit;

    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_zero(res, ctx);
        if (NFLOAT_IS_POS_INF(x))
            return nfloat_pos_inf(res, ctx);
        return nfloat_nan(res, ctx);
    }

    n = NFLOAT_CTX_NLIMBS(ctx);
    e = NFLOAT_EXP(x);
    sgnbit = NFLOAT_SGNBIT(x);

    /* log1p(x) = x (1 + delta), |delta| < |x| */
    if (e < -FLINT_BITS * n)
        return _nfloat_set_tiny_rel_err(res, x, sgnbit ? 1 : -1, ctx);

    if (e <= -1)
    {
        /* |x| < 1/2: s = |x| directly */
        ulong y[NFLOAT_ELEM_MAX_LIMBS + 2];
        ulong s[NFLOAT_ELEM_MAX_LIMBS + 1];
        ulong err;
        slong w, ue;
        int trunc;

        w = _nfloat_elem_limbs(ctx, -e + 1);
        trunc = _nfloat_get_fixed(s, w, w, NFLOAT_D(x), n, e);
        ue = _nfloat_log1pm_fixed(y, &err, s, trunc, w, sgnbit);
        return _nfloat_set_mpn_err(res, y, w + 1, ue, err, err, sgnbit, ctx);
    }

    if (sgnbit)
    {
        ulong t[NFLOAT_MAX_LIMBS + 1];

        /* x <= -1 */
        if (e >= 1)
        {
            if (e == 1 && NFLOAT_D(x)[n - 1] == (UWORD(1) << (FLINT_BITS - 1))
                && flint_mpn_zero_p(NFLOAT_D(x), n - 1))
                return nfloat_neg_inf(res, ctx);
            return nfloat_nan(res, ctx);
        }

        /* x in (-1, -1/2]: u = 1 - |x| in (0, 1/2] exactly */
        _nfloat_copy_limbs(t, NFLOAT_D(x), n);
        mpn_neg(t, t, n);
        {
            ulong u[NFLOAT_MAX_ALLOC];
            int status = nfloat_set_mpn_2exp(u, t, n, 0, 0, ctx);
            FLINT_ASSERT(status == GR_SUCCESS);
            (void) status;
            return _nfloat_log_mpn(res, NULL, NFLOAT_D(u), n, NFLOAT_EXP(u), 0, 0, 0, ctx);
        }
    }
    else
    {
        slong w = _nfloat_elem_limbs(ctx, 0);

        /* log1p(x) = log(x) + log1p(1/x) with 0 < log1p(1/x) < 2^(1-e),
           below one unit of the fixed-point logarithm */
        if (e > FLINT_BITS * w + 2)
            return _nfloat_log_mpn(res, NULL, NFLOAT_D(x), n, e, 0, 0, 1, ctx);

        /* x >= 1/2: u = 1 + x exactly, as an integer with fn fraction
           limbs (x's unit is 2^(e - FLINT_BITS n), e >= 0) */
        {
            ulong t[NFLOAT_ELEM_MAX_LIMBS + 3];
            slong s = e - FLINT_BITS * n;
            slong fn = (s >= 0) ? 0 : (-s + FLINT_BITS - 1) / FLINT_BITS;
            slong len = fn + e / FLINT_BITS + 2;
            unsigned int c;

            _nfloat_get_fixed(t, len, fn, NFLOAT_D(x), n, e);
            mpn_add_1(t + fn, t + fn, len - fn, 1);
            while (t[len - 1] == 0)
                len--;
            c = flint_clz(t[len - 1]);
            if (c != 0)
                mpn_lshift(t, t, len, c);
            return _nfloat_log_mpn(res, NULL, t, len, FLINT_BITS * (len - fn) - c, 0, 0, 0, ctx);
        }
    }
}

void
_nfloat_log_fix(_nfloat_fix_t fix, nn_srcptr d, slong dn, slong E, slong extra, gr_ctx_t ctx)
{
    _nfloat_log_mpn(NULL, fix, d, dn, E, 0, extra, 0, ctx);
}
