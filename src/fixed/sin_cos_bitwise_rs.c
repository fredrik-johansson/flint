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

/* sin and cos on [0, 1) by the rotation analog of the bitwise exp
   algorithm, following the identity of [Joh2022]: with
   alpha = 2 atan(b/a),
   exp(i x) = exp(i (x - k alpha)) ((a + b i)/(a - b i))^k, whose
   rotation unit is a ratio of conjugates and hence EXACTLY
   unimodular -- no normalizing square root arises.  With a = 2^i,
   b = 1 the rotation angles are 2 atan(2^-i).

   The TABLE, however, stores the UNDOUBLED angles
   A_i = atan(2^-i) (thread-locally cached like the exp
   logarithms), and phase 1 reduces x/2 rather than x.  This is not
   merely an economy of one table shared with
   fixed_atan_bitwise_rs: the windowed single-limb decision model of
   _fixed_bitwise_reduce requires tab_i < 2^-i, which A_i satisfies
   (as log(1 + 2^-i) does) but 2 atan(2^-i) ~ 2^(1-i) does NOT --
   with doubled entries the residual can keep one bit above the
   window limb after a boundary step, desynchronizing the model.
   Concavity of atan gives A_{i-1} < 2 A_i, so each index is used at
   most once, and sum_{i>=1} A_i ~ 0.8980 > 1/2 covers x/2.

   Writing x/2 = sum_{used} A_i + t', the used factors account for
   the rotation by sum 2 A_i and the residual angle is 2 t'.  The
   reduction therefore runs to index r + 1, so that
   2 t' < 2 A_{r+1} < 2^-r meets the series contract.

   Phase 2 evaluates sin and cos of the residual t < 2^-r: for short
   series both at once by fixed_sin_cos_rs (shared powers of t^2);
   once the series gets long, only the odd sine series, with
   cos t = sqrt(1 - sin^2 t) by one squaring and one integer square
   root -- the same tradeoff and threshold shape as exp's sinh path
   (for small n the square root costs more than the extra terms).

   Phase 3 reconstructs: with W = prod (1 + i 2^-i) over the used
   factors (the real 2^i denominators cancel between W and its
   conjugate), each factor is two shifts and two add/subtracts on
   the old components, and

       (cos x, sin x) = components of (c + i s) W^2 / |W|^2.

   Two squarings give Re W^2 = wx^2 - wy^2 and |W|^2 = wx^2 + wy^2
   from the same pair, one product gives Im W^2 = 2 wx wy, a complex
   multiplication by (c + i s) and two divisions by the real |W|^2
   finish.  The self-normalization is structural: the ratio of the
   SAME truncated W against its own conjugate is exactly unimodular
   up to the final roundings, so the reconstruction truncations
   perturb only the angle of the correction (O(r) ulp), never its
   modulus.  |W| < prod sqrt(1 + 4^-i) = 1.16444, so wx in
   (0.72, 1.17) and wy in [0, 0.92) fit n + 1 limbs, and
   |W|^2 in [1, 1.3559).

   Error accounting (output ulp 2^(-FLINT_BITS n)): reduction table
   and creep < 2 (r + 2) ulp on the angle; series <= 15 ulp per
   component (plus squaring and square root floors on the sqrt
   path); reconstruction truncations < 2 ulp of angle per used
   factor; the final squarings, products and division floors are
   constants.  Total below
   FIXED_SIN_COS_BITWISE_RS_MAX_ERR(n, r) = 6 r + 128. */

FLINT_TLS_PREFIX nn_ptr _fixed_atans = NULL;
FLINT_TLS_PREFIX slong _fixed_atans_n = 0;
FLINT_TLS_PREFIX slong _fixed_atans_r = 0;
static FLINT_TLS_PREFIX int _fixed_atans_cleanup_registered = 0;

static void
_fixed_atans_cleanup(void)
{
    flint_free(_fixed_atans);
    _fixed_atans = NULL;
    _fixed_atans_n = 0;
    _fixed_atans_r = 0;
    _fixed_atans_cleanup_registered = 0;
}

#define ATANS_A(i) (_fixed_atans + (i) * _fixed_atans_n)

void
_fixed_atans_ensure(slong nc, slong rc)
{
    arb_t v, t;
    arf_t lb;
    fmpz_t f;
    slong i, prec;

    if (nc <= _fixed_atans_n && rc <= _fixed_atans_r)
        return;

    nc = FLINT_MAX(nc, _fixed_atans_n);
    rc = FLINT_MAX(rc, _fixed_atans_r);

    flint_free(_fixed_atans);
    _fixed_atans = flint_malloc((rc + 1) * nc * sizeof(ulong));
    _fixed_atans_n = nc;
    _fixed_atans_r = rc;

    prec = FLINT_BITS * nc + 64;

    arb_init(v);
    arb_init(t);
    arf_init(lb);
    fmpz_init(f);

    /* entry 0 is unused (the reductions start at i = 1) but is
       tabulated anyway: A_0 = atan(1) = pi/4 < 1 fits the fraction
       format */
    for (i = 0; i <= rc; i++)
    {
        /* A_i = atan(2^-i), truncated (floor) to nc fraction
           limbs */
        arb_one(v);
        arb_mul_2exp_si(v, v, -i);
        arb_atan(t, v, prec);

        arb_get_lbound_arf(lb, t, prec);
        arf_mul_2exp_si(lb, lb, FLINT_BITS * nc);
        arf_get_fmpz(f, lb, ARF_RND_FLOOR);
        FLINT_ASSERT(fmpz_sgn(f) >= 0);
        fmpz_get_ui_array(ATANS_A(i), nc, f);
    }

    arb_clear(v);
    arb_clear(t);
    arf_clear(lb);
    fmpz_clear(f);

    if (!_fixed_atans_cleanup_registered)
    {
        flint_register_cleanup_function(_fixed_atans_cleanup);
        _fixed_atans_cleanup_registered = 1;
    }
}

/* use the sine series plus a squaring and a square root instead of
   the combined sin_cos series once the series gets long: the
   combined evaluation shares the powers of t^2 and costs about
   1.35x a single series, so the crossover sits later than exp's
   sinh threshold (initial placement; margins are flat) */
#define SINCOS_USE_SQRT(wn, r) \
    (FLINT_BITS * (wn) >= (((r) >= 64) ? 180 : 64) * (slong) (r))

/* Reconstruction shared by the generic code and the small paths:
   (cos x, sin x) = components of (c + i s) W^2 / |W|^2, with
   (wx, wy) = W in wn = n + 1 limbs (wx carries the units limb) and
   (s, cc) = sin and cos of the reduced angle, also wn limbs.  Either
   output may be NULL. */
static void
_fixed_sin_cos_recon(nn_ptr ysin, nn_ptr ycos, nn_srcptr s,
    nn_srcptr cc, nn_srcptr wx, nn_srcptr wy, slong n)
{
    slong wn = n + 1, i;
    nn_ptr A, B, D, N, q;
    TMP_INIT;

    TMP_START;
    A = TMP_ALLOC(((2 * wn) + (2 * wn) + (2 * wn)
        + (2 * (2 * wn + 2)) + (3 * n + 4)) * sizeof(ulong));
    B = A + 2 * wn;
    D = B + 2 * wn;
    N = D + 2 * wn;
    q = N + 2 * (2 * wn + 2);

    /* tail: (cos x, sin x) = components of (c + i s) W^2 / |W|^2.

       Integer views (b = FLINT_BITS): WX = wx 2^(b n) etc., all
       positive.  A = WX^2 - WY^2 = Re W^2 2^(2 b n) and
       D = WX^2 + WY^2 = |W|^2 2^(2 b n) come from the same two
       squarings; B = 2 WX WY = Im W^2 2^(2 b n).  The used angles
       sum below 1 < pi/2, so Re W^2 >= |W|^2 cos 1 > 0.54 |W|^2 and
       everything stays positive, including the cos numerator
       (cos x >= cos 1).  Each 2 wn-limb value has its fraction
       point below limb 2 n (integer parts <= 1.36 in limbs
       2n..2n+1). */
    mpn_sqr(A, wx, wn);
    mpn_sqr(B, wy, wn);
    flint_mpn_copyi(D, A, 2 * wn);
    mpn_add_n(D, D, B, 2 * wn);         /* |W|^2 in [1, 1.3559) */
    mpn_sub_n(A, A, B, 2 * wn);         /* Re W^2 > 0 */
    mpn_mul_n(B, wx, wy, wn);
    mpn_lshift(B, B, 2 * wn, 1);        /* Im W^2 in [0, 1.15) */

    for (i = 0; i < 2; i++)
    {
        /* cos: (c A - s B) / D;  sin: (c B + s A) / D */
        nn_ptr out = (i == 0) ? ycos : ysin;
        nn_srcptr P1 = (i == 0) ? A : B;
        nn_srcptr P2 = (i == 0) ? B : A;
        nn_ptr T1 = N, T2 = N + 2 * wn + 1;

        if (out == NULL)
            continue;

        /* keep the top wn + 1 limbs of P (limbs n..2n+1: value
           P 2^(b n) truncated, the dropped low limbs contribute
           less than one output ulp after the multiplication);
           T = trig * Ptop has (wn + 1) + wn = 2 wn + 1 limbs and
           represents trig * P * 2^(2 b n), point below limb 2 n */
        mpn_mul(T1, P1 + n, wn + 1, cc, wn);
        mpn_mul(T2, P2 + n, wn + 1, s, wn);
        if (i == 0)
            mpn_sub_n(T1, T1, T2, 2 * wn + 1);
        else
            mpn_add_n(T1, T1, T2, 2 * wn + 1);

        /* out = floor(T1 2^(b n) / D'): D' = D's significant limbs
           0..2n (top limb D[2n] >= 1); the numerator gains n zero
           limbs below (3 n + 3 in total), giving an (n + 3)-limb
           quotient with zero top limbs and out = 2^(b n) exactly
           when cos x = 1 */
        flint_mpn_zero(q, n);
        flint_mpn_copyi(q + n, T1, 2 * wn + 1);
        mpn_tdiv_qr(N, T2, 0, q, 3 * n + 3, D, 2 * n + 1);
        flint_mpn_copyi(out, N, n + 1);
    }

    TMP_END;
}

#if FLINT_BITS == 64
#include "sin_cos_bitwise_rs_small.inc"
#endif

/* smallest n at which each reduction parameter becomes optimal
   (generated by src/fixed/tune/tune-bitwise-r.c; x86-64 defaults) */
static const int _fixed_sin_cos_bitwise_rs_r_tab[] =
    {16, 32, 64, 128, 192, 256};
static const short _fixed_sin_cos_bitwise_rs_n_tab[] =
    {1, 5, 11, 38, 87, 128};

void
fixed_sin_cos_bitwise_rs(nn_ptr ysin, nn_ptr ycos, nn_srcptr x,
    slong n, int r)
{
#if FLINT_BITS == 64
    int r0;
#endif
    slong wn, nc, num, j;
    nn_ptr t, s, cc, wx, wy, va, vb, A;
    slong * used;
    TMP_INIT;

    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT(r == 0 || r >= 16);
    /* on 32-bit limbs, n = 1 would clamp r below the series
       contract */
    FLINT_ASSERT(FLINT_BITS == 64 || n >= 2);

#if FLINT_BITS == 64
    r0 = r;
#endif

    if (r == 0)
    {
        for (j = 0; j + 1 < (slong) (sizeof(_fixed_sin_cos_bitwise_rs_n_tab)
                / sizeof(short))
                && n >= _fixed_sin_cos_bitwise_rs_n_tab[j + 1]; j++)
            ;
        r = _fixed_sin_cos_bitwise_rs_r_tab[j];
    }

    r = FLINT_MAX(r, 16);
    r = FLINT_MIN((slong) r, FLINT_BITS * n - 16);

#if FLINT_BITS == 64
    /* r = 0: the specialized sizes run with a compile-time constant r
       and the hand-written series built for it, in which sin and cos
       share their single squaring */
    /* the tangent half-angle reconstruction: one division for
       t = tan(x/2), then sin and cos from it, sharing a single
       squaring.  No |W|^2, no normalization -- see tan_bitwise_rs.c.
       Measured a third to a half faster than the conjugate-ratio tail
       across these sizes. */
    if (r0 == 0 && _fixed_tan_halfangle(ysin, ycos, NULL, x, n))
        return;

    if (n <= 4 && r0 == 0
            && _fixed_sin_cos_bitwise_rs_opt(ysin, ycos, x, n))
        return;

    if (n <= 7)
    {
        _fixed_sin_cos_bitwise_rs_small_tab[n](ysin, ycos, x, r);
        return;
    }
#endif

    wn = n + 1;

    /* the reduction runs on x/2 against the undoubled table and one
       index further (see above), so that doubling the residual
       still meets the 2^-r series contract */
    _fixed_atans_ensure(n, r + 1);
    nc = _fixed_atans_n;

    TMP_START;
    t = TMP_ALLOC((n + wn + wn + wn + wn + n + n + 2 * n)
        * sizeof(ulong));
    s = t + n;              /* sin of the residual, wn limbs */
    cc = s + wn;            /* cos of the residual, wn limbs */
    wx = cc + wn;
    wy = wx + wn;
    va = wy + wn;
    vb = va + n;
    A = vb + n;             /* scratch for the sqrt path (2 n) */
    used = TMP_ALLOC(FIXED_BITWISE_REDUCE_USED_ALLOC(r + 1)
        * sizeof(slong));

    /* phase 1: reduce x/2 (to r + 1, see above); the bit shifted
       out of x costs at most one ulp of the angle */
    mpn_rshift(t, x, n, 1);
    num = _fixed_bitwise_reduce(t, n, r + 1, 1, _fixed_atans, nc,
        used);

    /* the residual angle is twice the reduced remainder */
    mpn_lshift(t, t, n, 1);

    /* phase 2: sin and cos of the residual t < 2^-r */
    if (r < 32)
    {
        /* the wider-range series (hardcoded to n = 6, generic
           beyond); the sqrt path is unavailable here, as fixed_sin_rs
           requires t < 2^-32 */
        _fixed_sin_cos_rs16(s, cc, t, n);
    }
    else if (!SINCOS_USE_SQRT(n, r))
    {
        fixed_sin_cos_rs(s, cc, t, n);
    }
    else
    {
        /* sine only (odd series, half the terms), then
           cos = sqrt(1 - sin^2): the mirror of exp's sinh path;
           for small n the square root costs more than the extra
           series terms, hence the threshold */
        fixed_sin_rs(s, t, n);

        if (flint_mpn_zero_p(s, n))
        {
            /* sin = 0: cos = 1 exactly */
            flint_mpn_zero(cc, n);
            cc[n] = 1;
        }
        else
        {
            /* M = (1 - s^2) 2^(2 b n) = 2^(2 b n) - S^2 with
               S = s 2^(b n) (b = FLINT_BITS): since 0 < S^2 <
               2^(2 b n), this is exactly the 2n-limb negation of
               S^2; cos = floor(sqrt(M)) < 2^(b n) has n limbs */
            mpn_sqr(A, s, n);
            mpn_neg(A, A, 2 * n);
            mpn_sqrtrem(cc, NULL, A, 2 * n);
            cc[n] = 0;
        }
    }

    /* phase 3: W = prod (1 + i 2^-i) over the used factors */
    flint_mpn_zero(wx, n);
    wx[n] = 1;
    flint_mpn_zero(wy, wn);

    for (j = 0; j < num; j++)
    {
        slong ii = used[j];
        slong qq = ii / FLINT_BITS, b = ii - qq * FLINT_BITS;

        /* va = trunc(wx >> ii), vb = trunc(wy >> ii), n limbs */
        flint_mpn_zero(va, n);
        flint_mpn_zero(vb, n);
        if (b != 0)
        {
            mpn_rshift(va, wx + qq, wn - qq, (int) b);
            mpn_rshift(vb, wy + qq, wn - qq, (int) b);
        }
        else
        {
            flint_mpn_copyi(va, wx + qq, wn - qq);
            flint_mpn_copyi(vb, wy + qq, wn - qq);
        }
        /* (wx, wy) <- (wx - wy 2^-ii, wy + wx 2^-ii), old values */
        mpn_sub(wx, wx, wn, vb, n);
        mpn_add(wy, wy, wn, va, n);
    }

    /* phase 3 tail (shared with the small paths) */
    _fixed_sin_cos_recon(ysin, ycos, s, cc, wx, wy, n);

    TMP_END;
}
