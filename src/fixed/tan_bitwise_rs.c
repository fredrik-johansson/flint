/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "longlong.h"
#include "mpn_extras.h"
#include "fixed.h"

#include "tan_rotate.inc"

/* the sizes for which the rotation is emitted in registers */
#define FIXED_TAN_ROTATE_NMAX 6

/* The tangent half-angle reconstruction.

   The angle is reduced exactly as for fixed_sin_cos_bitwise_rs: with the
   cached angles A_i = atan(2^-i),

       x/2 = sum_{i in used} A_i + t',        t' < 2^-r

   and the rotation W = prod (1 + i 2^-i) over the used indices is built
   in the same way, two shifts and two add/subtracts per factor.

   The difference is what is done with it.  Since arg(W) = sum A_i,

       tan(sum A_i) = Im W / Re W = wy / wx,

   and because TAN IS A RATIO, the growth of |W| cancels: no |W|^2, no
   normalizing square root, no unimodular correction -- the whole reason
   for the conjugate-ratio identity used by sin/cos simply evaporates.
   With u = tan(t') the addition formula gives the half-angle tangent in
   ONE division,

       t = tan(x/2) = (wy + wx u) / (wx - wy u),

   whose denominator is close to wx in (0.72, 1.17): no cancellation.
   Everything else follows from t by the half-angle formulas

       sin x = 2t / (1 + t^2),
       cos x = (1 - t^2) / (1 + t^2),
       tan x = 2t / (1 - t^2),

   sharing the single squaring t^2.  Note 0 <= t < tan(1/2) = 0.5464, so
   2t may EXCEED one and the outputs carry a unit limb.  The divisions
   are carried out with numerator and denominator halved, sin = t/D and
   cos = ((1 - t^2)/2)/D with D = (1 + t^2)/2 in [0.5, 0.6492), and
   tan as twice t/(1 - t^2): every divisor is then a NORMALIZED n-limb
   value (top bit set), one limb shorter than 1 + t^2 itself, and every
   division is well conditioned.

   Measured against the conjugate-ratio reconstruction this is a third to
   a half faster across n = 1..8, and it yields sin, cos or tan
   separately rather than computing two functions and discarding one. */

#if FLINT_BITS == 64

/* the reduction parameter each series was built for, re-swept after
   the register rotation and the normalized-divisor reconstruction:
   the tail getting cheaper moved every optimum down, dramatically so
   at n = 4..6 (27, 27, 30 previously) */
static const int _fixed_tan_opt_r[] = { 0, 4, 6, 9, 14, 15, 16 };

static void (* const _fixed_tan_opt_tab[])(nn_ptr, nn_srcptr) = {
    NULL,
    _fixed_tan_rs_opt_1, _fixed_tan_rs_opt_2, _fixed_tan_rs_opt_3,
    _fixed_tan_rs_opt_4, _fixed_tan_rs_opt_5, _fixed_tan_rs_opt_6
};

#define FIXED_TAN_NMAX 6

#else

/* no hand-written tangent series off 64-bit limbs: every size then
   takes tan(t') from the sine and cosine series with one division.
   The register rotation is portable and is used regardless. */
#define FIXED_TAN_NMAX 0

#endif

/* divide (num, nn) by (den, dn), whose top limb may be zero -- the
   denominators here sit near one, so mpn_tdiv_qr's normalization
   requirement has to be met by trimming the length rather than assumed */
static void
_fixed_divide(nn_ptr q, nn_ptr num, slong nn, nn_srcptr den, slong dn,
    nn_ptr scratch)
{
    while (dn > 1 && den[dn - 1] == 0)
        dn--;

    mpn_tdiv_qr(q, scratch, 0, num, nn, den, dn);
}

#if FLINT_BITS == 64

/* n = 1 without a single mpn call past the shared reduction: the
   rotation and series were already in registers, and here the whole
   reconstruction follows -- the two products against wx and wy are
   umul_ppmm chains, and each division collapses to one udiv_qrnnd
   against a divisor rounded to 64 normalized bits (the numerator and
   denominator of t are conditionally halved together when the
   denominator reaches one, dropping its unit limb; the sin/cos and
   tan divisors are normalized by the same halving as in the generic
   tail).  The roundings cost at most ~3 ulp against budgets of
   6r + 128 and 8r + 256.  Measured at roughly a THIRD of the cost of
   the generic path at this size: the mpn dispatch and buffer traffic
   dominated the call. */
static void
_fixed_tan_halfangle_1(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan,
    nn_srcptr x)
{
    slong num, nc;
    int r = _fixed_tan_opt_r[1];
    ulong t[1], wx[2], wy[2];
    ulong U, X0, X1, Y0, ph, pl, qh, ql, n1, n0, T1, T0, D1, D0;
    ulong dh, nh, nl, Q, R, rr, T, s;
    slong * used;
    TMP_INIT;

    _fixed_atans_ensure(1, r);
    nc = _fixed_atans_n;

    TMP_START;
    used = TMP_ALLOC(FIXED_BITWISE_REDUCE_USED_ALLOC(r) * sizeof(slong));

    t[0] = x[0] >> 1;
    num = _fixed_bitwise_reduce(t, 1, r, 1, _fixed_atans, nc, used);
    _fixed_tan_rotate_tab[1](wx, wy, used, num);

    _fixed_tan_rs_opt_1(t, t);          /* u = tan(t') */
    U = t[0];

    X0 = wx[0]; X1 = wx[1]; Y0 = wy[0];

    /* numerator wy + wx u and denominator wx - wy u, each a fraction
       limb plus a unit limb that is 0 or 1 */
    umul_ppmm(ph, pl, X0, U);
    add_ssaaaa(n1, n0, UWORD(0), ph, UWORD(0), X1 ? U : UWORD(0));
    add_ssaaaa(T1, T0, n1, n0, UWORD(0), Y0);
    umul_ppmm(qh, ql, Y0, U);
    sub_ddmmss(D1, D0, X1, X0, UWORD(0), qh);
    (void) pl; (void) ql; (void) T1;

    /* t: the numerator t (wx - wy u) < 0.64 never uses its unit limb;
       when the denominator reaches one, halve both */
    if (D1)
    {
        dh = (UWORD(1) << (FLINT_BITS - 1)) | (D0 >> 1);
        nh = T0 >> 1;
        nl = T0 << (FLINT_BITS - 1);
    }
    else
    {
        dh = D0;
        nh = T0;
        nl = 0;
    }
    udiv_qrnnd(Q, rr, nh, nl, dh);

    T = n_mulhi(Q, Q);                  /* t^2 */

    if (ysin != NULL || ycos != NULL)
    {
        if (T == 0)
        {
            if (ysin != NULL)
            {
                ysin[1] = Q >> (FLINT_BITS - 1);
                ysin[0] = Q << 1;       /* sin = 2t to within one ulp */
            }
            if (ycos != NULL)
            {
                ycos[0] = 0;
                ycos[1] = 1;            /* cos = 1 exactly */
            }
        }
        else
        {
            /* R = 1/(1 + t^2) = (2^127 - 1) / ((1 + t^2)/2) */
            udiv_qrnnd(R, rr, ~UWORD(0) >> 1, ~UWORD(0),
                (UWORD(1) << (FLINT_BITS - 1)) | (T >> 1));
            if (ysin != NULL)
            {
                s = n_mulhi(Q, R);      /* sin = 2 t R */
                ysin[1] = s >> (FLINT_BITS - 1);
                ysin[0] = s << 1;
            }
            if (ycos != NULL)
            {
                ycos[0] = n_mulhi(-T, R);   /* cos = (1 - t^2) R */
                ycos[1] = 0;
            }
        }
    }

    if (ytan != NULL)
    {
        if (T == 0)
        {
            ytan[1] = Q >> (FLINT_BITS - 1);
            ytan[0] = Q << 1;           /* tan = 2t to within one ulp */
        }
        else
        {
            /* tan/2 = t/(1 - t^2), then one doubling bit shift */
            udiv_qrnnd(s, rr, Q, UWORD(0), -T);
            ytan[1] = s >> (FLINT_BITS - 1);
            ytan[0] = s << 1;
        }
    }

    TMP_END;
}

#endif

int
_fixed_tan_halfangle(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan,
    nn_srcptr x, slong n, int r)
{
    slong wn = n + 1, nc, num, i, j;
    nn_ptr t, u, wx, wy, va, vb, T, D, N, Q, sc, ss, cc;
    slong * used;
    TMP_INIT;

    if (n < 1)
        return 0;

#if FLINT_BITS == 64
    if (n == 1 && r == 0)
    {
        _fixed_tan_halfangle_1(ysin, ycos, ytan, x);
        return 1;
    }
#endif

    /* the sizes with a tangent series of their own take the parameter
       it was built for; the rest take the caller's (at least 32, the
       contract of fixed_sin_cos_rs, which supplies tan(t') there) */
#if FLINT_BITS == 64
    if (n <= FIXED_TAN_NMAX && r == 0)
        r = _fixed_tan_opt_r[n];
    else
#endif
    {
        r = FLINT_MAX(r, 32);
        r = FLINT_MIN((slong) r, FLINT_BITS * n - 16);
    }

    _fixed_atans_ensure(n, r);
    nc = _fixed_atans_n;

    TMP_START;
    t = TMP_ALLOC((n + n + wn + wn + n + n + wn + wn
        + (2 * n + 2) + (2 * n + 2) + (2 * n + 4)
        + (n + 2) + (n + 2)) * sizeof(ulong));
    u = t + n;
    wx = u + n;
    wy = wx + wn;
    va = wy + wn;
    vb = va + n;
    T = vb + n;                 /* wn: t^2 / products */
    D = T + wn;                 /* wn: the divisor */
    N = D + wn;                 /* 2n+2: division numerator */
    Q = N + (2 * n + 2);        /* 2n+2: quotient */
    sc = Q + (2 * n + 2);       /* 2n+4: division scratch */
    ss = sc + (2 * n + 4);      /* n+2 each: sin, cos of the residual */
    cc = ss + (n + 2);
    used = TMP_ALLOC(FIXED_BITWISE_REDUCE_USED_ALLOC(r) * sizeof(slong));

    /* reduce x/2; the residual t' is wanted as it stands, so unlike the
       sin/cos reduction there is no doubling and no extra index */
    mpn_rshift(t, x, n, 1);
    num = _fixed_bitwise_reduce(t, n, r, 1, _fixed_atans, nc, used);

    /* W = prod (1 + i 2^-i).  For the small sizes this runs in
       registers (tan_rotate.inc): the generic loop below materializes
       wx >> i and wy >> i into scratch each time round -- zeroing them,
       filling them with a shift or a copy, then an mpn_sub and an
       mpn_add -- none of which is necessary when the components fit in
       registers and the shift is a pair of funnel shifts. */
    if (n <= FIXED_TAN_ROTATE_NMAX)
    {
        _fixed_tan_rotate_tab[n](wx, wy, used, num);
    }
    else
    {
        flint_mpn_zero(wx, n);
        wx[n] = 1;
        flint_mpn_zero(wy, wn);

        for (j = 0; j < num; j++)
        {
            slong ii = used[j], q = ii / FLINT_BITS;
            int b = (int) (ii - q * FLINT_BITS);

            flint_mpn_zero(va + (wn - q), q);
            flint_mpn_zero(vb + (n - q), q);
            if (b != 0)
            {
                mpn_rshift(va, wx + q, wn - q, b);
                if (n - q > 0)
                    mpn_rshift(vb, wy + q, n - q, b);
            }
            else
            {
                flint_mpn_copyi(va, wx + q, wn - q);
                if (n - q > 0)
                    flint_mpn_copyi(vb, wy + q, n - q);
            }
            mpn_sub(wx, wx, wn, vb, n);
            mpn_add(wy, wy, wn, va, n);
        }
    }

    /* u = tan(t').  Where a tangent series exists it is used directly;
       beyond that the sine and cosine series of the residual serve, at
       the price of one division -- the tangent numbers have no tidy
       recurrence, so tabulating them indefinitely is unattractive,
       while this costs a single extra division and reuses code that is
       there anyway. */
#if FLINT_BITS == 64
    if (n <= FIXED_TAN_NMAX)
    {
        _fixed_tan_opt_tab[n](u, t);
    }
    else
#endif
    {
        fixed_sin_cos_rs(ss, cc, t, n);
        flint_mpn_zero(N, n);
        flint_mpn_copyi(N + n, ss, wn);
        _fixed_divide(Q, N, wn + n, cc, wn, sc);
        flint_mpn_copyi(u, Q, n);
    }

    /* t = (wy + wx u) / (wx - wy u).  wx and wy are (n+1)-limb values
       scaled by 2^-64n, u an n-limb fraction, so the product at that
       scale is the top n+1 limbs of the (2n+1)-limb integer product. */
    flint_mpn_mul(N, wx, wn, u, n);
    mpn_add_n(T, wy, N + n, wn);            /* numerator   */
    flint_mpn_mul(N, wy, wn, u, n);
    mpn_sub_n(D, wx, N + n, wn);            /* denominator */

    /* The denominator sits in (0.72, 1.17) and the numerator, being
       t < tan(1/2) = 0.5464 times it, below 0.64: when the denominator
       reaches one, halving both leaves the quotient alone and drops
       the unit limbs, so the divisor is ALWAYS an n-limb value with
       its top bit set -- no trimming, and mpn_tdiv_qr skips its
       internal normalization shifts.  (Costs at most 2 ulp of floor
       against the unhalved division.) */
    if (D[n])
    {
        mpn_rshift(D, D, wn, 1);
        mpn_rshift(T, T, n, 1);
    }
    flint_mpn_zero(N, n);
    flint_mpn_copyi(N + n, T, n);
    mpn_tdiv_qr(Q, sc, 0, N, 2 * n, D, n);
    /* t = (Q, n): below tan(1/2) < 0.5464, so the higher limbs vanish */

    /* the one squaring, shared by all three outputs */
    flint_mpn_sqrhigh(T, Q, n);             /* t^2 < 0.2984 */

    /* The half-angle divisions below keep numerator and denominator
       HALVED: 1 + t^2 in [1, 1.2984) is the worst possible divisor,
       spending a whole limb to hold one bit above the point, while

           D = (1 + t^2)/2 in [0.5, 0.6492)

       is one limb shorter AND has its top bit set, which is what
       mpn_tdiv_qr wants; 1 - t^2 in (0.7015, 1) is normalized as it
       stands.  Halving the numerators alongside (sin = t/D rather
       than 2t/(1 + t^2), tan/2 = t/(1 - t^2) with a final one-bit
       shift) also keeps every 2n-limb-over-n-limb quotient within its
       documented n + 1 limbs -- the old (2n+1)-limb numerators over a
       trimmed n-limb divisor made mpn_tdiv_qr write n + 2 limbs, one
       past the contract. */
    if (ysin != NULL && ycos != NULL)
    {
        if (flint_mpn_zero_p(T, n))
        {
            /* t^2 vanishes at working precision: sin = 2t and cos = 1
               to within one ulp each, no division at all */
            flint_mpn_copyi(ysin, Q, n);
            ysin[n] = 0;
            mpn_lshift(ysin, ysin, n + 1, 1);
            flint_mpn_zero(ycos, n);
            ycos[n] = 1;
        }
        else
        {
            /* Both outputs from ONE division: the reciprocal
               R = 1/(1 + t^2) = 1/(2D) in (0.770, 1), then each output
               is a single mulhigh, sin = 2 t R and cos = (1 - t^2) R.
               Dividing 2^(128n-1) - 1 rather than 2^(128n-1) keeps R
               below 2^(64n) even when D = 2^(64n-1) exactly (t^2 one
               ulp).  Measured at a third to a half the cost of two
               divisions across n = 1..12, at most 2 ulp apart from
               them.  (A precomputed-inverse division, flint_mpn_preinvn
               + divrem_preinvn, was also measured and loses even to the
               two plain divisions at these sizes: the Newton setup only
               amortizes over many divisions by the same divisor.) */
            nn_ptr R = Q + n;   /* n + 1 quotient limbs, past t */

            mpn_rshift(D, T, n, 1);             /* D = (1 + t^2)/2 */
            D[n - 1] |= (UWORD(1) << (FLINT_BITS - 1));
            flint_mpn_store(N, 2 * n, ~UWORD(0));
            N[2 * n - 1] = ~UWORD(0) >> 1;
            mpn_tdiv_qr(R, sc, 0, N, 2 * n, D, n);

            flint_mpn_mulhigh_n(ysin, Q, R, n);
            ysin[n] = 0;
            mpn_lshift(ysin, ysin, n + 1, 1);

            mpn_neg(D, T, n);                   /* 1 - t^2 */
            flint_mpn_mulhigh_n(ycos, D, R, n);
            ycos[n] = 0;
        }
    }
    else if (ysin != NULL || ycos != NULL)
    {
        mpn_rshift(D, T, n, 1);             /* D = (1 + t^2)/2 */
        D[n - 1] |= (UWORD(1) << (FLINT_BITS - 1));

        if (ysin != NULL)
        {
            /* sin = 2t/(1 + t^2) = t/D */
            flint_mpn_zero(N, n);
            flint_mpn_copyi(N + n, Q, n);
            mpn_tdiv_qr(ysin, sc, 0, N, 2 * n, D, n);
        }
        if (ycos != NULL)
        {
            /* cos = (1 - t^2)/(1 + t^2) = ((1 - t^2)/2)/D, with
               (1 - t^2)/2 in (0.3507, 1/2] an n-limb fraction; the
               endpoint is t = 0, where cos comes out exactly 1 */
            flint_mpn_zero(N, n);
            if (flint_mpn_zero_p(T, n))
            {
                flint_mpn_zero(N + n, n);
                N[2 * n - 1] = UWORD(1) << (FLINT_BITS - 1);
            }
            else
            {
                mpn_neg(N + n, T, n);
                mpn_rshift(N + n, N + n, n, 1);
            }
            mpn_tdiv_qr(ycos, sc, 0, N, 2 * n, D, n);
        }
    }

    if (ytan != NULL)
    {
        if (flint_mpn_zero_p(T, n))
        {
            /* t^2 vanishes at working precision, so tan = 2t/(1 - t^2)
               is 2t to within one ulp */
            flint_mpn_copyi(ytan, Q, n);
            ytan[n] = 0;
            mpn_lshift(ytan, ytan, n + 1, 1);
        }
        else
        {
            /* tan/2 = t/(1 - t^2) < 0.78 stays below one, and the
               doubling is a one-bit shift of the (n+1)-limb quotient
               into the unit limb */
            mpn_neg(D, T, n);               /* normalized: > 0.7015 */
            flint_mpn_zero(N, n);
            flint_mpn_copyi(N + n, Q, n);
            mpn_tdiv_qr(ytan, sc, 0, N, 2 * n, D, n);
            mpn_lshift(ytan, ytan, n + 1, 1);
        }
    }

    (void) i;
    TMP_END;
    return 1;
}

void
fixed_tan_bitwise_rs(nn_ptr res, nn_srcptr x, slong n, int r)
{
    FLINT_ASSERT(n >= 1);

    if (_fixed_tan_halfangle(NULL, NULL, res, x, n, r))
        return;

    /* beyond the sizes with a tangent series of their own: sin and cos
       by their own reduction, then one division */
    {
        nn_ptr s, c, N, sc;
        slong i;
        TMP_INIT;

        TMP_START;
        s = TMP_ALLOC(((n + 1) + (n + 1) + (2 * n + 2) + (2 * n + 4))
            * sizeof(ulong));
        c = s + (n + 1);
        N = c + (n + 1);
        sc = N + (2 * n + 2);

        fixed_sin_cos_bitwise_rs(s, c, x, n, r);

        /* tan = sin / cos, bounded by tan(1) < 1.56.  The numerator
           sin < 0.842 never uses its unit limb, so a 2n-limb numerator
           over the trimmed divisor writes at most n + 1 quotient limbs
           (a (2n+1)-limb one would write n + 2, past the contract).
           cos >= cos(1) > 0.54 keeps the divisor at n limbs -- except
           when it comes out exactly 1 at working precision (x tiny),
           where the unit limb survives, the quotient is only n limbs,
           and res[n] is pre-zeroed to cover it. */
        flint_mpn_zero(N, n);
        flint_mpn_copyi(N + n, s, n);

        i = n + 1;
        while (i > 1 && c[i - 1] == 0)
            i--;
        res[n] = 0;
        mpn_tdiv_qr(res, sc, 0, N, 2 * n, c, i);

        TMP_END;
    }
}
