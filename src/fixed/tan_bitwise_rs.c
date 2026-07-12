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
   2t may EXCEED one and the outputs carry a unit limb; and 1 + t^2 and
   1 - t^2 lie in [1, 1.2984) and (0.7015, 1], so every division here is
   well conditioned.

   Measured against the conjugate-ratio reconstruction this is a third to
   a half faster across n = 1..8, and it yields sin, cos or tan
   separately rather than computing two functions and discarding one. */

#if FLINT_BITS == 64

/* the reduction parameter each series was built for */
static const int _fixed_tan_opt_r[] = { 0, 5, 6, 15, 27, 27, 30 };

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
    mpn_mul(N, wx, wn, u, n);
    mpn_add_n(T, wy, N + n, wn);            /* numerator   */
    mpn_mul(N, wy, wn, u, n);
    mpn_sub_n(D, wx, N + n, wn);            /* denominator */

    flint_mpn_zero(N, n);
    flint_mpn_copyi(N + n, T, wn);
    _fixed_divide(Q, N, wn + n, D, wn, sc);
    /* t = (Q, n): below tan(1/2) < 0.5464, so the higher limbs vanish */

    /* the one squaring, shared by all three outputs */
    flint_mpn_sqrhigh(T, Q, n);             /* t^2 < 0.2984 */

    if (ysin != NULL || ycos != NULL)
    {
        flint_mpn_copyi(D, T, n);
        D[n] = 1;                           /* 1 + t^2 */

        if (ysin != NULL)
        {
            flint_mpn_zero(N, 2 * n + 1);
            flint_mpn_copyi(N + n, Q, n);
            mpn_lshift(N, N, 2 * n + 1, 1); /* 2t, which may exceed one */
            _fixed_divide(ysin, N, 2 * n + 1, D, wn, sc);
        }
        if (ycos != NULL)
        {
            flint_mpn_zero(N, 2 * n + 1);
            if (flint_mpn_zero_p(T, n))
                N[2 * n] = 1;               /* t = 0: cos = 1 exactly */
            else
                mpn_neg(N + n, T, n);       /* 1 - t^2 */
            _fixed_divide(ycos, N, 2 * n + 1, D, wn, sc);
        }
    }

    if (ytan != NULL)
    {
        if (flint_mpn_zero_p(T, n))
        {
            D[n] = 1;                       /* 1 - t^2 = 1 */
            flint_mpn_zero(D, n);
        }
        else
        {
            mpn_neg(D, T, n);               /* 1 - t^2 in (0.7015, 1) */
            D[n] = 0;
        }
        flint_mpn_zero(N, 2 * n + 1);
        flint_mpn_copyi(N + n, Q, n);
        mpn_lshift(N, N, 2 * n + 1, 1);     /* 2t */
        _fixed_divide(ytan, N, 2 * n + 1, D, wn, sc);
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

        /* tan = sin / cos; cos >= cos(1) > 0.54, so the quotient is
           bounded by tan(1) < 1.56 and needs the unit limb */
        flint_mpn_zero(N, n);
        flint_mpn_copyi(N + n, s, n + 1);

        i = n + 1;
        while (i > 1 && c[i - 1] == 0)
            i--;
        mpn_tdiv_qr(res, sc, 0, N, 2 * n + 1, c, i);

        TMP_END;
    }
}
