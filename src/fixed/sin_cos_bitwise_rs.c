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

void
fixed_sin_cos_bitwise_rs(nn_ptr ysin, nn_ptr ycos, nn_srcptr x,
    slong n, int r)
{
    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT(r == 0 || r >= 32);

    /* Everything is done by the tangent half-angle reconstruction (see
       tan_bitwise_rs.c), which measured faster at every size tried,
       from n = 1 to 192.  The conjugate-ratio reconstruction that used
       to live here -- with its |W|^2, its two squarings and its two
       wide divisions -- is gone, and with it the doubled-angle
       reduction it needed. */
    _fixed_tan_halfangle(ysin, ycos, NULL, x, n, r);
}
