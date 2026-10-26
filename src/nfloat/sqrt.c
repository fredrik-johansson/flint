/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "longlong.h"
#include "mpn_extras.h"
#include "gr.h"
#include "nfloat.h"
#include "impl.h"
#include "elem.h"

/*
    Square root: for x = A 2^e with A = a / B^n in [1/2, 1) and t = e mod 2,
    sqrt(x) = sqrt(A / 2^t) 2^((e + t) / 2) where sqrt(A / 2^t) is in
    [1/2, 1). Its truncated mantissa floor(sqrt(a B^n / 2^t)) is computed
    exactly (flint_mpn_sqrthigh_exact, i.e. flint_mpn_sqrtrem), which also
    tells whether the root is exact, so all rounding modes are handled
    exactly; at one and two limbs the register kernels of
    mpn_extras/fp_small.h give the truncation directly when its guard
    bits certify it.
*/
int
nfloat_sqrt(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    nn_srcptr a;
    slong n, e, exp;
    ulong t;
    int inexact;

    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_NEG_INF(x))
            return nfloat_nan(res, ctx);
        else
            return nfloat_set(res, x, ctx);
    }

    if (NFLOAT_SGNBIT(x))
        return nfloat_nan(res, ctx);

    n = NFLOAT_CTX_NLIMBS(ctx);
    e = NFLOAT_EXP(x);
    t = e & 1;
    a = NFLOAT_D(x);

#if FLINT_BITS == 64
    if (n == 1)
    {
        /* the exact floor of sqrt(a 2^(FLINT_BITS - t)), in registers */
        ulong s1;

        inexact = _nfloat_isqrt_2(&s1, a[0] >> t, (a[0] << (FLINT_BITS - 1)) & (-t));
        exp = (e + (slong) t) / 2;
        NFLOAT_D(res)[0] = s1;
        if (inexact && nfloat_should_round_up(0, ctx))
            NFLOAT_MANT_INCREMENT(NFLOAT_D(res), n, exp);
        NFLOAT_EXP(res) = exp;
        NFLOAT_SGNBIT(res) = 0;
        return GR_SUCCESS;
    }
    else if (n == 2)
    {
        /* sqrt(A 2^-t) at three limbs within NFLOAT_SMALL_SQRT_ERR_3
           units: the truncation to two limbs is certified, and inexact,
           when the guard limb is not within that many units of 0 or B */
        ulong s2, s1, s0;

        _nfloat_small_sqrt_3(&s2, &s1, &s0, a[1], a[0], 0, t);
        if (s0 > NFLOAT_SMALL_SQRT_ERR_3 && s0 < UWORD_MAX - NFLOAT_SMALL_SQRT_ERR_3)
        {
            exp = (e + (slong) t) / 2;
            NFLOAT_D(res)[1] = s2;
            NFLOAT_D(res)[0] = s1;
            if (nfloat_should_round_up(0, ctx))
                NFLOAT_MANT_INCREMENT(NFLOAT_D(res), n, exp);
            NFLOAT_EXP(res) = exp;
            NFLOAT_SGNBIT(res) = 0;
            return GR_SUCCESS;
        }
    }
#endif

    /* the exact truncated root of the mantissa padded to 2n limbs, with
       the exactness flag (flint_mpn_sqrthigh_exact) */
    inexact = flint_mpn_sqrthigh_exact(NFLOAT_D(res), a, n, t, n - 1);
    exp = (e + (slong) t) / 2;

    if (inexact && nfloat_should_round_up(0, ctx))
        NFLOAT_MANT_INCREMENT(NFLOAT_D(res), n, exp);

    NFLOAT_EXP(res) = exp;
    NFLOAT_SGNBIT(res) = 0;
    return GR_SUCCESS;
}

/*
    Reciprocal square root: with x = A 2^e, t = e mod 2 and
    M = (4A)^(-1/2) in (1/2, 2^(-1/2)] (t = 0) or M = (2A)^(-1/2) in
    (2^(-1/2), 1) (t = 1, A != 1/2), we have
    1/sqrt(x) = M 2^((1 - t) - (e - t) / 2).

    The shared flint_mpn_rsqrthigh computes (Y, n + 1) within
    FLINT_MPN_RSQRTHIGH_GUARD_ERR_BOUND units, so that the n-limb
    truncation is certified (and inexact) when the guard limb D = Y[0] is
    more than that bound away from 0 and B - 1; flint_mpn_rsqrthigh_exact
    gives floor(M B^n) and the exactness otherwise. At one and two limbs
    the register kernels of mpn_extras/fp_small.h are used directly.
*/

int
nfloat_rsqrt(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    slong n, e, exp;
    ulong t;
    nn_srcptr a;
    int inexact;

    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_pos_inf(res, ctx);
        else if (NFLOAT_IS_POS_INF(x))
            return nfloat_zero(res, ctx);
        else
            return nfloat_nan(res, ctx);
    }

    if (NFLOAT_SGNBIT(x))
        return nfloat_nan(res, ctx);

    n = NFLOAT_CTX_NLIMBS(ctx);
    e = NFLOAT_EXP(x);
    t = e & 1;
    a = NFLOAT_D(x);

    /* x = 2^(e - 1) with e - 1 even */
    if (t && a[n - 1] == (UWORD(1) << (FLINT_BITS - 1)) && flint_mpn_zero_p(a, n - 1))
    {
        if (res != x)
            _nfloat_copy_limbs(NFLOAT_D(res), a, n);
        NFLOAT_EXP(res) = (1 - e) / 2 + 1;
        NFLOAT_SGNBIT(res) = 0;
        return GR_SUCCESS;
    }

    exp = (1 - (slong) t) - (e - (slong) t) / 2;

    if (n == 1)
    {
#if FLINT_BITS == 64
        ulong m;
        int directed = NFLOAT_CTX_HAS_DIRECTED_ROUNDING(ctx);

        if (_nfloat_rsqrt_1_approx(&m, a[0], t, directed))
        {
            NFLOAT_D(res)[0] = m;
            inexact = 1;
        }
        else
#endif
        {
            inexact = flint_mpn_rsqrthigh_exact(NFLOAT_D(res), a, 1, t, 0);
        }
    }
#if FLINT_BITS == 64
    else if (n == 2)
    {
        /* three limbs in registers within NFLOAT_SMALL_RSQRT_ERR_3 units */
        ulong z2, z1, z0;

        _nfloat_small_rsqrt_3(&z2, &z1, &z0, a[1], a[0], t);

        if (!NFLOAT_CTX_HAS_DIRECTED_ROUNDING(ctx) ||
            (z0 > NFLOAT_SMALL_RSQRT_ERR_3 && z0 < UWORD_MAX - NFLOAT_SMALL_RSQRT_ERR_3))
        {
            NFLOAT_D(res)[1] = z2;
            NFLOAT_D(res)[0] = z1;
            inexact = NFLOAT_CTX_HAS_DIRECTED_ROUNDING(ctx);
        }
        else
        {
            inexact = flint_mpn_rsqrthigh_exact(NFLOAT_D(res), a, n, t, n - 1);
        }
    }
#endif
    else
    {
        ulong Y[NFLOAT_MAX_LIMBS + 1];

        flint_mpn_rsqrthigh(Y, a, n, t, n);

        if (!NFLOAT_CTX_HAS_DIRECTED_ROUNDING(ctx))
        {
            _nfloat_copy_limbs(NFLOAT_D(res), Y + 1, n);
            inexact = 0;
        }
        else if (Y[0] > FLINT_MPN_RSQRTHIGH_GUARD_ERR_BOUND &&
                 Y[0] < UWORD_MAX - FLINT_MPN_RSQRTHIGH_GUARD_ERR_BOUND)
        {
            _nfloat_copy_limbs(NFLOAT_D(res), Y + 1, n);
            inexact = 1;
        }
        else
        {
            inexact = flint_mpn_rsqrthigh_exact(NFLOAT_D(res), a, n, t, n - 1);
        }
    }

    if (inexact && nfloat_should_round_up(0, ctx))
        NFLOAT_MANT_INCREMENT(NFLOAT_D(res), n, exp);

    NFLOAT_EXP(res) = exp;
    NFLOAT_SGNBIT(res) = 0;
    return GR_SUCCESS;
}
