/*
    Copyright (C) 2014, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "mpn_extras/fp_small.h"
#include "arf.h"

/* Roots with at most this many limbs plus guard limb are computed on the
   stack. */
#define ARF_SQRT_STACK_LIMBS 64

/*
    sqrt(x) = S 2^((e + t) / 2) for x = A 2^e, t = e mod 2 and
    S = sqrt(A / 2^t) in [1/2, 1), from flint_mpn_sqrthigh (n limbs and a
    guard limb, within FLINT_MPN_SQRTHIGH_GUARD_ERR_BOUND units) when the
    guard bits certify the rounding (always for a relaxed mode, which
    rounds to nearest), else from the exact truncated root with its
    remainder (flint_mpn_sqrthigh_exact).
*/
int
arf_sqrt(arf_ptr z, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    slong xn, n, fix;
    nn_srcptr xptr;
    ulong q_stack[ARF_SQRT_STACK_LIMBS + 1];
    nn_ptr q;
    int t, inexact, relaxed;
    fmpz_t e;

    if (arf_is_special(x))
    {
        if (arf_is_neg_inf(x))
            arf_nan(z);
        else
            arf_set(z, x);
        return 0;
    }

    if (ARF_SGNBIT(x))
    {
        arf_nan(z);
        return 0;
    }

    if (ARF_IS_POW2(x) && fmpz_is_odd(ARF_EXPREF(x)))
    {
        arf_set(z, x);
        fmpz_cdiv_q_2exp(ARF_EXPREF(z), ARF_EXPREF(z), 1);
        return 0;
    }

    if (prec == ARF_PREC_EXACT)
        flint_throw(FLINT_ERROR, "arf_sqrt: exact precision not supported\n");

    ARF_GET_MPN_READONLY(xptr, xn, x);
    t = fmpz_is_odd(ARF_EXPREF(x));
    relaxed = (rnd >= ARF_RND_FAST);

#if FLINT_BITS == 64
    if (xn == 1 && prec <= FLINT_BITS && !COEFF_IS_MPZ(ARF_EXP(x)))
    {
        /* the exact floor s of sqrt(a B / 2^t) (one limb, top bit set) with
           its remainder r = a B / 2^t - s^2 <= 2 s, in registers: s is
           rounded to prec bits with [r != 0] as the sticky bit, and for
           prec = FLINT_BITS with the half bit [r > s] (sqrt >= s + 1/2
           iff r >= s + 1/4; a tie is impossible) */
        ulong s, r1, r0, low, half, mask;
        slong ee = (ARF_EXP(x) + t) >> 1;
        int up, tie_even;

        _flint_mpn_isqrt_2_rem(&s, &r1, &r0, xptr[0] >> t, (xptr[0] << (FLINT_BITS - 1)) & (-(ulong) t));
        inexact = ((r1 | r0) != 0);

        /* the bits of s below prec, and whether the value is above the
           midpoint of the rounding interval */
        if (prec == FLINT_BITS)
        {
            mask = 0;
            low = 0;
            up = (r1 != 0 || r0 > s);
            tie_even = 0;
        }
        else
        {
            mask = (UWORD(1) << (FLINT_BITS - prec)) - 1;
            half = UWORD(1) << (FLINT_BITS - prec - 1);
            low = s & mask;
            inexact |= (low != 0);
            up = (low > half) || (low == half && (r1 | r0) != 0);
            tie_even = (low == half && (r1 | r0) == 0 && ((s >> (FLINT_BITS - prec)) & 1));
        }

        if (relaxed)
            rnd = ARF_RND_NEAR;
        switch (rnd)
        {
            case ARF_RND_NEAR: up = up || tie_even; break;
            case ARF_RND_UP:
            case ARF_RND_CEIL: up = inexact; break;
            default: up = 0;
        }

        s = (s & ~mask) + (up ? (mask + 1) : 0);
        if (s == 0)
        {
            s = UWORD(1) << (FLINT_BITS - 1);
            ee++;
        }

        /* z = s B^-1 2^ee */
        ARF_DEMOTE(z);
        ARF_NOPTR_D(z)[0] = s;
        ARF_XSIZE(z) = ARF_MAKE_XSIZE(1, 0);
        fmpz_set_si(ARF_EXPREF(z), ee);
        return inexact;
    }
#endif

    n = (prec + FLINT_BITS - 1) / FLINT_BITS;
    q = (n <= ARF_SQRT_STACK_LIMBS) ? q_stack : flint_malloc((n + 1) * sizeof(ulong));

    /* exponent (e + t) / 2 of the root (e + t is even) */
    fmpz_init(e);
    if (!COEFF_IS_MPZ(ARF_EXP(x)))
        *e = (ARF_EXP(x) + t) >> 1;
    else
    {
        fmpz_add_ui(e, ARF_EXPREF(x), t);
        fmpz_tdiv_q_2exp(e, e, 1);
    }

    flint_mpn_sqrthigh(q, xptr, xn, t, n);

    if (relaxed)
    {
        /* rounded to nearest; the inexact flag is certain when the guard
           window would certify a strict rounding, and otherwise decided
           from the input (which is needed intact for that) */
        if (_arf_set_round_approx(z, q, n, FLINT_BIT_COUNT(FLINT_MPN_SQRTHIGH_GUARD_ERR_BOUND), e, 0, prec, ARF_RND_NEAR))
        {
            inexact = 1;
        }
        else if (z == x)
        {
            arf_t u;
            arf_init(u);
            _arf_set_round_approx(u, q, n, FLINT_BIT_COUNT(FLINT_MPN_SQRTHIGH_GUARD_ERR_BOUND), e, 0, prec, rnd);
            inexact = _arf_sqrt_inexact(u, x);
            arf_swap(z, u);
            arf_clear(u);
        }
        else
        {
            _arf_set_round_approx(z, q, n, FLINT_BIT_COUNT(FLINT_MPN_SQRTHIGH_GUARD_ERR_BOUND), e, 0, prec, rnd);
            inexact = _arf_sqrt_inexact(z, x);
        }
    }
    else if (_arf_set_round_approx(z, q, n, FLINT_BIT_COUNT(FLINT_MPN_SQRTHIGH_GUARD_ERR_BOUND), e, 0, prec, rnd))
    {
        inexact = 1;
    }
    else
    {
        /* the floor of the exact root at n + 1 limbs, with a sticky bit for
           the remainder (the guard limb lies entirely below the rounding
           point) */
        inexact = flint_mpn_sqrthigh_exact(q, xptr, xn, t, n);
        q[0] |= inexact;
        inexact = _arf_set_round_mpn(z, &fix, q, n + 1, 0, prec, rnd);
        _fmpz_add_fast(ARF_EXPREF(z), e, fix);
    }

    fmpz_clear(e);
    if (q != q_stack)
        flint_free(q);
    return inexact;
}

int
arf_sqrt_ui(arf_t z, ulong x, slong prec, arf_rnd_t rnd)
{
    arf_t t;
    arf_init_set_ui(t, x); /* no need to free */
    return arf_sqrt(z, t, prec, rnd);
}

int
arf_sqrt_fmpz(arf_t z, const fmpz_t x, slong prec, arf_rnd_t rnd)
{
    arf_t t;
    int inexact;
    arf_init(t);
    arf_set_fmpz(t, x);
    inexact = arf_sqrt(z, t, prec, rnd);
    arf_clear(t);
    return inexact;
}
