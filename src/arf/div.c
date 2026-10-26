/*
    Copyright (C) 2014, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "arf.h"

static void
arf_div_special(arf_t z, const arf_t x, const arf_t y)
{
    if ((arf_is_zero(x) && !arf_is_zero(y) && !arf_is_nan(y)) ||
        (arf_is_inf(y) && !arf_is_special(x)))
    {
        arf_zero(z);
    }
    else if (arf_is_zero(y) || (arf_is_special(x) && arf_is_special(y)) ||
        arf_is_nan(x) || arf_is_nan(y))
    {
        arf_nan(z);
    }
    else if (arf_sgn(x) == arf_sgn(y))
        arf_pos_inf(z);
    else
        arf_neg_inf(z);
}

/*
    Sets z to the n + 1 limbs q (top bit set, a guard limb at the bottom)
    with exponent exp (the value being q B^-(n+1) 2^exp), which lie within
    2^e units of the exact result T, rounded to prec bits: for a relaxed
    mode by rounding to nearest, which gives an error of at most 1/2 ulp
    plus 2^e units of the guard limb (below 2^(e - FLINT_BITS) ulp, since
    prec <= n FLINT_BITS); for a strict mode only when the guard limb
    certifies that T rounds like q: at least FLINT_BITS bits of q lie
    below the rounding point, so the bits e .. FLINT_BITS - 2 of the
    guard limb form a window which, when neither all zero nor all ones,
    guarantees that T is inexact and that no carry or borrow of the
    error reaches the rounding bits (cf. _arf_mul_rounding_certain).
    Returns 1 if z was set and 0 if the check failed (then the exact
    result must be computed); the inexact flag of the rounding is not
    that of T and is not returned.
*/
int
_arf_set_round_approx(arf_ptr z, nn_ptr q, slong n, slong e, const fmpz_t exp, int sgnbit, slong prec, arf_rnd_t rnd)
{
    slong fix;

    if (rnd >= ARF_RND_FAST)
    {
        rnd = ARF_RND_NEAR;
    }
    else
    {
        ulong mask = (UWORD(1) << (FLINT_BITS - 1 - e)) - 1, v = (q[0] >> e) & mask;
        if (v == 0 || v == mask)
            return 0;
    }

    fix = _arf_set_round_small_certified(z, q + 1, n, q[0], sgnbit, prec, rnd);
    _fmpz_add_fast(ARF_EXPREF(z), exp, fix);
    return 1;
}

/* whether q != x / y, for nonzero finite q, x, y: X = Y Q for the odd
   parts implies bits(X) = bits(Y) + bits(Q) - {0, 1} and equal lowest
   limbs, in which case one exact multiplication decides */
static ulong
_arf_odd_low_limb(const arf_t x)
{
    nn_srcptr xp;
    slong xn;
    unsigned int c;

    ARF_GET_MPN_READONLY(xp, xn, x);
    c = flint_ctz(xp[0]);
    if (c == 0)
        return xp[0];
    return (xp[0] >> c) | (xn > 1 ? (xp[1] << (FLINT_BITS - c)) : 0);
}

int
_arf_div_inexact(const arf_t q, const arf_t x, const arf_t y)
{
    slong bq, bx, by;
    arf_t t;
    int inexact;

    bq = arf_bits(q);
    bx = arf_bits(x);
    by = arf_bits(y);
    if (bx < bq + by - 1 || bx > bq + by)
        return 1;
    if (_arf_odd_low_limb(x) != _arf_odd_low_limb(y) * _arf_odd_low_limb(q))
        return 1;

    arf_init(t);
    arf_mul(t, q, y, ARF_PREC_EXACT, ARF_RND_DOWN);
    inexact = !arf_equal(t, x);
    arf_clear(t);
    return inexact;
}

int
_arf_sqrt_inexact(const arf_t q, const arf_t x)
{
    slong bq, bx;
    ulong ql;
    arf_t t;
    int inexact;

    bq = arf_bits(q);
    bx = arf_bits(x);
    if (bx < 2 * bq - 1 || bx > 2 * bq)
        return 1;
    ql = _arf_odd_low_limb(q);
    if (_arf_odd_low_limb(x) != ql * ql)
        return 1;

    arf_init(t);
    arf_mul(t, q, q, ARF_PREC_EXACT, ARF_RND_DOWN);
    inexact = !arf_equal(t, x);
    arf_clear(t);
    return inexact;
}

/* Quotients with at most this many limbs plus guard limb are computed on
   the stack. */
#define ARF_DIV_STACK_LIMBS 64

/* the exact quotient with remainder (flint_mpn_divhigh_exact), rounded */
static int
_arf_div_exact(arf_ptr z, nn_srcptr xptr, slong xn, nn_srcptr yptr, slong yn,
    int sgnbit, const fmpz_t xexp, const fmpz_t yexp, slong n, slong prec, arf_rnd_t rnd)
{
    ulong q_stack[ARF_DIV_STACK_LIMBS + 1], s;
    nn_ptr q;
    slong fix;
    int inexact;

    q = (n <= ARF_DIV_STACK_LIMBS) ? q_stack : flint_malloc((n + 1) * sizeof(ulong));

    /* the floor of the exact quotient at n + 1 limbs, with a sticky bit
       for the remainder (the guard limb lies entirely below the rounding
       point) */
    inexact = flint_mpn_divhigh_exact(q, &s, xptr, xn, yptr, yn, n);
    q[0] |= inexact;

    inexact = _arf_set_round_mpn(z, &fix, q, n + 1, sgnbit, prec, rnd);
    _fmpz_sub2_fast(ARF_EXPREF(z), xexp, yexp, fix + (slong) s);

    if (q != q_stack)
        flint_free(q);
    return inexact;
}

int
arf_div(arf_ptr z, arf_srcptr x, arf_srcptr y, slong prec, arf_rnd_t rnd)
{
    slong xn, yn, n;
    nn_srcptr xptr, yptr;
    ulong q_stack[ARF_DIV_STACK_LIMBS + 1], s;
    nn_ptr q;
    int sgnbit, ok;
    fmpz_t e;

    if (arf_is_special(x) || arf_is_special(y))
    {
        arf_div_special(z, x, y);
        return 0;
    }

    ARF_GET_MPN_READONLY(xptr, xn, x);
    ARF_GET_MPN_READONLY(yptr, yn, y);
    sgnbit = ARF_SGNBIT(x) ^ ARF_SGNBIT(y);

    /* Division by a power of two */
    if (yn == 1 && yptr[0] == LIMB_TOP)
    {
        int inexact;
        fmpz_t t;
        fmpz_init_set(t, ARF_EXPREF(y));

        if (rnd >= ARF_RND_FAST)
            rnd = arf_rnd_relaxed_to_strict(rnd);

        if (ARF_SGNBIT(y))
            inexact = arf_neg_round(z, x, prec, rnd);
        else
            inexact = arf_set_round(z, x, prec, rnd);

        _fmpz_sub2_fast(ARF_EXPREF(z), ARF_EXPREF(z), t, 1);
        fmpz_clear(t);
        return inexact;
    }

    if (prec == ARF_PREC_EXACT)
        flint_throw(FLINT_ERROR, "arf_div: exact precision not supported\n");

    n = (prec + FLINT_BITS - 1) / FLINT_BITS;
    q = (n <= ARF_DIV_STACK_LIMBS) ? q_stack : flint_malloc((n + 1) * sizeof(ulong));

    /* the quotient with a guard limb: x / y = 2^s (q / B^(n+1)) 2^(xexp - yexp) */
    s = flint_mpn_divhigh(q, xptr, xn, yptr, yn, n);

    fmpz_init(e);
    _fmpz_sub2_fast(e, ARF_EXPREF(x), ARF_EXPREF(y), (slong) s);

    if (rnd >= ARF_RND_FAST)
    {
        /* rounded to nearest; the inexact flag is certain when the guard
           window would certify a strict rounding, and otherwise decided
           from the inputs (which are needed intact for that) */
        int inexact;

        if (_arf_set_round_approx(z, q, n, FLINT_BIT_COUNT(FLINT_MPN_DIVHIGH_GUARD_ERR_BOUND), e, sgnbit, prec, ARF_RND_NEAR))
        {
            inexact = 1;
        }
        else if (z == x || z == y)
        {
            arf_t t;
            arf_init(t);
            _arf_set_round_approx(t, q, n, FLINT_BIT_COUNT(FLINT_MPN_DIVHIGH_GUARD_ERR_BOUND), e, sgnbit, prec, rnd);
            inexact = _arf_div_inexact(t, x, y);
            arf_swap(z, t);
            arf_clear(t);
        }
        else
        {
            _arf_set_round_approx(z, q, n, FLINT_BIT_COUNT(FLINT_MPN_DIVHIGH_GUARD_ERR_BOUND), e, sgnbit, prec, rnd);
            inexact = _arf_div_inexact(z, x, y);
        }

        fmpz_clear(e);
        if (q != q_stack)
            flint_free(q);
        return inexact;
    }

    ok = _arf_set_round_approx(z, q, n, FLINT_BIT_COUNT(FLINT_MPN_DIVHIGH_GUARD_ERR_BOUND), e, sgnbit, prec, rnd);
    fmpz_clear(e);
    if (q != q_stack)
        flint_free(q);

    if (ok)
    {
        /* certified inexact (an exact quotient has a zero guard limb) */
        return 1;
    }

    return _arf_div_exact(z, xptr, xn, yptr, yn, sgnbit, ARF_EXPREF(x), ARF_EXPREF(y), n, prec, rnd);
}
