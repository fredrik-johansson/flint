/*
    Copyright (C) 2025 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ARB_IMPL_H
#define ARB_IMPL_H

#include <stdint.h>
#include "acb_types.h"
#include "mpn_extras.h"
#include "arb.h"

ulong euler_mod_p_powsum_noredc(ulong n, ulong p, const unsigned int * divtab);
ulong euler_mod_p_powsum(ulong n, ulong p, const unsigned int * divtab);
ulong euler_mod_p_powsum_1(ulong n, ulong p);
void divisor_table_odd(unsigned int * tab, slong len);
void arb_gamma_stirling_eval(arb_t s, const arb_t z, slong nterms, int digamma, slong prec);
int _arb_log_ui_smooth(arb_t res, ulong n, slong prec);
void arb_gamma_stirling_coeff(arb_t b, ulong k, int digamma, slong prec);
/* Scratch limbs required by _arb_dot_addmul_generic for an sn-limb
   accumulator: the product of two operands of up to sn + 2 limbs plus one
   limb for shifting (a high product with guard limb is shorter). */
#define ARB_DOT_TMP_LIMBS(sn) (2 * ((sn) + 2) + 1)
/* The guard limbs of the one-pass ball rounding (the lowest limb of a
   high product, quotient or root, whose error is charged as one unit of
   2^(e - prec - MAG_BITS)) are laid out for 64-bit limbs: with 32-bit
   limbs the guarded computations are asked for one more limb. */
#if FLINT_BITS == 64
# define ARB_GUARD_EXTRA 0
#else
# define ARB_GUARD_EXTRA 1
#endif

/* arb_set_round_mpn_2exp_si with err (NULL for none) added to the radius
   and, with extra = 1, one more unit of 2^(exp - prec - MAG_BITS) (which
   covers a lower approximation x whose deficit is below that unit, i.e.
   below 2^(d - MAG_BITS) units of x[0] with d the number of discarded
   bits, at least FLINT_BITS - 1 for a product with a guard limb) */
int _arb_mul_add_rad(arb_t res, const arb_t x, const arb_t y, int neg, slong prec, mag_srcptr zr);
int _arb_set_round_arf_err(arb_t res, arf_srcptr x, int neg, slong prec, mag_srcptr err);
void _arb_set_round_mpn_2exp_si_err(arb_t res, nn_srcptr x, slong xn, int sgnbit, slong exp, slong prec, int extra, mag_srcptr err);
void _arb_set_round_mpn_2exp_si_inplace_err(arb_t res, slong xn, int sgnbit, slong exp, slong prec, int extra, mag_srcptr err);
/* {x, xn} 2^exp +/- eps for an unknown eps > 0 more than FLINT_BITS bits
   below 2^(exp + xn FLINT_BITS - prec) (see set_round_mpn_2exp_si.c) */
void _arb_set_round_mpn_2exp_si_eps(arb_t res, nn_srcptr x, slong xn, int sgnbit, slong exp, int sub, slong prec, mag_srcptr err);

/* The previous implementations of the ball arithmetic (midpoint rounded
   by arf, one ulp added to the radius), the fallbacks of the one-pass
   versions for special values and non-LAGOM exponents. */
void _arb_mul_arf_rounding(arb_t z, const arb_t x, const arb_t y, slong prec);
void _arb_add_rounding(arb_t z, const arb_t x, const arb_t y, int ysub, slong prec);
void _arb_addmul_rounding(arb_t z, const arb_t x, const arb_t y, int neg, slong prec);
/* res = x + (-1)^ysub y, res = z + (-1)^neg x y, in one pass */
void _arb_add_sub(arb_t res, const arb_t x, const arb_t y, int ysub, slong prec);
void _arb_addmul(arb_t res, const arb_t z, const arb_t x, const arb_t y, int neg, slong prec);
/* res = x / y for finite nonzero LAGOM x, y, rounded to nearest at prec
   >= 2 bits with the rounding error as radius, plus err (NULL for none) */
void _arb_div_mid(arb_t res, arf_srcptr xm, arf_srcptr ym, slong prec, mag_srcptr err);
/* err = |xm| yr + |ym| xr + xr yr (+ zr if not NULL) rounded up, for
   LAGOM operands with at least one of xr, yr nonzero */
void _arb_mul_error(mag_t err, arf_srcptr xm, mag_srcptr xr, arf_srcptr ym, mag_srcptr yr, mag_srcptr zr);
void _arb_dot_addmul_generic(nn_ptr sum, nn_ptr serr, nn_ptr tmp, slong sn, nn_srcptr xptr, slong xn, nn_srcptr yptr, slong yn, int negative, flint_bitcnt_t shift);
/* The final steps of the dot products, shared by arb_dot, acb_dot and the
   approx variants (see arb/dot.c): the errors besides the final rounding
   as a mag, the accumulator rounded to nearest with the radius, the
   midpoint alone in the direction rnd (ARF_RND_NEAR: the same midpoint). */
void _arb_dot_round(arb_t res, nn_ptr sum, slong sn, int negative, slong sum_exp, slong prec, uint64_t serr, uint64_t srad, slong srad_exp);
void _arb_dot_round_mid(arf_t res, nn_ptr sum, slong sn, int negative, slong sum_exp, slong prec, arf_rnd_t rnd);
/* the midpoints of arb_dot / acb_dot for len <= 1, for the approx variants */
void _arb_dot_mid_short(arf_t res, arf_srcptr initial, int subtract, arf_srcptr x, arf_srcptr y, slong len, slong prec);
void _acb_dot_mid_short(arf_t rre, arf_t rim, arf_srcptr ire, arf_srcptr iim, int subtract, arf_srcptr xre, arf_srcptr xim, arf_srcptr yre, arf_srcptr yim, slong len, slong prec);
/* the midpoints of arb_dot_simple / acb_dot_simple (the fallbacks) of the
   exact operands; the strides and the real-to-imaginary offsets in arf units */
void _arb_dot_mid_simple(arf_t res, arf_srcptr initial, int subtract, arf_srcptr x, slong xstep, arf_srcptr y, slong ystep, slong len, slong prec);
void _acb_dot_mid_simple(arf_t rre, arf_t rim, arf_srcptr ire, arf_srcptr iim, int subtract, arf_srcptr xre, slong xstep, slong ximoff, arf_srcptr yre, slong ystep, slong yimoff, slong len, slong prec);
void _arb_dot_add_generic(nn_ptr sum, nn_ptr serr, nn_ptr tmp, slong sn, nn_srcptr xptr, slong xn, int negative, flint_bitcnt_t shift);
void mag_borwein_error(mag_t err, slong n);
void arb_exp_arf_huge(arb_t z, const arf_t x, slong mag, slong prec, int minus_one);
slong _arb_compute_bs_exponents(slong * tab, slong n);
slong _arb_get_exp_pos(const slong * tab, slong step);
void mag_set_ui_2exp_small(mag_t z, ulong x, slong e);
void mag_agm(mag_t res, const mag_t x, const mag_t y);
void arb_gamma_stirling_bound(mag_ptr err, const arb_t x, slong k0, slong knum, slong n);
void acb_gamma_stirling_bound(mag_ptr err, const acb_t z, slong k0, slong knum, slong n);
double d_lambertw_branch1(double x);
void arb_exp_taylor_sum_rs_generic(arb_t res, const arb_t x, slong N, slong prec);
void arb_sin_cos_fmpz_div_2exp_bsplit(arb_t wsin, arb_t wcos, const fmpz_t x, flint_bitcnt_t r, slong prec);
void arb_sin_cos_taylor_sum_rs(arb_t s, const arb_t x, slong N, int cosine, slong prec);

/*
    Sets res to (hi, lo, lo2) 2^exp3 (as a 3-limb integer, hi != 0) with
    further nonzero bits below if sticky != 0, rounded to nearest at
    prec <= 2 FLINT_BITS bits with the rounding error as radius, entirely
    in registers (the analogue of arb_set_round_mpn_2exp_si for products
    of operands of at most two limbs). The bits below lo2 only enter the
    sticky part, which is exact as long as hi has at most
    FLINT_BITS - MAG_BITS leading zero bits (then the normalisation cannot
    move them into the MAG_BITS bits of measured error).
*/
/*
    rad = man 2^e (0 <= man < 2^MAG_BITS, a rounding error in units of
    2^e) plus err (NULL for none), rounded up: the sum of two
    MAG_BITS-bit mantissas aligned to the larger exponent. The exponent of
    err must be small (LAGOM).
*/
FLINT_FORCE_INLINE void
_arb_set_rad_man_2exp_add(mag_ptr rad, ulong man, slong e, mag_srcptr err)
{
    slong bc, e1;

    if (man == 0)
    {
        if (err == NULL)
            mag_zero(rad);
        else
            mag_set(rad, err);
        return;
    }

    bc = FLINT_BIT_COUNT(man);
    man = man << (MAG_BITS - bc);
    e1 = e + bc;
    _fmpz_demote(MAG_EXPREF(rad));

    if (err != NULL)
    {
        ulong m2 = MAG_MAN(err);
        slong e2 = MAG_EXP(err), d;

        if (e1 >= e2)
        {
            d = e1 - e2;
            man = man + ((d >= MAG_BITS) ? UWORD(1) : (m2 >> d) + 1);
        }
        else
        {
            d = e2 - e1;
            man = m2 + ((d >= MAG_BITS) ? UWORD(1) : (man >> d) + 1);
            e1 = e2;
        }
        MAG_MAN(rad) = man;
        MAG_EXP(rad) = e1;
        MAG_FAST_ADJUST_ONE_TOO_LARGE(rad);
    }
    else
    {
        MAG_MAN(rad) = man;
        MAG_EXP(rad) = e1;
    }
}

/*
    Sets res to {x, xn} 2^exp (x nonzero, sign given by sgnbit) rounded to
    nearest at prec bits, with the radius set to an upper bound for the
    rounding error, in one pass. With

        |A: prec bits|B: MAG_BITS bits|C: tail|

    the midpoint is A, incremented if the top bit of B is set; the error is
    then (2^MAG_BITS - B) or (B + [C != 0]) in units of 2^-MAG_BITS ulp,
    which is stored directly as the radius (tight to about
    2^(-prec-MAG_BITS) relative to the value). The exponents must be small
    (LAGOM) and prec >= 2.

    With rerr != NULL the radius is left alone and the rounding error (in
    units of 2^(e - prec - MAG_BITS), e the exponent of the midpoint; 0
    when exact) is stored in *rerr instead, for callers combining it with
    other error terms themselves; err is then ignored.
*/
FLINT_FORCE_INLINE void
_arb_set_round_mpn_2exp_si_core(arb_t res, nn_srcptr x, slong xn, int sgnbit, slong exp, slong prec, int in_place, int extra, mag_srcptr err, ulong * rerr)
{
    slong lead, L, d, zn, sb, i, top_bits, pos, carry = 0;
    ulong B, cy, mask;
    int C, up;
    nn_ptr zptr;
    arf_ptr mid = arb_midref(res);

    /* leading zero limbs (fixed-point values with a zero integer part) */
    while (x[xn - 1] == 0)
    {
        xn--;
        if (xn == 0)
        {
            if (rerr != NULL)
            {
                arf_zero(mid);
                *rerr = 0;
            }
            else
                arb_zero(res);
            return;
        }
    }

    /* trailing zero limbs (frames with padding below the value): they
       would count as discarded bits and hide an exact result. Not for a
       lower approximation (extra != 0), which is never exact. */
    if (extra == 0)
    {
        while (x[0] == 0)
        {
            x++;
            xn--;
            exp += FLINT_BITS;
        }
    }

    lead = flint_clz(x[xn - 1]);
    L = xn * FLINT_BITS - lead;      /* bit length */
    d = L - prec;                    /* discarded bits */

    if (d <= 0)
    {
        /* exact: the xn limbs normalised, sized first (the lowest limb
           may shift out entirely; x[0] != 0 after the strip, so no more
           than that one) */
        slong low_out = (lead != 0) && ((x[0] << lead) == 0);
        nn_srcptr xs = x + low_out;
        ulong low = low_out ? (x[0] >> (FLINT_BITS - lead)) : 0, prev, cur;

        FLINT_ASSERT(extra == 0);
        exp = exp + L;
        zn = xn - low_out;

        if (in_place)
            zptr = ARF_PTR_D(mid);
        else
            ARF_GET_MPN_WRITE(zptr, zn, mid);

        /* forwards (low limbs first), safe for zptr <= xs */
        if (lead == 0)
        {
            if (zptr != xs)
                _flint_mpn_copy_small(zptr, xs, zn);
        }
        else
        {
            prev = 0;
            for (i = 0; i < zn; i++)
            {
                cur = xs[i];
                zptr[i] = (cur << lead) | (prev >> (FLINT_BITS - lead));
                prev = cur;
            }
            zptr[0] |= low;
        }

        if (in_place && zn <= ARF_NOPTR_LIMBS)
        {
            ulong t0 = zptr[0], t1 = (zn == 2) ? zptr[1] : 0;
            ARF_DEMOTE(mid);
            ARF_NOPTR_D(mid)[0] = t0;
            ARF_NOPTR_D(mid)[1] = t1;
        }
        ARF_XSIZE(mid) = ARF_MAKE_XSIZE(zn, sgnbit);
        _fmpz_demote(ARF_EXPREF(mid));
        *ARF_EXPREF(mid) = exp;
        if (rerr != NULL)
            *rerr = 0;
        else if (err == NULL)
            mag_zero(arb_radref(res));
        else
            mag_set(arb_radref(res), err);
        return;
    }

    /* B: the top MAG_BITS of the discarded part (left-aligned if there are
       fewer), C: whether anything below is nonzero. Bit positions count
       from bit 0 of x[0]; the discarded part is bits [0, d). */
    top_bits = FLINT_MIN(d, MAG_BITS);
    pos = d - top_bits;              /* lowest bit of B */
    {
        slong q = (ulong) pos / FLINT_BITS, r = (ulong) pos % FLINT_BITS;
        B = x[q] >> r;
        if (r != 0 && q + 1 < xn)
            B |= x[q + 1] << (FLINT_BITS - r);
        B &= ((UWORD(1) << top_bits) - 1);
        B <<= (MAG_BITS - top_bits);

        /* C: bits [0, pos); usually settled by the first limb looked at */
        C = (r != 0) && ((x[q] & ((UWORD(1) << r) - 1)) != 0);
        for (i = q - 1; i >= 0 && !C; i--)
            C = (x[i] != 0);
    }

    up = (B >> (MAG_BITS - 1)) & 1;

    /* an exact tie (the discarded part exactly one half) rounds to even,
       as the arf functions do with ARF_RND_NEAR, so that the midpoint is
       the correctly rounded value; a lower approximation (extra != 0) is
       not exact and keeps rounding up */
    if (up && B == MAG_ONE_HALF && !C && extra == 0)
        up = (x[(ulong) d / FLINT_BITS] >> ((ulong) d % FLINT_BITS)) & 1;

    /* The kept mantissa: the top prec bits of x, top-aligned in zn limbs
       (zn <= xn since prec < L), i.e. the top zn limbs of x 2^lead with
       the low sb bits cleared. */
    zn = (prec + FLINT_BITS - 1) / FLINT_BITS;
    sb = zn * FLINT_BITS - prec;
    exp = exp + L;     /* arf exponent of the value */

    if (in_place)
    {
        /* x is the output's own buffer (with xn limbs allocated): shift in
           place, then move the top zn limbs down */
        nn_ptr xw = (nn_ptr) x;
        ulong borrow = (xn > zn && lead != 0) ? (x[xn - zn - 1] >> (FLINT_BITS - lead)) : 0;

        if (lead != 0)
            mpn_lshift(xw + xn - zn, xw + xn - zn, zn, lead);
        xw[xn - zn] |= borrow;

        zptr = ARF_PTR_D(mid);
        _flint_mpn_copy_small(zptr, xw + xn - zn, zn);
    }
    else
    {
        /* Produce the rounded mantissa directly in the output buffer. */
        ARF_GET_MPN_WRITE(zptr, zn, mid);

        if (lead == 0)
        {
            _flint_mpn_copy_small(zptr, x + xn - zn, zn);
        }
        else
        {
            mpn_lshift(zptr, x + xn - zn, zn, lead);
            if (xn > zn)
                zptr[0] |= x[xn - zn - 1] >> (FLINT_BITS - lead);
        }
    }

    mask = (sb == 0) ? ~UWORD(0) : (~UWORD(0) << sb);
    zptr[0] &= mask;

    if (up)
    {
        cy = (sb == 0) ? 1 : (UWORD(1) << sb);
        for (i = 0; i < zn; i++)
        {
            zptr[i] += cy;
            cy = (zptr[i] < cy);
            if (!cy)
                break;
        }
        if (cy)
        {
            zptr[zn - 1] = UWORD(1) << (FLINT_BITS - 1);
            _flint_mpn_zero_small(zptr, zn - 1);
            exp++;
            carry = 1;
        }
    }

    /* trailing zero limbs (rare): shift down */
    if (zptr[0] == 0)
    {
        i = 1;
        while (zptr[i] == 0)
            i++;
        flint_mpn_copyi(zptr, zptr + i, zn - i);
        zn -= i;
        if (zn <= ARF_NOPTR_LIMBS)
        {
            /* demote to inline limbs */
            ulong t0 = zptr[0], t1 = (zn == 2) ? zptr[1] : 0;
            ARF_DEMOTE(mid);
            ARF_NOPTR_D(mid)[0] = t0;
            ARF_NOPTR_D(mid)[1] = t1;
        }
    }

    ARF_XSIZE(mid) = ARF_MAKE_XSIZE(zn, sgnbit);
    _fmpz_demote(ARF_EXPREF(mid));
    *ARF_EXPREF(mid) = exp;

    /* radius: error in units of 2^(e - prec - MAG_BITS) for the exponent
       e of the value (before a carry of the rounding), plus extra units for
       a lower approximation, plus err; for rerr in units of the exponent
       of the result (halved after a carry) */
    {
        ulong man = (up ? ((UWORD(1) << MAG_BITS) - B) : (B + C)) + extra;
        FLINT_ASSERT(extra == 0 || d >= MAG_BITS + 1);
        if (rerr != NULL)
            *rerr = carry ? (man + 1) >> 1 : man;
        else
            _arb_set_rad_man_2exp_add(arb_radref(res), man, exp - carry - prec - MAG_BITS, err);
    }
}

/*
    The core function for a value held in registers: res = (-1)^sgnbit
    (hi, lo, lo2) 2^exp3 + (unknown bits below if sticky != 0) rounded to
    nearest at prec <= 2 FLINT_BITS bits, hi != 0, with the rounding error
    plus err as radius (or in *rerr, as above).
*/
FLINT_FORCE_INLINE void
_arb_set_round_3limbs(arb_t res, ulong hi, ulong lo, ulong lo2, ulong sticky,
    int sgnbit, slong exp3, slong prec, mag_srcptr err, ulong * rerr)
{
    slong lead, exp, d, carry = 0;
    ulong B, man, cy;
    int C, up;
    nn_ptr zptr;
    arf_ptr mid = arb_midref(res);

    lead = flint_clz(hi);
    /* the sticky bits lie below lo2; after the normalisation they must
       still lie below the MAG_BITS bits of rounding error */
    FLINT_ASSERT(lead <= FLINT_BITS - MAG_BITS || sticky == 0);
    exp = exp3 + 3 * FLINT_BITS - lead;   /* arf exponent of the value */

    /* normalise the top two limbs; the bits shifted out of lo2 at the top
       enter lo, the rest of lo2 stays below */
    if (lead != 0)
    {
        hi = (hi << lead) | (lo >> (FLINT_BITS - lead));
        lo = (lo << lead) | (lo2 >> (FLINT_BITS - lead));
        lo2 = lo2 << lead;
    }

    /* discarded bits of (hi, lo) */
    d = 2 * FLINT_BITS - prec;

    /* B: bits [d - top_bits, d) of (hi, lo), left-aligned to MAG_BITS;
       C: bits below, plus lo2 and sticky. If d < MAG_BITS the remaining
       MAG_BITS - d bits of B come from lo2. */
    if (d >= MAG_BITS)
    {
        slong pos = d - MAG_BITS;
        if (pos >= FLINT_BITS)
            B = hi >> (pos - FLINT_BITS);
        else if (pos == 0)
            B = lo;
        else
            B = (lo >> pos) | (hi << (FLINT_BITS - pos));
        B &= (UWORD(1) << MAG_BITS) - 1;
        if (pos >= FLINT_BITS)
            C = (lo != 0) || ((hi & ((UWORD(1) << (pos - FLINT_BITS)) - 1)) != 0);
        else
            C = (pos != 0) && ((lo & ((UWORD(1) << pos) - 1)) != 0);
        C |= (lo2 != 0) || (sticky != 0);
    }
    else
    {
        /* d < MAG_BITS: the low d bits of lo, then MAG_BITS - d bits of lo2 */
        B = ((lo & ((UWORD(1) << d) - 1)) << (MAG_BITS - d)) | (lo2 >> (FLINT_BITS - (MAG_BITS - d)));
        C = ((lo2 << (MAG_BITS - d)) != 0) || (sticky != 0);
    }

    /* clear the discarded bits */
    if (d >= FLINT_BITS)
    {
        lo = 0;
        hi &= ~UWORD(0) << (d - FLINT_BITS);
    }
    else
    {
        lo &= ~UWORD(0) << d;
    }

    up = B >> (MAG_BITS - 1);

    /* an exact tie rounds to even (as _arb_set_round_mpn_2exp_si_core) */
    if (up && B == MAG_ONE_HALF && !C)
        up = (d >= FLINT_BITS) ? (hi >> (d - FLINT_BITS)) & 1 : (lo >> d) & 1;

    if (up)
    {
        if (d >= FLINT_BITS)
        {
            cy = UWORD(1) << (d - FLINT_BITS);
            hi += cy;
            if (hi < cy)
            {
                hi = UWORD(1) << (FLINT_BITS - 1);
                exp++;
                carry = 1;
            }
        }
        else
        {
            cy = UWORD(1) << d;
            lo += cy;
            hi += (lo < cy);
            if (hi == 0)
            {
                hi = UWORD(1) << (FLINT_BITS - 1);
                lo = 0;
                exp++;
                carry = 1;
            }
        }
    }

    ARF_DEMOTE(mid);
    zptr = ARF_NOPTR_D(mid);
    if (lo == 0)
    {
        ARF_XSIZE(mid) = ARF_MAKE_XSIZE(1, sgnbit);
        zptr[0] = hi;
    }
    else
    {
        ARF_XSIZE(mid) = ARF_MAKE_XSIZE(2, sgnbit);
        zptr[0] = lo;
        zptr[1] = hi;
    }
    fmpz_set_si(ARF_EXPREF(mid), exp);

    /* radius: the rounding error man 2^(exp - prec - MAG_BITS), plus err
       (or the error alone for the caller, as in the core function) */
    man = up ? ((UWORD(1) << MAG_BITS) - B) : (B + C);
    if (rerr != NULL)
        *rerr = carry ? (man + 1) >> 1 : man;
    else
        _arb_set_rad_man_2exp_add(arb_radref(res), man, exp - carry - prec - MAG_BITS, err);
}

/*
    res = (-1)^sgnbit {q, n + 1} 2^exp for n + 1 limbs with the top bit
    set and n >= 1, rounded to nearest at prec <= n FLINT_BITS bits with
    the rounding error as radius, plus extra units of
    2^(e - prec - MAG_BITS) (e the exponent of the value; one unit covers
    an approximation within 2^(FLINT_BITS - MAG_BITS) units of the lowest
    limb, in either direction) and err (NULL for none). The specialisation
    of _arb_set_round_mpn_2exp_si_core to the quotients and roots with a
    guard limb: at least FLINT_BITS bits lie below the rounding point, so
    the error bits are found without scanning. rerr as in the core
    function (the error, including extra, returned instead of the radius).
*/
FLINT_FORCE_INLINE void
_arb_set_round_guarded(arb_t res, nn_srcptr q, slong n, int sgnbit, slong exp, slong prec, int extra, mag_srcptr err, ulong * rerr)
{
    slong zn, sb, E, d, pos, qi, r, i, carry = 0;
    ulong B, C, mask, cy, man;
    int up;
    nn_ptr zptr;
    arf_ptr mid = arb_midref(res);

    FLINT_ASSERT(q[n] >> (FLINT_BITS - 1));
    FLINT_ASSERT(prec >= 2 && prec <= n * FLINT_BITS);

    zn = (prec + FLINT_BITS - 1) / FLINT_BITS;
    sb = zn * FLINT_BITS - prec;
    E = exp + (n + 1) * FLINT_BITS;        /* arf exponent of the value */
    d = (n + 1) * FLINT_BITS - prec;       /* discarded bits, >= FLINT_BITS */
    pos = d - MAG_BITS;                    /* lowest bit of B */

    qi = (ulong) pos / FLINT_BITS;
    r = (ulong) pos % FLINT_BITS;
    B = q[qi] >> r;
    if (r > FLINT_BITS - MAG_BITS)
        B |= q[qi + 1] << (FLINT_BITS - r);
    B &= (UWORD(1) << MAG_BITS) - 1;
    C = (r != 0) ? (q[qi] & ((UWORD(1) << r) - 1)) : 0;
    for (i = 0; i < qi; i++)
        C |= q[i];
    C = (C != 0);
    up = B >> (MAG_BITS - 1);

    /* an exact tie (extra == 0 and the discarded part exactly one half)
       rounds to even, as the core rounding; approximations round up */
    if (up && B == MAG_ONE_HALF && !C && extra == 0)
        up = (q[(ulong) d / FLINT_BITS] >> ((ulong) d % FLINT_BITS)) & 1;

    ARF_GET_MPN_WRITE(zptr, zn, mid);
    _flint_mpn_copy_small(zptr, q + n + 1 - zn, zn);

    mask = (sb == 0) ? ~UWORD(0) : (~UWORD(0) << sb);
    zptr[0] &= mask;

    if (up)
    {
        cy = (sb == 0) ? 1 : (UWORD(1) << sb);
        for (i = 0; i < zn; i++)
        {
            zptr[i] += cy;
            cy = (zptr[i] < cy);
            if (!cy)
                break;
        }
        if (cy)
        {
            zptr[zn - 1] = UWORD(1) << (FLINT_BITS - 1);
            _flint_mpn_zero_small(zptr, zn - 1);
            E++;
            carry = 1;
        }
    }

    /* trailing zero limbs (rare): shift down */
    if (zptr[0] == 0)
    {
        i = 1;
        while (zptr[i] == 0)
            i++;
        flint_mpn_copyi(zptr, zptr + i, zn - i);
        zn -= i;
        if (zn <= ARF_NOPTR_LIMBS)
        {
            ulong t0 = zptr[0], t1 = (zn == 2) ? zptr[1] : 0;
            ARF_DEMOTE(mid);
            ARF_NOPTR_D(mid)[0] = t0;
            ARF_NOPTR_D(mid)[1] = t1;
        }
    }

    ARF_XSIZE(mid) = ARF_MAKE_XSIZE(zn, sgnbit);
    _fmpz_demote(ARF_EXPREF(mid));
    *ARF_EXPREF(mid) = E;

    man = (up ? ((UWORD(1) << MAG_BITS) - B) : (B + C)) + extra;
    if (rerr != NULL)
        *rerr = carry ? (man + 1) >> 1 : man;
    else
        _arb_set_rad_man_2exp_add(arb_radref(res), man, E - carry - prec - MAG_BITS, err);
}

/* the rounding window below is formed on the stack up to this many limbs
   of precision */
#define ARB_ROUND_WINDOW_MAX 64

/*
    _arb_set_round_mpn_2exp_si_core (not in place, exact x) for long inputs:
    when x has at least one limb more than the precision needs, the top
    zn + 1 limbs are normalised into a window and rounded by
    _arb_set_round_guarded, the limbs below entering only as a sticky
    flag (one unit of the error granularity when nonzero, which is what
    the core charges for a nonzero tail as well), instead of scanning x.
*/
FLINT_FORCE_INLINE void
_arb_set_round_mpn_2exp_si_long(arb_t res, nn_srcptr x, slong xn, int sgnbit, slong exp, slong prec, mag_srcptr err, ulong * rerr)
{
    slong zn, top, lead, i;
    ulong q[ARB_ROUND_WINDOW_MAX + 1], sticky;

    /* leading zero limbs (cancellation in a sum) */
    while (xn > 0 && x[xn - 1] == 0)
        xn--;

    /* (prec may be ARF_PREC_EXACT: compare before counting limbs, which
       would overflow) */
    if (prec > ARB_ROUND_WINDOW_MAX * FLINT_BITS)
    {
        _arb_set_round_mpn_2exp_si_core(res, x, xn, sgnbit, exp, prec, 0, 0, err, rerr);
        return;
    }

    zn = (prec + FLINT_BITS - 1) / FLINT_BITS;

    if (xn <= zn)
    {
        _arb_set_round_mpn_2exp_si_core(res, x, xn, sgnbit, exp, prec, 0, 0, err, rerr);
        return;
    }

    top = xn - (zn + 1);    /* the lowest limb of the window */
    lead = flint_clz(x[xn - 1]);

    if (lead == 0)
    {
        _flint_mpn_copy_small(q, x + top, zn + 1);
        sticky = (top > 0) ? x[top - 1] : 0;
    }
    else
    {
        for (i = zn; i >= 1; i--)
            q[i] = (x[top + i] << lead) | (x[top + i - 1] >> (FLINT_BITS - lead));
        q[0] = (x[top] << lead) | ((top > 0) ? (x[top - 1] >> (FLINT_BITS - lead)) : 0);
        sticky = (top > 0) ? (x[top - 1] << lead) : 0;
    }

    for (i = top - 2; i >= 0 && sticky == 0; i--)
        sticky = x[i];

    _arb_set_round_guarded(res, q, zn, sgnbit, exp + top * FLINT_BITS - lead, prec, sticky != 0, err, rerr);
}

/* res = (hi, lo) 2^exp exactly (hi != 0), radius err (NULL for zero), in
   registers */
FLINT_FORCE_INLINE void
_arb_set_exact_2limbs(arb_t res, ulong hi, ulong lo, int sgnbit, slong exp, mag_srcptr err)
{
    arf_ptr mid = arb_midref(res);
    slong lead = flint_clz(hi);
    nn_ptr zptr;

    if (lead != 0)
    {
        hi = (hi << lead) | (lo >> (FLINT_BITS - lead));
        lo = lo << lead;
    }

    ARF_DEMOTE(mid);
    zptr = ARF_NOPTR_D(mid);
    if (lo == 0)
    {
        ARF_XSIZE(mid) = ARF_MAKE_XSIZE(1, sgnbit);
        zptr[0] = hi;
    }
    else
    {
        ARF_XSIZE(mid) = ARF_MAKE_XSIZE(2, sgnbit);
        zptr[0] = lo;
        zptr[1] = hi;
    }
    _fmpz_demote(ARF_EXPREF(mid));
    *ARF_EXPREF(mid) = exp + 2 * FLINT_BITS - lead;

    if (err == NULL)
        mag_zero(arb_radref(res));
    else
        mag_set(arb_radref(res), err);
}

/*
    Small exact sums: X + (-1)^(xs != ys) Y in a frame of F limbs, where
    X = {xp, xn} occupies the top xn limbs of the frame and Y = {yp, yn}
    is placed shift bits further down. Requires xn <= F and
    shift / FLINT_BITS + yn + (shift % FLINT_BITS != 0) <= F, i.e. Y lies
    entirely within the frame. Writes F + 1 limbs to t (the last one the
    carry), returns the number of limbs of the result (0 if it is zero)
    and sets *sgnbit to its sign. The value is {t, zn} 2^(e - F FLINT_BITS)
    where e is the exponent of the top of the frame (that of X).
*/
#define ARB_ADD_FRAME_MAX 8

FLINT_FORCE_INLINE slong
_arb_add_frame(nn_ptr t, slong F, nn_srcptr xp, slong xn, int xs,
    nn_srcptr yp, slong yn, int ys, slong shift, int * sgnbit)
{
    ulong yw[ARB_ADD_FRAME_MAX], cy;
    slong sl = (ulong) shift / FLINT_BITS, sb = (ulong) shift % FLINT_BITS, base, i, zn;

    /* the copies and fills by explicit limb stores (a loop becomes a
       call to memset / memcpy or a rep string instruction, whose start-up
       cost dominates at these lengths) */
    _flint_mpn_zero_small(t, F);
    _flint_mpn_zero_small(yw, F);
    _flint_mpn_copy_small(t + F - xn, xp, xn);

    base = F - sl - yn - (sb != 0);
    if (sb == 0)
    {
        _flint_mpn_copy_small(yw + base, yp, yn);
    }
    else
    {
        ulong lo = 0;
        for (i = 0; i < yn; i++)
        {
            yw[base + i] = lo | (yp[i] << (FLINT_BITS - sb));
            lo = yp[i] >> sb;
        }
        yw[base + yn] = lo;
    }

    if (xs == ys)
    {
        cy = 0;
        for (i = 0; i < F; i++)
        {
            ulong a = t[i], u = a + yw[i], c1 = (u < a), v = u + cy;
            cy = c1 | (v < u);
            t[i] = v;
        }
        t[F] = cy;
        *sgnbit = xs;
        return F + cy;
    }
    else
    {
        cy = 0;
        for (i = 0; i < F; i++)
        {
            ulong a = t[i], b = yw[i], u = a - b, c1 = (a < b), v = u - cy;
            cy = c1 | (u < cy);
            t[i] = v;
        }
        if (cy)
        {
            cy = 0;
            for (i = 0; i < F; i++)
            {
                ulong a = t[i];
                t[i] = -a - cy;
                cy = (a != 0) | cy;
            }
            *sgnbit = ys;
        }
        else
            *sgnbit = xs;
        t[F] = 0;
        zn = F;
        while (zn > 0 && t[zn - 1] == 0)
            zn--;
        return zn;
    }
}

/*
    res = (-1)^sgnbit {t, zn} 2^fexp (zn >= 1, top limb nonzero), with
    further nonzero bits below if sticky != 0, rounded to nearest at
    prec bits with the rounding error as radius, plus err: in registers
    for prec <= 2 FLINT_BITS (the limbs below the top three only
    contribute to the sticky bit), otherwise by the general rounding
    (which also handles an exact result; sticky must then be 0).
*/
FLINT_FORCE_INLINE void
_arb_set_round_frame(arb_t res, nn_srcptr t, slong zn, int sgnbit, slong fexp, slong prec, mag_srcptr err, ulong sticky)
{
    if (prec <= 2 * FLINT_BITS)
    {
        /* the top three limbs of the normalised value, so that the
           MAG_BITS bits of rounding error below the precision are
           exact (a top limb with leading zeros would otherwise push bits
           of the fourth limb into the sticky part) */
        ulong hi = t[zn - 1], lo = (zn >= 2) ? t[zn - 2] : 0, lo2 = (zn >= 3) ? t[zn - 3] : 0, lo3 = (zn >= 4) ? t[zn - 4] : 0;
        slong i, lead = flint_clz(hi);
        for (i = zn - 5; i >= 0; i--)
            sticky |= t[i];
        if (lead != 0)
        {
            hi = (hi << lead) | (lo >> (FLINT_BITS - lead));
            lo = (lo << lead) | (lo2 >> (FLINT_BITS - lead));
            lo2 = (lo2 << lead) | (lo3 >> (FLINT_BITS - lead));
            lo3 = lo3 << lead;
        }
        sticky |= lo3;
        _arb_set_round_3limbs(res, hi, lo, lo2, sticky, sgnbit, fexp + (zn - 3) * FLINT_BITS - lead, prec, err, NULL);
    }
    else
    {
        FLINT_ASSERT(sticky == 0);
        _arb_set_round_mpn_2exp_si_err(res, t, zn, sgnbit, fexp, prec, 0, err);
    }
}

/*
    Small sums restricted to a window: the top W <= ARB_ADD_WINDOW_MAX
    limbs of X + (-1)^(xs != ys) Y, where X = {xp, xn} (xn <= W) occupies
    the top xn limbs of the window and Y = {yp, yn} is placed shift >= 0
    bits further down, extending below the window if shift is large. The
    bits of Y below the window only enter the result as a sticky flag
    (and as a borrow at the bottom of the window when subtracting, the
    remainder being again sticky), which is all the rounding to at most
    2 FLINT_BITS bits needs. Writes W + 1 limbs to t (the last the
    carry), returns the number of limbs (0 if the result is zero; the
    value is then exactly zero, since bits below the window imply
    |Y| < |X|), sets *sgnbit and *sticky. The value is {t, zn} 2^(e - W
    FLINT_BITS) plus the sticky part, e the exponent of the top of the
    window.
*/
#define ARB_ADD_WINDOW_MAX 4

FLINT_FORCE_INLINE slong
_arb_add_window(nn_ptr t, slong W, nn_srcptr xp, slong xn, int xs,
    nn_srcptr yp, slong yn, int ys, slong shift, int * sgnbit, ulong * sticky)
{
    ulong yw[ARB_ADD_WINDOW_MAX], rest = 0, cy;
    slong sl = (ulong) shift / FLINT_BITS, sb = (ulong) shift % FLINT_BITS, i, k, zn;

    _flint_mpn_zero_small(t, W);
    _flint_mpn_zero_small(yw, W);
    _flint_mpn_copy_small(t + W - xn, xp, xn);

    /* limb i of Y (lowest first) lands at window index k, shifted right
       by sb bits, with the shifted-out bits at index k - 1 */
    for (i = 0; i < yn; i++)
    {
        ulong hi, lo;

        k = W - yn + i - sl;
        if (sb == 0)
        {
            hi = yp[i];
            lo = 0;
        }
        else
        {
            hi = yp[i] >> sb;
            lo = yp[i] << (FLINT_BITS - sb);
        }

        if (k >= 0)
            yw[k] |= hi;
        else
            rest |= hi;
        if (k - 1 >= 0)
            yw[k - 1] |= lo;
        else
            rest |= lo;
    }

    *sticky = rest;

    if (xs == ys)
    {
        cy = 0;
        for (i = 0; i < W; i++)
        {
            ulong a = t[i], s = a + yw[i], c1 = (s < a), s2 = s + cy;
            cy = c1 | (s2 < s);
            t[i] = s2;
        }
        t[W] = cy;
        *sgnbit = xs;
        return W + cy;
    }
    else
    {
        /* the part of Y below the window borrows one unit from the
           bottom of the window, leaving a sticky remainder */
        cy = (rest != 0);
        for (i = 0; i < W; i++)
        {
            ulong a = t[i], b = yw[i], d = a - b, c1 = (a < b), d2 = d - cy;
            cy = c1 | (d < cy);
            t[i] = d2;
        }
        t[W] = 0;

        if (cy)
        {
            /* |Y| > |X| (then Y lies within the window: no sticky part) */
            cy = 0;
            for (i = 0; i < W; i++)
            {
                ulong a = t[i], n = -a - cy;
                cy = (a != 0) | cy;
                t[i] = n;
            }
            *sgnbit = ys;
        }
        else
            *sgnbit = xs;

        zn = W;
        while (zn > 0 && t[zn - 1] == 0)
            zn--;
        return zn;
    }
}

#endif
