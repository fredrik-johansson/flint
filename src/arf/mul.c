/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arf.h"
#include "mpn_extras.h"


/*
    Rounds the normalised 2-limb value (hi, lo) (top bit set), with stop
    the limb below lo and sbot a sticky flag for anything further below,
    to prec <= 2 * FLINT_BITS bits in
    the given strict rounding mode, in registers. Sets the mantissa of z
    and its exponent to ex + ey + expfix (+1 on overflow to a power of
    two). Returns the inexact flag.
*/
FLINT_FORCE_INLINE int
_arf_mul_round_2(arf_ptr z, ulong hi, ulong lo, ulong stop, ulong sbot, int sgnbit,
    slong prec, arf_rnd_t rnd, slong expfix, arf_srcptr x, arf_srcptr y)
{
    slong k = 2 * FLINT_BITS - prec;   /* discarded bits */
    ulong half, rest, kept_lsb, mask, sticky = stop | sbot;
    int inexact, increment;
    nn_ptr zptr;

    if (rnd == ARF_RND_DOWN)
    {
        /* truncation only: no rounding decision, no carry */
        if (k == 0)
        {
            inexact = (sticky != 0);
        }
        else if (k < FLINT_BITS)
        {
            mask = (UWORD(1) << k) - 1;
            inexact = ((lo & mask) | sticky) != 0;
            lo &= ~mask;
        }
        else
        {
            mask = (UWORD(1) << (k - FLINT_BITS)) - 1;
            inexact = ((hi & mask) | lo | sticky) != 0;
            hi &= ~mask;
            lo = 0;
        }

        goto store;
    }

    if (k == 0)
    {
        /* the discarded part is the limb below lo (and further bits) */
        half = stop >> (FLINT_BITS - 1);
        rest = (stop << 1) | sbot;
        kept_lsb = lo & 1;
    }
    else if (k <= FLINT_BITS)
    {
        mask = (k == FLINT_BITS) ? ~UWORD(0) : ((UWORD(1) << k) - 1);
        half = (lo >> (k - 1)) & 1;
        rest = (lo & (mask >> 1)) | sticky;
        lo &= ~mask;
        kept_lsb = (k == FLINT_BITS) ? (hi & 1) : ((lo >> k) & 1);
    }
    else
    {
        slong kk = k - FLINT_BITS;   /* 1 .. FLINT_BITS - 1 */
        mask = (UWORD(1) << kk) - 1;
        half = (hi >> (kk - 1)) & 1;
        rest = (hi & (mask >> 1)) | lo | sticky;
        lo = 0;
        hi &= ~mask;
        kept_lsb = (hi >> kk) & 1;
    }

    inexact = (half | rest) != 0;

    if (rnd == ARF_RND_DOWN)
        increment = 0;
    else if (rnd == ARF_RND_NEAR)
        increment = half & ((rest != 0) | kept_lsb);
    else
        increment = inexact & (rnd == ARF_RND_UP || (rnd == ARF_RND_FLOOR ? sgnbit : !sgnbit));

    if (increment)
    {
        if (k < FLINT_BITS)
        {
            lo += UWORD(1) << k;
            hi += (lo < (UWORD(1) << k));
        }
        else
        {
            hi += UWORD(1) << (k - FLINT_BITS);
        }

        if (hi == 0)
        {
            /* overflow to a power of two */
            hi = UWORD(1) << (FLINT_BITS - 1);
            lo = 0;
            expfix++;
        }
    }

store:
    ARF_DEMOTE(z);
    _fmpz_add2_fast(ARF_EXPREF(z), ARF_EXPREF(x), ARF_EXPREF(y), expfix);
    zptr = ARF_NOPTR_D(z);

    if (lo == 0)
    {
        ARF_XSIZE(z) = ARF_MAKE_XSIZE(1, sgnbit);
        zptr[0] = hi;
    }
    else
    {
        ARF_XSIZE(z) = ARF_MAKE_XSIZE(2, sgnbit);
        zptr[0] = lo;
        zptr[1] = hi;
    }

    return inexact;
}

int
arf_mul(arf_ptr z, arf_srcptr x, arf_srcptr y, slong prec, arf_rnd_t rnd)
{
    slong xn, yn, fix;
    int sgnbit;

    xn = ARF_XSIZE(x);
    yn = ARF_XSIZE(y);
    sgnbit = (xn ^ yn) & 1;
    xn >>= 1;
    yn >>= 1;

    if (yn > xn)
    {
        FLINT_SWAP(arf_srcptr, x, y);
        FLINT_SWAP(slong, xn, yn);
    }

    /* Either operand is a special value. */
    if (yn == 0)
    {
        arf_mul_special(z, x, y);
        return 0;
    }

    if (prec < 1)
        prec = 1;

    /* Round-down (and FAST) products of 1- and 2-limb operands: the
       original register code. */
    if (xn <= 2 && (rnd == ARF_RND_DOWN || rnd == ARF_RND_FAST))
    {
        slong expfix, zn;
        ulong hi, lo;
        int ret;
        nn_ptr zptr;

        if (xn == 1)
        {
            lo = ARF_NOPTR_D(x)[0];
            hi = ARF_NOPTR_D(y)[0];

            umul_ppmm(hi, lo, hi, lo);

            /* Shift so that the top bit is set (branch free version). */
            fix = !(hi >> (FLINT_BITS - 1));
            hi = (hi << fix) | ((lo >> (FLINT_BITS - 1)) & fix);
            lo = (lo << fix);

            ARF_DEMOTE(z);

            if (lo == 0)
            {
                zn = 1;

                if (prec >= FLINT_BITS)
                {
                    lo = hi;
                    ret = 0;
                }
                else
                {
                    lo = MASK_LIMB(hi, FLINT_BITS - prec);
                    ret = (lo != hi);
                }
            }
            else
            {
                zn = 2;

                if (prec <= FLINT_BITS)  /* Must be inexact. */
                {
                    lo = MASK_LIMB(hi, FLINT_BITS - prec);
                    zn = ret = 1;
                }
                else if (prec >= 2 * FLINT_BITS) /* Must be exact. */
                {
                    ret = 0;
                }
                else  /* prec < FLINT_BITS < 2 * FLINT_BITS */
                {
                    ret = MASK_LIMB(lo, 2 * FLINT_BITS - prec) != lo;
                    lo = MASK_LIMB(lo, 2 * FLINT_BITS - prec);

                    if (lo == 0)
                    {
                        zn = 1;
                        lo = hi;
                    }
                }
            }

            _fmpz_add2_fast(ARF_EXPREF(z), ARF_EXPREF(x), ARF_EXPREF(y), -fix);
            ARF_XSIZE(z) = ARF_MAKE_XSIZE(zn, sgnbit);

            zptr = ARF_NOPTR_D(z);
            zptr[0] = lo;
            zptr[1] = hi;

            return ret;
        }
        else if (xn == 2)
        {
            ulong zz[4];
            ulong x1, x0, y1, y0;

            x0 = ARF_NOPTR_D(x)[0];
            x1 = ARF_NOPTR_D(x)[1];

            if (yn == 2)
            {
                y0 = ARF_NOPTR_D(y)[0];
                y1 = ARF_NOPTR_D(y)[1];

                FLINT_MPN_MUL_2X2(zz[3], zz[2], zz[1], zz[0], x1, x0, y1, y0);

                /* Likely case, must be inexact */
                if (prec <= 2 * FLINT_BITS)
                {
                    ARF_DEMOTE(z);

                    fix = !(zz[3] >> (FLINT_BITS - 1));
                    zz[3] = (zz[3] << fix) | ((zz[2] >> (FLINT_BITS - 1)) & fix);
                    zz[2] = (zz[2] << fix) | ((zz[1] >> (FLINT_BITS - 1)) & fix);

                    _fmpz_add2_fast(ARF_EXPREF(z), ARF_EXPREF(x), ARF_EXPREF(y), -fix);

                    /* Rounding */
                    if (prec != 2 * FLINT_BITS)
                    {
                        if (prec > FLINT_BITS)
                        {
                            zz[2] &= (LIMB_ONES << (2 * FLINT_BITS - prec));
                        }
                        else if (prec == FLINT_BITS)
                        {
                            zz[2] = 0;
                        }
                        else
                        {
                            zz[3] &= (LIMB_ONES << (FLINT_BITS - prec));
                            zz[2] = 0;
                        }
                    }

                    if (zz[2] == 0)
                    {
                        ARF_XSIZE(z) = ARF_MAKE_XSIZE(1, sgnbit);
                        ARF_NOPTR_D(z)[0] = zz[3];
                    }
                    else
                    {
                        ARF_XSIZE(z) = ARF_MAKE_XSIZE(2, sgnbit);
                        ARF_NOPTR_D(z)[0] = zz[2];
                        ARF_NOPTR_D(z)[1] = zz[3];
                    }

                    return 1;
                }
            }
            else
            {
                y0 = ARF_NOPTR_D(y)[0];
                FLINT_MPN_MUL_2X1(zz[2], zz[1], zz[0], x1, x0, y0);
            }

            zn = xn + yn;
            ret = _arf_set_round_mpn(z, &expfix, zz, zn, sgnbit, prec, ARF_RND_DOWN);
            _fmpz_add2_fast(ARF_EXPREF(z), ARF_EXPREF(x), ARF_EXPREF(y), expfix);
            return ret;
        }
    }

    if (xn <= 2 && prec <= 2 * FLINT_BITS)
    {
        /* other rounding modes: products of 2, 3 or 4 limbs rounded in
           registers; stop and sbot are the limbs below lo */
        ulong hi, lo, stop, sbot;

        if (rnd == ARF_RND_ACCURATE)
            rnd = ARF_RND_NEAR;

        if (xn == 1)
        {
            umul_ppmm(hi, lo, ARF_NOPTR_D(x)[0], ARF_NOPTR_D(y)[0]);
            stop = sbot = 0;
        }
        else
        {
            ulong zz[4];
            ulong x1 = ARF_NOPTR_D(x)[1], x0 = ARF_NOPTR_D(x)[0];

            if (yn == 2)
            {
                FLINT_MPN_MUL_2X2(zz[3], zz[2], zz[1], zz[0], x1, x0, ARF_NOPTR_D(y)[1], ARF_NOPTR_D(y)[0]);
                hi = zz[3];
                lo = zz[2];
                stop = zz[1];
                sbot = zz[0];
            }
            else
            {
                FLINT_MPN_MUL_2X1(zz[2], zz[1], zz[0], x1, x0, ARF_NOPTR_D(y)[0]);
                hi = zz[2];
                lo = zz[1];
                stop = zz[0];
                sbot = 0;
            }
        }

        /* normalise (branch free); the top bit of the limbs below moves
           into lo, the rest is only needed as a sticky flag */
        fix = !(hi >> (FLINT_BITS - 1));
        hi = (hi << fix) | ((lo >> (FLINT_BITS - 1)) & fix);
        lo = (lo << fix) | ((stop >> (FLINT_BITS - 1)) & fix);
        stop = (stop << fix) | ((sbot >> (FLINT_BITS - 1)) & fix);
        sbot = sbot << fix;

        return _arf_mul_round_2(z, hi, lo, stop, sbot, sgnbit, prec, rnd, -fix, x, y);
    }

    {
        /* outputs of up to 16 limbs from operands at least as long (or one
           of them longer): normalised high product of the (truncated)
           operands with inline certified rounding */
        slong outn = (prec + FLINT_BITS - 1) / FLINT_BITS;

        if (outn <= ARF_MUL_HIGH_SMALL_MAX && (yn >= outn || xn > outn))
        {
            int r = _arf_mul_high_small(z, x, y, outn, sgnbit, prec,
                (rnd >= ARF_RND_FAST) ? ARF_RND_NEAR : rnd);
            if (r != -1)
                return r;
        }
    }

    if (xn + yn <= ARF_MUL_SMALL_LIMBS)
    {
        /* small product: compute in full (a truncated product would only
           pay for a much smaller precision, which is rare and cheap here) */
        slong zn, alloc;
        nn_srcptr xptr, yptr;
        nn_ptr tmp;
        int inexact;
        ARF_MUL_TMP_DECL

        ARF_GET_MPN_READONLY(xptr, xn, x);
        ARF_GET_MPN_READONLY(yptr, yn, y);

        alloc = zn = xn + yn;
        ARF_MUL_TMP_ALLOC(tmp, alloc)

        FLINT_MPN_MUL_WITH_SPECIAL_CASES(tmp, xptr, xn, yptr, yn);

        inexact = _arf_set_round_mpn(z, &fix, tmp, zn, sgnbit, prec,
            (rnd >= ARF_RND_FAST) ? arf_rnd_relaxed_to_strict(rnd) : rnd);
        _fmpz_add2_fast(ARF_EXPREF(z), ARF_EXPREF(x), ARF_EXPREF(y), fix);
        ARF_MUL_TMP_FREE(tmp, alloc)

        return inexact;
    }

    return _arf_mul_mpn(z, x, y, xn, yn, sgnbit, prec, rnd);
}

int
arf_mul_mpz(arf_ptr z, arf_srcptr x, const mpz_t y, slong prec, arf_rnd_t rnd)
{
    slong xn, yn;
    slong fix, shift;
    int sgnbit, inexact;

    yn = FLINT_ABS(y->_mp_size);
    xn = ARF_XSIZE(x);
    xn >>= 1;
    sgnbit = ARF_SGNBIT(x) ^ (y->_mp_size < 0);

    /* Either operand is a special value. */
    if (xn == 0 || yn == 0)
    {
        if (arf_is_finite(x))
        {
            arf_zero(z);
        }
        else
        {
            arf_t t;
            arf_init_set_si(t, mpz_sgn(y));
            arf_mul(z, x, t, prec, rnd);
            arf_clear(t);
        }
        return 0;
    }
    else
    {
        slong zn, alloc;
        nn_srcptr xptr, yptr;
        nn_ptr tmp;
        ARF_MUL_TMP_DECL

        ARF_GET_MPN_READONLY(xptr, xn, x);
        yptr = y->_mp_d;

        alloc = zn = xn + yn;
        ARF_MUL_TMP_ALLOC(tmp, alloc)

        FLINT_MPN_MUL_WITH_SPECIAL_CASES(tmp, xptr, xn, yptr, yn);

        shift = yn * FLINT_BITS - (tmp[zn - 1] == 0) * FLINT_BITS;
        zn -= (tmp[zn - 1] == 0);

        inexact = _arf_set_round_mpn(z, &fix, tmp, zn, sgnbit, prec, rnd);

        _fmpz_add_fast(ARF_EXPREF(z), ARF_EXPREF(x), fix + shift);
        ARF_MUL_TMP_FREE(tmp, alloc)

        return inexact;
    }
}
