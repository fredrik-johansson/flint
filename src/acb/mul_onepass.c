/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "arb.h"
#include "acb.h"
#include "arb/impl.h"
#include "acb/impl.h"
#include "arb/terms.h"

/*
    One-pass complex multiplication of balls with LAGOM midpoints and
    radii: (a + b i)(c + d i) = (ac - bd) + (ad + bc) i, each part being
    the sum of two products (and, for acb_addmul, the old midpoint) added
    exactly in a limb frame and rounded to nearest once, the measured
    rounding error going into the radius as in arb_add and arb_mul. The
    products are exact when short and otherwise truncated to the working
    precision plus a guard limb by the normalised high products, the
    deficit being charged to the radius.

    The propagated error of a part is bounded componentwise: for the real
    part |a| cr + |b| dr + |c| ar + |d| br + ar cr + br dr, for the
    imaginary part likewise. The precision of a part is lowered to what
    its error leaves meaningful, as in arb_mul.
*/

/* scratch for the four products (and the accumulator) on the stack */
#define ACB_MUL_TMP_STACK 192

/*
    res = z + (-1)^neg (a + b i)(c + d i) for LAGOM balls with nonzero
    midpoints (z may be NULL: the plain product), in one pass per part.
    Returns 0 if the shape is not handled (the caller then falls back).
*/
int
_acb_mul_onepass(acb_t res, const acb_t z, const acb_t x, const acb_t y, int neg, slong prec)
{
    arf_srcptr a = arb_midref(acb_realref(x)), b = arb_midref(acb_imagref(x));
    arf_srcptr c = arb_midref(acb_realref(y)), d = arb_midref(acb_imagref(y));
    mag_srcptr ar = arb_radref(acb_realref(x)), br = arb_radref(acb_imagref(x));
    mag_srcptr cr = arb_radref(acb_realref(y)), dr = arb_radref(acb_imagref(y));
    mag_t am, bm, cm, dm, er, fr;
    mag_srcptr erp, frp;
    _arb_term_struct T[3], U[3];
    nn_srcptr ap, bp, cp, dp;
    slong nT, nU, outn, an, bn, cn, dn, maxn, slot, alloc;
    slong ea, eb, ec, ed;
    ulong tmp_stack[ACB_MUL_TMP_STACK];
    nn_ptr tmp;
    int sqr, inexact_input, sa, sb, sc, sd;

    /* exponents within +/- 2^(FLINT_BITS - 6), well inside the LAGOM range
       (so the infinities and NaN, whose exponents lie far outside, are
       excluded), in one combined check: the OR of the offset exponents
       is below 2^(FLINT_BITS - 5) exactly when all of them are; nonzero
       midpoints */
    {
        const ulong M = UWORD(1) << (FLINT_BITS - 6);
        ulong range = (ulong) (ARF_EXP(a) + M) | (ulong) (ARF_EXP(b) + M)
                    | (ulong) (ARF_EXP(c) + M) | (ulong) (ARF_EXP(d) + M)
                    | (ulong) (MAG_EXP(ar) + M) | (ulong) (MAG_EXP(br) + M)
                    | (ulong) (MAG_EXP(cr) + M) | (ulong) (MAG_EXP(dr) + M);

        if (z != NULL)
            range |= (ulong) (ARF_EXP(arb_midref(acb_realref(z))) + M) | (ulong) (ARF_EXP(arb_midref(acb_imagref(z))) + M)
                   | (ulong) (MAG_EXP(arb_radref(acb_realref(z))) + M) | (ulong) (MAG_EXP(arb_radref(acb_imagref(z))) + M);

        if (range >= 2 * M || ARF_IS_SPECIAL(a) || ARF_IS_SPECIAL(b) || ARF_IS_SPECIAL(c) || ARF_IS_SPECIAL(d)
            || prec < 2 || prec == ARF_PREC_EXACT)
            return 0;
    }

    sqr = (x == y);
    inexact_input = (MAG_MAN(ar) | MAG_MAN(br) | MAG_MAN(cr) | MAG_MAN(dr)) != 0;

    /* propagated errors: for the real part |a| cr + |b| dr + |c| ar +
       |d| br + ar cr + br dr, for the imaginary part likewise, each as one
       aligned sum of the six products */
    if (inexact_input)
    {
        mag_t zero;
        MAG_MAN(zero) = 0;
        MAG_EXP(zero) = 0;

        mag_fast_init_set_arf(am, a);
        mag_fast_init_set_arf(bm, b);

        if (sqr)
        {
            /* 2|a| ar + ar^2 + 2|b| br + br^2 and 2 (|a| br + |b| ar + ar br) */
            *cm = *am; MAG_EXP(cm)++;
            *dm = *bm; MAG_EXP(dm)++;
            _arb_mag_sum6(er, cm, ar, ar, ar, dm, br, br, br, zero, zero, zero, zero);
            _arb_mag_sum6(fr, am, br, bm, ar, ar, br, zero, zero, zero, zero, zero, zero);
            MAG_EXP(fr)++;
        }
        else
        {
            mag_fast_init_set_arf(cm, c);
            mag_fast_init_set_arf(dm, d);
            _arb_mag_sum6(er, am, cr, bm, dr, cm, ar, dm, br, ar, cr, br, dr);
            _arb_mag_sum6(fr, am, dr, bm, cr, cm, br, dm, ar, ar, dr, br, cr);
        }

        if (z != NULL)
        {
            mag_fast_add(er, er, arb_radref(acb_realref(z)));
            mag_fast_add(fr, fr, arb_radref(acb_imagref(z)));
        }

        erp = (MAG_MAN(er) != 0) ? er : NULL;
        frp = (MAG_MAN(fr) != 0) ? fr : NULL;
    }
    else if (z != NULL)
    {
        /* copies: z may be the output */
        *er = *arb_radref(acb_realref(z));
        *fr = *arb_radref(acb_imagref(z));
        erp = (MAG_MAN(er) != 0) ? er : NULL;
        frp = (MAG_MAN(fr) != 0) ? fr : NULL;
    }
    else
        erp = frp = NULL;

    ARF_GET_MPN_READONLY(ap, an, a);
    ARF_GET_MPN_READONLY(bp, bn, b);
    ARF_GET_MPN_READONLY(cp, cn, c);
    ARF_GET_MPN_READONLY(dp, dn, d);
    ea = ARF_EXP(a);
    eb = ARF_EXP(b);
    ec = ARF_EXP(c);
    ed = ARF_EXP(d);
    /* signs of the four products ac, -bd, ad, bc (with the overall sign) */
    sa = ARF_SGNBIT(a);
    sb = ARF_SGNBIT(b);
    sc = ARF_SGNBIT(c) ^ neg;
    sd = ARF_SGNBIT(d) ^ neg;

#if FLINT_BITS == 64
    /* small operands and precision: exact products, sums in registers
       (the frames and the error windows are laid out for 64-bit limbs) */
    if (z == NULL && prec <= 2 * FLINT_BITS && ((an - 1) | (bn - 1) | (cn - 1) | (dn - 1)) == 0)
    {
        /* all products before the output (which may alias the inputs)
           is written */
        ulong h0, l0, h1, l1, h2, l2, h3, l3;

        if (sqr)
        {
            umul_ppmm(h0, l0, ap[0], ap[0]);
            umul_ppmm(h1, l1, bp[0], bp[0]);
            umul_ppmm(h2, l2, ap[0], bp[0]);
            _arb_round_sum_1limb(acb_realref(res), h0, l0, 2 * ea, neg, h1, l1, 2 * eb, neg ^ 1, prec, erp);
            if (h2 >> (FLINT_BITS - 1))
                _arb_set_round_3limbs(acb_imagref(res), h2, l2, 0, 0, sa ^ sb ^ neg, ea + eb + 1 - 3 * FLINT_BITS, prec, frp, NULL);
            else
                _arb_set_round_3limbs(acb_imagref(res), (h2 << 1) | (l2 >> (FLINT_BITS - 1)), l2 << 1, 0, 0, sa ^ sb ^ neg, ea + eb - 3 * FLINT_BITS, prec, frp, NULL);
        }
        else
        {
            umul_ppmm(h0, l0, ap[0], cp[0]);
            umul_ppmm(h1, l1, bp[0], dp[0]);
            umul_ppmm(h2, l2, ap[0], dp[0]);
            umul_ppmm(h3, l3, bp[0], cp[0]);
            _arb_round_sum_1limb(acb_realref(res), h0, l0, ea + ec, sa ^ sc, h1, l1, eb + ed, sb ^ sd ^ 1, prec, erp);
            _arb_round_sum_1limb(acb_imagref(res), h2, l2, ea + ed, sa ^ sd, h3, l3, eb + ec, sb ^ sc, prec, frp);
        }
        return 1;
    }

    if (prec <= 2 * FLINT_BITS && ((an - 1) | (bn - 1) | (cn - 1) | (dn - 1)) <= 1
        && (z == NULL || (ARF_SIZE(arb_midref(acb_realref(z))) <= 2 && ARF_SIZE(arb_midref(acb_imagref(z))) <= 2)))
    {
        /* operands of at most two limbs, with a short accumulator: the
           terms of a part (exact products of up to four limbs) added in a
           four-limb frame in registers */
        ulong pt[4][4], zt[2][2];
        arf_srcptr zm;
        int short_terms = ((an - 1) | (bn - 1) | (cn - 1) | (dn - 1)) == 0;

        if (sqr)
        {
            _arb_term_exact(T + 0, pt[0], ap, an, ap, an, 2 * ea, neg);
            _arb_term_exact(T + 1, pt[1], bp, bn, bp, bn, 2 * eb, neg ^ 1);
            _arb_term_exact(U + 0, pt[2], ap, an, bp, bn, ea + eb + 1, sa ^ sb ^ neg);
            nT = 2;
            nU = 1;
        }
        else
        {
            _arb_term_exact(T + 0, pt[0], ap, an, cp, cn, ea + ec, sa ^ sc);
            _arb_term_exact(T + 1, pt[1], bp, bn, dp, dn, eb + ed, sb ^ sd ^ 1);
            _arb_term_exact(U + 0, pt[2], ap, an, dp, dn, ea + ed, sa ^ sd);
            _arb_term_exact(U + 1, pt[3], bp, bn, cp, cn, eb + ec, sb ^ sc);
            nT = 2;
            nU = 2;
        }

        if (z != NULL)
        {
            /* the midpoints of z (read before the output is written) */
            zm = arb_midref(acb_realref(z));
            if (!arf_is_zero(zm))
            {
                nn_srcptr zp;
                slong zn;
                ARF_GET_MPN_READONLY(zp, zn, zm);
                zt[0][0] = zp[0];
                zt[0][1] = zp[zn - 1];
                T[nT].p = zt[0]; T[nT].n = zn; T[nT].exp = ARF_EXP(zm); T[nT].e = -1; T[nT].sgn = ARF_SGNBIT(zm);
                nT++;
                short_terms &= (zn == 1);
            }

            zm = arb_midref(acb_imagref(z));
            if (!arf_is_zero(zm))
            {
                nn_srcptr zp;
                slong zn;
                ARF_GET_MPN_READONLY(zp, zn, zm);
                zt[1][0] = zp[0];
                zt[1][1] = zp[zn - 1];
                U[nU].p = zt[1]; U[nU].n = zn; U[nU].exp = ARF_EXP(zm); U[nU].e = -1; U[nU].sgn = ARF_SGNBIT(zm);
                nU++;
                short_terms &= (zn == 1);
            }
        }

        if (short_terms)
        {
            _arb_round_acc4(acb_realref(res), T, nT, prec, erp, 1);
            if (nU == 1)
                _arb_set_round_frame(acb_imagref(res), U[0].p, U[0].n, U[0].sgn, U[0].exp - U[0].n * FLINT_BITS, prec, frp, 0);
            else
                _arb_round_acc4(acb_imagref(res), U, nU, prec, frp, 1);
        }
        else
        {
            _arb_round_acc4(acb_realref(res), T, nT, prec, erp, 0);
            if (nU == 1)
                _arb_set_round_frame(acb_imagref(res), U[0].p, U[0].n, U[0].sgn, U[0].exp - U[0].n * FLINT_BITS, prec, frp, 0);
            else
                _arb_round_acc4(acb_imagref(res), U, nU, prec, frp, 0);
        }
        return 1;
    }
#endif

    /* the products, in slots of one size: a full product of two operands
       or the high product with its guard limb */
    outn = (prec + FLINT_BITS - 1) / FLINT_BITS + 1;
    maxn = FLINT_MAX(FLINT_MAX(an, bn), FLINT_MAX(cn, dn));
    slot = (outn >= 2 * maxn) ? 2 * maxn : 2 * maxn + 2;
    alloc = (sqr ? 3 : 4) * slot;
    if (z != NULL)
        alloc += ARF_SIZE(arb_midref(acb_realref(z))) + ARF_SIZE(arb_midref(acb_imagref(z)));
    tmp = (alloc <= ACB_MUL_TMP_STACK) ? tmp_stack
        : (alloc <= ARF_MUL_TLS_ALLOC) ? _arf_mul_tls_tmp(alloc) : flint_malloc(sizeof(ulong) * alloc);

    if (sqr)
    {
        /* real part a^2 - b^2, imaginary part 2ab */
        int full = _arb_product_full(an, bn, outn) && _arb_product_full(an, an, outn) && _arb_product_full(bn, bn, outn);

        if (full || _arb_products_cancel(ap + an - 1, ap + an - 1, 2 * ea, bp + bn - 1, bp + bn - 1, 2 * eb))
        {
            _arb_term_exact(T + 0, tmp, ap, an, ap, an, 2 * ea, neg);
            _arb_term_exact(T + 1, tmp + slot, bp, bn, bp, bn, 2 * eb, neg ^ 1);
        }
        else
        {
            _arb_term_high(T + 0, tmp, ap, an, ap, an, outn, 2 * ea, neg);
            _arb_term_high(T + 1, tmp + slot, bp, bn, bp, bn, outn, 2 * eb, neg ^ 1);
        }

        if (full)
            _arb_term_exact(U + 0, tmp + 2 * slot, ap, an, bp, bn, ea + eb + 1, sa ^ sb ^ neg);
        else
            _arb_term_high(U + 0, tmp + 2 * slot, ap, an, bp, bn, outn, ea + eb + 1, sa ^ sb ^ neg);
        nT = 2;
        nU = 1;
    }
    else
    {
        int full;

        if (an == bn && cn == dn)
            full = _arb_product_full(an, cn, outn);
        else
            full = _arb_product_full(an, cn, outn) && _arb_product_full(bn, dn, outn)
                && _arb_product_full(an, dn, outn) && _arb_product_full(bn, cn, outn);

        /* exact products for a part whose products nearly cancel */
        if (full || (((sa ^ sc) == (sb ^ sd)) && _arb_products_cancel(ap + an - 1, cp + cn - 1, ea + ec, bp + bn - 1, dp + dn - 1, eb + ed)))
        {
            _arb_term_exact(T + 0, tmp, ap, an, cp, cn, ea + ec, sa ^ sc);
            _arb_term_exact(T + 1, tmp + slot, bp, bn, dp, dn, eb + ed, sb ^ sd ^ 1);
        }
        else
        {
            _arb_term_high(T + 0, tmp, ap, an, cp, cn, outn, ea + ec, sa ^ sc);
            _arb_term_high(T + 1, tmp + slot, bp, bn, dp, dn, outn, eb + ed, sb ^ sd ^ 1);
        }

        if (full || (((sa ^ sd) != (sb ^ sc)) && _arb_products_cancel(ap + an - 1, dp + dn - 1, ea + ed, bp + bn - 1, cp + cn - 1, eb + ec)))
        {
            _arb_term_exact(U + 0, tmp + 2 * slot, ap, an, dp, dn, ea + ed, sa ^ sd);
            _arb_term_exact(U + 1, tmp + 3 * slot, bp, bn, cp, cn, eb + ec, sb ^ sc);
        }
        else
        {
            _arb_term_high(U + 0, tmp + 2 * slot, ap, an, dp, dn, outn, ea + ed, sa ^ sd);
            _arb_term_high(U + 1, tmp + 3 * slot, bp, bn, cp, cn, outn, eb + ec, sb ^ sc);
        }

        nT = 2;
        nU = 2;
    }

    if (z != NULL)
    {
        /* the midpoints of z as exact terms, copied since z may be the
           output (a zero part contributes nothing; its radius is already
           in the error) */
        nn_ptr zt = tmp + (sqr ? 3 : 4) * slot;
        arf_srcptr zm;
        nn_srcptr zp;
        slong zn;

        zm = arb_midref(acb_realref(z));
        if (!arf_is_zero(zm))
        {
            ARF_GET_MPN_READONLY(zp, zn, zm);
            _flint_mpn_copy_small(zt, zp, zn);
            T[nT].p = zt; T[nT].n = zn; T[nT].exp = ARF_EXP(zm); T[nT].e = -1; T[nT].sgn = ARF_SGNBIT(zm);
            nT++;
            zt += zn;
        }

        zm = arb_midref(acb_imagref(z));
        if (!arf_is_zero(zm))
        {
            ARF_GET_MPN_READONLY(zp, zn, zm);
            _flint_mpn_copy_small(zt, zp, zn);
            U[nU].p = zt; U[nU].n = zn; U[nU].exp = ARF_EXP(zm); U[nU].e = -1; U[nU].sgn = ARF_SGNBIT(zm);
            nU++;
        }
    }

    /* all terms are formed (the inputs may alias the output) */
    if (nT == 2)
        _arb_round_sum2(acb_realref(res), T, T + 1, prec, erp);
    else
        _arb_round_sum(acb_realref(res), T, nT, prec, erp);

    if (nU == 2)
        _arb_round_sum2(acb_imagref(res), U, U + 1, prec, frp);
    else if (nU == 1)
        _arb_round_term1(acb_imagref(res), U, prec, frp);
    else
        _arb_round_sum(acb_imagref(res), U, nU, prec, frp);

    if (alloc > ARF_MUL_TLS_ALLOC)
        flint_free(tmp);

    return 1;
}
