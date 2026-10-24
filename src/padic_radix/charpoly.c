/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr.h"
#include "gr_vec.h"
#include "gr_mat.h"
#include "gr_poly.h"
#include "padic_radix.h"
#include "impl.h"

/*
    Matrix-vector products with a matrix packed once in the limb-aligned
    form (impl.h), for algorithms that multiply the same matrix by many
    vectors, such as the Berkowitz algorithm. Each output is computed with
    the semantics of padic_radix_dot_strided_delayed.
*/

typedef struct
{
    slong r;
    slong c;
    slong L;            /* limbs per residue */
    slong * vv;
    slong * ve;
    slong * NN;
    slong * ww;
    slong * sz;
    nn_ptr U;
}
_packed_mat_struct;

/* pack A (all entries) with residues of L limbs; returns 0 if an exponent
   is not moderate */
static int
_packed_mat_init(_packed_mat_struct * P, const gr_mat_t A, slong L, gr_ctx_t ctx)
{
    const radix_struct * radix = PADIC_RADIX_CTX_RADIX(ctx);
    padic_radix_pack_meta_struct meta;
    slong i, j, q, len = A->r * A->c;

    P->r = A->r;
    P->c = A->c;
    P->L = L;
    P->vv = flint_malloc(sizeof(slong) * 5 * FLINT_MAX(len, 1));
    P->ve = P->vv + len;
    P->NN = P->ve + len;
    P->ww = P->NN + len;
    P->sz = P->ww + len;
    P->U = flint_malloc(sizeof(ulong) * L * FLINT_MAX(len, 1));

    for (i = 0; i < A->r; i++)
        for (j = 0; j < A->c; j++)
        {
            const padic_radix_struct * x = gr_mat_entry_srcptr(A, i, j, ctx);
            q = i * A->c + j;
            if (!_padic_radix_pack_meta(&meta, x, radix))
                return 0;
            P->vv[q] = meta.vv; P->ve[q] = meta.ve; P->NN[q] = meta.NN; P->ww[q] = meta.ww;
            P->sz[q] = (meta.vv >= PADIC_RADIX_PACK_NONE) ? 0 :
                _padic_radix_pack_unit(P->U + L * q, x, meta.ww, L, radix);
            if (P->sz[q] == 0)
                flint_mpn_zero(P->U + L * q, L);
        }

    return 1;
}

static void
_packed_mat_clear(_packed_mat_struct * P)
{
    flint_free(P->vv);
    flint_free(P->U);
}

/*
    res[i - r0] = sum_{j < len} P_ij vec[j] for r0 <= i < r1. res must not
    alias vec. Returns GR_UNABLE (res untouched) if some output needs more
    than P->L limbs or an exponent of vec is not moderate.
*/
static int
_packed_mat_vec(padic_radix_struct * res, const _packed_mat_struct * P, slong r0, slong r1,
    const padic_radix_struct * vec, slong len, gr_ctx_t ctx)
{
    const radix_struct * radix = PADIC_RADIX_CTX_RADIX(ctx);
    slong e = radix->exp, L = P->L, i, j;
    slong * xv, * xe, * xN, * xw, * xs, * oN, * oW, * os;
    nn_ptr XU, acc;
    padic_radix_pack_meta_struct meta;
    int status = GR_SUCCESS;
    TMP_INIT;

    TMP_START;
    xv = TMP_ALLOC(sizeof(slong) * (5 * len + 3 * (r1 - r0)));
    xe = xv + len; xN = xe + len; xw = xN + len; xs = xw + len;
    oN = xs + len; oW = oN + (r1 - r0); os = oW + (r1 - r0);

    for (j = 0; j < len; j++)
    {
        if (!_padic_radix_pack_meta(&meta, vec + j, radix))
            goto unable;
        xv[j] = meta.vv; xe[j] = meta.ve; xN[j] = meta.NN; xw[j] = meta.ww;
    }

    for (i = r0; i < r1; i++)
    {
        const slong * avv = P->vv + i * P->c, * ave = P->ve + i * P->c;
        const slong * aNN = P->NN + i * P->c, * aww = P->ww + i * P->c;
        slong vmin = 2 * PADIC_RADIX_PACK_BIG, Nt = 2 * PADIC_RADIX_PACK_BIG;
        slong W = 2 * PADIC_RADIX_PACK_BIG, N, s;

        for (j = 0; j < len; j++)
        {
            slong x1 = avv[j] + xv[j];
            slong e1 = aNN[j] + xe[j];
            slong e2 = ave[j] + xN[j];
            slong f = aww[j] + xw[j];
            vmin = FLINT_MIN(vmin, x1);
            e1 = FLINT_MIN(e1, e2);
            Nt = FLINT_MIN(Nt, e1);
            W = FLINT_MIN(W, f);
        }

        if (vmin >= PADIC_RADIX_PACK_NONE)
        {
            N = (Nt >= PADIC_RADIX_PACK_NONE) ? PADIC_RADIX_EXACT : Nt;
            s = 0;
        }
        else
        {
            N = _padic_radix_pack_prec(vmin, Nt, ctx);
            s = (N <= vmin) ? 0 : (slong) radix_div_exp((ulong) (N - e * W) + e - 1, radix);
            if (s > L)
                goto unable;
        }

        oN[i - r0] = N;
        oW[i - r0] = W;
        os[i - r0] = s;
    }

    XU = TMP_ALLOC(sizeof(ulong) * L * FLINT_MAX(len, 1));
    for (j = 0; j < len; j++)
    {
        xs[j] = (xv[j] >= PADIC_RADIX_PACK_NONE) ? 0 :
            _padic_radix_pack_unit(XU + L * j, vec + j, xw[j], L, radix);
        if (xs[j] == 0)
            flint_mpn_zero(XU + L * j, L);
    }

    acc = TMP_ALLOC(sizeof(ulong) * 3 * L);

    for (i = r0; i < r1; i++)
    {
        const slong * aww = P->ww + i * P->c, * asz = P->sz + i * P->c;
        nn_srcptr ua = P->U + L * i * P->c;
        slong N = oN[i - r0], W = oW[i - r0], nslots = os[i - r0];
        padic_radix_struct * c = res + (i - r0);
        nn_ptr d;

        if (nslots == 0)
        {
            radix_integer_zero(&c->u, radix);
            c->v = 0;
            c->N = N;
            status |= _padic_radix_finalize(c, ctx);
            continue;
        }

        flint_mpn_zero(acc, 3 * nslots);

        if (L == 1)
        {
            ulong c0 = 0, c1 = 0, c2 = 0, ph, pl;
            for (j = 0; j < len; j++)
            {
                if (aww[j] + xw[j] == W)
                {
                    umul_ppmm(ph, pl, ua[j], XU[j]);
                    add_sssaaaaaa(c2, c1, c0, c2, c1, c0, 0, ph, pl);
                }
            }
            acc[0] = c0; acc[1] = c1; acc[2] = c2;
        }
        else if (L == 2)
        {
            ulong c0 = 0, c1 = 0, c2 = 0, d0 = 0, d1 = 0, d2 = 0, ph, pl;
            for (j = 0; j < len; j++)
            {
                slong o = aww[j] + xw[j] - W;
                ulong a0 = ua[2 * j], a1 = ua[2 * j + 1], b0 = XU[2 * j], b1 = XU[2 * j + 1];

                if (o == 0)
                {
                    umul_ppmm(ph, pl, a0, b0);
                    add_sssaaaaaa(c2, c1, c0, c2, c1, c0, 0, ph, pl);
                    if (nslots == 2)
                    {
                        umul_ppmm(ph, pl, a0, b1);
                        add_sssaaaaaa(d2, d1, d0, d2, d1, d0, 0, ph, pl);
                        umul_ppmm(ph, pl, a1, b0);
                        add_sssaaaaaa(d2, d1, d0, d2, d1, d0, 0, ph, pl);
                    }
                }
                else if (o == 1 && nslots == 2)
                {
                    umul_ppmm(ph, pl, a0, b0);
                    add_sssaaaaaa(d2, d1, d0, d2, d1, d0, 0, ph, pl);
                }
            }
            acc[0] = c0; acc[1] = c1; acc[2] = c2;
            acc[3] = d0; acc[4] = d1; acc[5] = d2;
        }
        else
        {
            for (j = 0; j < len; j++)
            {
                slong an = asz[j], bn = xs[j], o;

                o = aww[j] + xw[j] - W;
                if (o >= nslots || an == 0 || bn == 0)
                    continue;

                _padic_radix_pack_addmul(acc, ua + L * j, an, XU + L * j, bn, o, nslots);
            }
        }

        d = radix_integer_fit_limbs(&c->u, nslots + 3, radix);
        c->u.size = _padic_radix_pack_normalise(d, acc, nslots, radix);
        status |= _padic_radix_pack_finish(c, W, N, ctx);
    }

    TMP_END;
    return status;

unable:
    TMP_END;
    return GR_UNABLE;
}

/*
    The Berkowitz algorithm (as _gr_mat_charpoly_berkowitz), with the
    matrix-vector products a_k = M_t a_(k-1) done by _packed_mat_vec on a
    matrix packed once. Falls back to the generic algorithm for contexts
    without a relative precision (which bounds the windows a priori), for
    matrices with exact nonzero entries, and for extreme exponents.
*/
int
_padic_radix_mat_charpoly_berkowitz(padic_radix_struct * cp, const gr_mat_t mat, gr_ctx_t ctx)
{
    const radix_struct * radix = PADIC_RADIX_CTX_RADIX(ctx);
    slong n = mat->r, e = radix->exp, r = PADIC_RADIX_CTX_PREC_REL(ctx), L, Lmax;
    slong i, k, t;
    padic_radix_struct * a, * A;
    _packed_mat_struct P;
    int status = GR_SUCCESS;

    if (n <= 5 || r == PADIC_RADIX_PREC_INF || r > PADIC_RADIX_PACK_MODERATE)
        return _gr_mat_charpoly_berkowitz(cp, mat, ctx);

    /* exact nonzero entries can give dot products with only exact terms,
       which keep exact sums exact: leave those to the generic algorithm */
    for (i = 0; i < n; i++)
        for (k = 0; k < n; k++)
        {
            const padic_radix_struct * x = gr_mat_entry_srcptr(mat, i, k, ctx);
            if (x->u.size != 0 && x->N == PADIC_RADIX_EXACT)
                return _gr_mat_charpoly_berkowitz(cp, mat, ctx);
        }

    /* residues of ceil(r / e) limbs suffice when the lowest term of an
       output is aligned with its frame (e.g. all valuations 0); in general
       an output window is at most
       r + 2e - 2 digits (see impl.h), and the matrix is repacked with that
       many limbs the first time a product needs them */
    L = FLINT_MAX(1, (slong) radix_div_exp((ulong) (r + e - 1), radix));
    Lmax = (slong) radix_div_exp((ulong) (r + 3 * e - 3), radix);

    if (!_packed_mat_init(&P, mat, L, ctx))
    {
        _packed_mat_clear(&P);
        return _gr_mat_charpoly_berkowitz(cp, mat, ctx);
    }

#define MAT(ii, jj) gr_mat_entry_srcptr(mat, ii, jj, ctx)

    a = gr_heap_init_vec(n * n, ctx);
    A = a + (n - 1) * n;

    status |= _gr_vec_zero(cp, n + 1, ctx);
    status |= padic_radix_neg(cp, MAT(0, 0), ctx);

    for (t = 1; t < n && status == GR_SUCCESS; t++)
    {
        for (i = 0; i <= t; i++)
            status |= padic_radix_set(a + i, MAT(i, t), ctx);

        status |= padic_radix_set(A, MAT(t, t), ctx);

        for (k = 1; k < t; k++)
        {
            int s = _packed_mat_vec(a + k * n, &P, 0, t + 1, a + (k - 1) * n, t + 1, ctx);
            if (s == GR_UNABLE && P.L < Lmax)
            {
                _packed_mat_clear(&P);
                if (!_packed_mat_init(&P, mat, Lmax, ctx))
                    flint_throw(FLINT_ERROR, "(%s)\n", __func__);
                s = _packed_mat_vec(a + k * n, &P, 0, t + 1, a + (k - 1) * n, t + 1, ctx);
            }
            if (s == GR_UNABLE)
            {
                s = GR_SUCCESS;
                for (i = 0; i <= t; i++)
                    s |= _gr_vec_dot(a + k * n + i, NULL, 0, MAT(i, 0), a + (k - 1) * n, t + 1, ctx);
            }
            status |= s;
            status |= padic_radix_set(A + k, a + k * n + t, ctx);
        }

        {
            int s = _packed_mat_vec(A + t, &P, t, t + 1, a + (t - 1) * n, t + 1, ctx);
            if (s == GR_UNABLE && P.L < Lmax)
            {
                _packed_mat_clear(&P);
                if (!_packed_mat_init(&P, mat, Lmax, ctx))
                    flint_throw(FLINT_ERROR, "(%s)\n", __func__);
                s = _packed_mat_vec(A + t, &P, t, t + 1, a + (t - 1) * n, t + 1, ctx);
            }
            if (s == GR_UNABLE)
                s = _gr_vec_dot(A + t, NULL, 0, MAT(t, 0), a + (t - 1) * n, t + 1, ctx);
            status |= s;
        }

        for (k = 0; k <= t; k++)
        {
            status |= _gr_vec_dot_rev(cp + k, cp + k, 1, A, cp, k, ctx);
            status |= padic_radix_sub(cp + k, cp + k, A + k, ctx);
        }
    }

    for (i = n; i > 0; i--)
        padic_radix_swap(cp + i, cp + i - 1, ctx);

    status |= padic_radix_one(cp, ctx);
    status |= _gr_poly_reverse(cp, cp, n + 1, n + 1, ctx);

    gr_heap_clear_vec(a, n * n, ctx);
    _packed_mat_clear(&P);

#undef MAT

    return status;
}

int
padic_radix_mat_charpoly(padic_radix_struct * cp, const gr_mat_t mat, gr_ctx_t ctx)
{
    return _padic_radix_mat_charpoly_berkowitz(cp, mat, ctx);
}
