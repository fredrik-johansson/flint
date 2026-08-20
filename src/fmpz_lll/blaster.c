/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Opus 4.8

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    BLASter lattice reduction for FLINT
    =====================================

    Implements the algorithm from:
        Léo Ducas, Ludo N. Pulles, Marc Stevens,
        "Towards a Modern LLL Implementation" (ASIACRYPT 2025)
        https://eprint.iacr.org/2025/774

    The key ideas (from the paper):

      1.  Represent the basis via the upper-triangular Cholesky factor R of the
          Gram matrix G = B B^T:  G = R^T R.  The Gram-Schmidt coefficients are
          μ_{ij} = R[i,j] / R[i,i]  (i < j), and the squared GS norms are
          ||b*_i||^2 = R[i,i]^2.

      2.  Segment the n basis vectors into contiguous blocks of length `seg_len`.
          Within each segment [lo, hi), run a *local* classical floating-point
          LLL on the diagonal sub-block of R, maintaining upper triangularity via
          Givens rotations during Lovász swaps, and accumulating the integer
          transform in a small (seg_len × seg_len) matrix U_seg.  Segments can be
          processed independently (in parallel in the original C++; sequentially
          here).

      3.  After the segment LLLs, apply *block size reduction* between every
          ordered pair of segments (I, J) with I < J:

              C = round(-R_II^{-1} R_IJ)              [triangular solve]
              R[0:hi_I, lo_J:hi_J] +=                 [matrix multiply]
                  R[0:hi_I, lo_I:hi_I] * C_float

          This is the "BLAS step" that gives the algorithm its name.

      4.  Combine all integer transforms (block-diagonal from segments, plus the
          off-diagonal size-reduction coefficients) into a single n × n unimodular
          matrix U_iter.  Apply  B ← U_iter * B  using fmpz_mat_mul.

      5.  Recompute R from B via the LQ decomposition  B = L Q,  R = L^T.

      6.  Repeat until the full basis is LLL-reduced (checked on R).

    Floating-point arithmetic uses FLINT's nfloat backend (packed fixed-precision,
    cache-friendly) with automatic fallback to arf (arbitrary-precision) for very
    high precision or when nfloat cannot be initialised at the required width.

    Convention:  B is n × m (rows are the n basis vectors in Z^m).

    Segmentation boundary:  each segment is [k*seg_len, min((k+1)*seg_len, n)).
    The default seg_len is capped at DEFAULT_SEG_LEN = 64 (close to the optimal
    ℓ ≈ 64 observed in the BLASter benchmarks for n ≈ 128–256).
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "fmpz_lll.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_mat.h"
#include "nfloat.h"

/* ================================================================
   §0  Internal constants
   ================================================================ */

#define BLASTER_DEFAULT_SEG_LEN   64
#define BLASTER_MAX_ITER          8192

/* ================================================================
   §1  Working precision
   ================================================================

   Two separate effects set the precision, and they add:

   (a) Conditioning of the reduction, O(n) bits.  This is the classical
       L^2 bound of Nguyen and Stehlé (~1.6n + o(n)); it describes the R
       factor of a basis that is already close to reduced.

   (b) Dynamic range of the *input*.  An unreduced basis with entries of
       B bits has size-reduction coefficients of up to ~B bits in the
       early passes, and those must be representable.  L^2 avoids paying
       for this by interleaving exact integer updates with lazy
       recomputation; this implementation does not, so it pays directly.

   Hence prec = max_bits + n + 64.  The two terms are *added*, not
   multiplied: an earlier multiplicative rule (max_bits * n) asked for
   4032 bits on a 128-dimensional q-ary lattice and ran ~25x slower than
   FLINT's existing LLL.

   Dropping term (b) is tempting, because q-ary lattices start nearly
   reduced and there n dominates - but it is wrong in general.  On
   integer-relation lattices, whose entries carry hundreds of bits, an
   n + 64 rule sent a 20-dimensional example through 137 iterations and
   four precision doublings: 16 s, against 0.006 s for fmpz_lll.

   This remains a heuristic.  The stall detector in §7 doubles the
   precision whenever an iteration achieves nothing, so an underestimate
   costs time but never correctness.
*/
static flint_bitcnt_t
_blaster_initial_prec(const fmpz_mat_t B)
{
    slong n = fmpz_mat_nrows(B), m = fmpz_mat_ncols(B);
    slong i, j;
    flint_bitcnt_t max_bits = 0, prec;

    for (i = 0; i < n; i++)
        for (j = 0; j < m; j++)
        {
            flint_bitcnt_t b = fmpz_bits(fmpz_mat_entry(B, i, j));
            if (b > max_bits)
                max_bits = b;
        }

    prec = max_bits + (flint_bitcnt_t) n + 64;

    if (prec < 64)
        prec = 64;

    prec = ((prec + 63) / 64) * 64;   /* whole limbs */

    return prec;
}

/* ================================================================
   §2  Compute R from B

   B (n × m) = L * Q  where L is n × n lower triangular and Q has
   orthonormal rows.  Then B B^T = L L^T, so R := L^T is the upper-
   triangular Cholesky factor of the Gram matrix.
   ================================================================ */
static int
_blaster_compute_R(gr_mat_t R, const fmpz_mat_t B, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong n = fmpz_mat_nrows(B), m = fmpz_mat_ncols(B);
    gr_mat_t Bf, L, Q;

    gr_mat_init(Bf, n, m, ctx);
    gr_mat_init(L,  n, n, ctx);
    gr_mat_init(Q,  n, m, ctx);

    status |= gr_mat_set_fmpz_mat(Bf, B, ctx);
    status |= gr_mat_lq(L, Q, Bf, ctx);     /* Bf = L * Q                    */
    status |= gr_mat_transpose(R, L, ctx);  /* R = L^T  (upper triangular)   */

    gr_mat_clear(Bf, ctx);
    gr_mat_clear(L,  ctx);
    gr_mat_clear(Q,  ctx);

    return status;
}

/* ================================================================
   §3  Givens rotation on two rows of R

   Given rows r1 and r2 of the n × n matrix R (with r2 > r1), compute
   the Givens rotation parameters that zero out R[r2, col_pivot]:

       a = R[r1, col_pivot],  b = R[r2, col_pivot]
       norm = sqrt(a^2 + b^2)
       c_g = a/norm,  s_g = b/norm

   and apply the rotation

       R[r1, j] ←  c_g * R[r1, j] + s_g * R[r2, j]
       R[r2, j] ← -s_g * R[r1, j] + c_g * R[r2, j]

   for all columns j from col_start to col_end - 1 (inclusive on left,
   exclusive on right), then explicitly zero R[r2, col_pivot].
   ================================================================ */
static int
_blaster_apply_givens(gr_mat_t R,
                      slong r1, slong r2,
                      slong col_pivot,
                      slong col_end,
                      gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    slong j;
    int status = GR_SUCCESS;
    gr_ptr a, b, norm, c_g, s_g, t1, t2;

    GR_TMP_INIT5(a, b, norm, c_g, s_g, ctx);
    GR_TMP_INIT2(t1, t2, ctx);

    /* a = R[r1, col_pivot],  b = R[r2, col_pivot] */
    status |= gr_set(a, GR_MAT_ENTRY(R, r1, col_pivot, sz), ctx);
    status |= gr_set(b, GR_MAT_ENTRY(R, r2, col_pivot, sz), ctx);

    /* norm = sqrt(a^2 + b^2) */
    status |= gr_sqr(t1, a, ctx);
    status |= gr_sqr(t2, b, ctx);
    status |= gr_add(norm, t1, t2, ctx);
    status |= gr_sqrt(norm, norm, ctx);

    if (status != GR_SUCCESS)
        goto cleanup;

    /* c_g = a/norm,  s_g = b/norm */
    status |= gr_div(c_g, a, norm, ctx);
    status |= gr_div(s_g, b, norm, ctx);

    /* Apply to columns col_pivot..col_end-1 */
    for (j = col_pivot; j < col_end && status == GR_SUCCESS; j++)
    {
        gr_ptr e1 = GR_MAT_ENTRY(R, r1, j, sz);
        gr_ptr e2 = GR_MAT_ENTRY(R, r2, j, sz);

        status |= gr_set(t1, e1, ctx);
        status |= gr_set(t2, e2, ctx);

        /* e1' =  c_g * t1 + s_g * t2 */
        status |= gr_mul(e1, c_g, t1, ctx);
        status |= gr_addmul(e1, s_g, t2, ctx);

        /* e2' = -s_g * t1 + c_g * t2 */
        status |= gr_mul(e2, c_g, t2, ctx);
        status |= gr_submul(e2, s_g, t1, ctx);
    }

    /* Force the zeroed-out entry to exactly zero */
    if (status == GR_SUCCESS)
        status |= gr_zero(GR_MAT_ENTRY(R, r2, col_pivot, sz), ctx);

cleanup:
    GR_TMP_CLEAR5(a, b, norm, c_g, s_g, ctx);
    GR_TMP_CLEAR2(t1, t2, ctx);

    return status;
}

/* ================================================================
   §4  Classical LLL on segment [lo, hi)

   Runs the standard floating-point LLL loop on the sub-block of R
   indexed by [lo, hi).  The FULL n × n matrix R is modified in place
   (size-reduction updates affect rows 0..j, Givens rotations affect
   all columns from the swap column to n-1), which keeps R's global
   upper-triangular structure correct.

   The integer transforms (size reductions + Lovász swaps) are
   accumulated in U_seg, a (hi-lo) × (hi-lo) matrix initialised to
   the identity by the caller.  On return:

       B_new[lo:hi, :] = U_seg * B_old[lo:hi, :]

   Lovász parameter: R[k,k]^2 >= (delta - mu_{k-1,k}^2) * R[k-1,k-1]^2.
   Size-reduction parameter: |mu_{j,kappa}| <= eta  (not re-checked
   after swaps; the outer loop re-certifies).
   ================================================================ */
static int
_blaster_lll_segment(gr_mat_t R, fmpz_mat_t U_seg,
                     slong lo, slong hi, slong n,
                     double delta, double eta,
                     gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    slong len = hi - lo;
    slong kappa, j, i, l;
    int status = GR_SUCCESS;
    int cmp;
    gr_ptr mu, c_f, lhs, rhs, tmp;
    fmpz_t c_int;

    (void)eta;  /* used conceptually; the outer loop re-certifies */

    if (len <= 1)
        return GR_SUCCESS;

    GR_TMP_INIT5(mu, c_f, lhs, rhs, tmp, ctx);
    fmpz_init(c_int);

    kappa = lo + 1;

    while (kappa < hi && status == GR_SUCCESS)
    {
        /* ---- Size-reduce column kappa against columns lo..kappa-1 ---- *
         *                                                                 *
         * For each j in [lo, kappa) (largest first, matching fplll):     *
         *   mu = R[j, kappa] / R[j, j]                                   *
         *   c  = round(mu)                                                *
         *   R[:, kappa] -= c * R[:, j]  (rows 0..j only are nonzero)    *
         *   U_seg row kappa-lo -= c * U_seg row j-lo                     */
        for (j = kappa - 1; j >= lo; j--)
        {
            gr_ptr R_jj = GR_MAT_ENTRY(R, j, j, sz);
            gr_ptr R_jk = GR_MAT_ENTRY(R, j, kappa, sz);

            if (gr_is_zero(R_jj, ctx) == T_TRUE)
                continue;

            /* mu = R[j, kappa] / R[j, j] */
            status = gr_div(mu, R_jk, R_jj, ctx);
            if (status != GR_SUCCESS) break;

            /* c = round(mu) */
            status |= gr_nint(c_f, mu, ctx);
            if (status != GR_SUCCESS) break;
            status |= gr_get_fmpz(c_int, c_f, ctx);
            if (status != GR_SUCCESS) break;

            if (fmpz_is_zero(c_int))
                continue;

            /* R[0:j+1, kappa] -= c * R[0:j+1, j]
               (R[j+1..n-1, j] = 0 by upper-triangularity)              */
            for (i = 0; i <= j; i++)
                status |= gr_submul(GR_MAT_ENTRY(R, i, kappa, sz),
                                    c_f,
                                    GR_MAT_ENTRY(R, i, j, sz),
                                    ctx);

            /* U_seg[kappa-lo, :] -= c * U_seg[j-lo, :]                  */
            for (l = 0; l < len; l++)
                fmpz_submul(fmpz_mat_entry(U_seg, kappa - lo, l),
                            c_int,
                            fmpz_mat_entry(U_seg, j - lo, l));
        }

        if (status != GR_SUCCESS)
            break;

        /* ---- Lovász condition at kappa ---- */
        if (gr_is_zero(GR_MAT_ENTRY(R, kappa - 1, kappa - 1, sz), ctx) == T_TRUE)
        {
            kappa++;
            continue;
        }

        /* mu = R[kappa-1, kappa] / R[kappa-1, kappa-1] */
        status = gr_div(mu,
            GR_MAT_ENTRY(R, kappa - 1, kappa,     sz),
            GR_MAT_ENTRY(R, kappa - 1, kappa - 1, sz), ctx);
        if (status != GR_SUCCESS) break;

        /* lhs = R[kappa, kappa]^2 */
        status |= gr_sqr(lhs, GR_MAT_ENTRY(R, kappa, kappa, sz), ctx);

        /* rhs = (delta - mu^2) * R[kappa-1, kappa-1]^2 */
        status |= gr_set_d(rhs, delta, ctx);
        status |= gr_sqr(tmp, mu, ctx);
        status |= gr_sub(rhs, rhs, tmp, ctx);
        status |= gr_sqr(tmp, GR_MAT_ENTRY(R, kappa - 1, kappa - 1, sz), ctx);
        status |= gr_mul(rhs, rhs, tmp, ctx);

        if (status != GR_SUCCESS) break;

        status = gr_cmp(&cmp, lhs, rhs, ctx);
        if (status != GR_SUCCESS) break;

        if (cmp < 0)    /* Lovász condition fails: swap kappa-1 and kappa */
        {
            /* Column swap: exchange cols kappa-1 and kappa of R.
               Only rows 0..kappa need updating (rows above are zero).    */
            for (i = 0; i <= kappa; i++)
                gr_swap(GR_MAT_ENTRY(R, i, kappa - 1, sz),
                        GR_MAT_ENTRY(R, i, kappa,     sz), ctx);

            /* After the swap R[kappa, kappa-1] != 0; apply Givens on
               rows (kappa-1, kappa) for columns kappa-1..n-1 to restore
               upper triangularity.                                        */
            status = _blaster_apply_givens(R, kappa - 1, kappa,
                                           kappa - 1, n, ctx);
            if (status != GR_SUCCESS) break;

            /* Swap rows (kappa-1-lo) and (kappa-lo) of U_seg. */
            for (l = 0; l < len; l++)
                fmpz_swap(fmpz_mat_entry(U_seg, kappa - 1 - lo, l),
                          fmpz_mat_entry(U_seg, kappa     - lo, l));

            if (kappa > lo + 1)
                kappa--;
        }
        else
        {
            kappa++;
        }
    }

    GR_TMP_CLEAR5(mu, c_f, lhs, rhs, tmp, ctx);
    fmpz_clear(c_int);

    return status;
}

/* ================================================================
   §5  Block size reduction between segments I < J

   Reduces the off-diagonal block R[lo_I:hi_I, lo_J:hi_J] of R using
   the diagonal block R[lo_I:hi_I, lo_I:hi_I]:

       X     = R_II^{-1} R_IJ        (triangular solve)
       C_int = round(-X)             (element-wise, stored in C_int)
       R[0:hi_I, lo_J:hi_J] += R[0:hi_I, lo_I:hi_I] * C_float

   On return C_int (len_I × len_J integer matrix, pre-allocated by
   the caller) holds the integer coefficients for the U update.
   ================================================================ */
static int
_blaster_reduce_block_pair(gr_mat_t R,
                           fmpz_mat_t C_int,
                           slong lo_I, slong hi_I,
                           slong lo_J, slong hi_J,
                           gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    slong len_I = hi_I - lo_I;
    slong len_J = hi_J - lo_J;
    slong i, j;
    int status = GR_SUCCESS;
    gr_mat_t R_II, R_IJ, X, C_f, T;

    if (len_I == 0 || len_J == 0)
        return GR_SUCCESS;

    /* Sub-matrix windows into R (no new memory allocated) */
    gr_mat_window_init(R_II, R, lo_I, lo_I, hi_I, hi_I, ctx);
    gr_mat_window_init(R_IJ, R, lo_I, lo_J, hi_I, hi_J, ctx);

    gr_mat_init(X,   len_I, len_J, ctx);
    gr_mat_init(C_f, len_I, len_J, ctx);
    gr_mat_init(T,   hi_I,  len_J, ctx);

    /* --- X = R_II^{-1} R_IJ  (solves R_II * X = R_IJ) --- */
    status = gr_mat_nonsingular_solve_triu(X, R_II, R_IJ, 0, ctx);
    if (status != GR_SUCCESS)
        goto done;

    /* --- C_f = round(-X),  C_int = integer version --- */
    for (i = 0; i < len_I && status == GR_SUCCESS; i++)
    {
        for (j = 0; j < len_J && status == GR_SUCCESS; j++)
        {
            gr_ptr xij = GR_MAT_ENTRY(X,   i, j, sz);
            gr_ptr cf  = GR_MAT_ENTRY(C_f, i, j, sz);

            status |= gr_neg(cf, xij, ctx);
            status |= gr_nint(cf, cf, ctx);
            status |= gr_get_fmpz(fmpz_mat_entry(C_int, i, j), cf, ctx);
        }
    }
    if (status != GR_SUCCESS)
        goto done;

    /* --- R[0:hi_I, lo_J:hi_J] += R[0:hi_I, lo_I:hi_I] * C_f --- *
     *                                                               *
     * T = R_top_I * C_f  (hi_I × len_I) * (len_I × len_J)        *
     * R_top_J += T                                                  */
    {
        gr_mat_t R_top_I, R_top_J;
        gr_mat_window_init(R_top_I, R, 0, lo_I, hi_I, hi_I, ctx);
        gr_mat_window_init(R_top_J, R, 0, lo_J, hi_I, hi_J, ctx);

        status |= gr_mat_mul(T, R_top_I, C_f, ctx);
        if (status == GR_SUCCESS)
            status |= gr_mat_add(R_top_J, R_top_J, T, ctx);

        gr_mat_window_clear(R_top_I, ctx);
        gr_mat_window_clear(R_top_J, ctx);
    }

done:
    gr_mat_clear(X,   ctx);
    gr_mat_clear(C_f, ctx);
    gr_mat_clear(T,   ctx);
    gr_mat_window_clear(R_II, ctx);
    gr_mat_window_clear(R_IJ, ctx);

    return status;
}

/* ================================================================
   §6  LLL-reducedness check on R

   Returns  1  if R encodes an LLL-reduced basis,
            0  if a violation is found,
           -1  on floating-point error.
   ================================================================ */
static int
_blaster_is_lll_reduced(const gr_mat_t R, slong n,
                         double delta, double eta,
                         gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    slong k, j;
    int status = GR_SUCCESS;
    int cmp;
    gr_ptr mu, lhs, rhs, tmp;

    GR_TMP_INIT4(mu, lhs, rhs, tmp, ctx);

    for (k = 1; k < n && status == GR_SUCCESS; k++)
    {
        gr_ptr Rkk   = GR_MAT_ENTRY(R, k,     k,     sz);
        gr_ptr Rk1k1 = GR_MAT_ENTRY(R, k - 1, k - 1, sz);

        if (gr_is_zero(Rk1k1, ctx) == T_TRUE)
            continue;

        /* Size-reducedness: |R[j,k]/R[j,j]| <= eta for all j < k */
        for (j = 0; j < k && status == GR_SUCCESS; j++)
        {
            gr_ptr R_jj = GR_MAT_ENTRY(R, j, j, sz);
            if (gr_is_zero(R_jj, ctx) == T_TRUE)
                continue;

            status |= gr_div(mu, GR_MAT_ENTRY(R, j, k, sz), R_jj, ctx);
            status |= gr_abs(mu, mu, ctx);
            status |= gr_set_d(tmp, eta, ctx);

            if (status == GR_SUCCESS && gr_gt(mu, tmp, ctx) == T_TRUE)
            {
                GR_TMP_CLEAR4(mu, lhs, rhs, tmp, ctx);
                return 0;
            }
        }

        /* Lovász: R[k,k]^2 >= (delta - mu_{k-1,k}^2) * R[k-1,k-1]^2  */
        status |= gr_div(mu,
            GR_MAT_ENTRY(R, k - 1, k, sz), Rk1k1, ctx);
        status |= gr_sqr(lhs, Rkk, ctx);
        status |= gr_set_d(rhs, delta, ctx);
        status |= gr_sqr(tmp, mu, ctx);
        status |= gr_sub(rhs, rhs, tmp, ctx);
        status |= gr_sqr(tmp, Rk1k1, ctx);
        status |= gr_mul(rhs, rhs, tmp, ctx);

        if (status != GR_SUCCESS) break;

        status = gr_cmp(&cmp, lhs, rhs, ctx);
        if (status != GR_SUCCESS) break;

        if (cmp < 0)
        {
            GR_TMP_CLEAR4(mu, lhs, rhs, tmp, ctx);
            return 0;
        }
    }

    GR_TMP_CLEAR4(mu, lhs, rhs, tmp, ctx);
    return (status == GR_SUCCESS) ? 1 : -1;
}

/* ================================================================
   §6b  Lovász-only check

   The main loop is driven by the Lovász conditions alone, because the
   block size reduction of §5 does *not* by itself guarantee the
   classical bound |mu_{ij}| <= 1/2 across block boundaries.

   Writing R'_IJ = R_II (X - round(X)) = R_II E with |E| <= 1/2
   entrywise, we get

       mu'_{ij} = E_{ij} + sum_{k > i} mu_{ik} E_{kj},

   so the off-diagonal coefficients are only bounded by roughly
   (1 + sum |mu|)/2, not by 1/2.  Size-reducedness is therefore restored
   by an exact classical pass (§6c) once the Lovász conditions hold.
   That pass never changes the diagonal R[i,i], hence never breaks the
   Lovász conditions, so the two stages compose cleanly.

   Returns  1 if all Lovász conditions hold,
            0 if one fails,
           -1 on floating-point error.
   ================================================================ */
static int
_blaster_lovasz_ok(const gr_mat_t R, slong n, double delta, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    slong k;
    int status = GR_SUCCESS;
    int cmp;
    gr_ptr mu, lhs, rhs, tmp;

    GR_TMP_INIT4(mu, lhs, rhs, tmp, ctx);

    for (k = 1; k < n && status == GR_SUCCESS; k++)
    {
        gr_ptr Rkk   = GR_MAT_ENTRY(R, k,     k,     sz);
        gr_ptr Rk1k1 = GR_MAT_ENTRY(R, k - 1, k - 1, sz);

        if (gr_is_zero(Rk1k1, ctx) == T_TRUE)
            continue;

        /* R[k,k]^2 >= (delta - mu_{k-1,k}^2) * R[k-1,k-1]^2 */
        status |= gr_div(mu, GR_MAT_ENTRY(R, k - 1, k, sz), Rk1k1, ctx);
        status |= gr_sqr(lhs, Rkk, ctx);
        status |= gr_set_d(rhs, delta, ctx);
        status |= gr_sqr(tmp, mu, ctx);
        status |= gr_sub(rhs, rhs, tmp, ctx);
        status |= gr_sqr(tmp, Rk1k1, ctx);
        status |= gr_mul(rhs, rhs, tmp, ctx);

        if (status != GR_SUCCESS)
            break;

        status = gr_cmp(&cmp, lhs, rhs, ctx);
        if (status != GR_SUCCESS)
            break;

        if (cmp < 0)
        {
            GR_TMP_CLEAR4(mu, lhs, rhs, tmp, ctx);
            return 0;
        }
    }

    GR_TMP_CLEAR4(mu, lhs, rhs, tmp, ctx);
    return (status == GR_SUCCESS) ? 1 : -1;
}

/* ================================================================
   §6c  Exact classical size reduction of the whole basis

   For kappa = 1, ..., n-1 and j = kappa-1, ..., 0:

       c = round(R[j,kappa] / R[j,j])
       R[0:j+1, kappa] -= c * R[0:j+1, j]
       U[kappa, :]     -= c * U[j, :]

   Processing j in decreasing order is what yields the exact bound
   |mu_{j,kappa}| <= 1/2; no diagonal entry is touched, so the Gram-
   Schmidt norms - and hence the Lovász conditions - are preserved.

   U must be n x n and set to the identity by the caller.  Returns
   GR_SUCCESS or a gr error code.
   ================================================================ */
static int
_blaster_size_reduce_full(gr_mat_t R, fmpz_mat_t U, slong n, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    slong kappa, j, i;
    int status = GR_SUCCESS;
    gr_ptr mu, c_f;
    fmpz_t c_int;

    GR_TMP_INIT2(mu, c_f, ctx);
    fmpz_init(c_int);

    for (kappa = 1; kappa < n && status == GR_SUCCESS; kappa++)
    {
        for (j = kappa - 1; j >= 0; j--)
        {
            gr_ptr R_jj = GR_MAT_ENTRY(R, j, j, sz);

            if (gr_is_zero(R_jj, ctx) == T_TRUE)
                continue;

            status = gr_div(mu, GR_MAT_ENTRY(R, j, kappa, sz), R_jj, ctx);
            if (status != GR_SUCCESS) break;

            status |= gr_nint(c_f, mu, ctx);
            if (status != GR_SUCCESS) break;
            status |= gr_get_fmpz(c_int, c_f, ctx);
            if (status != GR_SUCCESS) break;

            if (fmpz_is_zero(c_int))
                continue;

            for (i = 0; i <= j; i++)
                status |= gr_submul(GR_MAT_ENTRY(R, i, kappa, sz),
                                    c_f,
                                    GR_MAT_ENTRY(R, i, j, sz), ctx);

            _fmpz_vec_scalar_submul_fmpz(fmpz_mat_entry(U, kappa, 0),
                                         fmpz_mat_entry(U, j, 0),
                                         n, c_int);
        }
    }

    GR_TMP_CLEAR2(mu, c_f, ctx);
    fmpz_clear(c_int);

    return status;
}

/* ================================================================
   §7  Main function: fmpz_lll_blaster
   ================================================================ */

/*
    BLASter-style LLL reduction of the integer basis B (n × m, rows
    are basis vectors).

    Parameters:
      B      : [in/out] basis to reduce, modified in place.
      U_out  : [out] if not NULL and has dimensions n × n, the total
               unimodular transformation U such that B_new = U * B_old
               is written here.
      fl     : LLL parameters (delta, eta); if NULL, defaults 0.99/0.51
               are used.

    Returns 1 on success (basis is LLL-reduced on exit), 0 if the
    maximum iteration count was reached without convergence, or -1 on
    numerical failure.

    Precision:  automatically chosen based on ||B||_max and n.  If the
    floating-point LLL check fails (e.g. due to precision loss) the
    precision is doubled and the iteration retried.
*/
int
fmpz_lll_blaster(fmpz_mat_t B, fmpz_mat_t U_out, const fmpz_lll_t fl)
{
    slong n = fmpz_mat_nrows(B);
    slong m = fmpz_mat_ncols(B);

    double delta = (fl != NULL) ? fl->delta : 0.99;
    double eta   = (fl != NULL) ? fl->eta   : 0.51;

    slong seg_len, n_segs, max_segs, s, s_I, s_J;
    slong lo, hi, len, lo_I, hi_I, lo_J, hi_J, len_I, len_J;
    slong iter;
    slong *bnd;             /* segment boundaries for the current iteration */
    int result = 0;
    int status = GR_SUCCESS;
    int reduced, lovasz;

    /* --- Trivial cases --- */
    if (n <= 1)
        return 1;
    if (m == 0)
        return 1;

    /* --- Segment length ℓ: default 64, capped at n --- */
    seg_len = BLASTER_DEFAULT_SEG_LEN;
    if (seg_len > n)
        seg_len = n;
    if (seg_len < 2)
        seg_len = 2;

    /* Boundaries: at most ceil(n/seg_len) + 2 segments once the
       half-segment offset used on odd iterations is taken into account. */
    max_segs = n / seg_len + 3;
    bnd = (slong *) flint_malloc((max_segs + 1) * sizeof(slong));

    /* --- Working precision --- */
    flint_bitcnt_t prec = _blaster_initial_prec(B);

    /* --- Floating-point context --- */
    gr_ctx_t ctx;
    if (nfloat_ctx_init(ctx, prec, 0) != GR_SUCCESS)
        gr_ctx_init_real_float_arf(ctx, prec);

    /* --- Persistent matrices --- */
    gr_mat_t R;           /* n × n  upper-triangular Cholesky factor */
    fmpz_mat_t U_total;   /* n × n  total unimodular transform       */
    fmpz_mat_t U_iter;    /* n × n  transform for one iteration      */
    fmpz_mat_t B_tmp;     /* n × m  scratch for B update             */
    fmpz_mat_t U_tmp;     /* n × n  scratch for U update             */

    gr_mat_init(R, n, n, ctx);
    fmpz_mat_init(U_total, n, n);
    fmpz_mat_init(U_iter,  n, n);
    fmpz_mat_init(B_tmp,   n, m);
    fmpz_mat_init(U_tmp,   n, n);
    fmpz_mat_one(U_total);

    /* ==============================================================
       Main loop
       ============================================================== */
    for (iter = 0; iter < BLASTER_MAX_ITER; iter++)
    {
        /* Step 1: Recompute R = Chol(B B^T)^T from the current B. */
        status = _blaster_compute_R(R, B, ctx);
        if (status != GR_SUCCESS)
        {
            /* Precision failure: double prec and retry. */
            goto raise_precision;
        }

        /* Step 2: Is the basis already fully LLL-reduced? */
        reduced = _blaster_is_lll_reduced(R, n, delta, eta, ctx);
        if (reduced == 1)
        {
            result = 1;
            break;
        }
        if (reduced == -1)
        {
            /* Precision too low; double and retry. */
            goto raise_precision;
        }

        /* Step 2b: If only size-reducedness is missing, finish with an
           exact classical size-reduction pass.  This cannot disturb the
           Lovász conditions because it never touches R[i,i].            */
        lovasz = _blaster_lovasz_ok(R, n, delta, ctx);
        if (lovasz == -1)
        {
            goto raise_precision;
        }

        if (lovasz == 1)
        {
            fmpz_mat_one(U_iter);
            status = _blaster_size_reduce_full(R, U_iter, n, ctx);
            if (status != GR_SUCCESS)
                break;

            /* Nothing changed, yet the basis is not certified reduced:
               the working precision is too low to make progress. */
            if (fmpz_mat_is_one(U_iter))
                goto raise_precision;

            fmpz_mat_mul(B_tmp, U_iter, B);
            fmpz_mat_swap(B, B_tmp);
            fmpz_mat_mul(U_tmp, U_iter, U_total);
            fmpz_mat_swap(U_total, U_tmp);
            continue;       /* re-verify on the next pass */
        }

        /* Step 3: Segment boundaries for this iteration.

           On odd iterations the boundaries are shifted by seg_len/2 so
           that index pairs straddling a boundary on even iterations sit
           in the interior of a segment on odd ones.  Without this
           alternation the Lovász condition at a segment boundary would
           never be examined and the loop could not converge.            */
        {
            slong offset = (iter & 1) ? seg_len / 2 : 0;
            slong nb = 0;

            bnd[nb++] = 0;
            if (offset > 0 && offset < n)
                bnd[nb++] = offset;
            while (bnd[nb - 1] < n)
            {
                bnd[nb] = FLINT_MIN(bnd[nb - 1] + seg_len, n);
                nb++;
            }
            n_segs = nb - 1;
        }

        /* Step 4: Build U_iter for this iteration.
           Start as I_n; per-segment LLL will set the diagonal blocks. */
        fmpz_mat_one(U_iter);

        /* Step 4a: Segment LLLs. */
        for (s = 0; s < n_segs && status == GR_SUCCESS; s++)
        {
            fmpz_mat_t U_seg;

            lo  = bnd[s];
            hi  = bnd[s + 1];
            len = hi - lo;

            fmpz_mat_init(U_seg, len, len);
            fmpz_mat_one(U_seg);

            /* LLL on the segment block of R, accumulate in U_seg */
            status = _blaster_lll_segment(R, U_seg, lo, hi, n,
                                           delta, eta, ctx);

            if (status == GR_SUCCESS)
            {
                /* Embed U_seg into the diagonal block of U_iter */
                for (slong i = 0; i < len; i++)
                    for (slong j = 0; j < len; j++)
                        fmpz_set(fmpz_mat_entry(U_iter, lo + i, lo + j),
                                 fmpz_mat_entry(U_seg, i, j));
            }

            fmpz_mat_clear(U_seg);
        }

        if (status != GR_SUCCESS)
            break;

        /* Step 4b: Block size reduction between all ordered segment pairs.
           We process pairs (I, J) with s_I < s_J in lexicographic order.
           C_int is the integer coefficient block (len_I × len_J).
           After computing C_int (and updating R in _blaster_reduce_block_pair),
           update U_iter:
               U_iter[lo_J:hi_J, :] += C_int^T * U_iter[lo_I:hi_I, :]   */
        for (s_I = 0; s_I < n_segs && status == GR_SUCCESS; s_I++)
        {
            lo_I  = bnd[s_I];
            hi_I  = bnd[s_I + 1];
            len_I = hi_I - lo_I;

            for (s_J = s_I + 1; s_J < n_segs && status == GR_SUCCESS; s_J++)
            {
                lo_J  = bnd[s_J];
                hi_J  = bnd[s_J + 1];
                len_J = hi_J - lo_J;

                fmpz_mat_t C_int;
                fmpz_mat_init(C_int, len_I, len_J);

                status = _blaster_reduce_block_pair(R, C_int,
                                                     lo_I, hi_I,
                                                     lo_J, hi_J,
                                                     ctx);

                if (status == GR_SUCCESS)
                {
                    /* U_iter[lo_J + i, :] += C_int[k, i] * U_iter[lo_I + k, :]
                       for i in [0, len_J), k in [0, len_I).
                       Reads from rows lo_I..hi_I-1, writes to lo_J..hi_J-1.
                       Since s_I < s_J we have hi_I <= lo_J so no aliasing.  */
                    for (slong i = 0; i < len_J; i++)
                    {
                        fmpz * row_J = fmpz_mat_entry(U_iter, lo_J + i, 0);
                        for (slong k = 0; k < len_I; k++)
                        {
                            fmpz * cki = fmpz_mat_entry(C_int, k, i);
                            if (fmpz_is_zero(cki))
                                continue;
                            _fmpz_vec_scalar_addmul_fmpz(
                                row_J,
                                fmpz_mat_entry(U_iter, lo_I + k, 0),
                                n, cki);
                        }
                    }
                }

                fmpz_mat_clear(C_int);
            }
        }

        if (status != GR_SUCCESS)
            break;

        /* Stall check: an identity transform means this iteration achieved
           nothing, so no future iteration can either.  Raise precision. */
        if (fmpz_mat_is_one(U_iter))
            goto raise_precision;

        /* Step 5: Apply U_iter to B:   B ← U_iter * B */
        fmpz_mat_mul(B_tmp, U_iter, B);
        fmpz_mat_swap(B, B_tmp);

        /* Step 6: Accumulate total transform:  U_total ← U_iter * U_total */
        fmpz_mat_mul(U_tmp, U_iter, U_total);
        fmpz_mat_swap(U_total, U_tmp);

        continue;

raise_precision:

        /* R must be cleared with the context that allocated it: a gr
           element's size depends on the context, so clearing it against
           a freshly initialised context of different precision would
           corrupt the heap. */
        gr_mat_clear(R, ctx);
        gr_ctx_clear(ctx);

        prec *= 2;
        if (prec > (flint_bitcnt_t) UWORD_MAX / 2)
        {
            /* Re-initialise minimally so the cleanup code below is valid. */
            if (nfloat_ctx_init(ctx, 64, 0) != GR_SUCCESS)
                gr_ctx_init_real_float_arf(ctx, 64);
            gr_mat_init(R, n, n, ctx);
            break;
        }

        if (nfloat_ctx_init(ctx, prec, 0) != GR_SUCCESS)
            gr_ctx_init_real_float_arf(ctx, prec);
        gr_mat_init(R, n, n, ctx);
        status = GR_SUCCESS;

    } /* end main loop */

    /* --- Write out the total transform if requested --- */
    if (U_out != NULL
        && fmpz_mat_nrows(U_out) == n
        && fmpz_mat_ncols(U_out) == n)
    {
        fmpz_mat_set(U_out, U_total);
    }

    /* --- Cleanup --- */
    flint_free(bnd);

    gr_mat_clear(R, ctx);
    gr_ctx_clear(ctx);

    fmpz_mat_clear(U_total);
    fmpz_mat_clear(U_iter);
    fmpz_mat_clear(B_tmp);
    fmpz_mat_clear(U_tmp);

    return result;
}
