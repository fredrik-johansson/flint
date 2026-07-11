/* Hand-optimized replacements for the n = 1, 2 and 3 specialized exp
   series (see WINNERS_OPT in dev/gen_fixed_exp_rs_hard.py, which skips
   these two sizes and defers to this file).  At these sizes the RS
   code generator's machinery -- the shift/split framework, the calls
   to flint_mpn_sqrhigh / flint_mpn_mulhigh_n / mpn_addmul_1, the
   fac20_inv rescaling -- dwarfs the handful of word multiplies the
   series actually needs, so the evaluators are written directly.

   Both compute exp(x) for a reduced argument x with x < 2^-r (r = 12
   for n = 1, r = 16 for n = 2), returning n fraction limbs and the
   unit limb, which is always 1 (exp(x) < 1 + 2x < 2 here).

   The work is done at MIXED precision.  Writing the series as

       exp(x) = 1 + x + x^2 P,   P = C2 + x Q,   Q = C3 + x R,   ...

   with C_k = 1/k!, an error in an inner Horner value is attenuated by
   the powers of x multiplying it: err(P) enters the result scaled by
   x^2 < 2^-2r, and each level below that by a further x < 2^-r.  So at
   n = 2 (output ulp 2^-128, x < 2^-16), holding the series error near
   2^-122,

       P needs 90 bits, Q 74, R 58, S 42, U 26

   -- only P and Q need two limbs, while R, S and U (three of the five
   Horner steps) fit in ONE and are evaluated on x's top limb alone.
   The same argument at n = 1 leaves P needing just 34 bits, and at
   n = 3 (output ulp 2^-192, x < 2^-16) the eleven-term Horner splits

       V_2 154 bits, V_3 138, V_4 122, V_5 106, V_6 90, V_7 74,
       V_8 58, V_9 42, V_10 26

   so only V_2 and V_3 need three limbs, four need two and three need
   one -- leaving just four genuine 3-limb multiplies (x V_4, x V_3,
   x^2 and x^2 V_2), for which flint_mpn_mulhigh_n and
   flint_mpn_sqrhigh are used; everything above them stays in
   registers.  (An earlier attempt that ran the whole n = 3 Horner at
   three limbs was no faster than the generated code; the mixed form
   is.)

   The products keep only the top limbs: the 2x2 form drops the a0*b0
   partial entirely, multiplying by a one-limb value drops a column,
   and the squaring uses the symmetry of its two cross products.  The
   error budget -- FIXED_EXP_BITWISE_RS_MAX_ERR is thousands of ulp
   against an output of hundreds of bits -- absorbs all of this many
   times over; the measured error is a handful of ulp. */


#include "hand_mulhi.inc"

/* c_k = floor(2^64 / k!), one limb */
#define FIXED_EXP_E2 UWORD(0x8000000000000000)
#define FIXED_EXP_E3 UWORD(0x2aaaaaaaaaaaaaaa)
#define FIXED_EXP_E4 UWORD(0x0aaaaaaaaaaaaaaa)
#define FIXED_EXP_E5 UWORD(0x0222222222222222)
#define FIXED_EXP_E6 UWORD(0x005b05b05b05b05b)
#define FIXED_EXP_E7 UWORD(0x000d00d00d00d00d)

/* c_k = floor(2^128 / k!), two limbs */
#define FIXED_EXP_D2_HI UWORD(0x8000000000000000)
#define FIXED_EXP_D2_LO UWORD(0x0000000000000000)
#define FIXED_EXP_D3_HI UWORD(0x2aaaaaaaaaaaaaaa)
#define FIXED_EXP_D3_LO UWORD(0xaaaaaaaaaaaaaaaa)


/* exp(x), x < 2^-12, one fraction limb.  N = 4 terms: the dropped
   x^5/5! is about 2^-67, below the 2^-64 output ulp.  Written as
   1 + x + x^2 (C2 + x C3 + x^2 C4), reusing x^2 so that the two small
   terms are formed independently rather than in a serial chain. */
void
_fixed_exp_rs_opt_1(nn_ptr res, nn_srcptr x)
{
    ulong X = x[0], x2, p, h, l;

    x2 = n_mulhi(X, X);
    p = FIXED_EXP_E2 + n_mulhi(X, FIXED_EXP_E3)
                     + n_mulhi(x2, FIXED_EXP_E4);
    p = n_mulhi(x2, p);

    add_ssaaaa(h, l, UWORD(0), X, UWORD(0), p);
    res[0] = l;
    res[1] = 1 + h;
}

/* exp(x), x < 2^-16, two fraction limbs.  N = 7 terms: the dropped
   x^8/8! is about 2^-143, below the 2^-128 output ulp.

   R = C4 + x (C5 + x (C6 + x C7)) needs only 58 bits, so it and the
   two Horner steps below it run at a single limb, on x's top limb
   alone; Q = C3 + x R and P = C2 + x Q then need two. */
void
_fixed_exp_rs_opt_2(nn_ptr res, nn_srcptr x)
{
    ulong X1 = x[1], X0 = x[0];
    ulong u, s, r, q1, q0, p1, p0, h, l, x2h, x2l;

    /* inner tail at one limb */
    u = FIXED_EXP_E6 + n_mulhi(X1, FIXED_EXP_E7);
    s = FIXED_EXP_E5 + n_mulhi(X1, u);
    r = FIXED_EXP_E4 + n_mulhi(X1, s);

    /* Q = C3 + x R (R is a single limb) */
    _fixed_mulhi_2x1(&h, &l, X1, X0, r);
    add_ssaaaa(q1, q0, FIXED_EXP_D3_HI, FIXED_EXP_D3_LO, h, l);

    /* P = C2 + x Q */
    _fixed_mulhi_2x2_sloppy(&h, &l, X1, X0, q1, q0);
    add_ssaaaa(p1, p0, FIXED_EXP_D2_HI, FIXED_EXP_D2_LO, h, l);

    /* x^2 P, then the fraction x + x^2 P */
    _fixed_sqrhi_2x2(&x2h, &x2l, X1, X0);
    _fixed_mulhi_2x2_sloppy(&h, &l, x2h, x2l, p1, p0);
    add_ssaaaa(h, l, h, l, X1, X0);

    res[0] = l;
    res[1] = h;
    res[2] = 1;
}

/* 1-limb 1/k! for the innermost n = 3 Horner steps */
#define FIXED_EXP_G8  UWORD(0x0001a01a01a01a01)
#define FIXED_EXP_G9  UWORD(0x00002e3bc74aad8e)
#define FIXED_EXP_G10 UWORD(0x0000049f93edde27)

/* 2-limb 1/k! for the middle n = 3 Horner steps */
#define FIXED_EXP_H4_HI UWORD(0x0aaaaaaaaaaaaaaa)
#define FIXED_EXP_H4_LO UWORD(0xaaaaaaaaaaaaaaaa)
#define FIXED_EXP_H5_HI UWORD(0x0222222222222222)
#define FIXED_EXP_H5_LO UWORD(0x2222222222222222)
#define FIXED_EXP_H6_HI UWORD(0x005b05b05b05b05b)
#define FIXED_EXP_H6_LO UWORD(0x05b05b05b05b05b0)
#define FIXED_EXP_H7_HI UWORD(0x000d00d00d00d00d)
#define FIXED_EXP_H7_LO UWORD(0x00d00d00d00d00d0)

/* 1/3! to three limbs; 1/2! is just the top bit, added in directly */
static const ulong _fixed_exp_K3[3] = {
    UWORD(0xaaaaaaaaaaaaaaaa),
    UWORD(0xaaaaaaaaaaaaaaaa),
    UWORD(0x2aaaaaaaaaaaaaaa)
};

/* Top 3 limbs of a (3 limbs) times a 2-limb fraction b, which embeds
   as (0, b0, b1): the whole b-limb-0 column of the product vanishes,
   so three full products and two high-only ones suffice, against the
   six of a general 3x3. */
FLINT_FORCE_INLINE void
_fixed_mulhi_3x2(ulong * r2, ulong * r1, ulong * r0,
    ulong a2, ulong a1, ulong a0, ulong b1, ulong b0)
{
    ulong t2, t1, t0, ph, pl;

    umul_ppmm(t2, t1, a2, b1);          /* weight 4 */
    t0 = 0;
    umul_ppmm(ph, pl, a1, b1);          /* weight 3 */
    add_sssaaaaaa(t2, t1, t0, t2, t1, t0, UWORD(0), ph, pl);
    umul_ppmm(ph, pl, a2, b0);          /* weight 3 */
    add_sssaaaaaa(t2, t1, t0, t2, t1, t0, UWORD(0), ph, pl);
    ph = n_mulhi(a0, b1);               /* weight 2: high word only */
    add_sssaaaaaa(t2, t1, t0, t2, t1, t0, UWORD(0), UWORD(0), ph);
    ph = n_mulhi(a1, b0);               /* weight 2: high word only */
    add_sssaaaaaa(t2, t1, t0, t2, t1, t0, UWORD(0), UWORD(0), ph);

    *r0 = t0;
    *r1 = t1;
    *r2 = t2;
}

/* exp(x), x < 2^-16, three fraction limbs.  N = 11 terms: the dropped
   x^11/11! is about 2^-201, below the 2^-192 output ulp.  The Horner
   chain runs at one limb (V_10, V_9, V_8), then two (V_7 .. V_4), then
   three (V_3, V_2), as the error analysis above allows. */
void
_fixed_exp_rs_opt_3(nn_ptr res, nn_srcptr x)
{
    ulong X2 = x[2], X1 = x[1];
    ulong v, w1, w0, h, l, cy;
    ulong V3[3], V2[3], T[3], x2[3];

    /* one limb */
    v = FIXED_EXP_G10;
    v = FIXED_EXP_G9 + n_mulhi(X2, v);
    v = FIXED_EXP_G8 + n_mulhi(X2, v);

    /* two limbs (the first step multiplies by the one-limb V_8) */
    _fixed_mulhi_2x1(&h, &l, X2, X1, v);
    add_ssaaaa(w1, w0, FIXED_EXP_H7_HI, FIXED_EXP_H7_LO, h, l);
    _fixed_mulhi_2x2_sloppy(&h, &l, X2, X1, w1, w0);
    add_ssaaaa(w1, w0, FIXED_EXP_H6_HI, FIXED_EXP_H6_LO, h, l);
    _fixed_mulhi_2x2_sloppy(&h, &l, X2, X1, w1, w0);
    add_ssaaaa(w1, w0, FIXED_EXP_H5_HI, FIXED_EXP_H5_LO, h, l);
    _fixed_mulhi_2x2_sloppy(&h, &l, X2, X1, w1, w0);
    add_ssaaaa(w1, w0, FIXED_EXP_H4_HI, FIXED_EXP_H4_LO, h, l);

    /* three limbs: V_3 = C3 + x V_4 (V_4 has no low limb), then
       V_2 = 1/2 + x V_3 */
    _fixed_mulhi_3x2(&V3[2], &V3[1], &V3[0], x[2], x[1], x[0], w1, w0);
    mpn_add_n(V3, V3, _fixed_exp_K3, 3);
    flint_mpn_mulhigh_n(V2, x, V3, 3);
    V2[2] += UWORD(0x8000000000000000);

    /* x^2 V_2, then the fraction x + x^2 V_2 */
    flint_mpn_sqrhigh(x2, x, 3);
    flint_mpn_mulhigh_n(T, x2, V2, 3);
    cy = mpn_add_n(res, x, T, 3);
    res[3] = 1 + cy;
}

/* exp(x), x < 2^-16, four fraction limbs.  N = 14 terms.  The level
   widths from the analysis above are

       V_2 4, V_3 4, V_4..V_7 3, V_8..V_11 2, V_12 V_13 1

   so only the two outermost Horner steps run at the full width.  With
   the series this cheap the reduction becomes the dominant cost, and
   r = 16 measured best -- where the generated series wanted r = 32. */
void
_fixed_exp_rs_opt_4(nn_ptr res, nn_srcptr x)
{
    ulong v, w1, w0, h, l, cy;
    ulong A[4], B[4], T[4], x2[4];
    ulong X2 = x[3], X1 = x[2];

    v = UWORD(0x00000000b092309d);

    /* V_12 = C_12 + x V_13  (1 x 1) */
    v = UWORD(0x00000008f76c77fc) + n_mulhi(X2, v);

    /* V_11 = C_11 + x V_12  (2 x 1) */
    _fixed_mulhi_2x1(&h, &l, X2, X1, v);
    add_ssaaaa(w1, w0, UWORD(0x0000006b99159fd5), UWORD(0x138e3f9d1f92e0df), h, l);

    /* V_10 = C_10 + x V_11  (2 x 2) */
    _fixed_mulhi_2x2_sloppy(&h, &l, X2, X1, w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x0000049f93edde27), UWORD(0xd71cbbc05b4fa999), h, l);

    /* V_9 = C_9 + x V_10  (2 x 2) */
    _fixed_mulhi_2x2_sloppy(&h, &l, X2, X1, w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x00002e3bc74aad8e), UWORD(0x671f5583911ca002), h, l);

    /* V_8 = C_8 + x V_9  (2 x 2) */
    _fixed_mulhi_2x2_sloppy(&h, &l, X2, X1, w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x0001a01a01a01a01), UWORD(0xa01a01a01a01a01a), h, l);

    /* V_7 = C_7 + x V_8  (3 x 2) */
    _fixed_mulhi_3x2(&A[2], &A[1], &A[0], x[3], x[2], x[1],
        w1, w0);
    B[0] = UWORD(0x0d00d00d00d00d00);
    B[1] = UWORD(0x00d00d00d00d00d0);
    B[2] = UWORD(0x000d00d00d00d00d);
    mpn_add_n(A, A, B, 3);

    /* V_6 = C_6 + x V_7  (3 x 3) */
    flint_mpn_mulhigh_n(T, x + 1, A, 3);
    B[0] = UWORD(0x5b05b05b05b05b05);
    B[1] = UWORD(0x05b05b05b05b05b0);
    B[2] = UWORD(0x005b05b05b05b05b);
    mpn_add_n(A, T, B, 3);

    /* V_5 = C_5 + x V_6  (3 x 3) */
    flint_mpn_mulhigh_n(T, x + 1, A, 3);
    B[0] = UWORD(0x2222222222222222);
    B[1] = UWORD(0x2222222222222222);
    B[2] = UWORD(0x0222222222222222);
    mpn_add_n(A, T, B, 3);

    /* V_4 = C_4 + x V_5  (3 x 3) */
    flint_mpn_mulhigh_n(T, x + 1, A, 3);
    B[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    B[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    B[2] = UWORD(0x0aaaaaaaaaaaaaaa);
    mpn_add_n(A, T, B, 3);

    /* V_3 = C_3 + x V_4  (4 x 3) */
    B[0] = UWORD(0);
    B[1] = A[0];
    B[2] = A[1];
    B[3] = A[2];
    flint_mpn_mulhigh_n(T, x, B, 4);
    B[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    B[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    B[2] = UWORD(0xaaaaaaaaaaaaaaaa);
    B[3] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_add_n(A, T, B, 4);

    /* V_2 = 1/2 + x V_3 (C_2 is just the top bit) */
    flint_mpn_mulhigh_n(A, x, A, 4);
    A[3] += UWORD(0x8000000000000000);

    /* x^2 V_2 */
    flint_mpn_sqrhigh(x2, x, 4);
    flint_mpn_mulhigh_n(T, x2, A, 4);
    cy = mpn_add_n(res, x, T, 4);
    res[4] = 1 + cy;
}
