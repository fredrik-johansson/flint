/* atanh, n = 1, r = 8: 3 levels, widths [1, 1, 1] */
static void
hs_atanh_1_8(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x2492492492492492);

    /* V_1 = c + z V_2   (1 x 1 limbs) */
    v = n_mulhi(z[0], v);
    v = UWORD(0x3333333333333333) + v;

    /* V_0 = c + z V_1   (1 x 1 limbs) */
    v = n_mulhi(z[0], v);
    v = UWORD(0x5555555555555555) + v;

    /* atanh = x + x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_add_n(res, x, V, 1);
}

/* atanh, n = 1, r = 9: 2 levels, widths [1, 1] */
static void
hs_atanh_1_9(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x3333333333333333);

    /* V_0 = c + z V_1   (1 x 1 limbs) */
    v = n_mulhi(z[0], v);
    v = UWORD(0x5555555555555555) + v;

    /* atanh = x + x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_add_n(res, x, V, 1);
}

/* atanh, n = 1, r = 10: 2 levels, widths [1, 1] */
static void
hs_atanh_1_10(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x3333333333333333);

    /* V_0 = c + z V_1   (1 x 1 limbs) */
    v = n_mulhi(z[0], v);
    v = UWORD(0x5555555555555555) + v;

    /* atanh = x + x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_add_n(res, x, V, 1);
}

/* atanh, n = 1, r = 11: 2 levels, widths [1, 1] */
static void
hs_atanh_1_11(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x3333333333333333);

    /* V_0 = c + z V_1   (1 x 1 limbs) */
    v = n_mulhi(z[0], v);
    v = UWORD(0x5555555555555555) + v;

    /* atanh = x + x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_add_n(res, x, V, 1);
}

/* atanh, n = 1, r = 12: 2 levels, widths [1, 1] */
static void
hs_atanh_1_12(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x3333333333333333);

    /* V_0 = c + z V_1   (1 x 1 limbs) */
    v = n_mulhi(z[0], v);
    v = UWORD(0x5555555555555555) + v;

    /* atanh = x + x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_add_n(res, x, V, 1);
}

/* atanh, n = 1, r = 13: 1 levels, widths [1] */
static void
hs_atanh_1_13(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_add_n(res, x, V, 1);
}

/* atanh, n = 1, r = 14: 1 levels, widths [1] */
static void
hs_atanh_1_14(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_add_n(res, x, V, 1);
}

/* atanh, n = 1, r = 15: 1 levels, widths [1] */
static void
hs_atanh_1_15(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_add_n(res, x, V, 1);
}

/* atanh, n = 1, r = 16: 1 levels, widths [1] */
static void
hs_atanh_1_16(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_add_n(res, x, V, 1);
}

/* atanh, n = 1, r = 17: 1 levels, widths [1] */
static void
hs_atanh_1_17(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_add_n(res, x, V, 1);
}

/* atanh, n = 1, r = 18: 1 levels, widths [1] */
static void
hs_atanh_1_18(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_add_n(res, x, V, 1);
}

/* atanh, n = 1, r = 19: 1 levels, widths [1] */
static void
hs_atanh_1_19(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_add_n(res, x, V, 1);
}

/* atanh, n = 1, r = 20: 1 levels, widths [1] */
static void
hs_atanh_1_20(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_add_n(res, x, V, 1);
}

/* atanh, n = 1, r = 21: 1 levels, widths [1] */
static void
hs_atanh_1_21(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_add_n(res, x, V, 1);
}

/* atanh, n = 1, r = 22: 1 levels, widths [1] */
static void
hs_atanh_1_22(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_add_n(res, x, V, 1);
}

/* atanh, n = 1, r = 23: 1 levels, widths [1] */
static void
hs_atanh_1_23(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_add_n(res, x, V, 1);
}

/* atanh, n = 1, r = 24: 1 levels, widths [1] */
static void
hs_atanh_1_24(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_add_n(res, x, V, 1);
}

/* atanh, n = 1, r = 25: 1 levels, widths [1] */
static void
hs_atanh_1_25(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_add_n(res, x, V, 1);
}

/* atanh, n = 1, r = 26: 1 levels, widths [1] */
static void
hs_atanh_1_26(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_add_n(res, x, V, 1);
}

/* atanh, n = 1, r = 27: 1 levels, widths [1] */
static void
hs_atanh_1_27(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_add_n(res, x, V, 1);
}

/* atanh, n = 1, r = 28: 1 levels, widths [1] */
static void
hs_atanh_1_28(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_add_n(res, x, V, 1);
}

/* atanh, n = 1, r = 29: 1 levels, widths [1] */
static void
hs_atanh_1_29(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_add_n(res, x, V, 1);
}

/* atanh, n = 1, r = 30: 1 levels, widths [1] */
static void
hs_atanh_1_30(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_add_n(res, x, V, 1);
}

/* atanh, n = 1, r = 31: 1 levels, widths [1] */
static void
hs_atanh_1_31(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_add_n(res, x, V, 1);
}

/* atanh, n = 1, r = 32: 1 levels, widths [1] */
static void
hs_atanh_1_32(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_add_n(res, x, V, 1);
}

/* atanh, n = 1, r = 33: 1 levels, widths [1] */
static void
hs_atanh_1_33(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_add_n(res, x, V, 1);
}

/* atanh, n = 1, r = 34: 1 levels, widths [1] */
static void
hs_atanh_1_34(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_add_n(res, x, V, 1);
}

/* atanh, n = 1, r = 35: 1 levels, widths [1] */
static void
hs_atanh_1_35(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_add_n(res, x, V, 1);
}

/* atanh, n = 1, r = 36: 1 levels, widths [1] */
static void
hs_atanh_1_36(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_add_n(res, x, V, 1);
}

/* atanh, n = 1, r = 37: 1 levels, widths [1] */
static void
hs_atanh_1_37(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_add_n(res, x, V, 1);
}

/* atanh, n = 1, r = 38: 1 levels, widths [1] */
static void
hs_atanh_1_38(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_add_n(res, x, V, 1);
}

/* atanh, n = 1, r = 39: 1 levels, widths [1] */
static void
hs_atanh_1_39(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_add_n(res, x, V, 1);
}

/* atanh, n = 1, r = 40: 1 levels, widths [1] */
static void
hs_atanh_1_40(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_add_n(res, x, V, 1);
}

/* atanh, n = 1, r = 41: 1 levels, widths [1] */
static void
hs_atanh_1_41(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_add_n(res, x, V, 1);
}

/* atanh, n = 1, r = 42: 1 levels, widths [1] */
static void
hs_atanh_1_42(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_add_n(res, x, V, 1);
}

/* atanh, n = 1, r = 43: 1 levels, widths [1] */
static void
hs_atanh_1_43(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_add_n(res, x, V, 1);
}

/* atanh, n = 1, r = 44: 1 levels, widths [1] */
static void
hs_atanh_1_44(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_add_n(res, x, V, 1);
}

/* atanh, n = 1, r = 45: 1 levels, widths [1] */
static void
hs_atanh_1_45(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_add_n(res, x, V, 1);
}

/* atanh, n = 1, r = 46: 1 levels, widths [1] */
static void
hs_atanh_1_46(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_add_n(res, x, V, 1);
}

/* atanh, n = 1, r = 47: 1 levels, widths [1] */
static void
hs_atanh_1_47(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_add_n(res, x, V, 1);
}

/* atanh, n = 1, r = 48: 1 levels, widths [1] */
static void
hs_atanh_1_48(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_add_n(res, x, V, 1);
}

/* atanh, n = 2, r = 8: 7 levels, widths [2, 2, 2, 1, 1, 1, 1] */
static void
hs_atanh_2_8(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x1111111111111111);

    /* V_5 = c + z V_6   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x13b13b13b13b13b1) + v;

    /* V_4 = c + z V_5   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x1745d1745d1745d1) + v;

    /* V_3 = c + z V_4   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x1c71c71c71c71c71) + v;

    /* V_2 = c + z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    add_ssaaaa(w1, w0, UWORD(0x2492492492492492), UWORD(0x4924924924924924), h, l);

    /* V_1 = c + z V_2   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[1], z[0], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c + z V_1   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[1], z[0], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atanh = x + x^3 V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_add_n(res, x, V, 2);
}

/* atanh, n = 2, r = 9: 6 levels, widths [2, 2, 1, 1, 1, 1] */
static void
hs_atanh_2_9(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x13b13b13b13b13b1);

    /* V_4 = c + z V_5   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x1745d1745d1745d1) + v;

    /* V_3 = c + z V_4   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x1c71c71c71c71c71) + v;

    /* V_2 = c + z V_3   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x2492492492492492) + v;

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    add_ssaaaa(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c + z V_1   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[1], z[0], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atanh = x + x^3 V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_add_n(res, x, V, 2);
}

/* atanh, n = 2, r = 10: 5 levels, widths [2, 2, 1, 1, 1] */
static void
hs_atanh_2_10(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x1745d1745d1745d1);

    /* V_3 = c + z V_4   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x1c71c71c71c71c71) + v;

    /* V_2 = c + z V_3   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x2492492492492492) + v;

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    add_ssaaaa(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c + z V_1   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[1], z[0], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atanh = x + x^3 V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_add_n(res, x, V, 2);
}

/* atanh, n = 2, r = 11: 5 levels, widths [2, 2, 1, 1, 1] */
static void
hs_atanh_2_11(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x1745d1745d1745d1);

    /* V_3 = c + z V_4   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x1c71c71c71c71c71) + v;

    /* V_2 = c + z V_3   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x2492492492492492) + v;

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    add_ssaaaa(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c + z V_1   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[1], z[0], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atanh = x + x^3 V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_add_n(res, x, V, 2);
}

/* atanh, n = 2, r = 12: 4 levels, widths [2, 1, 1, 1] */
static void
hs_atanh_2_12(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x1c71c71c71c71c71);

    /* V_2 = c + z V_3   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x2492492492492492) + v;

    /* V_1 = c + z V_2   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x3333333333333333) + v;

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    add_ssaaaa(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atanh = x + x^3 V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_add_n(res, x, V, 2);
}

/* atanh, n = 2, r = 13: 4 levels, widths [2, 1, 1, 1] */
static void
hs_atanh_2_13(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x1c71c71c71c71c71);

    /* V_2 = c + z V_3   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x2492492492492492) + v;

    /* V_1 = c + z V_2   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x3333333333333333) + v;

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    add_ssaaaa(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atanh = x + x^3 V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_add_n(res, x, V, 2);
}

/* atanh, n = 2, r = 14: 3 levels, widths [2, 1, 1] */
static void
hs_atanh_2_14(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x2492492492492492);

    /* V_1 = c + z V_2   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x3333333333333333) + v;

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    add_ssaaaa(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atanh = x + x^3 V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_add_n(res, x, V, 2);
}

/* atanh, n = 2, r = 15: 3 levels, widths [2, 1, 1] */
static void
hs_atanh_2_15(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x2492492492492492);

    /* V_1 = c + z V_2   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x3333333333333333) + v;

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    add_ssaaaa(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atanh = x + x^3 V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_add_n(res, x, V, 2);
}

/* atanh, n = 2, r = 16: 3 levels, widths [2, 1, 1] */
static void
hs_atanh_2_16(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x2492492492492492);

    /* V_1 = c + z V_2   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x3333333333333333) + v;

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    add_ssaaaa(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atanh = x + x^3 V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_add_n(res, x, V, 2);
}

/* atanh, n = 2, r = 17: 3 levels, widths [2, 1, 1] */
static void
hs_atanh_2_17(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x2492492492492492);

    /* V_1 = c + z V_2   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x3333333333333333) + v;

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    add_ssaaaa(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atanh = x + x^3 V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_add_n(res, x, V, 2);
}

/* atanh, n = 2, r = 18: 2 levels, widths [2, 1] */
static void
hs_atanh_2_18(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x3333333333333333);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    add_ssaaaa(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atanh = x + x^3 V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_add_n(res, x, V, 2);
}

/* atanh, n = 2, r = 19: 2 levels, widths [2, 1] */
static void
hs_atanh_2_19(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x3333333333333333);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    add_ssaaaa(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atanh = x + x^3 V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_add_n(res, x, V, 2);
}

/* atanh, n = 2, r = 20: 2 levels, widths [1, 1] */
static void
hs_atanh_2_20(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x3333333333333333);

    /* V_0 = c + z V_1   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x5555555555555555) + v;

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_add_n(res, x, V, 2);
}

/* atanh, n = 2, r = 21: 2 levels, widths [1, 1] */
static void
hs_atanh_2_21(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x3333333333333333);

    /* V_0 = c + z V_1   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x5555555555555555) + v;

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_add_n(res, x, V, 2);
}

/* atanh, n = 2, r = 22: 2 levels, widths [1, 1] */
static void
hs_atanh_2_22(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x3333333333333333);

    /* V_0 = c + z V_1   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x5555555555555555) + v;

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_add_n(res, x, V, 2);
}

/* atanh, n = 2, r = 23: 2 levels, widths [1, 1] */
static void
hs_atanh_2_23(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x3333333333333333);

    /* V_0 = c + z V_1   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x5555555555555555) + v;

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_add_n(res, x, V, 2);
}

/* atanh, n = 2, r = 24: 2 levels, widths [1, 1] */
static void
hs_atanh_2_24(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x3333333333333333);

    /* V_0 = c + z V_1   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x5555555555555555) + v;

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_add_n(res, x, V, 2);
}

/* atanh, n = 2, r = 25: 2 levels, widths [1, 1] */
static void
hs_atanh_2_25(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x3333333333333333);

    /* V_0 = c + z V_1   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x5555555555555555) + v;

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_add_n(res, x, V, 2);
}

/* atanh, n = 2, r = 26: 1 levels, widths [1] */
static void
hs_atanh_2_26(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_add_n(res, x, V, 2);
}

/* atanh, n = 2, r = 27: 1 levels, widths [1] */
static void
hs_atanh_2_27(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_add_n(res, x, V, 2);
}

/* atanh, n = 2, r = 28: 1 levels, widths [1] */
static void
hs_atanh_2_28(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_add_n(res, x, V, 2);
}

/* atanh, n = 2, r = 29: 1 levels, widths [1] */
static void
hs_atanh_2_29(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_add_n(res, x, V, 2);
}

/* atanh, n = 2, r = 30: 1 levels, widths [1] */
static void
hs_atanh_2_30(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_add_n(res, x, V, 2);
}

/* atanh, n = 2, r = 31: 1 levels, widths [1] */
static void
hs_atanh_2_31(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_add_n(res, x, V, 2);
}

/* atanh, n = 2, r = 32: 1 levels, widths [1] */
static void
hs_atanh_2_32(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_add_n(res, x, V, 2);
}

/* atanh, n = 2, r = 33: 1 levels, widths [1] */
static void
hs_atanh_2_33(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_add_n(res, x, V, 2);
}

/* atanh, n = 2, r = 34: 1 levels, widths [1] */
static void
hs_atanh_2_34(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_add_n(res, x, V, 2);
}

/* atanh, n = 2, r = 35: 1 levels, widths [1] */
static void
hs_atanh_2_35(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_add_n(res, x, V, 2);
}

/* atanh, n = 2, r = 36: 1 levels, widths [1] */
static void
hs_atanh_2_36(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_add_n(res, x, V, 2);
}

/* atanh, n = 2, r = 37: 1 levels, widths [1] */
static void
hs_atanh_2_37(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_add_n(res, x, V, 2);
}

/* atanh, n = 2, r = 38: 1 levels, widths [1] */
static void
hs_atanh_2_38(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_add_n(res, x, V, 2);
}

/* atanh, n = 2, r = 39: 1 levels, widths [1] */
static void
hs_atanh_2_39(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_add_n(res, x, V, 2);
}

/* atanh, n = 2, r = 40: 1 levels, widths [1] */
static void
hs_atanh_2_40(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_add_n(res, x, V, 2);
}

/* atanh, n = 2, r = 41: 1 levels, widths [1] */
static void
hs_atanh_2_41(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_add_n(res, x, V, 2);
}

/* atanh, n = 2, r = 42: 1 levels, widths [1] */
static void
hs_atanh_2_42(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_add_n(res, x, V, 2);
}

/* atanh, n = 2, r = 43: 1 levels, widths [1] */
static void
hs_atanh_2_43(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_add_n(res, x, V, 2);
}

/* atanh, n = 2, r = 44: 1 levels, widths [1] */
static void
hs_atanh_2_44(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_add_n(res, x, V, 2);
}

/* atanh, n = 2, r = 45: 1 levels, widths [1] */
static void
hs_atanh_2_45(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_add_n(res, x, V, 2);
}

/* atanh, n = 2, r = 46: 1 levels, widths [1] */
static void
hs_atanh_2_46(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_add_n(res, x, V, 2);
}

/* atanh, n = 2, r = 47: 1 levels, widths [1] */
static void
hs_atanh_2_47(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_add_n(res, x, V, 2);
}

/* atanh, n = 2, r = 48: 1 levels, widths [1] */
static void
hs_atanh_2_48(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_add_n(res, x, V, 2);
}

/* atanh, n = 3, r = 8: 11 levels, widths [3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1] */
static void
hs_atanh_3_8(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0b21642c8590b216);

    /* V_9 = c + z V_10   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x0c30c30c30c30c30) + v;

    /* V_8 = c + z V_9   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x0d79435e50d79435) + v;

    /* V_7 = c + z V_8   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x0f0f0f0f0f0f0f0f) + v;

    /* V_6 = c + z V_7   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x1111111111111111), UWORD(0x1111111111111111), h, l);

    /* V_5 = c + z V_6   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x13b13b13b13b13b1), UWORD(0x3b13b13b13b13b13), h, l);

    /* V_4 = c + z V_5   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x1745d1745d1745d1), UWORD(0x745d1745d1745d17), h, l);

    /* V_3 = c + z V_4   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x1c71c71c71c71c71), UWORD(0xc71c71c71c71c71c), h, l);

    /* V_2 = c + z V_3   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x9249249249249249);
    P[1] = UWORD(0x4924924924924924);
    P[2] = UWORD(0x2492492492492492);
    mpn_add_n(V, P, T, 3);

    /* V_1 = c + z V_2   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x3333333333333333);
    P[1] = UWORD(0x3333333333333333);
    P[2] = UWORD(0x3333333333333333);
    mpn_add_n(V, P, T, 3);

    /* V_0 = c + z V_1   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 3);

    /* atanh = x + x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_add_n(res, x, V, 3);
}

/* atanh, n = 3, r = 9: 9 levels, widths [3, 3, 2, 2, 2, 2, 1, 1, 1] */
static void
hs_atanh_3_9(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0d79435e50d79435);

    /* V_7 = c + z V_8   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x0f0f0f0f0f0f0f0f) + v;

    /* V_6 = c + z V_7   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x1111111111111111) + v;

    /* V_5 = c + z V_6   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x13b13b13b13b13b1), UWORD(0x3b13b13b13b13b13), h, l);

    /* V_4 = c + z V_5   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x1745d1745d1745d1), UWORD(0x745d1745d1745d17), h, l);

    /* V_3 = c + z V_4   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x1c71c71c71c71c71), UWORD(0xc71c71c71c71c71c), h, l);

    /* V_2 = c + z V_3   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x2492492492492492), UWORD(0x4924924924924924), h, l);

    /* V_1 = c + z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x3333333333333333);
    P[1] = UWORD(0x3333333333333333);
    P[2] = UWORD(0x3333333333333333);
    mpn_add_n(V, P, T, 3);

    /* V_0 = c + z V_1   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 3);

    /* atanh = x + x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_add_n(res, x, V, 3);
}

/* atanh, n = 3, r = 10: 8 levels, widths [3, 3, 2, 2, 2, 1, 1, 1] */
static void
hs_atanh_3_10(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0f0f0f0f0f0f0f0f);

    /* V_6 = c + z V_7   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x1111111111111111) + v;

    /* V_5 = c + z V_6   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x13b13b13b13b13b1) + v;

    /* V_4 = c + z V_5   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x1745d1745d1745d1), UWORD(0x745d1745d1745d17), h, l);

    /* V_3 = c + z V_4   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x1c71c71c71c71c71), UWORD(0xc71c71c71c71c71c), h, l);

    /* V_2 = c + z V_3   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x2492492492492492), UWORD(0x4924924924924924), h, l);

    /* V_1 = c + z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x3333333333333333);
    P[1] = UWORD(0x3333333333333333);
    P[2] = UWORD(0x3333333333333333);
    mpn_add_n(V, P, T, 3);

    /* V_0 = c + z V_1   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 3);

    /* atanh = x + x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_add_n(res, x, V, 3);
}

/* atanh, n = 3, r = 11: 8 levels, widths [3, 3, 2, 2, 2, 1, 1, 1] */
static void
hs_atanh_3_11(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0f0f0f0f0f0f0f0f);

    /* V_6 = c + z V_7   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x1111111111111111) + v;

    /* V_5 = c + z V_6   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x13b13b13b13b13b1) + v;

    /* V_4 = c + z V_5   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x1745d1745d1745d1), UWORD(0x745d1745d1745d17), h, l);

    /* V_3 = c + z V_4   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x1c71c71c71c71c71), UWORD(0xc71c71c71c71c71c), h, l);

    /* V_2 = c + z V_3   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x2492492492492492), UWORD(0x4924924924924924), h, l);

    /* V_1 = c + z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x3333333333333333);
    P[1] = UWORD(0x3333333333333333);
    P[2] = UWORD(0x3333333333333333);
    mpn_add_n(V, P, T, 3);

    /* V_0 = c + z V_1   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 3);

    /* atanh = x + x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_add_n(res, x, V, 3);
}

/* atanh, n = 3, r = 12: 7 levels, widths [3, 2, 2, 2, 1, 1, 1] */
static void
hs_atanh_3_12(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x1111111111111111);

    /* V_5 = c + z V_6   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x13b13b13b13b13b1) + v;

    /* V_4 = c + z V_5   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x1745d1745d1745d1) + v;

    /* V_3 = c + z V_4   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x1c71c71c71c71c71), UWORD(0xc71c71c71c71c71c), h, l);

    /* V_2 = c + z V_3   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x2492492492492492), UWORD(0x4924924924924924), h, l);

    /* V_1 = c + z V_2   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 3);

    /* atanh = x + x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_add_n(res, x, V, 3);
}

/* atanh, n = 3, r = 13: 6 levels, widths [3, 2, 2, 2, 1, 1] */
static void
hs_atanh_3_13(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x13b13b13b13b13b1);

    /* V_4 = c + z V_5   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x1745d1745d1745d1) + v;

    /* V_3 = c + z V_4   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x1c71c71c71c71c71), UWORD(0xc71c71c71c71c71c), h, l);

    /* V_2 = c + z V_3   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x2492492492492492), UWORD(0x4924924924924924), h, l);

    /* V_1 = c + z V_2   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 3);

    /* atanh = x + x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_add_n(res, x, V, 3);
}

/* atanh, n = 3, r = 14: 6 levels, widths [3, 2, 2, 1, 1, 1] */
static void
hs_atanh_3_14(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x13b13b13b13b13b1);

    /* V_4 = c + z V_5   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x1745d1745d1745d1) + v;

    /* V_3 = c + z V_4   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x1c71c71c71c71c71) + v;

    /* V_2 = c + z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x2492492492492492), UWORD(0x4924924924924924), h, l);

    /* V_1 = c + z V_2   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 3);

    /* atanh = x + x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_add_n(res, x, V, 3);
}

/* atanh, n = 3, r = 15: 5 levels, widths [3, 2, 2, 1, 1] */
static void
hs_atanh_3_15(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x1745d1745d1745d1);

    /* V_3 = c + z V_4   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x1c71c71c71c71c71) + v;

    /* V_2 = c + z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x2492492492492492), UWORD(0x4924924924924924), h, l);

    /* V_1 = c + z V_2   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 3);

    /* atanh = x + x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_add_n(res, x, V, 3);
}

/* atanh, n = 3, r = 16: 5 levels, widths [3, 2, 2, 1, 1] */
static void
hs_atanh_3_16(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x1745d1745d1745d1);

    /* V_3 = c + z V_4   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x1c71c71c71c71c71) + v;

    /* V_2 = c + z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x2492492492492492), UWORD(0x4924924924924924), h, l);

    /* V_1 = c + z V_2   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 3);

    /* atanh = x + x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_add_n(res, x, V, 3);
}

/* atanh, n = 3, r = 17: 5 levels, widths [3, 2, 2, 1, 1] */
static void
hs_atanh_3_17(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x1745d1745d1745d1);

    /* V_3 = c + z V_4   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x1c71c71c71c71c71) + v;

    /* V_2 = c + z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x2492492492492492), UWORD(0x4924924924924924), h, l);

    /* V_1 = c + z V_2   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 3);

    /* atanh = x + x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_add_n(res, x, V, 3);
}

/* atanh, n = 3, r = 18: 4 levels, widths [3, 2, 1, 1] */
static void
hs_atanh_3_18(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x1c71c71c71c71c71);

    /* V_2 = c + z V_3   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x2492492492492492) + v;

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 3);

    /* atanh = x + x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_add_n(res, x, V, 3);
}

/* atanh, n = 3, r = 19: 4 levels, widths [3, 2, 1, 1] */
static void
hs_atanh_3_19(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x1c71c71c71c71c71);

    /* V_2 = c + z V_3   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x2492492492492492) + v;

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 3);

    /* atanh = x + x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_add_n(res, x, V, 3);
}

/* atanh, n = 3, r = 20: 4 levels, widths [2, 2, 1, 1] */
static void
hs_atanh_3_20(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x1c71c71c71c71c71);

    /* V_2 = c + z V_3   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x2492492492492492) + v;

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c + z V_1   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_add_n(res, x, V, 3);
}

/* atanh, n = 3, r = 21: 3 levels, widths [2, 2, 1] */
static void
hs_atanh_3_21(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x2492492492492492);

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c + z V_1   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_add_n(res, x, V, 3);
}

/* atanh, n = 3, r = 22: 3 levels, widths [2, 2, 1] */
static void
hs_atanh_3_22(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x2492492492492492);

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c + z V_1   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_add_n(res, x, V, 3);
}

/* atanh, n = 3, r = 23: 3 levels, widths [2, 2, 1] */
static void
hs_atanh_3_23(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x2492492492492492);

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c + z V_1   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_add_n(res, x, V, 3);
}

/* atanh, n = 3, r = 24: 3 levels, widths [2, 2, 1] */
static void
hs_atanh_3_24(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x2492492492492492);

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c + z V_1   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_add_n(res, x, V, 3);
}

/* atanh, n = 3, r = 25: 3 levels, widths [2, 1, 1] */
static void
hs_atanh_3_25(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x2492492492492492);

    /* V_1 = c + z V_2   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x3333333333333333) + v;

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_add_n(res, x, V, 3);
}

/* atanh, n = 3, r = 26: 3 levels, widths [2, 1, 1] */
static void
hs_atanh_3_26(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x2492492492492492);

    /* V_1 = c + z V_2   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x3333333333333333) + v;

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_add_n(res, x, V, 3);
}

/* atanh, n = 3, r = 27: 3 levels, widths [2, 1, 1] */
static void
hs_atanh_3_27(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x2492492492492492);

    /* V_1 = c + z V_2   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x3333333333333333) + v;

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_add_n(res, x, V, 3);
}

/* atanh, n = 3, r = 28: 2 levels, widths [2, 1] */
static void
hs_atanh_3_28(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x3333333333333333);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_add_n(res, x, V, 3);
}

/* atanh, n = 3, r = 29: 2 levels, widths [2, 1] */
static void
hs_atanh_3_29(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x3333333333333333);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_add_n(res, x, V, 3);
}

/* atanh, n = 3, r = 30: 2 levels, widths [2, 1] */
static void
hs_atanh_3_30(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x3333333333333333);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_add_n(res, x, V, 3);
}

/* atanh, n = 3, r = 31: 2 levels, widths [2, 1] */
static void
hs_atanh_3_31(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x3333333333333333);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_add_n(res, x, V, 3);
}

/* atanh, n = 3, r = 32: 2 levels, widths [2, 1] */
static void
hs_atanh_3_32(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x3333333333333333);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_add_n(res, x, V, 3);
}

/* atanh, n = 3, r = 33: 2 levels, widths [2, 1] */
static void
hs_atanh_3_33(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x3333333333333333);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_add_n(res, x, V, 3);
}

/* atanh, n = 3, r = 34: 2 levels, widths [2, 1] */
static void
hs_atanh_3_34(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x3333333333333333);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_add_n(res, x, V, 3);
}

/* atanh, n = 3, r = 35: 2 levels, widths [2, 1] */
static void
hs_atanh_3_35(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x3333333333333333);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_add_n(res, x, V, 3);
}

/* atanh, n = 3, r = 36: 2 levels, widths [2, 1] */
static void
hs_atanh_3_36(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x3333333333333333);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_add_n(res, x, V, 3);
}

/* atanh, n = 3, r = 37: 2 levels, widths [2, 1] */
static void
hs_atanh_3_37(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x3333333333333333);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_add_n(res, x, V, 3);
}

/* atanh, n = 3, r = 38: 1 levels, widths [2] */
static void
hs_atanh_3_38(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    w1 = UWORD(0x5555555555555555); w0 = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_add_n(res, x, V, 3);
}

/* atanh, n = 3, r = 39: 1 levels, widths [2] */
static void
hs_atanh_3_39(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    w1 = UWORD(0x5555555555555555); w0 = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_add_n(res, x, V, 3);
}

/* atanh, n = 3, r = 40: 1 levels, widths [2] */
static void
hs_atanh_3_40(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    w1 = UWORD(0x5555555555555555); w0 = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_add_n(res, x, V, 3);
}

/* atanh, n = 3, r = 41: 1 levels, widths [1] */
static void
hs_atanh_3_41(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = v;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_add_n(res, x, V, 3);
}

/* atanh, n = 3, r = 42: 1 levels, widths [1] */
static void
hs_atanh_3_42(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = v;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_add_n(res, x, V, 3);
}

/* atanh, n = 3, r = 43: 1 levels, widths [1] */
static void
hs_atanh_3_43(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = v;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_add_n(res, x, V, 3);
}

/* atanh, n = 3, r = 44: 1 levels, widths [1] */
static void
hs_atanh_3_44(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = v;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_add_n(res, x, V, 3);
}

/* atanh, n = 3, r = 45: 1 levels, widths [1] */
static void
hs_atanh_3_45(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = v;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_add_n(res, x, V, 3);
}

/* atanh, n = 3, r = 46: 1 levels, widths [1] */
static void
hs_atanh_3_46(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = v;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_add_n(res, x, V, 3);
}

/* atanh, n = 3, r = 47: 1 levels, widths [1] */
static void
hs_atanh_3_47(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = v;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_add_n(res, x, V, 3);
}

/* atanh, n = 3, r = 48: 1 levels, widths [1] */
static void
hs_atanh_3_48(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x5555555555555555);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = v;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_add_n(res, x, V, 3);
}

/* atanh, n = 4, r = 8: 15 levels, widths [4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1] */
static void
hs_atanh_4_8(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0842108421084210);

    /* V_13 = c + z V_14   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x08d3dcb08d3dcb08) + v;

    /* V_12 = c + z V_13   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x097b425ed097b425) + v;

    /* V_11 = c + z V_12   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x0a3d70a3d70a3d70) + v;

    /* V_10 = c + z V_11   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x0b21642c8590b216), UWORD(0x42c8590b21642c85), h, l);

    /* V_9 = c + z V_10   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x0c30c30c30c30c30), UWORD(0xc30c30c30c30c30c), h, l);

    /* V_8 = c + z V_9   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x0d79435e50d79435), UWORD(0xe50d79435e50d794), h, l);

    /* V_7 = c + z V_8   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x0f0f0f0f0f0f0f0f), UWORD(0x0f0f0f0f0f0f0f0f), h, l);

    /* V_6 = c + z V_7   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x1111111111111111);
    P[1] = UWORD(0x1111111111111111);
    P[2] = UWORD(0x1111111111111111);
    mpn_add_n(V, P, T, 3);

    /* V_5 = c + z V_6   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xb13b13b13b13b13b);
    P[1] = UWORD(0x3b13b13b13b13b13);
    P[2] = UWORD(0x13b13b13b13b13b1);
    mpn_add_n(V, P, T, 3);

    /* V_4 = c + z V_5   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x45d1745d1745d174);
    P[1] = UWORD(0x745d1745d1745d17);
    P[2] = UWORD(0x1745d1745d1745d1);
    mpn_add_n(V, P, T, 3);

    /* V_3 = c + z V_4   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x71c71c71c71c71c7);
    P[1] = UWORD(0xc71c71c71c71c71c);
    P[2] = UWORD(0x1c71c71c71c71c71);
    mpn_add_n(V, P, T, 3);

    /* V_2 = c + z V_3   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x2492492492492492);
    P[1] = UWORD(0x9249249249249249);
    P[2] = UWORD(0x4924924924924924);
    P[3] = UWORD(0x2492492492492492);
    mpn_add_n(V, P, T, 4);

    /* V_1 = c + z V_2   (4 x 4 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x3333333333333333);
    P[1] = UWORD(0x3333333333333333);
    P[2] = UWORD(0x3333333333333333);
    P[3] = UWORD(0x3333333333333333);
    mpn_add_n(V, P, T, 4);

    /* V_0 = c + z V_1   (4 x 4 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    P[3] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 4);

    /* atanh = x + x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_add_n(res, x, V, 4);
}

/* atanh, n = 4, r = 9: 13 levels, widths [4, 4, 3, 3, 3, 3, 2, 2, 2, 1, 1, 1, 1] */
static void
hs_atanh_4_9(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x097b425ed097b425);

    /* V_11 = c + z V_12   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x0a3d70a3d70a3d70) + v;

    /* V_10 = c + z V_11   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x0b21642c8590b216) + v;

    /* V_9 = c + z V_10   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x0c30c30c30c30c30) + v;

    /* V_8 = c + z V_9   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x0d79435e50d79435), UWORD(0xe50d79435e50d794), h, l);

    /* V_7 = c + z V_8   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x0f0f0f0f0f0f0f0f), UWORD(0x0f0f0f0f0f0f0f0f), h, l);

    /* V_6 = c + z V_7   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x1111111111111111), UWORD(0x1111111111111111), h, l);

    /* V_5 = c + z V_6   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xb13b13b13b13b13b);
    P[1] = UWORD(0x3b13b13b13b13b13);
    P[2] = UWORD(0x13b13b13b13b13b1);
    mpn_add_n(V, P, T, 3);

    /* V_4 = c + z V_5   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x45d1745d1745d174);
    P[1] = UWORD(0x745d1745d1745d17);
    P[2] = UWORD(0x1745d1745d1745d1);
    mpn_add_n(V, P, T, 3);

    /* V_3 = c + z V_4   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x71c71c71c71c71c7);
    P[1] = UWORD(0xc71c71c71c71c71c);
    P[2] = UWORD(0x1c71c71c71c71c71);
    mpn_add_n(V, P, T, 3);

    /* V_2 = c + z V_3   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x9249249249249249);
    P[1] = UWORD(0x4924924924924924);
    P[2] = UWORD(0x2492492492492492);
    mpn_add_n(V, P, T, 3);

    /* V_1 = c + z V_2   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x3333333333333333);
    P[1] = UWORD(0x3333333333333333);
    P[2] = UWORD(0x3333333333333333);
    P[3] = UWORD(0x3333333333333333);
    mpn_add_n(V, P, T, 4);

    /* V_0 = c + z V_1   (4 x 4 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    P[3] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 4);

    /* atanh = x + x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_add_n(res, x, V, 4);
}

/* atanh, n = 4, r = 10: 12 levels, widths [4, 4, 3, 3, 3, 2, 2, 2, 1, 1, 1, 1] */
static void
hs_atanh_4_10(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0a3d70a3d70a3d70);

    /* V_10 = c + z V_11   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x0b21642c8590b216) + v;

    /* V_9 = c + z V_10   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x0c30c30c30c30c30) + v;

    /* V_8 = c + z V_9   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x0d79435e50d79435) + v;

    /* V_7 = c + z V_8   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x0f0f0f0f0f0f0f0f), UWORD(0x0f0f0f0f0f0f0f0f), h, l);

    /* V_6 = c + z V_7   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x1111111111111111), UWORD(0x1111111111111111), h, l);

    /* V_5 = c + z V_6   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x13b13b13b13b13b1), UWORD(0x3b13b13b13b13b13), h, l);

    /* V_4 = c + z V_5   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x45d1745d1745d174);
    P[1] = UWORD(0x745d1745d1745d17);
    P[2] = UWORD(0x1745d1745d1745d1);
    mpn_add_n(V, P, T, 3);

    /* V_3 = c + z V_4   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x71c71c71c71c71c7);
    P[1] = UWORD(0xc71c71c71c71c71c);
    P[2] = UWORD(0x1c71c71c71c71c71);
    mpn_add_n(V, P, T, 3);

    /* V_2 = c + z V_3   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x9249249249249249);
    P[1] = UWORD(0x4924924924924924);
    P[2] = UWORD(0x2492492492492492);
    mpn_add_n(V, P, T, 3);

    /* V_1 = c + z V_2   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x3333333333333333);
    P[1] = UWORD(0x3333333333333333);
    P[2] = UWORD(0x3333333333333333);
    P[3] = UWORD(0x3333333333333333);
    mpn_add_n(V, P, T, 4);

    /* V_0 = c + z V_1   (4 x 4 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    P[3] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 4);

    /* atanh = x + x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_add_n(res, x, V, 4);
}

/* atanh, n = 4, r = 11: 10 levels, widths [4, 4, 3, 3, 3, 2, 2, 1, 1, 1] */
static void
hs_atanh_4_11(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0c30c30c30c30c30);

    /* V_8 = c + z V_9   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x0d79435e50d79435) + v;

    /* V_7 = c + z V_8   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x0f0f0f0f0f0f0f0f) + v;

    /* V_6 = c + z V_7   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x1111111111111111), UWORD(0x1111111111111111), h, l);

    /* V_5 = c + z V_6   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x13b13b13b13b13b1), UWORD(0x3b13b13b13b13b13), h, l);

    /* V_4 = c + z V_5   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x45d1745d1745d174);
    P[1] = UWORD(0x745d1745d1745d17);
    P[2] = UWORD(0x1745d1745d1745d1);
    mpn_add_n(V, P, T, 3);

    /* V_3 = c + z V_4   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x71c71c71c71c71c7);
    P[1] = UWORD(0xc71c71c71c71c71c);
    P[2] = UWORD(0x1c71c71c71c71c71);
    mpn_add_n(V, P, T, 3);

    /* V_2 = c + z V_3   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x9249249249249249);
    P[1] = UWORD(0x4924924924924924);
    P[2] = UWORD(0x2492492492492492);
    mpn_add_n(V, P, T, 3);

    /* V_1 = c + z V_2   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x3333333333333333);
    P[1] = UWORD(0x3333333333333333);
    P[2] = UWORD(0x3333333333333333);
    P[3] = UWORD(0x3333333333333333);
    mpn_add_n(V, P, T, 4);

    /* V_0 = c + z V_1   (4 x 4 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    P[3] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 4);

    /* atanh = x + x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_add_n(res, x, V, 4);
}

/* atanh, n = 4, r = 12: 9 levels, widths [4, 3, 3, 3, 2, 2, 2, 1, 1] */
static void
hs_atanh_4_12(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0d79435e50d79435);

    /* V_7 = c + z V_8   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x0f0f0f0f0f0f0f0f) + v;

    /* V_6 = c + z V_7   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x1111111111111111), UWORD(0x1111111111111111), h, l);

    /* V_5 = c + z V_6   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x13b13b13b13b13b1), UWORD(0x3b13b13b13b13b13), h, l);

    /* V_4 = c + z V_5   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x1745d1745d1745d1), UWORD(0x745d1745d1745d17), h, l);

    /* V_3 = c + z V_4   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x71c71c71c71c71c7);
    P[1] = UWORD(0xc71c71c71c71c71c);
    P[2] = UWORD(0x1c71c71c71c71c71);
    mpn_add_n(V, P, T, 3);

    /* V_2 = c + z V_3   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x9249249249249249);
    P[1] = UWORD(0x4924924924924924);
    P[2] = UWORD(0x2492492492492492);
    mpn_add_n(V, P, T, 3);

    /* V_1 = c + z V_2   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x3333333333333333);
    P[1] = UWORD(0x3333333333333333);
    P[2] = UWORD(0x3333333333333333);
    mpn_add_n(V, P, T, 3);

    /* V_0 = c + z V_1   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    P[3] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 4);

    /* atanh = x + x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_add_n(res, x, V, 4);
}

/* atanh, n = 4, r = 13: 9 levels, widths [4, 3, 3, 3, 2, 2, 1, 1, 1] */
static void
hs_atanh_4_13(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0d79435e50d79435);

    /* V_7 = c + z V_8   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x0f0f0f0f0f0f0f0f) + v;

    /* V_6 = c + z V_7   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x1111111111111111) + v;

    /* V_5 = c + z V_6   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x13b13b13b13b13b1), UWORD(0x3b13b13b13b13b13), h, l);

    /* V_4 = c + z V_5   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x1745d1745d1745d1), UWORD(0x745d1745d1745d17), h, l);

    /* V_3 = c + z V_4   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x71c71c71c71c71c7);
    P[1] = UWORD(0xc71c71c71c71c71c);
    P[2] = UWORD(0x1c71c71c71c71c71);
    mpn_add_n(V, P, T, 3);

    /* V_2 = c + z V_3   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x9249249249249249);
    P[1] = UWORD(0x4924924924924924);
    P[2] = UWORD(0x2492492492492492);
    mpn_add_n(V, P, T, 3);

    /* V_1 = c + z V_2   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x3333333333333333);
    P[1] = UWORD(0x3333333333333333);
    P[2] = UWORD(0x3333333333333333);
    mpn_add_n(V, P, T, 3);

    /* V_0 = c + z V_1   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    P[3] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 4);

    /* atanh = x + x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_add_n(res, x, V, 4);
}

/* atanh, n = 4, r = 14: 8 levels, widths [4, 3, 3, 2, 2, 2, 1, 1] */
static void
hs_atanh_4_14(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0f0f0f0f0f0f0f0f);

    /* V_6 = c + z V_7   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x1111111111111111) + v;

    /* V_5 = c + z V_6   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x13b13b13b13b13b1), UWORD(0x3b13b13b13b13b13), h, l);

    /* V_4 = c + z V_5   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x1745d1745d1745d1), UWORD(0x745d1745d1745d17), h, l);

    /* V_3 = c + z V_4   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x1c71c71c71c71c71), UWORD(0xc71c71c71c71c71c), h, l);

    /* V_2 = c + z V_3   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x9249249249249249);
    P[1] = UWORD(0x4924924924924924);
    P[2] = UWORD(0x2492492492492492);
    mpn_add_n(V, P, T, 3);

    /* V_1 = c + z V_2   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x3333333333333333);
    P[1] = UWORD(0x3333333333333333);
    P[2] = UWORD(0x3333333333333333);
    mpn_add_n(V, P, T, 3);

    /* V_0 = c + z V_1   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    P[3] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 4);

    /* atanh = x + x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_add_n(res, x, V, 4);
}

/* atanh, n = 4, r = 15: 7 levels, widths [4, 3, 3, 2, 2, 1, 1] */
static void
hs_atanh_4_15(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x1111111111111111);

    /* V_5 = c + z V_6   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x13b13b13b13b13b1) + v;

    /* V_4 = c + z V_5   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x1745d1745d1745d1), UWORD(0x745d1745d1745d17), h, l);

    /* V_3 = c + z V_4   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x1c71c71c71c71c71), UWORD(0xc71c71c71c71c71c), h, l);

    /* V_2 = c + z V_3   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x9249249249249249);
    P[1] = UWORD(0x4924924924924924);
    P[2] = UWORD(0x2492492492492492);
    mpn_add_n(V, P, T, 3);

    /* V_1 = c + z V_2   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x3333333333333333);
    P[1] = UWORD(0x3333333333333333);
    P[2] = UWORD(0x3333333333333333);
    mpn_add_n(V, P, T, 3);

    /* V_0 = c + z V_1   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    P[3] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 4);

    /* atanh = x + x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_add_n(res, x, V, 4);
}

/* atanh, n = 4, r = 16: 7 levels, widths [4, 3, 3, 2, 2, 1, 1] */
static void
hs_atanh_4_16(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x1111111111111111);

    /* V_5 = c + z V_6   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x13b13b13b13b13b1) + v;

    /* V_4 = c + z V_5   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x1745d1745d1745d1), UWORD(0x745d1745d1745d17), h, l);

    /* V_3 = c + z V_4   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x1c71c71c71c71c71), UWORD(0xc71c71c71c71c71c), h, l);

    /* V_2 = c + z V_3   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x9249249249249249);
    P[1] = UWORD(0x4924924924924924);
    P[2] = UWORD(0x2492492492492492);
    mpn_add_n(V, P, T, 3);

    /* V_1 = c + z V_2   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x3333333333333333);
    P[1] = UWORD(0x3333333333333333);
    P[2] = UWORD(0x3333333333333333);
    mpn_add_n(V, P, T, 3);

    /* V_0 = c + z V_1   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    P[3] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 4);

    /* atanh = x + x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_add_n(res, x, V, 4);
}

/* atanh, n = 4, r = 17: 6 levels, widths [4, 3, 3, 2, 1, 1] */
static void
hs_atanh_4_17(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x13b13b13b13b13b1);

    /* V_4 = c + z V_5   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x1745d1745d1745d1) + v;

    /* V_3 = c + z V_4   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x1c71c71c71c71c71), UWORD(0xc71c71c71c71c71c), h, l);

    /* V_2 = c + z V_3   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x9249249249249249);
    P[1] = UWORD(0x4924924924924924);
    P[2] = UWORD(0x2492492492492492);
    mpn_add_n(V, P, T, 3);

    /* V_1 = c + z V_2   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x3333333333333333);
    P[1] = UWORD(0x3333333333333333);
    P[2] = UWORD(0x3333333333333333);
    mpn_add_n(V, P, T, 3);

    /* V_0 = c + z V_1   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    P[3] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 4);

    /* atanh = x + x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_add_n(res, x, V, 4);
}

/* atanh, n = 4, r = 18: 6 levels, widths [4, 3, 2, 2, 1, 1] */
static void
hs_atanh_4_18(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x13b13b13b13b13b1);

    /* V_4 = c + z V_5   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x1745d1745d1745d1) + v;

    /* V_3 = c + z V_4   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x1c71c71c71c71c71), UWORD(0xc71c71c71c71c71c), h, l);

    /* V_2 = c + z V_3   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x2492492492492492), UWORD(0x4924924924924924), h, l);

    /* V_1 = c + z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x3333333333333333);
    P[1] = UWORD(0x3333333333333333);
    P[2] = UWORD(0x3333333333333333);
    mpn_add_n(V, P, T, 3);

    /* V_0 = c + z V_1   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    P[3] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 4);

    /* atanh = x + x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_add_n(res, x, V, 4);
}

/* atanh, n = 4, r = 19: 6 levels, widths [4, 3, 2, 2, 1, 1] */
static void
hs_atanh_4_19(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x13b13b13b13b13b1);

    /* V_4 = c + z V_5   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x1745d1745d1745d1) + v;

    /* V_3 = c + z V_4   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x1c71c71c71c71c71), UWORD(0xc71c71c71c71c71c), h, l);

    /* V_2 = c + z V_3   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x2492492492492492), UWORD(0x4924924924924924), h, l);

    /* V_1 = c + z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x3333333333333333);
    P[1] = UWORD(0x3333333333333333);
    P[2] = UWORD(0x3333333333333333);
    mpn_add_n(V, P, T, 3);

    /* V_0 = c + z V_1   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    P[3] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 4);

    /* atanh = x + x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_add_n(res, x, V, 4);
}

/* atanh, n = 4, r = 20: 5 levels, widths [3, 3, 2, 2, 1] */
static void
hs_atanh_4_20(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x1745d1745d1745d1);

    /* V_3 = c + z V_4   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x1c71c71c71c71c71), UWORD(0xc71c71c71c71c71c), h, l);

    /* V_2 = c + z V_3   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x2492492492492492), UWORD(0x4924924924924924), h, l);

    /* V_1 = c + z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x3333333333333333);
    P[1] = UWORD(0x3333333333333333);
    P[2] = UWORD(0x3333333333333333);
    mpn_add_n(V, P, T, 3);

    /* V_0 = c + z V_1   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 3);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_add_n(res, x, V, 4);
}

/* atanh, n = 4, r = 21: 5 levels, widths [3, 3, 2, 1, 1] */
static void
hs_atanh_4_21(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x1745d1745d1745d1);

    /* V_3 = c + z V_4   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x1c71c71c71c71c71) + v;

    /* V_2 = c + z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x2492492492492492), UWORD(0x4924924924924924), h, l);

    /* V_1 = c + z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x3333333333333333);
    P[1] = UWORD(0x3333333333333333);
    P[2] = UWORD(0x3333333333333333);
    mpn_add_n(V, P, T, 3);

    /* V_0 = c + z V_1   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 3);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_add_n(res, x, V, 4);
}

/* atanh, n = 4, r = 22: 5 levels, widths [3, 3, 2, 1, 1] */
static void
hs_atanh_4_22(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x1745d1745d1745d1);

    /* V_3 = c + z V_4   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x1c71c71c71c71c71) + v;

    /* V_2 = c + z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x2492492492492492), UWORD(0x4924924924924924), h, l);

    /* V_1 = c + z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x3333333333333333);
    P[1] = UWORD(0x3333333333333333);
    P[2] = UWORD(0x3333333333333333);
    mpn_add_n(V, P, T, 3);

    /* V_0 = c + z V_1   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 3);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_add_n(res, x, V, 4);
}

/* atanh, n = 4, r = 23: 4 levels, widths [3, 3, 2, 1] */
static void
hs_atanh_4_23(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x1c71c71c71c71c71);

    /* V_2 = c + z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x2492492492492492), UWORD(0x4924924924924924), h, l);

    /* V_1 = c + z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x3333333333333333);
    P[1] = UWORD(0x3333333333333333);
    P[2] = UWORD(0x3333333333333333);
    mpn_add_n(V, P, T, 3);

    /* V_0 = c + z V_1   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 3);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_add_n(res, x, V, 4);
}

/* atanh, n = 4, r = 24: 4 levels, widths [3, 3, 2, 1] */
static void
hs_atanh_4_24(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x1c71c71c71c71c71);

    /* V_2 = c + z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x2492492492492492), UWORD(0x4924924924924924), h, l);

    /* V_1 = c + z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x3333333333333333);
    P[1] = UWORD(0x3333333333333333);
    P[2] = UWORD(0x3333333333333333);
    mpn_add_n(V, P, T, 3);

    /* V_0 = c + z V_1   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 3);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_add_n(res, x, V, 4);
}

/* atanh, n = 4, r = 25: 4 levels, widths [3, 2, 2, 1] */
static void
hs_atanh_4_25(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x1c71c71c71c71c71);

    /* V_2 = c + z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x2492492492492492), UWORD(0x4924924924924924), h, l);

    /* V_1 = c + z V_2   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 3);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_add_n(res, x, V, 4);
}

/* atanh, n = 4, r = 26: 4 levels, widths [3, 2, 2, 1] */
static void
hs_atanh_4_26(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x1c71c71c71c71c71);

    /* V_2 = c + z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x2492492492492492), UWORD(0x4924924924924924), h, l);

    /* V_1 = c + z V_2   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 3);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_add_n(res, x, V, 4);
}

/* atanh, n = 4, r = 27: 4 levels, widths [3, 2, 1, 1] */
static void
hs_atanh_4_27(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x1c71c71c71c71c71);

    /* V_2 = c + z V_3   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x2492492492492492) + v;

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 3);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_add_n(res, x, V, 4);
}

/* atanh, n = 4, r = 28: 4 levels, widths [3, 2, 1, 1] */
static void
hs_atanh_4_28(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x1c71c71c71c71c71);

    /* V_2 = c + z V_3   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x2492492492492492) + v;

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 3);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_add_n(res, x, V, 4);
}

/* atanh, n = 4, r = 29: 3 levels, widths [3, 2, 1] */
static void
hs_atanh_4_29(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x2492492492492492);

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 3);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_add_n(res, x, V, 4);
}

/* atanh, n = 4, r = 30: 3 levels, widths [3, 2, 1] */
static void
hs_atanh_4_30(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x2492492492492492);

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 3);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_add_n(res, x, V, 4);
}

/* atanh, n = 4, r = 31: 3 levels, widths [3, 2, 1] */
static void
hs_atanh_4_31(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x2492492492492492);

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 3);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_add_n(res, x, V, 4);
}

/* atanh, n = 4, r = 32: 3 levels, widths [3, 2, 1] */
static void
hs_atanh_4_32(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x2492492492492492);

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 3);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_add_n(res, x, V, 4);
}

/* atanh, n = 4, r = 33: 3 levels, widths [3, 2, 1] */
static void
hs_atanh_4_33(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x2492492492492492);

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 3);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_add_n(res, x, V, 4);
}

/* atanh, n = 4, r = 34: 3 levels, widths [3, 2, 1] */
static void
hs_atanh_4_34(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x2492492492492492);

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 3);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_add_n(res, x, V, 4);
}

/* atanh, n = 4, r = 35: 3 levels, widths [3, 2, 1] */
static void
hs_atanh_4_35(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x2492492492492492);

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 3);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_add_n(res, x, V, 4);
}

/* atanh, n = 4, r = 36: 3 levels, widths [3, 2, 1] */
static void
hs_atanh_4_36(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x2492492492492492);

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 3);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_add_n(res, x, V, 4);
}

/* atanh, n = 4, r = 37: 2 levels, widths [3, 2] */
static void
hs_atanh_4_37(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    w1 = UWORD(0x3333333333333333); w0 = UWORD(0x3333333333333333);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 3);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_add_n(res, x, V, 4);
}

/* atanh, n = 4, r = 38: 2 levels, widths [3, 1] */
static void
hs_atanh_4_38(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x3333333333333333);

    /* V_0 = c + z V_1   (3 x 1 limbs) */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = v;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 3);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_add_n(res, x, V, 4);
}

/* atanh, n = 4, r = 39: 2 levels, widths [3, 1] */
static void
hs_atanh_4_39(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x3333333333333333);

    /* V_0 = c + z V_1   (3 x 1 limbs) */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = v;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 3);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_add_n(res, x, V, 4);
}

/* atanh, n = 4, r = 40: 2 levels, widths [3, 1] */
static void
hs_atanh_4_40(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x3333333333333333);

    /* V_0 = c + z V_1   (3 x 1 limbs) */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = v;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_add_n(V, P, T, 3);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_add_n(res, x, V, 4);
}

/* atanh, n = 4, r = 41: 2 levels, widths [2, 1] */
static void
hs_atanh_4_41(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x3333333333333333);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = w0; P[3] = w1;
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_add_n(res, x, V, 4);
}

/* atanh, n = 4, r = 42: 2 levels, widths [2, 1] */
static void
hs_atanh_4_42(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x3333333333333333);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = w0; P[3] = w1;
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_add_n(res, x, V, 4);
}

/* atanh, n = 4, r = 43: 2 levels, widths [2, 1] */
static void
hs_atanh_4_43(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x3333333333333333);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = w0; P[3] = w1;
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_add_n(res, x, V, 4);
}

/* atanh, n = 4, r = 44: 2 levels, widths [2, 1] */
static void
hs_atanh_4_44(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x3333333333333333);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = w0; P[3] = w1;
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_add_n(res, x, V, 4);
}

/* atanh, n = 4, r = 45: 2 levels, widths [2, 1] */
static void
hs_atanh_4_45(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x3333333333333333);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = w0; P[3] = w1;
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_add_n(res, x, V, 4);
}

/* atanh, n = 4, r = 46: 2 levels, widths [2, 1] */
static void
hs_atanh_4_46(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x3333333333333333);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = w0; P[3] = w1;
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_add_n(res, x, V, 4);
}

/* atanh, n = 4, r = 47: 2 levels, widths [2, 1] */
static void
hs_atanh_4_47(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x3333333333333333);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = w0; P[3] = w1;
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_add_n(res, x, V, 4);
}

/* atanh, n = 4, r = 48: 2 levels, widths [2, 1] */
static void
hs_atanh_4_48(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x3333333333333333);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atanh = x + x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = w0; P[3] = w1;
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_add_n(res, x, V, 4);
}

/* atan, n = 1, r = 8: 3 levels, widths [1, 1, 1] */
static void
hs_atan_1_8(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x2492492492492492);

    /* V_1 = c - z V_2   (1 x 1 limbs) */
    v = n_mulhi(z[0], v);
    v = UWORD(0x3333333333333333) - v;

    /* V_0 = c - z V_1   (1 x 1 limbs) */
    v = n_mulhi(z[0], v);
    v = UWORD(0x5555555555555555) - v;

    /* atan = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* atan, n = 1, r = 9: 2 levels, widths [1, 1] */
static void
hs_atan_1_9(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x3333333333333333);

    /* V_0 = c - z V_1   (1 x 1 limbs) */
    v = n_mulhi(z[0], v);
    v = UWORD(0x5555555555555555) - v;

    /* atan = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* atan, n = 1, r = 10: 2 levels, widths [1, 1] */
static void
hs_atan_1_10(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x3333333333333333);

    /* V_0 = c - z V_1   (1 x 1 limbs) */
    v = n_mulhi(z[0], v);
    v = UWORD(0x5555555555555555) - v;

    /* atan = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* atan, n = 1, r = 11: 2 levels, widths [1, 1] */
static void
hs_atan_1_11(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x3333333333333333);

    /* V_0 = c - z V_1   (1 x 1 limbs) */
    v = n_mulhi(z[0], v);
    v = UWORD(0x5555555555555555) - v;

    /* atan = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* atan, n = 1, r = 12: 2 levels, widths [1, 1] */
static void
hs_atan_1_12(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x3333333333333333);

    /* V_0 = c - z V_1   (1 x 1 limbs) */
    v = n_mulhi(z[0], v);
    v = UWORD(0x5555555555555555) - v;

    /* atan = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* atan, n = 1, r = 13: 1 levels, widths [1] */
static void
hs_atan_1_13(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* atan, n = 1, r = 14: 1 levels, widths [1] */
static void
hs_atan_1_14(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* atan, n = 1, r = 15: 1 levels, widths [1] */
static void
hs_atan_1_15(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* atan, n = 1, r = 16: 1 levels, widths [1] */
static void
hs_atan_1_16(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* atan, n = 1, r = 17: 1 levels, widths [1] */
static void
hs_atan_1_17(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* atan, n = 1, r = 18: 1 levels, widths [1] */
static void
hs_atan_1_18(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* atan, n = 1, r = 19: 1 levels, widths [1] */
static void
hs_atan_1_19(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* atan, n = 1, r = 20: 1 levels, widths [1] */
static void
hs_atan_1_20(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* atan, n = 1, r = 21: 1 levels, widths [1] */
static void
hs_atan_1_21(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* atan, n = 1, r = 22: 1 levels, widths [1] */
static void
hs_atan_1_22(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* atan, n = 1, r = 23: 1 levels, widths [1] */
static void
hs_atan_1_23(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* atan, n = 1, r = 24: 1 levels, widths [1] */
static void
hs_atan_1_24(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* atan, n = 1, r = 25: 1 levels, widths [1] */
static void
hs_atan_1_25(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* atan, n = 1, r = 26: 1 levels, widths [1] */
static void
hs_atan_1_26(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* atan, n = 1, r = 27: 1 levels, widths [1] */
static void
hs_atan_1_27(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* atan, n = 1, r = 28: 1 levels, widths [1] */
static void
hs_atan_1_28(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* atan, n = 1, r = 29: 1 levels, widths [1] */
static void
hs_atan_1_29(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* atan, n = 1, r = 30: 1 levels, widths [1] */
static void
hs_atan_1_30(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* atan, n = 1, r = 31: 1 levels, widths [1] */
static void
hs_atan_1_31(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* atan, n = 1, r = 32: 1 levels, widths [1] */
static void
hs_atan_1_32(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* atan, n = 1, r = 33: 1 levels, widths [1] */
static void
hs_atan_1_33(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* atan, n = 1, r = 34: 1 levels, widths [1] */
static void
hs_atan_1_34(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* atan, n = 1, r = 35: 1 levels, widths [1] */
static void
hs_atan_1_35(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* atan, n = 1, r = 36: 1 levels, widths [1] */
static void
hs_atan_1_36(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* atan, n = 1, r = 37: 1 levels, widths [1] */
static void
hs_atan_1_37(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* atan, n = 1, r = 38: 1 levels, widths [1] */
static void
hs_atan_1_38(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* atan, n = 1, r = 39: 1 levels, widths [1] */
static void
hs_atan_1_39(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* atan, n = 1, r = 40: 1 levels, widths [1] */
static void
hs_atan_1_40(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* atan, n = 1, r = 41: 1 levels, widths [1] */
static void
hs_atan_1_41(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* atan, n = 1, r = 42: 1 levels, widths [1] */
static void
hs_atan_1_42(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* atan, n = 1, r = 43: 1 levels, widths [1] */
static void
hs_atan_1_43(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* atan, n = 1, r = 44: 1 levels, widths [1] */
static void
hs_atan_1_44(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* atan, n = 1, r = 45: 1 levels, widths [1] */
static void
hs_atan_1_45(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* atan, n = 1, r = 46: 1 levels, widths [1] */
static void
hs_atan_1_46(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* atan, n = 1, r = 47: 1 levels, widths [1] */
static void
hs_atan_1_47(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* atan, n = 1, r = 48: 1 levels, widths [1] */
static void
hs_atan_1_48(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* atan, n = 2, r = 8: 7 levels, widths [2, 2, 2, 1, 1, 1, 1] */
static void
hs_atan_2_8(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x1111111111111111);

    /* V_5 = c - z V_6   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x13b13b13b13b13b1) - v;

    /* V_4 = c - z V_5   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x1745d1745d1745d1) - v;

    /* V_3 = c - z V_4   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x1c71c71c71c71c71) - v;

    /* V_2 = c - z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    sub_ddmmss(w1, w0, UWORD(0x2492492492492492), UWORD(0x4924924924924924), h, l);

    /* V_1 = c - z V_2   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[1], z[0], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c - z V_1   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[1], z[0], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atan = x - x^3 V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* atan, n = 2, r = 9: 6 levels, widths [2, 2, 1, 1, 1, 1] */
static void
hs_atan_2_9(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x13b13b13b13b13b1);

    /* V_4 = c - z V_5   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x1745d1745d1745d1) - v;

    /* V_3 = c - z V_4   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x1c71c71c71c71c71) - v;

    /* V_2 = c - z V_3   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x2492492492492492) - v;

    /* V_1 = c - z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    sub_ddmmss(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c - z V_1   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[1], z[0], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atan = x - x^3 V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* atan, n = 2, r = 10: 5 levels, widths [2, 2, 1, 1, 1] */
static void
hs_atan_2_10(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x1745d1745d1745d1);

    /* V_3 = c - z V_4   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x1c71c71c71c71c71) - v;

    /* V_2 = c - z V_3   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x2492492492492492) - v;

    /* V_1 = c - z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    sub_ddmmss(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c - z V_1   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[1], z[0], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atan = x - x^3 V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* atan, n = 2, r = 11: 5 levels, widths [2, 2, 1, 1, 1] */
static void
hs_atan_2_11(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x1745d1745d1745d1);

    /* V_3 = c - z V_4   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x1c71c71c71c71c71) - v;

    /* V_2 = c - z V_3   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x2492492492492492) - v;

    /* V_1 = c - z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    sub_ddmmss(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c - z V_1   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[1], z[0], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atan = x - x^3 V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* atan, n = 2, r = 12: 4 levels, widths [2, 1, 1, 1] */
static void
hs_atan_2_12(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x1c71c71c71c71c71);

    /* V_2 = c - z V_3   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x2492492492492492) - v;

    /* V_1 = c - z V_2   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x3333333333333333) - v;

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    sub_ddmmss(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atan = x - x^3 V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* atan, n = 2, r = 13: 4 levels, widths [2, 1, 1, 1] */
static void
hs_atan_2_13(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x1c71c71c71c71c71);

    /* V_2 = c - z V_3   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x2492492492492492) - v;

    /* V_1 = c - z V_2   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x3333333333333333) - v;

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    sub_ddmmss(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atan = x - x^3 V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* atan, n = 2, r = 14: 3 levels, widths [2, 1, 1] */
static void
hs_atan_2_14(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x2492492492492492);

    /* V_1 = c - z V_2   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x3333333333333333) - v;

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    sub_ddmmss(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atan = x - x^3 V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* atan, n = 2, r = 15: 3 levels, widths [2, 1, 1] */
static void
hs_atan_2_15(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x2492492492492492);

    /* V_1 = c - z V_2   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x3333333333333333) - v;

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    sub_ddmmss(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atan = x - x^3 V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* atan, n = 2, r = 16: 3 levels, widths [2, 1, 1] */
static void
hs_atan_2_16(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x2492492492492492);

    /* V_1 = c - z V_2   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x3333333333333333) - v;

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    sub_ddmmss(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atan = x - x^3 V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* atan, n = 2, r = 17: 3 levels, widths [2, 1, 1] */
static void
hs_atan_2_17(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x2492492492492492);

    /* V_1 = c - z V_2   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x3333333333333333) - v;

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    sub_ddmmss(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atan = x - x^3 V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* atan, n = 2, r = 18: 2 levels, widths [2, 1] */
static void
hs_atan_2_18(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x3333333333333333);

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    sub_ddmmss(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atan = x - x^3 V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* atan, n = 2, r = 19: 2 levels, widths [2, 1] */
static void
hs_atan_2_19(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x3333333333333333);

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    sub_ddmmss(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atan = x - x^3 V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* atan, n = 2, r = 20: 2 levels, widths [1, 1] */
static void
hs_atan_2_20(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x3333333333333333);

    /* V_0 = c - z V_1   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x5555555555555555) - v;

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* atan, n = 2, r = 21: 2 levels, widths [1, 1] */
static void
hs_atan_2_21(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x3333333333333333);

    /* V_0 = c - z V_1   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x5555555555555555) - v;

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* atan, n = 2, r = 22: 2 levels, widths [1, 1] */
static void
hs_atan_2_22(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x3333333333333333);

    /* V_0 = c - z V_1   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x5555555555555555) - v;

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* atan, n = 2, r = 23: 2 levels, widths [1, 1] */
static void
hs_atan_2_23(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x3333333333333333);

    /* V_0 = c - z V_1   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x5555555555555555) - v;

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* atan, n = 2, r = 24: 2 levels, widths [1, 1] */
static void
hs_atan_2_24(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x3333333333333333);

    /* V_0 = c - z V_1   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x5555555555555555) - v;

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* atan, n = 2, r = 25: 2 levels, widths [1, 1] */
static void
hs_atan_2_25(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x3333333333333333);

    /* V_0 = c - z V_1   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x5555555555555555) - v;

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* atan, n = 2, r = 26: 1 levels, widths [1] */
static void
hs_atan_2_26(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* atan, n = 2, r = 27: 1 levels, widths [1] */
static void
hs_atan_2_27(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* atan, n = 2, r = 28: 1 levels, widths [1] */
static void
hs_atan_2_28(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* atan, n = 2, r = 29: 1 levels, widths [1] */
static void
hs_atan_2_29(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* atan, n = 2, r = 30: 1 levels, widths [1] */
static void
hs_atan_2_30(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* atan, n = 2, r = 31: 1 levels, widths [1] */
static void
hs_atan_2_31(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* atan, n = 2, r = 32: 1 levels, widths [1] */
static void
hs_atan_2_32(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* atan, n = 2, r = 33: 1 levels, widths [1] */
static void
hs_atan_2_33(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* atan, n = 2, r = 34: 1 levels, widths [1] */
static void
hs_atan_2_34(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* atan, n = 2, r = 35: 1 levels, widths [1] */
static void
hs_atan_2_35(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* atan, n = 2, r = 36: 1 levels, widths [1] */
static void
hs_atan_2_36(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* atan, n = 2, r = 37: 1 levels, widths [1] */
static void
hs_atan_2_37(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* atan, n = 2, r = 38: 1 levels, widths [1] */
static void
hs_atan_2_38(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* atan, n = 2, r = 39: 1 levels, widths [1] */
static void
hs_atan_2_39(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* atan, n = 2, r = 40: 1 levels, widths [1] */
static void
hs_atan_2_40(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* atan, n = 2, r = 41: 1 levels, widths [1] */
static void
hs_atan_2_41(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* atan, n = 2, r = 42: 1 levels, widths [1] */
static void
hs_atan_2_42(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* atan, n = 2, r = 43: 1 levels, widths [1] */
static void
hs_atan_2_43(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* atan, n = 2, r = 44: 1 levels, widths [1] */
static void
hs_atan_2_44(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* atan, n = 2, r = 45: 1 levels, widths [1] */
static void
hs_atan_2_45(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* atan, n = 2, r = 46: 1 levels, widths [1] */
static void
hs_atan_2_46(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* atan, n = 2, r = 47: 1 levels, widths [1] */
static void
hs_atan_2_47(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* atan, n = 2, r = 48: 1 levels, widths [1] */
static void
hs_atan_2_48(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* atan, n = 3, r = 8: 11 levels, widths [3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1] */
static void
hs_atan_3_8(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0b21642c8590b216);

    /* V_9 = c - z V_10   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x0c30c30c30c30c30) - v;

    /* V_8 = c - z V_9   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x0d79435e50d79435) - v;

    /* V_7 = c - z V_8   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x0f0f0f0f0f0f0f0f) - v;

    /* V_6 = c - z V_7   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x1111111111111111), UWORD(0x1111111111111111), h, l);

    /* V_5 = c - z V_6   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x13b13b13b13b13b1), UWORD(0x3b13b13b13b13b13), h, l);

    /* V_4 = c - z V_5   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x1745d1745d1745d1), UWORD(0x745d1745d1745d17), h, l);

    /* V_3 = c - z V_4   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x1c71c71c71c71c71), UWORD(0xc71c71c71c71c71c), h, l);

    /* V_2 = c - z V_3   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x9249249249249249);
    P[1] = UWORD(0x4924924924924924);
    P[2] = UWORD(0x2492492492492492);
    mpn_sub_n(V, P, T, 3);

    /* V_1 = c - z V_2   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x3333333333333333);
    P[1] = UWORD(0x3333333333333333);
    P[2] = UWORD(0x3333333333333333);
    mpn_sub_n(V, P, T, 3);

    /* V_0 = c - z V_1   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 3);

    /* atan = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* atan, n = 3, r = 9: 9 levels, widths [3, 3, 2, 2, 2, 2, 1, 1, 1] */
static void
hs_atan_3_9(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0d79435e50d79435);

    /* V_7 = c - z V_8   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x0f0f0f0f0f0f0f0f) - v;

    /* V_6 = c - z V_7   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x1111111111111111) - v;

    /* V_5 = c - z V_6   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x13b13b13b13b13b1), UWORD(0x3b13b13b13b13b13), h, l);

    /* V_4 = c - z V_5   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x1745d1745d1745d1), UWORD(0x745d1745d1745d17), h, l);

    /* V_3 = c - z V_4   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x1c71c71c71c71c71), UWORD(0xc71c71c71c71c71c), h, l);

    /* V_2 = c - z V_3   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x2492492492492492), UWORD(0x4924924924924924), h, l);

    /* V_1 = c - z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x3333333333333333);
    P[1] = UWORD(0x3333333333333333);
    P[2] = UWORD(0x3333333333333333);
    mpn_sub_n(V, P, T, 3);

    /* V_0 = c - z V_1   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 3);

    /* atan = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* atan, n = 3, r = 10: 8 levels, widths [3, 3, 2, 2, 2, 1, 1, 1] */
static void
hs_atan_3_10(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0f0f0f0f0f0f0f0f);

    /* V_6 = c - z V_7   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x1111111111111111) - v;

    /* V_5 = c - z V_6   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x13b13b13b13b13b1) - v;

    /* V_4 = c - z V_5   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x1745d1745d1745d1), UWORD(0x745d1745d1745d17), h, l);

    /* V_3 = c - z V_4   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x1c71c71c71c71c71), UWORD(0xc71c71c71c71c71c), h, l);

    /* V_2 = c - z V_3   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x2492492492492492), UWORD(0x4924924924924924), h, l);

    /* V_1 = c - z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x3333333333333333);
    P[1] = UWORD(0x3333333333333333);
    P[2] = UWORD(0x3333333333333333);
    mpn_sub_n(V, P, T, 3);

    /* V_0 = c - z V_1   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 3);

    /* atan = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* atan, n = 3, r = 11: 8 levels, widths [3, 3, 2, 2, 2, 1, 1, 1] */
static void
hs_atan_3_11(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0f0f0f0f0f0f0f0f);

    /* V_6 = c - z V_7   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x1111111111111111) - v;

    /* V_5 = c - z V_6   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x13b13b13b13b13b1) - v;

    /* V_4 = c - z V_5   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x1745d1745d1745d1), UWORD(0x745d1745d1745d17), h, l);

    /* V_3 = c - z V_4   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x1c71c71c71c71c71), UWORD(0xc71c71c71c71c71c), h, l);

    /* V_2 = c - z V_3   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x2492492492492492), UWORD(0x4924924924924924), h, l);

    /* V_1 = c - z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x3333333333333333);
    P[1] = UWORD(0x3333333333333333);
    P[2] = UWORD(0x3333333333333333);
    mpn_sub_n(V, P, T, 3);

    /* V_0 = c - z V_1   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 3);

    /* atan = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* atan, n = 3, r = 12: 7 levels, widths [3, 2, 2, 2, 1, 1, 1] */
static void
hs_atan_3_12(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x1111111111111111);

    /* V_5 = c - z V_6   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x13b13b13b13b13b1) - v;

    /* V_4 = c - z V_5   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x1745d1745d1745d1) - v;

    /* V_3 = c - z V_4   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x1c71c71c71c71c71), UWORD(0xc71c71c71c71c71c), h, l);

    /* V_2 = c - z V_3   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x2492492492492492), UWORD(0x4924924924924924), h, l);

    /* V_1 = c - z V_2   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c - z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 3);

    /* atan = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* atan, n = 3, r = 13: 6 levels, widths [3, 2, 2, 2, 1, 1] */
static void
hs_atan_3_13(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x13b13b13b13b13b1);

    /* V_4 = c - z V_5   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x1745d1745d1745d1) - v;

    /* V_3 = c - z V_4   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x1c71c71c71c71c71), UWORD(0xc71c71c71c71c71c), h, l);

    /* V_2 = c - z V_3   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x2492492492492492), UWORD(0x4924924924924924), h, l);

    /* V_1 = c - z V_2   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c - z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 3);

    /* atan = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* atan, n = 3, r = 14: 6 levels, widths [3, 2, 2, 1, 1, 1] */
static void
hs_atan_3_14(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x13b13b13b13b13b1);

    /* V_4 = c - z V_5   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x1745d1745d1745d1) - v;

    /* V_3 = c - z V_4   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x1c71c71c71c71c71) - v;

    /* V_2 = c - z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x2492492492492492), UWORD(0x4924924924924924), h, l);

    /* V_1 = c - z V_2   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c - z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 3);

    /* atan = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* atan, n = 3, r = 15: 5 levels, widths [3, 2, 2, 1, 1] */
static void
hs_atan_3_15(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x1745d1745d1745d1);

    /* V_3 = c - z V_4   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x1c71c71c71c71c71) - v;

    /* V_2 = c - z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x2492492492492492), UWORD(0x4924924924924924), h, l);

    /* V_1 = c - z V_2   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c - z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 3);

    /* atan = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* atan, n = 3, r = 16: 5 levels, widths [3, 2, 2, 1, 1] */
static void
hs_atan_3_16(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x1745d1745d1745d1);

    /* V_3 = c - z V_4   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x1c71c71c71c71c71) - v;

    /* V_2 = c - z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x2492492492492492), UWORD(0x4924924924924924), h, l);

    /* V_1 = c - z V_2   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c - z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 3);

    /* atan = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* atan, n = 3, r = 17: 5 levels, widths [3, 2, 2, 1, 1] */
static void
hs_atan_3_17(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x1745d1745d1745d1);

    /* V_3 = c - z V_4   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x1c71c71c71c71c71) - v;

    /* V_2 = c - z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x2492492492492492), UWORD(0x4924924924924924), h, l);

    /* V_1 = c - z V_2   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c - z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 3);

    /* atan = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* atan, n = 3, r = 18: 4 levels, widths [3, 2, 1, 1] */
static void
hs_atan_3_18(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x1c71c71c71c71c71);

    /* V_2 = c - z V_3   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x2492492492492492) - v;

    /* V_1 = c - z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c - z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 3);

    /* atan = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* atan, n = 3, r = 19: 4 levels, widths [3, 2, 1, 1] */
static void
hs_atan_3_19(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x1c71c71c71c71c71);

    /* V_2 = c - z V_3   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x2492492492492492) - v;

    /* V_1 = c - z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c - z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 3);

    /* atan = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* atan, n = 3, r = 20: 4 levels, widths [2, 2, 1, 1] */
static void
hs_atan_3_20(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x1c71c71c71c71c71);

    /* V_2 = c - z V_3   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x2492492492492492) - v;

    /* V_1 = c - z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c - z V_1   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* atan, n = 3, r = 21: 3 levels, widths [2, 2, 1] */
static void
hs_atan_3_21(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x2492492492492492);

    /* V_1 = c - z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c - z V_1   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* atan, n = 3, r = 22: 3 levels, widths [2, 2, 1] */
static void
hs_atan_3_22(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x2492492492492492);

    /* V_1 = c - z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c - z V_1   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* atan, n = 3, r = 23: 3 levels, widths [2, 2, 1] */
static void
hs_atan_3_23(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x2492492492492492);

    /* V_1 = c - z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c - z V_1   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* atan, n = 3, r = 24: 3 levels, widths [2, 2, 1] */
static void
hs_atan_3_24(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x2492492492492492);

    /* V_1 = c - z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c - z V_1   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* atan, n = 3, r = 25: 3 levels, widths [2, 1, 1] */
static void
hs_atan_3_25(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x2492492492492492);

    /* V_1 = c - z V_2   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x3333333333333333) - v;

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* atan, n = 3, r = 26: 3 levels, widths [2, 1, 1] */
static void
hs_atan_3_26(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x2492492492492492);

    /* V_1 = c - z V_2   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x3333333333333333) - v;

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* atan, n = 3, r = 27: 3 levels, widths [2, 1, 1] */
static void
hs_atan_3_27(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x2492492492492492);

    /* V_1 = c - z V_2   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x3333333333333333) - v;

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* atan, n = 3, r = 28: 2 levels, widths [2, 1] */
static void
hs_atan_3_28(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x3333333333333333);

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* atan, n = 3, r = 29: 2 levels, widths [2, 1] */
static void
hs_atan_3_29(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x3333333333333333);

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* atan, n = 3, r = 30: 2 levels, widths [2, 1] */
static void
hs_atan_3_30(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x3333333333333333);

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* atan, n = 3, r = 31: 2 levels, widths [2, 1] */
static void
hs_atan_3_31(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x3333333333333333);

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* atan, n = 3, r = 32: 2 levels, widths [2, 1] */
static void
hs_atan_3_32(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x3333333333333333);

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* atan, n = 3, r = 33: 2 levels, widths [2, 1] */
static void
hs_atan_3_33(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x3333333333333333);

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* atan, n = 3, r = 34: 2 levels, widths [2, 1] */
static void
hs_atan_3_34(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x3333333333333333);

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* atan, n = 3, r = 35: 2 levels, widths [2, 1] */
static void
hs_atan_3_35(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x3333333333333333);

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* atan, n = 3, r = 36: 2 levels, widths [2, 1] */
static void
hs_atan_3_36(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x3333333333333333);

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* atan, n = 3, r = 37: 2 levels, widths [2, 1] */
static void
hs_atan_3_37(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x3333333333333333);

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* atan, n = 3, r = 38: 1 levels, widths [2] */
static void
hs_atan_3_38(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    w1 = UWORD(0x5555555555555555); w0 = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* atan, n = 3, r = 39: 1 levels, widths [2] */
static void
hs_atan_3_39(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    w1 = UWORD(0x5555555555555555); w0 = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* atan, n = 3, r = 40: 1 levels, widths [2] */
static void
hs_atan_3_40(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    w1 = UWORD(0x5555555555555555); w0 = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* atan, n = 3, r = 41: 1 levels, widths [1] */
static void
hs_atan_3_41(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = v;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* atan, n = 3, r = 42: 1 levels, widths [1] */
static void
hs_atan_3_42(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = v;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* atan, n = 3, r = 43: 1 levels, widths [1] */
static void
hs_atan_3_43(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = v;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* atan, n = 3, r = 44: 1 levels, widths [1] */
static void
hs_atan_3_44(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = v;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* atan, n = 3, r = 45: 1 levels, widths [1] */
static void
hs_atan_3_45(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = v;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* atan, n = 3, r = 46: 1 levels, widths [1] */
static void
hs_atan_3_46(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = v;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* atan, n = 3, r = 47: 1 levels, widths [1] */
static void
hs_atan_3_47(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = v;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* atan, n = 3, r = 48: 1 levels, widths [1] */
static void
hs_atan_3_48(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x5555555555555555);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = v;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* atan, n = 4, r = 8: 15 levels, widths [4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1] */
static void
hs_atan_4_8(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0842108421084210);

    /* V_13 = c - z V_14   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x08d3dcb08d3dcb08) - v;

    /* V_12 = c - z V_13   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x097b425ed097b425) - v;

    /* V_11 = c - z V_12   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x0a3d70a3d70a3d70) - v;

    /* V_10 = c - z V_11   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x0b21642c8590b216), UWORD(0x42c8590b21642c85), h, l);

    /* V_9 = c - z V_10   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x0c30c30c30c30c30), UWORD(0xc30c30c30c30c30c), h, l);

    /* V_8 = c - z V_9   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x0d79435e50d79435), UWORD(0xe50d79435e50d794), h, l);

    /* V_7 = c - z V_8   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x0f0f0f0f0f0f0f0f), UWORD(0x0f0f0f0f0f0f0f0f), h, l);

    /* V_6 = c - z V_7   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x1111111111111111);
    P[1] = UWORD(0x1111111111111111);
    P[2] = UWORD(0x1111111111111111);
    mpn_sub_n(V, P, T, 3);

    /* V_5 = c - z V_6   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xb13b13b13b13b13b);
    P[1] = UWORD(0x3b13b13b13b13b13);
    P[2] = UWORD(0x13b13b13b13b13b1);
    mpn_sub_n(V, P, T, 3);

    /* V_4 = c - z V_5   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x45d1745d1745d174);
    P[1] = UWORD(0x745d1745d1745d17);
    P[2] = UWORD(0x1745d1745d1745d1);
    mpn_sub_n(V, P, T, 3);

    /* V_3 = c - z V_4   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x71c71c71c71c71c7);
    P[1] = UWORD(0xc71c71c71c71c71c);
    P[2] = UWORD(0x1c71c71c71c71c71);
    mpn_sub_n(V, P, T, 3);

    /* V_2 = c - z V_3   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x2492492492492492);
    P[1] = UWORD(0x9249249249249249);
    P[2] = UWORD(0x4924924924924924);
    P[3] = UWORD(0x2492492492492492);
    mpn_sub_n(V, P, T, 4);

    /* V_1 = c - z V_2   (4 x 4 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x3333333333333333);
    P[1] = UWORD(0x3333333333333333);
    P[2] = UWORD(0x3333333333333333);
    P[3] = UWORD(0x3333333333333333);
    mpn_sub_n(V, P, T, 4);

    /* V_0 = c - z V_1   (4 x 4 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    P[3] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 4);

    /* atan = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* atan, n = 4, r = 9: 13 levels, widths [4, 4, 3, 3, 3, 3, 2, 2, 2, 1, 1, 1, 1] */
static void
hs_atan_4_9(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x097b425ed097b425);

    /* V_11 = c - z V_12   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x0a3d70a3d70a3d70) - v;

    /* V_10 = c - z V_11   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x0b21642c8590b216) - v;

    /* V_9 = c - z V_10   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x0c30c30c30c30c30) - v;

    /* V_8 = c - z V_9   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x0d79435e50d79435), UWORD(0xe50d79435e50d794), h, l);

    /* V_7 = c - z V_8   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x0f0f0f0f0f0f0f0f), UWORD(0x0f0f0f0f0f0f0f0f), h, l);

    /* V_6 = c - z V_7   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x1111111111111111), UWORD(0x1111111111111111), h, l);

    /* V_5 = c - z V_6   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xb13b13b13b13b13b);
    P[1] = UWORD(0x3b13b13b13b13b13);
    P[2] = UWORD(0x13b13b13b13b13b1);
    mpn_sub_n(V, P, T, 3);

    /* V_4 = c - z V_5   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x45d1745d1745d174);
    P[1] = UWORD(0x745d1745d1745d17);
    P[2] = UWORD(0x1745d1745d1745d1);
    mpn_sub_n(V, P, T, 3);

    /* V_3 = c - z V_4   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x71c71c71c71c71c7);
    P[1] = UWORD(0xc71c71c71c71c71c);
    P[2] = UWORD(0x1c71c71c71c71c71);
    mpn_sub_n(V, P, T, 3);

    /* V_2 = c - z V_3   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x9249249249249249);
    P[1] = UWORD(0x4924924924924924);
    P[2] = UWORD(0x2492492492492492);
    mpn_sub_n(V, P, T, 3);

    /* V_1 = c - z V_2   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x3333333333333333);
    P[1] = UWORD(0x3333333333333333);
    P[2] = UWORD(0x3333333333333333);
    P[3] = UWORD(0x3333333333333333);
    mpn_sub_n(V, P, T, 4);

    /* V_0 = c - z V_1   (4 x 4 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    P[3] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 4);

    /* atan = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* atan, n = 4, r = 10: 12 levels, widths [4, 4, 3, 3, 3, 2, 2, 2, 1, 1, 1, 1] */
static void
hs_atan_4_10(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0a3d70a3d70a3d70);

    /* V_10 = c - z V_11   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x0b21642c8590b216) - v;

    /* V_9 = c - z V_10   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x0c30c30c30c30c30) - v;

    /* V_8 = c - z V_9   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x0d79435e50d79435) - v;

    /* V_7 = c - z V_8   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x0f0f0f0f0f0f0f0f), UWORD(0x0f0f0f0f0f0f0f0f), h, l);

    /* V_6 = c - z V_7   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x1111111111111111), UWORD(0x1111111111111111), h, l);

    /* V_5 = c - z V_6   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x13b13b13b13b13b1), UWORD(0x3b13b13b13b13b13), h, l);

    /* V_4 = c - z V_5   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x45d1745d1745d174);
    P[1] = UWORD(0x745d1745d1745d17);
    P[2] = UWORD(0x1745d1745d1745d1);
    mpn_sub_n(V, P, T, 3);

    /* V_3 = c - z V_4   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x71c71c71c71c71c7);
    P[1] = UWORD(0xc71c71c71c71c71c);
    P[2] = UWORD(0x1c71c71c71c71c71);
    mpn_sub_n(V, P, T, 3);

    /* V_2 = c - z V_3   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x9249249249249249);
    P[1] = UWORD(0x4924924924924924);
    P[2] = UWORD(0x2492492492492492);
    mpn_sub_n(V, P, T, 3);

    /* V_1 = c - z V_2   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x3333333333333333);
    P[1] = UWORD(0x3333333333333333);
    P[2] = UWORD(0x3333333333333333);
    P[3] = UWORD(0x3333333333333333);
    mpn_sub_n(V, P, T, 4);

    /* V_0 = c - z V_1   (4 x 4 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    P[3] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 4);

    /* atan = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* atan, n = 4, r = 11: 10 levels, widths [4, 4, 3, 3, 3, 2, 2, 1, 1, 1] */
static void
hs_atan_4_11(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0c30c30c30c30c30);

    /* V_8 = c - z V_9   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x0d79435e50d79435) - v;

    /* V_7 = c - z V_8   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x0f0f0f0f0f0f0f0f) - v;

    /* V_6 = c - z V_7   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x1111111111111111), UWORD(0x1111111111111111), h, l);

    /* V_5 = c - z V_6   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x13b13b13b13b13b1), UWORD(0x3b13b13b13b13b13), h, l);

    /* V_4 = c - z V_5   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x45d1745d1745d174);
    P[1] = UWORD(0x745d1745d1745d17);
    P[2] = UWORD(0x1745d1745d1745d1);
    mpn_sub_n(V, P, T, 3);

    /* V_3 = c - z V_4   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x71c71c71c71c71c7);
    P[1] = UWORD(0xc71c71c71c71c71c);
    P[2] = UWORD(0x1c71c71c71c71c71);
    mpn_sub_n(V, P, T, 3);

    /* V_2 = c - z V_3   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x9249249249249249);
    P[1] = UWORD(0x4924924924924924);
    P[2] = UWORD(0x2492492492492492);
    mpn_sub_n(V, P, T, 3);

    /* V_1 = c - z V_2   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x3333333333333333);
    P[1] = UWORD(0x3333333333333333);
    P[2] = UWORD(0x3333333333333333);
    P[3] = UWORD(0x3333333333333333);
    mpn_sub_n(V, P, T, 4);

    /* V_0 = c - z V_1   (4 x 4 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    P[3] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 4);

    /* atan = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* atan, n = 4, r = 12: 9 levels, widths [4, 3, 3, 3, 2, 2, 2, 1, 1] */
static void
hs_atan_4_12(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0d79435e50d79435);

    /* V_7 = c - z V_8   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x0f0f0f0f0f0f0f0f) - v;

    /* V_6 = c - z V_7   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x1111111111111111), UWORD(0x1111111111111111), h, l);

    /* V_5 = c - z V_6   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x13b13b13b13b13b1), UWORD(0x3b13b13b13b13b13), h, l);

    /* V_4 = c - z V_5   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x1745d1745d1745d1), UWORD(0x745d1745d1745d17), h, l);

    /* V_3 = c - z V_4   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x71c71c71c71c71c7);
    P[1] = UWORD(0xc71c71c71c71c71c);
    P[2] = UWORD(0x1c71c71c71c71c71);
    mpn_sub_n(V, P, T, 3);

    /* V_2 = c - z V_3   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x9249249249249249);
    P[1] = UWORD(0x4924924924924924);
    P[2] = UWORD(0x2492492492492492);
    mpn_sub_n(V, P, T, 3);

    /* V_1 = c - z V_2   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x3333333333333333);
    P[1] = UWORD(0x3333333333333333);
    P[2] = UWORD(0x3333333333333333);
    mpn_sub_n(V, P, T, 3);

    /* V_0 = c - z V_1   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    P[3] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 4);

    /* atan = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* atan, n = 4, r = 13: 9 levels, widths [4, 3, 3, 3, 2, 2, 1, 1, 1] */
static void
hs_atan_4_13(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0d79435e50d79435);

    /* V_7 = c - z V_8   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x0f0f0f0f0f0f0f0f) - v;

    /* V_6 = c - z V_7   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x1111111111111111) - v;

    /* V_5 = c - z V_6   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x13b13b13b13b13b1), UWORD(0x3b13b13b13b13b13), h, l);

    /* V_4 = c - z V_5   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x1745d1745d1745d1), UWORD(0x745d1745d1745d17), h, l);

    /* V_3 = c - z V_4   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x71c71c71c71c71c7);
    P[1] = UWORD(0xc71c71c71c71c71c);
    P[2] = UWORD(0x1c71c71c71c71c71);
    mpn_sub_n(V, P, T, 3);

    /* V_2 = c - z V_3   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x9249249249249249);
    P[1] = UWORD(0x4924924924924924);
    P[2] = UWORD(0x2492492492492492);
    mpn_sub_n(V, P, T, 3);

    /* V_1 = c - z V_2   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x3333333333333333);
    P[1] = UWORD(0x3333333333333333);
    P[2] = UWORD(0x3333333333333333);
    mpn_sub_n(V, P, T, 3);

    /* V_0 = c - z V_1   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    P[3] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 4);

    /* atan = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* atan, n = 4, r = 14: 8 levels, widths [4, 3, 3, 2, 2, 2, 1, 1] */
static void
hs_atan_4_14(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0f0f0f0f0f0f0f0f);

    /* V_6 = c - z V_7   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x1111111111111111) - v;

    /* V_5 = c - z V_6   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x13b13b13b13b13b1), UWORD(0x3b13b13b13b13b13), h, l);

    /* V_4 = c - z V_5   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x1745d1745d1745d1), UWORD(0x745d1745d1745d17), h, l);

    /* V_3 = c - z V_4   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x1c71c71c71c71c71), UWORD(0xc71c71c71c71c71c), h, l);

    /* V_2 = c - z V_3   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x9249249249249249);
    P[1] = UWORD(0x4924924924924924);
    P[2] = UWORD(0x2492492492492492);
    mpn_sub_n(V, P, T, 3);

    /* V_1 = c - z V_2   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x3333333333333333);
    P[1] = UWORD(0x3333333333333333);
    P[2] = UWORD(0x3333333333333333);
    mpn_sub_n(V, P, T, 3);

    /* V_0 = c - z V_1   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    P[3] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 4);

    /* atan = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* atan, n = 4, r = 15: 7 levels, widths [4, 3, 3, 2, 2, 1, 1] */
static void
hs_atan_4_15(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x1111111111111111);

    /* V_5 = c - z V_6   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x13b13b13b13b13b1) - v;

    /* V_4 = c - z V_5   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x1745d1745d1745d1), UWORD(0x745d1745d1745d17), h, l);

    /* V_3 = c - z V_4   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x1c71c71c71c71c71), UWORD(0xc71c71c71c71c71c), h, l);

    /* V_2 = c - z V_3   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x9249249249249249);
    P[1] = UWORD(0x4924924924924924);
    P[2] = UWORD(0x2492492492492492);
    mpn_sub_n(V, P, T, 3);

    /* V_1 = c - z V_2   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x3333333333333333);
    P[1] = UWORD(0x3333333333333333);
    P[2] = UWORD(0x3333333333333333);
    mpn_sub_n(V, P, T, 3);

    /* V_0 = c - z V_1   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    P[3] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 4);

    /* atan = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* atan, n = 4, r = 16: 7 levels, widths [4, 3, 3, 2, 2, 1, 1] */
static void
hs_atan_4_16(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x1111111111111111);

    /* V_5 = c - z V_6   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x13b13b13b13b13b1) - v;

    /* V_4 = c - z V_5   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x1745d1745d1745d1), UWORD(0x745d1745d1745d17), h, l);

    /* V_3 = c - z V_4   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x1c71c71c71c71c71), UWORD(0xc71c71c71c71c71c), h, l);

    /* V_2 = c - z V_3   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x9249249249249249);
    P[1] = UWORD(0x4924924924924924);
    P[2] = UWORD(0x2492492492492492);
    mpn_sub_n(V, P, T, 3);

    /* V_1 = c - z V_2   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x3333333333333333);
    P[1] = UWORD(0x3333333333333333);
    P[2] = UWORD(0x3333333333333333);
    mpn_sub_n(V, P, T, 3);

    /* V_0 = c - z V_1   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    P[3] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 4);

    /* atan = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* atan, n = 4, r = 17: 6 levels, widths [4, 3, 3, 2, 1, 1] */
static void
hs_atan_4_17(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x13b13b13b13b13b1);

    /* V_4 = c - z V_5   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x1745d1745d1745d1) - v;

    /* V_3 = c - z V_4   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x1c71c71c71c71c71), UWORD(0xc71c71c71c71c71c), h, l);

    /* V_2 = c - z V_3   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x9249249249249249);
    P[1] = UWORD(0x4924924924924924);
    P[2] = UWORD(0x2492492492492492);
    mpn_sub_n(V, P, T, 3);

    /* V_1 = c - z V_2   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x3333333333333333);
    P[1] = UWORD(0x3333333333333333);
    P[2] = UWORD(0x3333333333333333);
    mpn_sub_n(V, P, T, 3);

    /* V_0 = c - z V_1   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    P[3] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 4);

    /* atan = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* atan, n = 4, r = 18: 6 levels, widths [4, 3, 2, 2, 1, 1] */
static void
hs_atan_4_18(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x13b13b13b13b13b1);

    /* V_4 = c - z V_5   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x1745d1745d1745d1) - v;

    /* V_3 = c - z V_4   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x1c71c71c71c71c71), UWORD(0xc71c71c71c71c71c), h, l);

    /* V_2 = c - z V_3   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x2492492492492492), UWORD(0x4924924924924924), h, l);

    /* V_1 = c - z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x3333333333333333);
    P[1] = UWORD(0x3333333333333333);
    P[2] = UWORD(0x3333333333333333);
    mpn_sub_n(V, P, T, 3);

    /* V_0 = c - z V_1   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    P[3] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 4);

    /* atan = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* atan, n = 4, r = 19: 6 levels, widths [4, 3, 2, 2, 1, 1] */
static void
hs_atan_4_19(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x13b13b13b13b13b1);

    /* V_4 = c - z V_5   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x1745d1745d1745d1) - v;

    /* V_3 = c - z V_4   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x1c71c71c71c71c71), UWORD(0xc71c71c71c71c71c), h, l);

    /* V_2 = c - z V_3   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x2492492492492492), UWORD(0x4924924924924924), h, l);

    /* V_1 = c - z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x3333333333333333);
    P[1] = UWORD(0x3333333333333333);
    P[2] = UWORD(0x3333333333333333);
    mpn_sub_n(V, P, T, 3);

    /* V_0 = c - z V_1   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    P[3] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 4);

    /* atan = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* atan, n = 4, r = 20: 5 levels, widths [3, 3, 2, 2, 1] */
static void
hs_atan_4_20(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x1745d1745d1745d1);

    /* V_3 = c - z V_4   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x1c71c71c71c71c71), UWORD(0xc71c71c71c71c71c), h, l);

    /* V_2 = c - z V_3   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x2492492492492492), UWORD(0x4924924924924924), h, l);

    /* V_1 = c - z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x3333333333333333);
    P[1] = UWORD(0x3333333333333333);
    P[2] = UWORD(0x3333333333333333);
    mpn_sub_n(V, P, T, 3);

    /* V_0 = c - z V_1   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 3);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* atan, n = 4, r = 21: 5 levels, widths [3, 3, 2, 1, 1] */
static void
hs_atan_4_21(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x1745d1745d1745d1);

    /* V_3 = c - z V_4   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x1c71c71c71c71c71) - v;

    /* V_2 = c - z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x2492492492492492), UWORD(0x4924924924924924), h, l);

    /* V_1 = c - z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x3333333333333333);
    P[1] = UWORD(0x3333333333333333);
    P[2] = UWORD(0x3333333333333333);
    mpn_sub_n(V, P, T, 3);

    /* V_0 = c - z V_1   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 3);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* atan, n = 4, r = 22: 5 levels, widths [3, 3, 2, 1, 1] */
static void
hs_atan_4_22(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x1745d1745d1745d1);

    /* V_3 = c - z V_4   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x1c71c71c71c71c71) - v;

    /* V_2 = c - z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x2492492492492492), UWORD(0x4924924924924924), h, l);

    /* V_1 = c - z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x3333333333333333);
    P[1] = UWORD(0x3333333333333333);
    P[2] = UWORD(0x3333333333333333);
    mpn_sub_n(V, P, T, 3);

    /* V_0 = c - z V_1   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 3);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* atan, n = 4, r = 23: 4 levels, widths [3, 3, 2, 1] */
static void
hs_atan_4_23(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x1c71c71c71c71c71);

    /* V_2 = c - z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x2492492492492492), UWORD(0x4924924924924924), h, l);

    /* V_1 = c - z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x3333333333333333);
    P[1] = UWORD(0x3333333333333333);
    P[2] = UWORD(0x3333333333333333);
    mpn_sub_n(V, P, T, 3);

    /* V_0 = c - z V_1   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 3);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* atan, n = 4, r = 24: 4 levels, widths [3, 3, 2, 1] */
static void
hs_atan_4_24(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x1c71c71c71c71c71);

    /* V_2 = c - z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x2492492492492492), UWORD(0x4924924924924924), h, l);

    /* V_1 = c - z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x3333333333333333);
    P[1] = UWORD(0x3333333333333333);
    P[2] = UWORD(0x3333333333333333);
    mpn_sub_n(V, P, T, 3);

    /* V_0 = c - z V_1   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 3);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* atan, n = 4, r = 25: 4 levels, widths [3, 2, 2, 1] */
static void
hs_atan_4_25(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x1c71c71c71c71c71);

    /* V_2 = c - z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x2492492492492492), UWORD(0x4924924924924924), h, l);

    /* V_1 = c - z V_2   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c - z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 3);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* atan, n = 4, r = 26: 4 levels, widths [3, 2, 2, 1] */
static void
hs_atan_4_26(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x1c71c71c71c71c71);

    /* V_2 = c - z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x2492492492492492), UWORD(0x4924924924924924), h, l);

    /* V_1 = c - z V_2   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c - z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 3);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* atan, n = 4, r = 27: 4 levels, widths [3, 2, 1, 1] */
static void
hs_atan_4_27(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x1c71c71c71c71c71);

    /* V_2 = c - z V_3   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x2492492492492492) - v;

    /* V_1 = c - z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c - z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 3);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* atan, n = 4, r = 28: 4 levels, widths [3, 2, 1, 1] */
static void
hs_atan_4_28(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x1c71c71c71c71c71);

    /* V_2 = c - z V_3   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x2492492492492492) - v;

    /* V_1 = c - z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c - z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 3);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* atan, n = 4, r = 29: 3 levels, widths [3, 2, 1] */
static void
hs_atan_4_29(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x2492492492492492);

    /* V_1 = c - z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c - z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 3);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* atan, n = 4, r = 30: 3 levels, widths [3, 2, 1] */
static void
hs_atan_4_30(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x2492492492492492);

    /* V_1 = c - z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c - z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 3);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* atan, n = 4, r = 31: 3 levels, widths [3, 2, 1] */
static void
hs_atan_4_31(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x2492492492492492);

    /* V_1 = c - z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c - z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 3);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* atan, n = 4, r = 32: 3 levels, widths [3, 2, 1] */
static void
hs_atan_4_32(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x2492492492492492);

    /* V_1 = c - z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c - z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 3);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* atan, n = 4, r = 33: 3 levels, widths [3, 2, 1] */
static void
hs_atan_4_33(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x2492492492492492);

    /* V_1 = c - z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c - z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 3);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* atan, n = 4, r = 34: 3 levels, widths [3, 2, 1] */
static void
hs_atan_4_34(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x2492492492492492);

    /* V_1 = c - z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c - z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 3);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* atan, n = 4, r = 35: 3 levels, widths [3, 2, 1] */
static void
hs_atan_4_35(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x2492492492492492);

    /* V_1 = c - z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c - z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 3);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* atan, n = 4, r = 36: 3 levels, widths [3, 2, 1] */
static void
hs_atan_4_36(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x2492492492492492);

    /* V_1 = c - z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x3333333333333333), UWORD(0x3333333333333333), h, l);

    /* V_0 = c - z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 3);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* atan, n = 4, r = 37: 2 levels, widths [3, 2] */
static void
hs_atan_4_37(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    w1 = UWORD(0x3333333333333333); w0 = UWORD(0x3333333333333333);

    /* V_0 = c - z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 3);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* atan, n = 4, r = 38: 2 levels, widths [3, 1] */
static void
hs_atan_4_38(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x3333333333333333);

    /* V_0 = c - z V_1   (3 x 1 limbs) */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = v;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 3);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* atan, n = 4, r = 39: 2 levels, widths [3, 1] */
static void
hs_atan_4_39(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x3333333333333333);

    /* V_0 = c - z V_1   (3 x 1 limbs) */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = v;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 3);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* atan, n = 4, r = 40: 2 levels, widths [3, 1] */
static void
hs_atan_4_40(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x3333333333333333);

    /* V_0 = c - z V_1   (3 x 1 limbs) */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = v;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5555555555555555);
    P[1] = UWORD(0x5555555555555555);
    P[2] = UWORD(0x5555555555555555);
    mpn_sub_n(V, P, T, 3);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* atan, n = 4, r = 41: 2 levels, widths [2, 1] */
static void
hs_atan_4_41(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x3333333333333333);

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = w0; P[3] = w1;
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* atan, n = 4, r = 42: 2 levels, widths [2, 1] */
static void
hs_atan_4_42(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x3333333333333333);

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = w0; P[3] = w1;
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* atan, n = 4, r = 43: 2 levels, widths [2, 1] */
static void
hs_atan_4_43(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x3333333333333333);

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = w0; P[3] = w1;
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* atan, n = 4, r = 44: 2 levels, widths [2, 1] */
static void
hs_atan_4_44(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x3333333333333333);

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = w0; P[3] = w1;
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* atan, n = 4, r = 45: 2 levels, widths [2, 1] */
static void
hs_atan_4_45(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x3333333333333333);

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = w0; P[3] = w1;
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* atan, n = 4, r = 46: 2 levels, widths [2, 1] */
static void
hs_atan_4_46(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x3333333333333333);

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = w0; P[3] = w1;
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* atan, n = 4, r = 47: 2 levels, widths [2, 1] */
static void
hs_atan_4_47(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x3333333333333333);

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = w0; P[3] = w1;
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* atan, n = 4, r = 48: 2 levels, widths [2, 1] */
static void
hs_atan_4_48(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x3333333333333333);

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x5555555555555555), UWORD(0x5555555555555555), h, l);

    /* atan = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = w0; P[3] = w1;
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* sin, n = 1, r = 8: 2 levels, widths [1, 1] */
static void
hs_sin_1_8(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x0222222222222222);

    /* V_0 = c - z V_1   (1 x 1 limbs) */
    v = n_mulhi(z[0], v);
    v = UWORD(0x2aaaaaaaaaaaaaaa) - v;

    /* sin = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* sin, n = 1, r = 9: 2 levels, widths [1, 1] */
static void
hs_sin_1_9(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x0222222222222222);

    /* V_0 = c - z V_1   (1 x 1 limbs) */
    v = n_mulhi(z[0], v);
    v = UWORD(0x2aaaaaaaaaaaaaaa) - v;

    /* sin = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* sin, n = 1, r = 10: 2 levels, widths [1, 1] */
static void
hs_sin_1_10(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x0222222222222222);

    /* V_0 = c - z V_1   (1 x 1 limbs) */
    v = n_mulhi(z[0], v);
    v = UWORD(0x2aaaaaaaaaaaaaaa) - v;

    /* sin = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* sin, n = 1, r = 11: 2 levels, widths [1, 1] */
static void
hs_sin_1_11(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x0222222222222222);

    /* V_0 = c - z V_1   (1 x 1 limbs) */
    v = n_mulhi(z[0], v);
    v = UWORD(0x2aaaaaaaaaaaaaaa) - v;

    /* sin = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* sin, n = 1, r = 12: 1 levels, widths [1] */
static void
hs_sin_1_12(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* sin, n = 1, r = 13: 1 levels, widths [1] */
static void
hs_sin_1_13(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* sin, n = 1, r = 14: 1 levels, widths [1] */
static void
hs_sin_1_14(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* sin, n = 1, r = 15: 1 levels, widths [1] */
static void
hs_sin_1_15(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* sin, n = 1, r = 16: 1 levels, widths [1] */
static void
hs_sin_1_16(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* sin, n = 1, r = 17: 1 levels, widths [1] */
static void
hs_sin_1_17(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* sin, n = 1, r = 18: 1 levels, widths [1] */
static void
hs_sin_1_18(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* sin, n = 1, r = 19: 1 levels, widths [1] */
static void
hs_sin_1_19(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* sin, n = 1, r = 20: 1 levels, widths [1] */
static void
hs_sin_1_20(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* sin, n = 1, r = 21: 1 levels, widths [1] */
static void
hs_sin_1_21(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* sin, n = 1, r = 22: 1 levels, widths [1] */
static void
hs_sin_1_22(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* sin, n = 1, r = 23: 1 levels, widths [1] */
static void
hs_sin_1_23(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* sin, n = 1, r = 24: 1 levels, widths [1] */
static void
hs_sin_1_24(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* sin, n = 1, r = 25: 1 levels, widths [1] */
static void
hs_sin_1_25(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* sin, n = 1, r = 26: 1 levels, widths [1] */
static void
hs_sin_1_26(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* sin, n = 1, r = 27: 1 levels, widths [1] */
static void
hs_sin_1_27(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* sin, n = 1, r = 28: 1 levels, widths [1] */
static void
hs_sin_1_28(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* sin, n = 1, r = 29: 1 levels, widths [1] */
static void
hs_sin_1_29(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* sin, n = 1, r = 30: 1 levels, widths [1] */
static void
hs_sin_1_30(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* sin, n = 1, r = 31: 1 levels, widths [1] */
static void
hs_sin_1_31(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* sin, n = 1, r = 32: 1 levels, widths [1] */
static void
hs_sin_1_32(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* sin, n = 1, r = 33: 1 levels, widths [1] */
static void
hs_sin_1_33(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* sin, n = 1, r = 34: 1 levels, widths [1] */
static void
hs_sin_1_34(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* sin, n = 1, r = 35: 1 levels, widths [1] */
static void
hs_sin_1_35(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* sin, n = 1, r = 36: 1 levels, widths [1] */
static void
hs_sin_1_36(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* sin, n = 1, r = 37: 1 levels, widths [1] */
static void
hs_sin_1_37(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* sin, n = 1, r = 38: 1 levels, widths [1] */
static void
hs_sin_1_38(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* sin, n = 1, r = 39: 1 levels, widths [1] */
static void
hs_sin_1_39(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* sin, n = 1, r = 40: 1 levels, widths [1] */
static void
hs_sin_1_40(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* sin, n = 1, r = 41: 1 levels, widths [1] */
static void
hs_sin_1_41(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* sin, n = 1, r = 42: 1 levels, widths [1] */
static void
hs_sin_1_42(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* sin, n = 1, r = 43: 1 levels, widths [1] */
static void
hs_sin_1_43(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* sin, n = 1, r = 44: 1 levels, widths [1] */
static void
hs_sin_1_44(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* sin, n = 1, r = 45: 1 levels, widths [1] */
static void
hs_sin_1_45(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* sin, n = 1, r = 46: 1 levels, widths [1] */
static void
hs_sin_1_46(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* sin, n = 1, r = 47: 1 levels, widths [1] */
static void
hs_sin_1_47(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* sin, n = 1, r = 48: 1 levels, widths [1] */
static void
hs_sin_1_48(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, x, 1);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 1);
    mpn_sub_n(res, x, V, 1);
}

/* sin, n = 2, r = 8: 5 levels, widths [2, 2, 2, 1, 1] */
static void
hs_sin_2_8(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x0000006b99159fd5);

    /* V_3 = c - z V_4   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x00002e3bc74aad8e) - v;

    /* V_2 = c - z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    sub_ddmmss(w1, w0, UWORD(0x000d00d00d00d00d), UWORD(0x00d00d00d00d00d0), h, l);

    /* V_1 = c - z V_2   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[1], z[0], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x0222222222222222), UWORD(0x2222222222222222), h, l);

    /* V_0 = c - z V_1   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[1], z[0], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x2aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* sin = x - x^3 V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* sin, n = 2, r = 9: 5 levels, widths [2, 2, 1, 1, 1] */
static void
hs_sin_2_9(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x0000006b99159fd5);

    /* V_3 = c - z V_4   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x00002e3bc74aad8e) - v;

    /* V_2 = c - z V_3   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x000d00d00d00d00d) - v;

    /* V_1 = c - z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    sub_ddmmss(w1, w0, UWORD(0x0222222222222222), UWORD(0x2222222222222222), h, l);

    /* V_0 = c - z V_1   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[1], z[0], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x2aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* sin = x - x^3 V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* sin, n = 2, r = 10: 4 levels, widths [2, 2, 1, 1] */
static void
hs_sin_2_10(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x00002e3bc74aad8e);

    /* V_2 = c - z V_3   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x000d00d00d00d00d) - v;

    /* V_1 = c - z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    sub_ddmmss(w1, w0, UWORD(0x0222222222222222), UWORD(0x2222222222222222), h, l);

    /* V_0 = c - z V_1   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[1], z[0], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x2aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* sin = x - x^3 V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* sin, n = 2, r = 11: 4 levels, widths [2, 2, 1, 1] */
static void
hs_sin_2_11(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x00002e3bc74aad8e);

    /* V_2 = c - z V_3   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x000d00d00d00d00d) - v;

    /* V_1 = c - z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    sub_ddmmss(w1, w0, UWORD(0x0222222222222222), UWORD(0x2222222222222222), h, l);

    /* V_0 = c - z V_1   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[1], z[0], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x2aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* sin = x - x^3 V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* sin, n = 2, r = 12: 4 levels, widths [2, 1, 1, 1] */
static void
hs_sin_2_12(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x00002e3bc74aad8e);

    /* V_2 = c - z V_3   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x000d00d00d00d00d) - v;

    /* V_1 = c - z V_2   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x0222222222222222) - v;

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    sub_ddmmss(w1, w0, UWORD(0x2aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* sin = x - x^3 V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* sin, n = 2, r = 13: 3 levels, widths [2, 1, 1] */
static void
hs_sin_2_13(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x000d00d00d00d00d);

    /* V_1 = c - z V_2   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x0222222222222222) - v;

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    sub_ddmmss(w1, w0, UWORD(0x2aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* sin = x - x^3 V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* sin, n = 2, r = 14: 3 levels, widths [2, 1, 1] */
static void
hs_sin_2_14(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x000d00d00d00d00d);

    /* V_1 = c - z V_2   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x0222222222222222) - v;

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    sub_ddmmss(w1, w0, UWORD(0x2aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* sin = x - x^3 V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* sin, n = 2, r = 15: 3 levels, widths [2, 1, 1] */
static void
hs_sin_2_15(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x000d00d00d00d00d);

    /* V_1 = c - z V_2   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x0222222222222222) - v;

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    sub_ddmmss(w1, w0, UWORD(0x2aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* sin = x - x^3 V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* sin, n = 2, r = 16: 3 levels, widths [2, 1, 1] */
static void
hs_sin_2_16(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x000d00d00d00d00d);

    /* V_1 = c - z V_2   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x0222222222222222) - v;

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    sub_ddmmss(w1, w0, UWORD(0x2aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* sin = x - x^3 V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* sin, n = 2, r = 17: 2 levels, widths [2, 1] */
static void
hs_sin_2_17(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x0222222222222222);

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    sub_ddmmss(w1, w0, UWORD(0x2aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* sin = x - x^3 V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* sin, n = 2, r = 18: 2 levels, widths [2, 1] */
static void
hs_sin_2_18(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x0222222222222222);

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    sub_ddmmss(w1, w0, UWORD(0x2aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* sin = x - x^3 V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* sin, n = 2, r = 19: 2 levels, widths [2, 1] */
static void
hs_sin_2_19(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x0222222222222222);

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    sub_ddmmss(w1, w0, UWORD(0x2aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* sin = x - x^3 V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* sin, n = 2, r = 20: 2 levels, widths [1, 1] */
static void
hs_sin_2_20(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x0222222222222222);

    /* V_0 = c - z V_1   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x2aaaaaaaaaaaaaaa) - v;

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* sin, n = 2, r = 21: 2 levels, widths [1, 1] */
static void
hs_sin_2_21(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x0222222222222222);

    /* V_0 = c - z V_1   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x2aaaaaaaaaaaaaaa) - v;

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* sin, n = 2, r = 22: 2 levels, widths [1, 1] */
static void
hs_sin_2_22(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x0222222222222222);

    /* V_0 = c - z V_1   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x2aaaaaaaaaaaaaaa) - v;

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* sin, n = 2, r = 23: 2 levels, widths [1, 1] */
static void
hs_sin_2_23(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x0222222222222222);

    /* V_0 = c - z V_1   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x2aaaaaaaaaaaaaaa) - v;

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* sin, n = 2, r = 24: 2 levels, widths [1, 1] */
static void
hs_sin_2_24(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x0222222222222222);

    /* V_0 = c - z V_1   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x2aaaaaaaaaaaaaaa) - v;

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* sin, n = 2, r = 25: 1 levels, widths [1] */
static void
hs_sin_2_25(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* sin, n = 2, r = 26: 1 levels, widths [1] */
static void
hs_sin_2_26(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* sin, n = 2, r = 27: 1 levels, widths [1] */
static void
hs_sin_2_27(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* sin, n = 2, r = 28: 1 levels, widths [1] */
static void
hs_sin_2_28(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* sin, n = 2, r = 29: 1 levels, widths [1] */
static void
hs_sin_2_29(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* sin, n = 2, r = 30: 1 levels, widths [1] */
static void
hs_sin_2_30(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* sin, n = 2, r = 31: 1 levels, widths [1] */
static void
hs_sin_2_31(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* sin, n = 2, r = 32: 1 levels, widths [1] */
static void
hs_sin_2_32(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* sin, n = 2, r = 33: 1 levels, widths [1] */
static void
hs_sin_2_33(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* sin, n = 2, r = 34: 1 levels, widths [1] */
static void
hs_sin_2_34(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* sin, n = 2, r = 35: 1 levels, widths [1] */
static void
hs_sin_2_35(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* sin, n = 2, r = 36: 1 levels, widths [1] */
static void
hs_sin_2_36(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* sin, n = 2, r = 37: 1 levels, widths [1] */
static void
hs_sin_2_37(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* sin, n = 2, r = 38: 1 levels, widths [1] */
static void
hs_sin_2_38(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* sin, n = 2, r = 39: 1 levels, widths [1] */
static void
hs_sin_2_39(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* sin, n = 2, r = 40: 1 levels, widths [1] */
static void
hs_sin_2_40(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* sin, n = 2, r = 41: 1 levels, widths [1] */
static void
hs_sin_2_41(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* sin, n = 2, r = 42: 1 levels, widths [1] */
static void
hs_sin_2_42(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* sin, n = 2, r = 43: 1 levels, widths [1] */
static void
hs_sin_2_43(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* sin, n = 2, r = 44: 1 levels, widths [1] */
static void
hs_sin_2_44(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* sin, n = 2, r = 45: 1 levels, widths [1] */
static void
hs_sin_2_45(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* sin, n = 2, r = 46: 1 levels, widths [1] */
static void
hs_sin_2_46(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* sin, n = 2, r = 47: 1 levels, widths [1] */
static void
hs_sin_2_47(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* sin, n = 2, r = 48: 1 levels, widths [1] */
static void
hs_sin_2_48(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, x, 2);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 2);
    mpn_sub_n(res, x, V, 2);
}

/* sin, n = 3, r = 8: 8 levels, widths [3, 3, 3, 2, 2, 2, 2, 1] */
static void
hs_sin_3_8(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x000000000000ca96);

    /* V_6 = c - z V_7   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x0000000000d73f9f), UWORD(0x399dc0f88ec32b58), h, l);

    /* V_5 = c - z V_6   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x00000000b092309d), UWORD(0x43684be51c198e91), h, l);

    /* V_4 = c - z V_5   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x0000006b99159fd5), UWORD(0x138e3f9d1f92e0df), h, l);

    /* V_3 = c - z V_4   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x00002e3bc74aad8e), UWORD(0x671f5583911ca002), h, l);

    /* V_2 = c - z V_3   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x0d00d00d00d00d00);
    P[1] = UWORD(0x00d00d00d00d00d0);
    P[2] = UWORD(0x000d00d00d00d00d);
    mpn_sub_n(V, P, T, 3);

    /* V_1 = c - z V_2   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x2222222222222222);
    P[1] = UWORD(0x2222222222222222);
    P[2] = UWORD(0x0222222222222222);
    mpn_sub_n(V, P, T, 3);

    /* V_0 = c - z V_1   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 3);

    /* sin = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* sin, n = 3, r = 9: 7 levels, widths [3, 3, 2, 2, 2, 2, 1] */
static void
hs_sin_3_9(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0000000000d73f9f);

    /* V_5 = c - z V_6   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x00000000b092309d), UWORD(0x43684be51c198e91), h, l);

    /* V_4 = c - z V_5   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x0000006b99159fd5), UWORD(0x138e3f9d1f92e0df), h, l);

    /* V_3 = c - z V_4   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x00002e3bc74aad8e), UWORD(0x671f5583911ca002), h, l);

    /* V_2 = c - z V_3   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x000d00d00d00d00d), UWORD(0x00d00d00d00d00d0), h, l);

    /* V_1 = c - z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x2222222222222222);
    P[1] = UWORD(0x2222222222222222);
    P[2] = UWORD(0x0222222222222222);
    mpn_sub_n(V, P, T, 3);

    /* V_0 = c - z V_1   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 3);

    /* sin = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* sin, n = 3, r = 10: 7 levels, widths [3, 3, 2, 2, 2, 1, 1] */
static void
hs_sin_3_10(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0000000000d73f9f);

    /* V_5 = c - z V_6   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x00000000b092309d) - v;

    /* V_4 = c - z V_5   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x0000006b99159fd5), UWORD(0x138e3f9d1f92e0df), h, l);

    /* V_3 = c - z V_4   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x00002e3bc74aad8e), UWORD(0x671f5583911ca002), h, l);

    /* V_2 = c - z V_3   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x000d00d00d00d00d), UWORD(0x00d00d00d00d00d0), h, l);

    /* V_1 = c - z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x2222222222222222);
    P[1] = UWORD(0x2222222222222222);
    P[2] = UWORD(0x0222222222222222);
    mpn_sub_n(V, P, T, 3);

    /* V_0 = c - z V_1   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 3);

    /* sin = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* sin, n = 3, r = 11: 6 levels, widths [3, 3, 2, 2, 2, 1] */
static void
hs_sin_3_11(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x00000000b092309d);

    /* V_4 = c - z V_5   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x0000006b99159fd5), UWORD(0x138e3f9d1f92e0df), h, l);

    /* V_3 = c - z V_4   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x00002e3bc74aad8e), UWORD(0x671f5583911ca002), h, l);

    /* V_2 = c - z V_3   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x000d00d00d00d00d), UWORD(0x00d00d00d00d00d0), h, l);

    /* V_1 = c - z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x2222222222222222);
    P[1] = UWORD(0x2222222222222222);
    P[2] = UWORD(0x0222222222222222);
    mpn_sub_n(V, P, T, 3);

    /* V_0 = c - z V_1   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 3);

    /* sin = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* sin, n = 3, r = 12: 6 levels, widths [3, 2, 2, 2, 1, 1] */
static void
hs_sin_3_12(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x00000000b092309d);

    /* V_4 = c - z V_5   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x0000006b99159fd5) - v;

    /* V_3 = c - z V_4   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x00002e3bc74aad8e), UWORD(0x671f5583911ca002), h, l);

    /* V_2 = c - z V_3   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x000d00d00d00d00d), UWORD(0x00d00d00d00d00d0), h, l);

    /* V_1 = c - z V_2   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x0222222222222222), UWORD(0x2222222222222222), h, l);

    /* V_0 = c - z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 3);

    /* sin = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* sin, n = 3, r = 13: 5 levels, widths [3, 2, 2, 2, 1] */
static void
hs_sin_3_13(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0000006b99159fd5);

    /* V_3 = c - z V_4   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x00002e3bc74aad8e), UWORD(0x671f5583911ca002), h, l);

    /* V_2 = c - z V_3   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x000d00d00d00d00d), UWORD(0x00d00d00d00d00d0), h, l);

    /* V_1 = c - z V_2   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x0222222222222222), UWORD(0x2222222222222222), h, l);

    /* V_0 = c - z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 3);

    /* sin = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* sin, n = 3, r = 14: 5 levels, widths [3, 2, 2, 1, 1] */
static void
hs_sin_3_14(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0000006b99159fd5);

    /* V_3 = c - z V_4   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x00002e3bc74aad8e) - v;

    /* V_2 = c - z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x000d00d00d00d00d), UWORD(0x00d00d00d00d00d0), h, l);

    /* V_1 = c - z V_2   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x0222222222222222), UWORD(0x2222222222222222), h, l);

    /* V_0 = c - z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 3);

    /* sin = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* sin, n = 3, r = 15: 5 levels, widths [3, 2, 2, 1, 1] */
static void
hs_sin_3_15(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0000006b99159fd5);

    /* V_3 = c - z V_4   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x00002e3bc74aad8e) - v;

    /* V_2 = c - z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x000d00d00d00d00d), UWORD(0x00d00d00d00d00d0), h, l);

    /* V_1 = c - z V_2   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x0222222222222222), UWORD(0x2222222222222222), h, l);

    /* V_0 = c - z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 3);

    /* sin = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* sin, n = 3, r = 16: 4 levels, widths [3, 2, 2, 1] */
static void
hs_sin_3_16(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x00002e3bc74aad8e);

    /* V_2 = c - z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x000d00d00d00d00d), UWORD(0x00d00d00d00d00d0), h, l);

    /* V_1 = c - z V_2   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x0222222222222222), UWORD(0x2222222222222222), h, l);

    /* V_0 = c - z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 3);

    /* sin = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* sin, n = 3, r = 17: 4 levels, widths [3, 2, 2, 1] */
static void
hs_sin_3_17(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x00002e3bc74aad8e);

    /* V_2 = c - z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x000d00d00d00d00d), UWORD(0x00d00d00d00d00d0), h, l);

    /* V_1 = c - z V_2   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x0222222222222222), UWORD(0x2222222222222222), h, l);

    /* V_0 = c - z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 3);

    /* sin = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* sin, n = 3, r = 18: 4 levels, widths [3, 2, 1, 1] */
static void
hs_sin_3_18(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x00002e3bc74aad8e);

    /* V_2 = c - z V_3   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x000d00d00d00d00d) - v;

    /* V_1 = c - z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x0222222222222222), UWORD(0x2222222222222222), h, l);

    /* V_0 = c - z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 3);

    /* sin = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* sin, n = 3, r = 19: 4 levels, widths [3, 2, 1, 1] */
static void
hs_sin_3_19(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x00002e3bc74aad8e);

    /* V_2 = c - z V_3   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x000d00d00d00d00d) - v;

    /* V_1 = c - z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x0222222222222222), UWORD(0x2222222222222222), h, l);

    /* V_0 = c - z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 3);

    /* sin = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* sin, n = 3, r = 20: 3 levels, widths [2, 2, 1] */
static void
hs_sin_3_20(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x000d00d00d00d00d);

    /* V_1 = c - z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x0222222222222222), UWORD(0x2222222222222222), h, l);

    /* V_0 = c - z V_1   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x2aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* sin, n = 3, r = 21: 3 levels, widths [2, 2, 1] */
static void
hs_sin_3_21(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x000d00d00d00d00d);

    /* V_1 = c - z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x0222222222222222), UWORD(0x2222222222222222), h, l);

    /* V_0 = c - z V_1   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x2aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* sin, n = 3, r = 22: 3 levels, widths [2, 2, 1] */
static void
hs_sin_3_22(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x000d00d00d00d00d);

    /* V_1 = c - z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x0222222222222222), UWORD(0x2222222222222222), h, l);

    /* V_0 = c - z V_1   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x2aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* sin, n = 3, r = 23: 3 levels, widths [2, 2, 1] */
static void
hs_sin_3_23(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x000d00d00d00d00d);

    /* V_1 = c - z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x0222222222222222), UWORD(0x2222222222222222), h, l);

    /* V_0 = c - z V_1   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x2aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* sin, n = 3, r = 24: 3 levels, widths [2, 2, 1] */
static void
hs_sin_3_24(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x000d00d00d00d00d);

    /* V_1 = c - z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x0222222222222222), UWORD(0x2222222222222222), h, l);

    /* V_0 = c - z V_1   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x2aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* sin, n = 3, r = 25: 3 levels, widths [2, 1, 1] */
static void
hs_sin_3_25(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x000d00d00d00d00d);

    /* V_1 = c - z V_2   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x0222222222222222) - v;

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x2aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* sin, n = 3, r = 26: 2 levels, widths [2, 1] */
static void
hs_sin_3_26(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0222222222222222);

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x2aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* sin, n = 3, r = 27: 2 levels, widths [2, 1] */
static void
hs_sin_3_27(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0222222222222222);

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x2aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* sin, n = 3, r = 28: 2 levels, widths [2, 1] */
static void
hs_sin_3_28(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0222222222222222);

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x2aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* sin, n = 3, r = 29: 2 levels, widths [2, 1] */
static void
hs_sin_3_29(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0222222222222222);

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x2aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* sin, n = 3, r = 30: 2 levels, widths [2, 1] */
static void
hs_sin_3_30(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0222222222222222);

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x2aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* sin, n = 3, r = 31: 2 levels, widths [2, 1] */
static void
hs_sin_3_31(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0222222222222222);

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x2aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* sin, n = 3, r = 32: 2 levels, widths [2, 1] */
static void
hs_sin_3_32(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0222222222222222);

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x2aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* sin, n = 3, r = 33: 2 levels, widths [2, 1] */
static void
hs_sin_3_33(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0222222222222222);

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x2aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* sin, n = 3, r = 34: 2 levels, widths [2, 1] */
static void
hs_sin_3_34(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0222222222222222);

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x2aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* sin, n = 3, r = 35: 2 levels, widths [2, 1] */
static void
hs_sin_3_35(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0222222222222222);

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x2aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* sin, n = 3, r = 36: 2 levels, widths [2, 1] */
static void
hs_sin_3_36(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0222222222222222);

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x2aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* sin, n = 3, r = 37: 2 levels, widths [2, 1] */
static void
hs_sin_3_37(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0222222222222222);

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    sub_ddmmss(w1, w0, UWORD(0x2aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* sin, n = 3, r = 38: 1 levels, widths [2] */
static void
hs_sin_3_38(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    w1 = UWORD(0x2aaaaaaaaaaaaaaa); w0 = UWORD(0xaaaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* sin, n = 3, r = 39: 1 levels, widths [2] */
static void
hs_sin_3_39(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    w1 = UWORD(0x2aaaaaaaaaaaaaaa); w0 = UWORD(0xaaaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* sin, n = 3, r = 40: 1 levels, widths [2] */
static void
hs_sin_3_40(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    w1 = UWORD(0x2aaaaaaaaaaaaaaa); w0 = UWORD(0xaaaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* sin, n = 3, r = 41: 1 levels, widths [1] */
static void
hs_sin_3_41(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = v;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* sin, n = 3, r = 42: 1 levels, widths [1] */
static void
hs_sin_3_42(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = v;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* sin, n = 3, r = 43: 1 levels, widths [1] */
static void
hs_sin_3_43(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = v;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* sin, n = 3, r = 44: 1 levels, widths [1] */
static void
hs_sin_3_44(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = v;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* sin, n = 3, r = 45: 1 levels, widths [1] */
static void
hs_sin_3_45(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = v;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* sin, n = 3, r = 46: 1 levels, widths [1] */
static void
hs_sin_3_46(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = v;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* sin, n = 3, r = 47: 1 levels, widths [1] */
static void
hs_sin_3_47(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = v;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* sin, n = 3, r = 48: 1 levels, widths [1] */
static void
hs_sin_3_48(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x2aaaaaaaaaaaaaaa);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = v;
    flint_mpn_mulhigh_n(T, z, x, 3);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 3);
    mpn_sub_n(res, x, V, 3);
}

/* sin, n = 4, r = 8: 10 levels, widths [4, 4, 4, 3, 3, 3, 3, 2, 2, 2] */
static void
hs_sin_4_8(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    w1 = UWORD(0x0000000000000000); w0 = UWORD(0x5c6e3bdb73d5c62f);

    /* V_8 = c - z V_9   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x0000000000000097), UWORD(0xa4da340a0ab92650), h, l);

    /* V_7 = c - z V_8   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x000000000000ca96), UWORD(0x3b81856a53593028), h, l);

    /* V_6 = c - z V_7   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x774657f48f5eaf63);
    P[1] = UWORD(0x399dc0f88ec32b58);
    P[2] = UWORD(0x0000000000d73f9f);
    mpn_sub_n(V, P, T, 3);

    /* V_5 = c - z V_6   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xd7b4269d9babdfa2);
    P[1] = UWORD(0x43684be51c198e91);
    P[2] = UWORD(0x00000000b092309d);
    mpn_sub_n(V, P, T, 3);

    /* V_4 = c - z V_5   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x71c7880adcbc46da);
    P[1] = UWORD(0x138e3f9d1f92e0df);
    P[2] = UWORD(0x0000006b99159fd5);
    mpn_sub_n(V, P, T, 3);

    /* V_3 = c - z V_4   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xe3bc74aad8e671f5);
    P[1] = UWORD(0x671f5583911ca002);
    P[2] = UWORD(0x00002e3bc74aad8e);
    mpn_sub_n(V, P, T, 3);

    /* V_2 = c - z V_3   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0xd00d00d00d00d00d);
    P[1] = UWORD(0x0d00d00d00d00d00);
    P[2] = UWORD(0x00d00d00d00d00d0);
    P[3] = UWORD(0x000d00d00d00d00d);
    mpn_sub_n(V, P, T, 4);

    /* V_1 = c - z V_2   (4 x 4 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x2222222222222222);
    P[1] = UWORD(0x2222222222222222);
    P[2] = UWORD(0x2222222222222222);
    P[3] = UWORD(0x0222222222222222);
    mpn_sub_n(V, P, T, 4);

    /* V_0 = c - z V_1   (4 x 4 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[3] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 4);

    /* sin = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* sin, n = 4, r = 9: 10 levels, widths [4, 4, 3, 3, 3, 3, 2, 2, 2, 1] */
static void
hs_sin_4_9(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0000000000000000);

    /* V_8 = c - z V_9   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x0000000000000097), UWORD(0xa4da340a0ab92650), h, l);

    /* V_7 = c - z V_8   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x000000000000ca96), UWORD(0x3b81856a53593028), h, l);

    /* V_6 = c - z V_7   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x0000000000d73f9f), UWORD(0x399dc0f88ec32b58), h, l);

    /* V_5 = c - z V_6   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xd7b4269d9babdfa2);
    P[1] = UWORD(0x43684be51c198e91);
    P[2] = UWORD(0x00000000b092309d);
    mpn_sub_n(V, P, T, 3);

    /* V_4 = c - z V_5   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x71c7880adcbc46da);
    P[1] = UWORD(0x138e3f9d1f92e0df);
    P[2] = UWORD(0x0000006b99159fd5);
    mpn_sub_n(V, P, T, 3);

    /* V_3 = c - z V_4   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xe3bc74aad8e671f5);
    P[1] = UWORD(0x671f5583911ca002);
    P[2] = UWORD(0x00002e3bc74aad8e);
    mpn_sub_n(V, P, T, 3);

    /* V_2 = c - z V_3   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x0d00d00d00d00d00);
    P[1] = UWORD(0x00d00d00d00d00d0);
    P[2] = UWORD(0x000d00d00d00d00d);
    mpn_sub_n(V, P, T, 3);

    /* V_1 = c - z V_2   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x2222222222222222);
    P[1] = UWORD(0x2222222222222222);
    P[2] = UWORD(0x2222222222222222);
    P[3] = UWORD(0x0222222222222222);
    mpn_sub_n(V, P, T, 4);

    /* V_0 = c - z V_1   (4 x 4 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[3] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 4);

    /* sin = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* sin, n = 4, r = 10: 9 levels, widths [4, 4, 3, 3, 3, 2, 2, 2, 1] */
static void
hs_sin_4_10(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0000000000000097);

    /* V_7 = c - z V_8   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x000000000000ca96), UWORD(0x3b81856a53593028), h, l);

    /* V_6 = c - z V_7   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x0000000000d73f9f), UWORD(0x399dc0f88ec32b58), h, l);

    /* V_5 = c - z V_6   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x00000000b092309d), UWORD(0x43684be51c198e91), h, l);

    /* V_4 = c - z V_5   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x71c7880adcbc46da);
    P[1] = UWORD(0x138e3f9d1f92e0df);
    P[2] = UWORD(0x0000006b99159fd5);
    mpn_sub_n(V, P, T, 3);

    /* V_3 = c - z V_4   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xe3bc74aad8e671f5);
    P[1] = UWORD(0x671f5583911ca002);
    P[2] = UWORD(0x00002e3bc74aad8e);
    mpn_sub_n(V, P, T, 3);

    /* V_2 = c - z V_3   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x0d00d00d00d00d00);
    P[1] = UWORD(0x00d00d00d00d00d0);
    P[2] = UWORD(0x000d00d00d00d00d);
    mpn_sub_n(V, P, T, 3);

    /* V_1 = c - z V_2   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x2222222222222222);
    P[1] = UWORD(0x2222222222222222);
    P[2] = UWORD(0x2222222222222222);
    P[3] = UWORD(0x0222222222222222);
    mpn_sub_n(V, P, T, 4);

    /* V_0 = c - z V_1   (4 x 4 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[3] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 4);

    /* sin = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* sin, n = 4, r = 11: 8 levels, widths [4, 4, 3, 3, 3, 2, 2, 1] */
static void
hs_sin_4_11(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x000000000000ca96);

    /* V_6 = c - z V_7   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x0000000000d73f9f), UWORD(0x399dc0f88ec32b58), h, l);

    /* V_5 = c - z V_6   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x00000000b092309d), UWORD(0x43684be51c198e91), h, l);

    /* V_4 = c - z V_5   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x71c7880adcbc46da);
    P[1] = UWORD(0x138e3f9d1f92e0df);
    P[2] = UWORD(0x0000006b99159fd5);
    mpn_sub_n(V, P, T, 3);

    /* V_3 = c - z V_4   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xe3bc74aad8e671f5);
    P[1] = UWORD(0x671f5583911ca002);
    P[2] = UWORD(0x00002e3bc74aad8e);
    mpn_sub_n(V, P, T, 3);

    /* V_2 = c - z V_3   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x0d00d00d00d00d00);
    P[1] = UWORD(0x00d00d00d00d00d0);
    P[2] = UWORD(0x000d00d00d00d00d);
    mpn_sub_n(V, P, T, 3);

    /* V_1 = c - z V_2   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x2222222222222222);
    P[1] = UWORD(0x2222222222222222);
    P[2] = UWORD(0x2222222222222222);
    P[3] = UWORD(0x0222222222222222);
    mpn_sub_n(V, P, T, 4);

    /* V_0 = c - z V_1   (4 x 4 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[3] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 4);

    /* sin = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* sin, n = 4, r = 12: 8 levels, widths [4, 3, 3, 3, 2, 2, 2, 1] */
static void
hs_sin_4_12(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x000000000000ca96);

    /* V_6 = c - z V_7   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x0000000000d73f9f), UWORD(0x399dc0f88ec32b58), h, l);

    /* V_5 = c - z V_6   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x00000000b092309d), UWORD(0x43684be51c198e91), h, l);

    /* V_4 = c - z V_5   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x0000006b99159fd5), UWORD(0x138e3f9d1f92e0df), h, l);

    /* V_3 = c - z V_4   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xe3bc74aad8e671f5);
    P[1] = UWORD(0x671f5583911ca002);
    P[2] = UWORD(0x00002e3bc74aad8e);
    mpn_sub_n(V, P, T, 3);

    /* V_2 = c - z V_3   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x0d00d00d00d00d00);
    P[1] = UWORD(0x00d00d00d00d00d0);
    P[2] = UWORD(0x000d00d00d00d00d);
    mpn_sub_n(V, P, T, 3);

    /* V_1 = c - z V_2   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x2222222222222222);
    P[1] = UWORD(0x2222222222222222);
    P[2] = UWORD(0x0222222222222222);
    mpn_sub_n(V, P, T, 3);

    /* V_0 = c - z V_1   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[3] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 4);

    /* sin = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* sin, n = 4, r = 13: 7 levels, widths [4, 3, 3, 3, 2, 2, 1] */
static void
hs_sin_4_13(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0000000000d73f9f);

    /* V_5 = c - z V_6   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x00000000b092309d), UWORD(0x43684be51c198e91), h, l);

    /* V_4 = c - z V_5   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x0000006b99159fd5), UWORD(0x138e3f9d1f92e0df), h, l);

    /* V_3 = c - z V_4   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xe3bc74aad8e671f5);
    P[1] = UWORD(0x671f5583911ca002);
    P[2] = UWORD(0x00002e3bc74aad8e);
    mpn_sub_n(V, P, T, 3);

    /* V_2 = c - z V_3   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x0d00d00d00d00d00);
    P[1] = UWORD(0x00d00d00d00d00d0);
    P[2] = UWORD(0x000d00d00d00d00d);
    mpn_sub_n(V, P, T, 3);

    /* V_1 = c - z V_2   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x2222222222222222);
    P[1] = UWORD(0x2222222222222222);
    P[2] = UWORD(0x0222222222222222);
    mpn_sub_n(V, P, T, 3);

    /* V_0 = c - z V_1   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[3] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 4);

    /* sin = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* sin, n = 4, r = 14: 7 levels, widths [4, 3, 3, 2, 2, 2, 1] */
static void
hs_sin_4_14(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0000000000d73f9f);

    /* V_5 = c - z V_6   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x00000000b092309d), UWORD(0x43684be51c198e91), h, l);

    /* V_4 = c - z V_5   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x0000006b99159fd5), UWORD(0x138e3f9d1f92e0df), h, l);

    /* V_3 = c - z V_4   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x00002e3bc74aad8e), UWORD(0x671f5583911ca002), h, l);

    /* V_2 = c - z V_3   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x0d00d00d00d00d00);
    P[1] = UWORD(0x00d00d00d00d00d0);
    P[2] = UWORD(0x000d00d00d00d00d);
    mpn_sub_n(V, P, T, 3);

    /* V_1 = c - z V_2   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x2222222222222222);
    P[1] = UWORD(0x2222222222222222);
    P[2] = UWORD(0x0222222222222222);
    mpn_sub_n(V, P, T, 3);

    /* V_0 = c - z V_1   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[3] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 4);

    /* sin = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* sin, n = 4, r = 15: 6 levels, widths [4, 3, 3, 2, 2, 1] */
static void
hs_sin_4_15(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x00000000b092309d);

    /* V_4 = c - z V_5   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x0000006b99159fd5), UWORD(0x138e3f9d1f92e0df), h, l);

    /* V_3 = c - z V_4   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x00002e3bc74aad8e), UWORD(0x671f5583911ca002), h, l);

    /* V_2 = c - z V_3   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x0d00d00d00d00d00);
    P[1] = UWORD(0x00d00d00d00d00d0);
    P[2] = UWORD(0x000d00d00d00d00d);
    mpn_sub_n(V, P, T, 3);

    /* V_1 = c - z V_2   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x2222222222222222);
    P[1] = UWORD(0x2222222222222222);
    P[2] = UWORD(0x0222222222222222);
    mpn_sub_n(V, P, T, 3);

    /* V_0 = c - z V_1   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[3] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 4);

    /* sin = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* sin, n = 4, r = 16: 6 levels, widths [4, 3, 3, 2, 2, 1] */
static void
hs_sin_4_16(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x00000000b092309d);

    /* V_4 = c - z V_5   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x0000006b99159fd5), UWORD(0x138e3f9d1f92e0df), h, l);

    /* V_3 = c - z V_4   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x00002e3bc74aad8e), UWORD(0x671f5583911ca002), h, l);

    /* V_2 = c - z V_3   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x0d00d00d00d00d00);
    P[1] = UWORD(0x00d00d00d00d00d0);
    P[2] = UWORD(0x000d00d00d00d00d);
    mpn_sub_n(V, P, T, 3);

    /* V_1 = c - z V_2   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x2222222222222222);
    P[1] = UWORD(0x2222222222222222);
    P[2] = UWORD(0x0222222222222222);
    mpn_sub_n(V, P, T, 3);

    /* V_0 = c - z V_1   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[3] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 4);

    /* sin = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* sin, n = 4, r = 17: 6 levels, widths [4, 3, 3, 2, 1, 1] */
static void
hs_sin_4_17(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x00000000b092309d);

    /* V_4 = c - z V_5   (1 x 1 limbs) */
    v = n_mulhi(z[3], v);
    v = UWORD(0x0000006b99159fd5) - v;

    /* V_3 = c - z V_4   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x00002e3bc74aad8e), UWORD(0x671f5583911ca002), h, l);

    /* V_2 = c - z V_3   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x0d00d00d00d00d00);
    P[1] = UWORD(0x00d00d00d00d00d0);
    P[2] = UWORD(0x000d00d00d00d00d);
    mpn_sub_n(V, P, T, 3);

    /* V_1 = c - z V_2   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x2222222222222222);
    P[1] = UWORD(0x2222222222222222);
    P[2] = UWORD(0x0222222222222222);
    mpn_sub_n(V, P, T, 3);

    /* V_0 = c - z V_1   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[3] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 4);

    /* sin = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* sin, n = 4, r = 18: 5 levels, widths [4, 3, 2, 2, 1] */
static void
hs_sin_4_18(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0000006b99159fd5);

    /* V_3 = c - z V_4   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x00002e3bc74aad8e), UWORD(0x671f5583911ca002), h, l);

    /* V_2 = c - z V_3   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x000d00d00d00d00d), UWORD(0x00d00d00d00d00d0), h, l);

    /* V_1 = c - z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x2222222222222222);
    P[1] = UWORD(0x2222222222222222);
    P[2] = UWORD(0x0222222222222222);
    mpn_sub_n(V, P, T, 3);

    /* V_0 = c - z V_1   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[3] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 4);

    /* sin = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* sin, n = 4, r = 19: 5 levels, widths [4, 3, 2, 2, 1] */
static void
hs_sin_4_19(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0000006b99159fd5);

    /* V_3 = c - z V_4   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x00002e3bc74aad8e), UWORD(0x671f5583911ca002), h, l);

    /* V_2 = c - z V_3   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x000d00d00d00d00d), UWORD(0x00d00d00d00d00d0), h, l);

    /* V_1 = c - z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x2222222222222222);
    P[1] = UWORD(0x2222222222222222);
    P[2] = UWORD(0x0222222222222222);
    mpn_sub_n(V, P, T, 3);

    /* V_0 = c - z V_1   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[3] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 4);

    /* sin = x - x^3 V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* sin, n = 4, r = 20: 5 levels, widths [3, 3, 2, 2, 1] */
static void
hs_sin_4_20(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0000006b99159fd5);

    /* V_3 = c - z V_4   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x00002e3bc74aad8e), UWORD(0x671f5583911ca002), h, l);

    /* V_2 = c - z V_3   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x000d00d00d00d00d), UWORD(0x00d00d00d00d00d0), h, l);

    /* V_1 = c - z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x2222222222222222);
    P[1] = UWORD(0x2222222222222222);
    P[2] = UWORD(0x0222222222222222);
    mpn_sub_n(V, P, T, 3);

    /* V_0 = c - z V_1   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 3);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* sin, n = 4, r = 21: 4 levels, widths [3, 3, 2, 1] */
static void
hs_sin_4_21(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x00002e3bc74aad8e);

    /* V_2 = c - z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x000d00d00d00d00d), UWORD(0x00d00d00d00d00d0), h, l);

    /* V_1 = c - z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x2222222222222222);
    P[1] = UWORD(0x2222222222222222);
    P[2] = UWORD(0x0222222222222222);
    mpn_sub_n(V, P, T, 3);

    /* V_0 = c - z V_1   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 3);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* sin, n = 4, r = 22: 4 levels, widths [3, 3, 2, 1] */
static void
hs_sin_4_22(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x00002e3bc74aad8e);

    /* V_2 = c - z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x000d00d00d00d00d), UWORD(0x00d00d00d00d00d0), h, l);

    /* V_1 = c - z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x2222222222222222);
    P[1] = UWORD(0x2222222222222222);
    P[2] = UWORD(0x0222222222222222);
    mpn_sub_n(V, P, T, 3);

    /* V_0 = c - z V_1   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 3);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* sin, n = 4, r = 23: 4 levels, widths [3, 3, 2, 1] */
static void
hs_sin_4_23(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x00002e3bc74aad8e);

    /* V_2 = c - z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x000d00d00d00d00d), UWORD(0x00d00d00d00d00d0), h, l);

    /* V_1 = c - z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x2222222222222222);
    P[1] = UWORD(0x2222222222222222);
    P[2] = UWORD(0x0222222222222222);
    mpn_sub_n(V, P, T, 3);

    /* V_0 = c - z V_1   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 3);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* sin, n = 4, r = 24: 4 levels, widths [3, 3, 2, 1] */
static void
hs_sin_4_24(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x00002e3bc74aad8e);

    /* V_2 = c - z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x000d00d00d00d00d), UWORD(0x00d00d00d00d00d0), h, l);

    /* V_1 = c - z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x2222222222222222);
    P[1] = UWORD(0x2222222222222222);
    P[2] = UWORD(0x0222222222222222);
    mpn_sub_n(V, P, T, 3);

    /* V_0 = c - z V_1   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 3);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* sin, n = 4, r = 25: 4 levels, widths [3, 2, 2, 1] */
static void
hs_sin_4_25(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x00002e3bc74aad8e);

    /* V_2 = c - z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x000d00d00d00d00d), UWORD(0x00d00d00d00d00d0), h, l);

    /* V_1 = c - z V_2   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x0222222222222222), UWORD(0x2222222222222222), h, l);

    /* V_0 = c - z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 3);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* sin, n = 4, r = 26: 4 levels, widths [3, 2, 2, 1] */
static void
hs_sin_4_26(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x00002e3bc74aad8e);

    /* V_2 = c - z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x000d00d00d00d00d), UWORD(0x00d00d00d00d00d0), h, l);

    /* V_1 = c - z V_2   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    sub_ddmmss(w1, w0, UWORD(0x0222222222222222), UWORD(0x2222222222222222), h, l);

    /* V_0 = c - z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 3);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* sin, n = 4, r = 27: 3 levels, widths [3, 2, 1] */
static void
hs_sin_4_27(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x000d00d00d00d00d);

    /* V_1 = c - z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x0222222222222222), UWORD(0x2222222222222222), h, l);

    /* V_0 = c - z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 3);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* sin, n = 4, r = 28: 3 levels, widths [3, 2, 1] */
static void
hs_sin_4_28(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x000d00d00d00d00d);

    /* V_1 = c - z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x0222222222222222), UWORD(0x2222222222222222), h, l);

    /* V_0 = c - z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 3);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* sin, n = 4, r = 29: 3 levels, widths [3, 2, 1] */
static void
hs_sin_4_29(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x000d00d00d00d00d);

    /* V_1 = c - z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x0222222222222222), UWORD(0x2222222222222222), h, l);

    /* V_0 = c - z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 3);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* sin, n = 4, r = 30: 3 levels, widths [3, 2, 1] */
static void
hs_sin_4_30(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x000d00d00d00d00d);

    /* V_1 = c - z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x0222222222222222), UWORD(0x2222222222222222), h, l);

    /* V_0 = c - z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 3);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* sin, n = 4, r = 31: 3 levels, widths [3, 2, 1] */
static void
hs_sin_4_31(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x000d00d00d00d00d);

    /* V_1 = c - z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x0222222222222222), UWORD(0x2222222222222222), h, l);

    /* V_0 = c - z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 3);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* sin, n = 4, r = 32: 3 levels, widths [3, 2, 1] */
static void
hs_sin_4_32(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x000d00d00d00d00d);

    /* V_1 = c - z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x0222222222222222), UWORD(0x2222222222222222), h, l);

    /* V_0 = c - z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 3);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* sin, n = 4, r = 33: 3 levels, widths [3, 2, 1] */
static void
hs_sin_4_33(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x000d00d00d00d00d);

    /* V_1 = c - z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x0222222222222222), UWORD(0x2222222222222222), h, l);

    /* V_0 = c - z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 3);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* sin, n = 4, r = 34: 3 levels, widths [3, 2, 1] */
static void
hs_sin_4_34(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x000d00d00d00d00d);

    /* V_1 = c - z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x0222222222222222), UWORD(0x2222222222222222), h, l);

    /* V_0 = c - z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 3);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* sin, n = 4, r = 35: 2 levels, widths [3, 2] */
static void
hs_sin_4_35(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    w1 = UWORD(0x0222222222222222); w0 = UWORD(0x2222222222222222);

    /* V_0 = c - z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 3);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* sin, n = 4, r = 36: 2 levels, widths [3, 2] */
static void
hs_sin_4_36(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    w1 = UWORD(0x0222222222222222); w0 = UWORD(0x2222222222222222);

    /* V_0 = c - z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 3);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* sin, n = 4, r = 37: 2 levels, widths [3, 2] */
static void
hs_sin_4_37(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    w1 = UWORD(0x0222222222222222); w0 = UWORD(0x2222222222222222);

    /* V_0 = c - z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 3);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* sin, n = 4, r = 38: 2 levels, widths [3, 1] */
static void
hs_sin_4_38(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0222222222222222);

    /* V_0 = c - z V_1   (3 x 1 limbs) */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = v;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 3);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* sin, n = 4, r = 39: 2 levels, widths [3, 1] */
static void
hs_sin_4_39(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0222222222222222);

    /* V_0 = c - z V_1   (3 x 1 limbs) */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = v;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 3);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* sin, n = 4, r = 40: 2 levels, widths [3, 1] */
static void
hs_sin_4_40(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0222222222222222);

    /* V_0 = c - z V_1   (3 x 1 limbs) */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = v;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x2aaaaaaaaaaaaaaa);
    mpn_sub_n(V, P, T, 3);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* sin, n = 4, r = 41: 2 levels, widths [2, 1] */
static void
hs_sin_4_41(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0222222222222222);

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x2aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = w0; P[3] = w1;
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* sin, n = 4, r = 42: 2 levels, widths [2, 1] */
static void
hs_sin_4_42(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0222222222222222);

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x2aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = w0; P[3] = w1;
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* sin, n = 4, r = 43: 2 levels, widths [2, 1] */
static void
hs_sin_4_43(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0222222222222222);

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x2aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = w0; P[3] = w1;
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* sin, n = 4, r = 44: 2 levels, widths [2, 1] */
static void
hs_sin_4_44(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0222222222222222);

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x2aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = w0; P[3] = w1;
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* sin, n = 4, r = 45: 2 levels, widths [2, 1] */
static void
hs_sin_4_45(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0222222222222222);

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x2aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = w0; P[3] = w1;
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* sin, n = 4, r = 46: 2 levels, widths [2, 1] */
static void
hs_sin_4_46(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0222222222222222);

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x2aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = w0; P[3] = w1;
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* sin, n = 4, r = 47: 2 levels, widths [2, 1] */
static void
hs_sin_4_47(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0222222222222222);

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x2aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = w0; P[3] = w1;
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* sin, n = 4, r = 48: 2 levels, widths [2, 1] */
static void
hs_sin_4_48(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0222222222222222);

    /* V_0 = c - z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    sub_ddmmss(w1, w0, UWORD(0x2aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* sin = x - x^3 V_0 */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = w0; P[3] = w1;
    flint_mpn_mulhigh_n(T, z, x, 4);      /* x^3 */
    flint_mpn_mulhigh_n(V, T, P, 4);
    mpn_sub_n(res, x, V, 4);
}

/* cos, n = 1, r = 8: 3 levels, widths [1, 1, 1] */
static void
hs_cos_1_8(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x005b05b05b05b05b);

    /* V_1 = c + z V_2   (1 x 1 limbs) */
    v = n_mulhi(z[0], v);
    v = UWORD(0x0aaaaaaaaaaaaaaa) + v;

    /* V_0 = c + z V_1   (1 x 1 limbs) */
    v = n_mulhi(z[0], v);
    v = UWORD(0x8000000000000000) + v;

    /* cos = 1 - z V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, P, 1);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 1))
    {
        flint_mpn_zero(res, 1);
        res[1] = 1;
    }
    else
    {
        mpn_neg(res, T, 1);
        res[1] = 0;
    }
}

/* cos, n = 1, r = 9: 3 levels, widths [1, 1, 1] */
static void
hs_cos_1_9(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x005b05b05b05b05b);

    /* V_1 = c + z V_2   (1 x 1 limbs) */
    v = n_mulhi(z[0], v);
    v = UWORD(0x0aaaaaaaaaaaaaaa) + v;

    /* V_0 = c + z V_1   (1 x 1 limbs) */
    v = n_mulhi(z[0], v);
    v = UWORD(0x8000000000000000) + v;

    /* cos = 1 - z V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, P, 1);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 1))
    {
        flint_mpn_zero(res, 1);
        res[1] = 1;
    }
    else
    {
        mpn_neg(res, T, 1);
        res[1] = 0;
    }
}

/* cos, n = 1, r = 10: 2 levels, widths [1, 1] */
static void
hs_cos_1_10(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x0aaaaaaaaaaaaaaa);

    /* V_0 = c + z V_1   (1 x 1 limbs) */
    v = n_mulhi(z[0], v);
    v = UWORD(0x8000000000000000) + v;

    /* cos = 1 - z V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, P, 1);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 1))
    {
        flint_mpn_zero(res, 1);
        res[1] = 1;
    }
    else
    {
        mpn_neg(res, T, 1);
        res[1] = 0;
    }
}

/* cos, n = 1, r = 11: 2 levels, widths [1, 1] */
static void
hs_cos_1_11(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x0aaaaaaaaaaaaaaa);

    /* V_0 = c + z V_1   (1 x 1 limbs) */
    v = n_mulhi(z[0], v);
    v = UWORD(0x8000000000000000) + v;

    /* cos = 1 - z V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, P, 1);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 1))
    {
        flint_mpn_zero(res, 1);
        res[1] = 1;
    }
    else
    {
        mpn_neg(res, T, 1);
        res[1] = 0;
    }
}

/* cos, n = 1, r = 12: 2 levels, widths [1, 1] */
static void
hs_cos_1_12(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x0aaaaaaaaaaaaaaa);

    /* V_0 = c + z V_1   (1 x 1 limbs) */
    v = n_mulhi(z[0], v);
    v = UWORD(0x8000000000000000) + v;

    /* cos = 1 - z V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, P, 1);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 1))
    {
        flint_mpn_zero(res, 1);
        res[1] = 1;
    }
    else
    {
        mpn_neg(res, T, 1);
        res[1] = 0;
    }
}

/* cos, n = 1, r = 13: 2 levels, widths [1, 1] */
static void
hs_cos_1_13(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x0aaaaaaaaaaaaaaa);

    /* V_0 = c + z V_1   (1 x 1 limbs) */
    v = n_mulhi(z[0], v);
    v = UWORD(0x8000000000000000) + v;

    /* cos = 1 - z V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, P, 1);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 1))
    {
        flint_mpn_zero(res, 1);
        res[1] = 1;
    }
    else
    {
        mpn_neg(res, T, 1);
        res[1] = 0;
    }
}

/* cos, n = 1, r = 14: 2 levels, widths [1, 1] */
static void
hs_cos_1_14(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x0aaaaaaaaaaaaaaa);

    /* V_0 = c + z V_1   (1 x 1 limbs) */
    v = n_mulhi(z[0], v);
    v = UWORD(0x8000000000000000) + v;

    /* cos = 1 - z V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, P, 1);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 1))
    {
        flint_mpn_zero(res, 1);
        res[1] = 1;
    }
    else
    {
        mpn_neg(res, T, 1);
        res[1] = 0;
    }
}

/* cos, n = 1, r = 15: 1 levels, widths [1] */
static void
hs_cos_1_15(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, P, 1);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 1))
    {
        flint_mpn_zero(res, 1);
        res[1] = 1;
    }
    else
    {
        mpn_neg(res, T, 1);
        res[1] = 0;
    }
}

/* cos, n = 1, r = 16: 1 levels, widths [1] */
static void
hs_cos_1_16(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, P, 1);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 1))
    {
        flint_mpn_zero(res, 1);
        res[1] = 1;
    }
    else
    {
        mpn_neg(res, T, 1);
        res[1] = 0;
    }
}

/* cos, n = 1, r = 17: 1 levels, widths [1] */
static void
hs_cos_1_17(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, P, 1);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 1))
    {
        flint_mpn_zero(res, 1);
        res[1] = 1;
    }
    else
    {
        mpn_neg(res, T, 1);
        res[1] = 0;
    }
}

/* cos, n = 1, r = 18: 1 levels, widths [1] */
static void
hs_cos_1_18(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, P, 1);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 1))
    {
        flint_mpn_zero(res, 1);
        res[1] = 1;
    }
    else
    {
        mpn_neg(res, T, 1);
        res[1] = 0;
    }
}

/* cos, n = 1, r = 19: 1 levels, widths [1] */
static void
hs_cos_1_19(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, P, 1);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 1))
    {
        flint_mpn_zero(res, 1);
        res[1] = 1;
    }
    else
    {
        mpn_neg(res, T, 1);
        res[1] = 0;
    }
}

/* cos, n = 1, r = 20: 1 levels, widths [1] */
static void
hs_cos_1_20(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, P, 1);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 1))
    {
        flint_mpn_zero(res, 1);
        res[1] = 1;
    }
    else
    {
        mpn_neg(res, T, 1);
        res[1] = 0;
    }
}

/* cos, n = 1, r = 21: 1 levels, widths [1] */
static void
hs_cos_1_21(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, P, 1);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 1))
    {
        flint_mpn_zero(res, 1);
        res[1] = 1;
    }
    else
    {
        mpn_neg(res, T, 1);
        res[1] = 0;
    }
}

/* cos, n = 1, r = 22: 1 levels, widths [1] */
static void
hs_cos_1_22(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, P, 1);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 1))
    {
        flint_mpn_zero(res, 1);
        res[1] = 1;
    }
    else
    {
        mpn_neg(res, T, 1);
        res[1] = 0;
    }
}

/* cos, n = 1, r = 23: 1 levels, widths [1] */
static void
hs_cos_1_23(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, P, 1);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 1))
    {
        flint_mpn_zero(res, 1);
        res[1] = 1;
    }
    else
    {
        mpn_neg(res, T, 1);
        res[1] = 0;
    }
}

/* cos, n = 1, r = 24: 1 levels, widths [1] */
static void
hs_cos_1_24(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, P, 1);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 1))
    {
        flint_mpn_zero(res, 1);
        res[1] = 1;
    }
    else
    {
        mpn_neg(res, T, 1);
        res[1] = 0;
    }
}

/* cos, n = 1, r = 25: 1 levels, widths [1] */
static void
hs_cos_1_25(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, P, 1);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 1))
    {
        flint_mpn_zero(res, 1);
        res[1] = 1;
    }
    else
    {
        mpn_neg(res, T, 1);
        res[1] = 0;
    }
}

/* cos, n = 1, r = 26: 1 levels, widths [1] */
static void
hs_cos_1_26(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, P, 1);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 1))
    {
        flint_mpn_zero(res, 1);
        res[1] = 1;
    }
    else
    {
        mpn_neg(res, T, 1);
        res[1] = 0;
    }
}

/* cos, n = 1, r = 27: 1 levels, widths [1] */
static void
hs_cos_1_27(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, P, 1);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 1))
    {
        flint_mpn_zero(res, 1);
        res[1] = 1;
    }
    else
    {
        mpn_neg(res, T, 1);
        res[1] = 0;
    }
}

/* cos, n = 1, r = 28: 1 levels, widths [1] */
static void
hs_cos_1_28(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, P, 1);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 1))
    {
        flint_mpn_zero(res, 1);
        res[1] = 1;
    }
    else
    {
        mpn_neg(res, T, 1);
        res[1] = 0;
    }
}

/* cos, n = 1, r = 29: 1 levels, widths [1] */
static void
hs_cos_1_29(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, P, 1);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 1))
    {
        flint_mpn_zero(res, 1);
        res[1] = 1;
    }
    else
    {
        mpn_neg(res, T, 1);
        res[1] = 0;
    }
}

/* cos, n = 1, r = 30: 1 levels, widths [1] */
static void
hs_cos_1_30(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, P, 1);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 1))
    {
        flint_mpn_zero(res, 1);
        res[1] = 1;
    }
    else
    {
        mpn_neg(res, T, 1);
        res[1] = 0;
    }
}

/* cos, n = 1, r = 31: 1 levels, widths [1] */
static void
hs_cos_1_31(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, P, 1);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 1))
    {
        flint_mpn_zero(res, 1);
        res[1] = 1;
    }
    else
    {
        mpn_neg(res, T, 1);
        res[1] = 0;
    }
}

/* cos, n = 1, r = 32: 1 levels, widths [1] */
static void
hs_cos_1_32(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, P, 1);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 1))
    {
        flint_mpn_zero(res, 1);
        res[1] = 1;
    }
    else
    {
        mpn_neg(res, T, 1);
        res[1] = 0;
    }
}

/* cos, n = 1, r = 33: 1 levels, widths [1] */
static void
hs_cos_1_33(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, P, 1);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 1))
    {
        flint_mpn_zero(res, 1);
        res[1] = 1;
    }
    else
    {
        mpn_neg(res, T, 1);
        res[1] = 0;
    }
}

/* cos, n = 1, r = 34: 1 levels, widths [1] */
static void
hs_cos_1_34(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, P, 1);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 1))
    {
        flint_mpn_zero(res, 1);
        res[1] = 1;
    }
    else
    {
        mpn_neg(res, T, 1);
        res[1] = 0;
    }
}

/* cos, n = 1, r = 35: 1 levels, widths [1] */
static void
hs_cos_1_35(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, P, 1);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 1))
    {
        flint_mpn_zero(res, 1);
        res[1] = 1;
    }
    else
    {
        mpn_neg(res, T, 1);
        res[1] = 0;
    }
}

/* cos, n = 1, r = 36: 1 levels, widths [1] */
static void
hs_cos_1_36(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, P, 1);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 1))
    {
        flint_mpn_zero(res, 1);
        res[1] = 1;
    }
    else
    {
        mpn_neg(res, T, 1);
        res[1] = 0;
    }
}

/* cos, n = 1, r = 37: 1 levels, widths [1] */
static void
hs_cos_1_37(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, P, 1);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 1))
    {
        flint_mpn_zero(res, 1);
        res[1] = 1;
    }
    else
    {
        mpn_neg(res, T, 1);
        res[1] = 0;
    }
}

/* cos, n = 1, r = 38: 1 levels, widths [1] */
static void
hs_cos_1_38(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, P, 1);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 1))
    {
        flint_mpn_zero(res, 1);
        res[1] = 1;
    }
    else
    {
        mpn_neg(res, T, 1);
        res[1] = 0;
    }
}

/* cos, n = 1, r = 39: 1 levels, widths [1] */
static void
hs_cos_1_39(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, P, 1);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 1))
    {
        flint_mpn_zero(res, 1);
        res[1] = 1;
    }
    else
    {
        mpn_neg(res, T, 1);
        res[1] = 0;
    }
}

/* cos, n = 1, r = 40: 1 levels, widths [1] */
static void
hs_cos_1_40(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, P, 1);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 1))
    {
        flint_mpn_zero(res, 1);
        res[1] = 1;
    }
    else
    {
        mpn_neg(res, T, 1);
        res[1] = 0;
    }
}

/* cos, n = 1, r = 41: 1 levels, widths [1] */
static void
hs_cos_1_41(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, P, 1);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 1))
    {
        flint_mpn_zero(res, 1);
        res[1] = 1;
    }
    else
    {
        mpn_neg(res, T, 1);
        res[1] = 0;
    }
}

/* cos, n = 1, r = 42: 1 levels, widths [1] */
static void
hs_cos_1_42(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, P, 1);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 1))
    {
        flint_mpn_zero(res, 1);
        res[1] = 1;
    }
    else
    {
        mpn_neg(res, T, 1);
        res[1] = 0;
    }
}

/* cos, n = 1, r = 43: 1 levels, widths [1] */
static void
hs_cos_1_43(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, P, 1);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 1))
    {
        flint_mpn_zero(res, 1);
        res[1] = 1;
    }
    else
    {
        mpn_neg(res, T, 1);
        res[1] = 0;
    }
}

/* cos, n = 1, r = 44: 1 levels, widths [1] */
static void
hs_cos_1_44(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, P, 1);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 1))
    {
        flint_mpn_zero(res, 1);
        res[1] = 1;
    }
    else
    {
        mpn_neg(res, T, 1);
        res[1] = 0;
    }
}

/* cos, n = 1, r = 45: 1 levels, widths [1] */
static void
hs_cos_1_45(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, P, 1);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 1))
    {
        flint_mpn_zero(res, 1);
        res[1] = 1;
    }
    else
    {
        mpn_neg(res, T, 1);
        res[1] = 0;
    }
}

/* cos, n = 1, r = 46: 1 levels, widths [1] */
static void
hs_cos_1_46(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, P, 1);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 1))
    {
        flint_mpn_zero(res, 1);
        res[1] = 1;
    }
    else
    {
        mpn_neg(res, T, 1);
        res[1] = 0;
    }
}

/* cos, n = 1, r = 47: 1 levels, widths [1] */
static void
hs_cos_1_47(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, P, 1);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 1))
    {
        flint_mpn_zero(res, 1);
        res[1] = 1;
    }
    else
    {
        mpn_neg(res, T, 1);
        res[1] = 0;
    }
}

/* cos, n = 1, r = 48: 1 levels, widths [1] */
static void
hs_cos_1_48(nn_ptr res, nn_srcptr x)
{
    ulong z[1], V[1], T[1], P[1];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 1);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = v;
    flint_mpn_mulhigh_n(T, z, P, 1);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 1))
    {
        flint_mpn_zero(res, 1);
        res[1] = 1;
    }
    else
    {
        mpn_neg(res, T, 1);
        res[1] = 0;
    }
}

/* cos, n = 2, r = 8: 6 levels, widths [2, 2, 2, 1, 1, 1] */
static void
hs_cos_2_8(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x00000008f76c77fc);

    /* V_4 = c + z V_5   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x0000049f93edde27) + v;

    /* V_3 = c + z V_4   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x0001a01a01a01a01) + v;

    /* V_2 = c + z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    add_ssaaaa(w1, w0, UWORD(0x005b05b05b05b05b), UWORD(0x05b05b05b05b05b0), h, l);

    /* V_1 = c + z V_2   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[1], z[0], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x0aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* V_0 = c + z V_1   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[1], z[0], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x8000000000000000), UWORD(0x0000000000000000), h, l);

    /* cos = 1 - z V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, P, 2);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 2))
    {
        flint_mpn_zero(res, 2);
        res[2] = 1;
    }
    else
    {
        mpn_neg(res, T, 2);
        res[2] = 0;
    }
}

/* cos, n = 2, r = 9: 5 levels, widths [2, 2, 2, 1, 1] */
static void
hs_cos_2_9(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x0000049f93edde27);

    /* V_3 = c + z V_4   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x0001a01a01a01a01) + v;

    /* V_2 = c + z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    add_ssaaaa(w1, w0, UWORD(0x005b05b05b05b05b), UWORD(0x05b05b05b05b05b0), h, l);

    /* V_1 = c + z V_2   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[1], z[0], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x0aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* V_0 = c + z V_1   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[1], z[0], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x8000000000000000), UWORD(0x0000000000000000), h, l);

    /* cos = 1 - z V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, P, 2);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 2))
    {
        flint_mpn_zero(res, 2);
        res[2] = 1;
    }
    else
    {
        mpn_neg(res, T, 2);
        res[2] = 0;
    }
}

/* cos, n = 2, r = 10: 5 levels, widths [2, 2, 1, 1, 1] */
static void
hs_cos_2_10(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x0000049f93edde27);

    /* V_3 = c + z V_4   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x0001a01a01a01a01) + v;

    /* V_2 = c + z V_3   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x005b05b05b05b05b) + v;

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    add_ssaaaa(w1, w0, UWORD(0x0aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* V_0 = c + z V_1   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[1], z[0], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x8000000000000000), UWORD(0x0000000000000000), h, l);

    /* cos = 1 - z V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, P, 2);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 2))
    {
        flint_mpn_zero(res, 2);
        res[2] = 1;
    }
    else
    {
        mpn_neg(res, T, 2);
        res[2] = 0;
    }
}

/* cos, n = 2, r = 11: 4 levels, widths [2, 2, 1, 1] */
static void
hs_cos_2_11(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x0001a01a01a01a01);

    /* V_2 = c + z V_3   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x005b05b05b05b05b) + v;

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    add_ssaaaa(w1, w0, UWORD(0x0aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* V_0 = c + z V_1   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[1], z[0], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x8000000000000000), UWORD(0x0000000000000000), h, l);

    /* cos = 1 - z V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, P, 2);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 2))
    {
        flint_mpn_zero(res, 2);
        res[2] = 1;
    }
    else
    {
        mpn_neg(res, T, 2);
        res[2] = 0;
    }
}

/* cos, n = 2, r = 12: 4 levels, widths [2, 2, 1, 1] */
static void
hs_cos_2_12(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x0001a01a01a01a01);

    /* V_2 = c + z V_3   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x005b05b05b05b05b) + v;

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    add_ssaaaa(w1, w0, UWORD(0x0aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* V_0 = c + z V_1   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[1], z[0], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x8000000000000000), UWORD(0x0000000000000000), h, l);

    /* cos = 1 - z V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, P, 2);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 2))
    {
        flint_mpn_zero(res, 2);
        res[2] = 1;
    }
    else
    {
        mpn_neg(res, T, 2);
        res[2] = 0;
    }
}

/* cos, n = 2, r = 13: 4 levels, widths [2, 2, 1, 1] */
static void
hs_cos_2_13(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x0001a01a01a01a01);

    /* V_2 = c + z V_3   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x005b05b05b05b05b) + v;

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    add_ssaaaa(w1, w0, UWORD(0x0aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* V_0 = c + z V_1   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[1], z[0], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x8000000000000000), UWORD(0x0000000000000000), h, l);

    /* cos = 1 - z V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, P, 2);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 2))
    {
        flint_mpn_zero(res, 2);
        res[2] = 1;
    }
    else
    {
        mpn_neg(res, T, 2);
        res[2] = 0;
    }
}

/* cos, n = 2, r = 14: 4 levels, widths [2, 2, 1, 1] */
static void
hs_cos_2_14(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x0001a01a01a01a01);

    /* V_2 = c + z V_3   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x005b05b05b05b05b) + v;

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    add_ssaaaa(w1, w0, UWORD(0x0aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* V_0 = c + z V_1   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[1], z[0], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x8000000000000000), UWORD(0x0000000000000000), h, l);

    /* cos = 1 - z V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, P, 2);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 2))
    {
        flint_mpn_zero(res, 2);
        res[2] = 1;
    }
    else
    {
        mpn_neg(res, T, 2);
        res[2] = 0;
    }
}

/* cos, n = 2, r = 15: 3 levels, widths [2, 1, 1] */
static void
hs_cos_2_15(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x005b05b05b05b05b);

    /* V_1 = c + z V_2   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x0aaaaaaaaaaaaaaa) + v;

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    add_ssaaaa(w1, w0, UWORD(0x8000000000000000), UWORD(0x0000000000000000), h, l);

    /* cos = 1 - z V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, P, 2);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 2))
    {
        flint_mpn_zero(res, 2);
        res[2] = 1;
    }
    else
    {
        mpn_neg(res, T, 2);
        res[2] = 0;
    }
}

/* cos, n = 2, r = 16: 3 levels, widths [2, 1, 1] */
static void
hs_cos_2_16(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x005b05b05b05b05b);

    /* V_1 = c + z V_2   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x0aaaaaaaaaaaaaaa) + v;

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    add_ssaaaa(w1, w0, UWORD(0x8000000000000000), UWORD(0x0000000000000000), h, l);

    /* cos = 1 - z V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, P, 2);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 2))
    {
        flint_mpn_zero(res, 2);
        res[2] = 1;
    }
    else
    {
        mpn_neg(res, T, 2);
        res[2] = 0;
    }
}

/* cos, n = 2, r = 17: 3 levels, widths [2, 1, 1] */
static void
hs_cos_2_17(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x005b05b05b05b05b);

    /* V_1 = c + z V_2   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x0aaaaaaaaaaaaaaa) + v;

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    add_ssaaaa(w1, w0, UWORD(0x8000000000000000), UWORD(0x0000000000000000), h, l);

    /* cos = 1 - z V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, P, 2);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 2))
    {
        flint_mpn_zero(res, 2);
        res[2] = 1;
    }
    else
    {
        mpn_neg(res, T, 2);
        res[2] = 0;
    }
}

/* cos, n = 2, r = 18: 3 levels, widths [2, 1, 1] */
static void
hs_cos_2_18(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x005b05b05b05b05b);

    /* V_1 = c + z V_2   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x0aaaaaaaaaaaaaaa) + v;

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    add_ssaaaa(w1, w0, UWORD(0x8000000000000000), UWORD(0x0000000000000000), h, l);

    /* cos = 1 - z V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, P, 2);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 2))
    {
        flint_mpn_zero(res, 2);
        res[2] = 1;
    }
    else
    {
        mpn_neg(res, T, 2);
        res[2] = 0;
    }
}

/* cos, n = 2, r = 19: 3 levels, widths [2, 1, 1] */
static void
hs_cos_2_19(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x005b05b05b05b05b);

    /* V_1 = c + z V_2   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x0aaaaaaaaaaaaaaa) + v;

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    add_ssaaaa(w1, w0, UWORD(0x8000000000000000), UWORD(0x0000000000000000), h, l);

    /* cos = 1 - z V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, P, 2);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 2))
    {
        flint_mpn_zero(res, 2);
        res[2] = 1;
    }
    else
    {
        mpn_neg(res, T, 2);
        res[2] = 0;
    }
}

/* cos, n = 2, r = 20: 2 levels, widths [2, 1] */
static void
hs_cos_2_20(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x0aaaaaaaaaaaaaaa);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    add_ssaaaa(w1, w0, UWORD(0x8000000000000000), UWORD(0x0000000000000000), h, l);

    /* cos = 1 - z V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, P, 2);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 2))
    {
        flint_mpn_zero(res, 2);
        res[2] = 1;
    }
    else
    {
        mpn_neg(res, T, 2);
        res[2] = 0;
    }
}

/* cos, n = 2, r = 21: 2 levels, widths [2, 1] */
static void
hs_cos_2_21(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x0aaaaaaaaaaaaaaa);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    add_ssaaaa(w1, w0, UWORD(0x8000000000000000), UWORD(0x0000000000000000), h, l);

    /* cos = 1 - z V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, P, 2);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 2))
    {
        flint_mpn_zero(res, 2);
        res[2] = 1;
    }
    else
    {
        mpn_neg(res, T, 2);
        res[2] = 0;
    }
}

/* cos, n = 2, r = 22: 2 levels, widths [2, 1] */
static void
hs_cos_2_22(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x0aaaaaaaaaaaaaaa);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    add_ssaaaa(w1, w0, UWORD(0x8000000000000000), UWORD(0x0000000000000000), h, l);

    /* cos = 1 - z V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, P, 2);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 2))
    {
        flint_mpn_zero(res, 2);
        res[2] = 1;
    }
    else
    {
        mpn_neg(res, T, 2);
        res[2] = 0;
    }
}

/* cos, n = 2, r = 23: 2 levels, widths [2, 1] */
static void
hs_cos_2_23(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x0aaaaaaaaaaaaaaa);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    add_ssaaaa(w1, w0, UWORD(0x8000000000000000), UWORD(0x0000000000000000), h, l);

    /* cos = 1 - z V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, P, 2);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 2))
    {
        flint_mpn_zero(res, 2);
        res[2] = 1;
    }
    else
    {
        mpn_neg(res, T, 2);
        res[2] = 0;
    }
}

/* cos, n = 2, r = 24: 2 levels, widths [2, 1] */
static void
hs_cos_2_24(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x0aaaaaaaaaaaaaaa);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    add_ssaaaa(w1, w0, UWORD(0x8000000000000000), UWORD(0x0000000000000000), h, l);

    /* cos = 1 - z V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, P, 2);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 2))
    {
        flint_mpn_zero(res, 2);
        res[2] = 1;
    }
    else
    {
        mpn_neg(res, T, 2);
        res[2] = 0;
    }
}

/* cos, n = 2, r = 25: 2 levels, widths [2, 1] */
static void
hs_cos_2_25(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x0aaaaaaaaaaaaaaa);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    add_ssaaaa(w1, w0, UWORD(0x8000000000000000), UWORD(0x0000000000000000), h, l);

    /* cos = 1 - z V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, P, 2);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 2))
    {
        flint_mpn_zero(res, 2);
        res[2] = 1;
    }
    else
    {
        mpn_neg(res, T, 2);
        res[2] = 0;
    }
}

/* cos, n = 2, r = 26: 2 levels, widths [2, 1] */
static void
hs_cos_2_26(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x0aaaaaaaaaaaaaaa);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    add_ssaaaa(w1, w0, UWORD(0x8000000000000000), UWORD(0x0000000000000000), h, l);

    /* cos = 1 - z V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, P, 2);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 2))
    {
        flint_mpn_zero(res, 2);
        res[2] = 1;
    }
    else
    {
        mpn_neg(res, T, 2);
        res[2] = 0;
    }
}

/* cos, n = 2, r = 27: 2 levels, widths [2, 1] */
static void
hs_cos_2_27(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x0aaaaaaaaaaaaaaa);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    add_ssaaaa(w1, w0, UWORD(0x8000000000000000), UWORD(0x0000000000000000), h, l);

    /* cos = 1 - z V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, P, 2);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 2))
    {
        flint_mpn_zero(res, 2);
        res[2] = 1;
    }
    else
    {
        mpn_neg(res, T, 2);
        res[2] = 0;
    }
}

/* cos, n = 2, r = 28: 2 levels, widths [2, 1] */
static void
hs_cos_2_28(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x0aaaaaaaaaaaaaaa);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[1], z[0], v);
    add_ssaaaa(w1, w0, UWORD(0x8000000000000000), UWORD(0x0000000000000000), h, l);

    /* cos = 1 - z V_0 */
    P[0] = w0; P[1] = w1;
    flint_mpn_mulhigh_n(T, z, P, 2);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 2))
    {
        flint_mpn_zero(res, 2);
        res[2] = 1;
    }
    else
    {
        mpn_neg(res, T, 2);
        res[2] = 0;
    }
}

/* cos, n = 2, r = 29: 2 levels, widths [1, 1] */
static void
hs_cos_2_29(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x0aaaaaaaaaaaaaaa);

    /* V_0 = c + z V_1   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x8000000000000000) + v;

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, P, 2);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 2))
    {
        flint_mpn_zero(res, 2);
        res[2] = 1;
    }
    else
    {
        mpn_neg(res, T, 2);
        res[2] = 0;
    }
}

/* cos, n = 2, r = 30: 2 levels, widths [1, 1] */
static void
hs_cos_2_30(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x0aaaaaaaaaaaaaaa);

    /* V_0 = c + z V_1   (1 x 1 limbs) */
    v = n_mulhi(z[1], v);
    v = UWORD(0x8000000000000000) + v;

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, P, 2);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 2))
    {
        flint_mpn_zero(res, 2);
        res[2] = 1;
    }
    else
    {
        mpn_neg(res, T, 2);
        res[2] = 0;
    }
}

/* cos, n = 2, r = 31: 1 levels, widths [1] */
static void
hs_cos_2_31(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, P, 2);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 2))
    {
        flint_mpn_zero(res, 2);
        res[2] = 1;
    }
    else
    {
        mpn_neg(res, T, 2);
        res[2] = 0;
    }
}

/* cos, n = 2, r = 32: 1 levels, widths [1] */
static void
hs_cos_2_32(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, P, 2);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 2))
    {
        flint_mpn_zero(res, 2);
        res[2] = 1;
    }
    else
    {
        mpn_neg(res, T, 2);
        res[2] = 0;
    }
}

/* cos, n = 2, r = 33: 1 levels, widths [1] */
static void
hs_cos_2_33(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, P, 2);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 2))
    {
        flint_mpn_zero(res, 2);
        res[2] = 1;
    }
    else
    {
        mpn_neg(res, T, 2);
        res[2] = 0;
    }
}

/* cos, n = 2, r = 34: 1 levels, widths [1] */
static void
hs_cos_2_34(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, P, 2);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 2))
    {
        flint_mpn_zero(res, 2);
        res[2] = 1;
    }
    else
    {
        mpn_neg(res, T, 2);
        res[2] = 0;
    }
}

/* cos, n = 2, r = 35: 1 levels, widths [1] */
static void
hs_cos_2_35(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, P, 2);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 2))
    {
        flint_mpn_zero(res, 2);
        res[2] = 1;
    }
    else
    {
        mpn_neg(res, T, 2);
        res[2] = 0;
    }
}

/* cos, n = 2, r = 36: 1 levels, widths [1] */
static void
hs_cos_2_36(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, P, 2);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 2))
    {
        flint_mpn_zero(res, 2);
        res[2] = 1;
    }
    else
    {
        mpn_neg(res, T, 2);
        res[2] = 0;
    }
}

/* cos, n = 2, r = 37: 1 levels, widths [1] */
static void
hs_cos_2_37(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, P, 2);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 2))
    {
        flint_mpn_zero(res, 2);
        res[2] = 1;
    }
    else
    {
        mpn_neg(res, T, 2);
        res[2] = 0;
    }
}

/* cos, n = 2, r = 38: 1 levels, widths [1] */
static void
hs_cos_2_38(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, P, 2);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 2))
    {
        flint_mpn_zero(res, 2);
        res[2] = 1;
    }
    else
    {
        mpn_neg(res, T, 2);
        res[2] = 0;
    }
}

/* cos, n = 2, r = 39: 1 levels, widths [1] */
static void
hs_cos_2_39(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, P, 2);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 2))
    {
        flint_mpn_zero(res, 2);
        res[2] = 1;
    }
    else
    {
        mpn_neg(res, T, 2);
        res[2] = 0;
    }
}

/* cos, n = 2, r = 40: 1 levels, widths [1] */
static void
hs_cos_2_40(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, P, 2);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 2))
    {
        flint_mpn_zero(res, 2);
        res[2] = 1;
    }
    else
    {
        mpn_neg(res, T, 2);
        res[2] = 0;
    }
}

/* cos, n = 2, r = 41: 1 levels, widths [1] */
static void
hs_cos_2_41(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, P, 2);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 2))
    {
        flint_mpn_zero(res, 2);
        res[2] = 1;
    }
    else
    {
        mpn_neg(res, T, 2);
        res[2] = 0;
    }
}

/* cos, n = 2, r = 42: 1 levels, widths [1] */
static void
hs_cos_2_42(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, P, 2);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 2))
    {
        flint_mpn_zero(res, 2);
        res[2] = 1;
    }
    else
    {
        mpn_neg(res, T, 2);
        res[2] = 0;
    }
}

/* cos, n = 2, r = 43: 1 levels, widths [1] */
static void
hs_cos_2_43(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, P, 2);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 2))
    {
        flint_mpn_zero(res, 2);
        res[2] = 1;
    }
    else
    {
        mpn_neg(res, T, 2);
        res[2] = 0;
    }
}

/* cos, n = 2, r = 44: 1 levels, widths [1] */
static void
hs_cos_2_44(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, P, 2);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 2))
    {
        flint_mpn_zero(res, 2);
        res[2] = 1;
    }
    else
    {
        mpn_neg(res, T, 2);
        res[2] = 0;
    }
}

/* cos, n = 2, r = 45: 1 levels, widths [1] */
static void
hs_cos_2_45(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, P, 2);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 2))
    {
        flint_mpn_zero(res, 2);
        res[2] = 1;
    }
    else
    {
        mpn_neg(res, T, 2);
        res[2] = 0;
    }
}

/* cos, n = 2, r = 46: 1 levels, widths [1] */
static void
hs_cos_2_46(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, P, 2);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 2))
    {
        flint_mpn_zero(res, 2);
        res[2] = 1;
    }
    else
    {
        mpn_neg(res, T, 2);
        res[2] = 0;
    }
}

/* cos, n = 2, r = 47: 1 levels, widths [1] */
static void
hs_cos_2_47(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, P, 2);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 2))
    {
        flint_mpn_zero(res, 2);
        res[2] = 1;
    }
    else
    {
        mpn_neg(res, T, 2);
        res[2] = 0;
    }
}

/* cos, n = 2, r = 48: 1 levels, widths [1] */
static void
hs_cos_2_48(nn_ptr res, nn_srcptr x)
{
    ulong z[2], V[2], T[2], P[2];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 2);

    v = UWORD(0x8000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = v;
    flint_mpn_mulhigh_n(T, z, P, 2);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 2))
    {
        flint_mpn_zero(res, 2);
        res[2] = 1;
    }
    else
    {
        mpn_neg(res, T, 2);
        res[2] = 0;
    }
}

/* cos, n = 3, r = 8: 8 levels, widths [3, 3, 3, 2, 2, 2, 2, 1] */
static void
hs_cos_3_8(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x00000000000d73f9);

    /* V_6 = c + z V_7   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x000000000c9cba54), UWORD(0x603e4e905d6f8a2e), h, l);

    /* V_5 = c + z V_6   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x00000008f76c77fc), UWORD(0x6c4bdaa26d4c3d67), h, l);

    /* V_4 = c + z V_5   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x0000049f93edde27), UWORD(0xd71cbbc05b4fa999), h, l);

    /* V_3 = c + z V_4   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x0001a01a01a01a01), UWORD(0xa01a01a01a01a01a), h, l);

    /* V_2 = c + z V_3   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x5b05b05b05b05b05);
    P[1] = UWORD(0x05b05b05b05b05b0);
    P[2] = UWORD(0x005b05b05b05b05b);
    mpn_add_n(V, P, T, 3);

    /* V_1 = c + z V_2   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x0aaaaaaaaaaaaaaa);
    mpn_add_n(V, P, T, 3);

    /* V_0 = c + z V_1   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 3);

    /* cos = 1 - z V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, P, 3);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 3))
    {
        flint_mpn_zero(res, 3);
        res[3] = 1;
    }
    else
    {
        mpn_neg(res, T, 3);
        res[3] = 0;
    }
}

/* cos, n = 3, r = 9: 8 levels, widths [3, 3, 3, 2, 2, 2, 1, 1] */
static void
hs_cos_3_9(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x00000000000d73f9);

    /* V_6 = c + z V_7   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x000000000c9cba54) + v;

    /* V_5 = c + z V_6   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x00000008f76c77fc), UWORD(0x6c4bdaa26d4c3d67), h, l);

    /* V_4 = c + z V_5   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x0000049f93edde27), UWORD(0xd71cbbc05b4fa999), h, l);

    /* V_3 = c + z V_4   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x0001a01a01a01a01), UWORD(0xa01a01a01a01a01a), h, l);

    /* V_2 = c + z V_3   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x5b05b05b05b05b05);
    P[1] = UWORD(0x05b05b05b05b05b0);
    P[2] = UWORD(0x005b05b05b05b05b);
    mpn_add_n(V, P, T, 3);

    /* V_1 = c + z V_2   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x0aaaaaaaaaaaaaaa);
    mpn_add_n(V, P, T, 3);

    /* V_0 = c + z V_1   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 3);

    /* cos = 1 - z V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, P, 3);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 3))
    {
        flint_mpn_zero(res, 3);
        res[3] = 1;
    }
    else
    {
        mpn_neg(res, T, 3);
        res[3] = 0;
    }
}

/* cos, n = 3, r = 10: 7 levels, widths [3, 3, 2, 2, 2, 2, 1] */
static void
hs_cos_3_10(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x000000000c9cba54);

    /* V_5 = c + z V_6   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x00000008f76c77fc), UWORD(0x6c4bdaa26d4c3d67), h, l);

    /* V_4 = c + z V_5   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x0000049f93edde27), UWORD(0xd71cbbc05b4fa999), h, l);

    /* V_3 = c + z V_4   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x0001a01a01a01a01), UWORD(0xa01a01a01a01a01a), h, l);

    /* V_2 = c + z V_3   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x005b05b05b05b05b), UWORD(0x05b05b05b05b05b0), h, l);

    /* V_1 = c + z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x0aaaaaaaaaaaaaaa);
    mpn_add_n(V, P, T, 3);

    /* V_0 = c + z V_1   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 3);

    /* cos = 1 - z V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, P, 3);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 3))
    {
        flint_mpn_zero(res, 3);
        res[3] = 1;
    }
    else
    {
        mpn_neg(res, T, 3);
        res[3] = 0;
    }
}

/* cos, n = 3, r = 11: 7 levels, widths [3, 3, 2, 2, 2, 1, 1] */
static void
hs_cos_3_11(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x000000000c9cba54);

    /* V_5 = c + z V_6   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x00000008f76c77fc) + v;

    /* V_4 = c + z V_5   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x0000049f93edde27), UWORD(0xd71cbbc05b4fa999), h, l);

    /* V_3 = c + z V_4   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x0001a01a01a01a01), UWORD(0xa01a01a01a01a01a), h, l);

    /* V_2 = c + z V_3   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x005b05b05b05b05b), UWORD(0x05b05b05b05b05b0), h, l);

    /* V_1 = c + z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x0aaaaaaaaaaaaaaa);
    mpn_add_n(V, P, T, 3);

    /* V_0 = c + z V_1   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 3);

    /* cos = 1 - z V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, P, 3);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 3))
    {
        flint_mpn_zero(res, 3);
        res[3] = 1;
    }
    else
    {
        mpn_neg(res, T, 3);
        res[3] = 0;
    }
}

/* cos, n = 3, r = 12: 6 levels, widths [3, 3, 2, 2, 2, 1] */
static void
hs_cos_3_12(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x00000008f76c77fc);

    /* V_4 = c + z V_5   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x0000049f93edde27), UWORD(0xd71cbbc05b4fa999), h, l);

    /* V_3 = c + z V_4   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x0001a01a01a01a01), UWORD(0xa01a01a01a01a01a), h, l);

    /* V_2 = c + z V_3   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x005b05b05b05b05b), UWORD(0x05b05b05b05b05b0), h, l);

    /* V_1 = c + z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x0aaaaaaaaaaaaaaa);
    mpn_add_n(V, P, T, 3);

    /* V_0 = c + z V_1   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 3);

    /* cos = 1 - z V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, P, 3);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 3))
    {
        flint_mpn_zero(res, 3);
        res[3] = 1;
    }
    else
    {
        mpn_neg(res, T, 3);
        res[3] = 0;
    }
}

/* cos, n = 3, r = 13: 6 levels, widths [3, 3, 2, 2, 1, 1] */
static void
hs_cos_3_13(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x00000008f76c77fc);

    /* V_4 = c + z V_5   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x0000049f93edde27) + v;

    /* V_3 = c + z V_4   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x0001a01a01a01a01), UWORD(0xa01a01a01a01a01a), h, l);

    /* V_2 = c + z V_3   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x005b05b05b05b05b), UWORD(0x05b05b05b05b05b0), h, l);

    /* V_1 = c + z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x0aaaaaaaaaaaaaaa);
    mpn_add_n(V, P, T, 3);

    /* V_0 = c + z V_1   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 3);

    /* cos = 1 - z V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, P, 3);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 3))
    {
        flint_mpn_zero(res, 3);
        res[3] = 1;
    }
    else
    {
        mpn_neg(res, T, 3);
        res[3] = 0;
    }
}

/* cos, n = 3, r = 14: 5 levels, widths [3, 3, 2, 2, 1] */
static void
hs_cos_3_14(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0000049f93edde27);

    /* V_3 = c + z V_4   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x0001a01a01a01a01), UWORD(0xa01a01a01a01a01a), h, l);

    /* V_2 = c + z V_3   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x005b05b05b05b05b), UWORD(0x05b05b05b05b05b0), h, l);

    /* V_1 = c + z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x0aaaaaaaaaaaaaaa);
    mpn_add_n(V, P, T, 3);

    /* V_0 = c + z V_1   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 3);

    /* cos = 1 - z V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, P, 3);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 3))
    {
        flint_mpn_zero(res, 3);
        res[3] = 1;
    }
    else
    {
        mpn_neg(res, T, 3);
        res[3] = 0;
    }
}

/* cos, n = 3, r = 15: 5 levels, widths [3, 2, 2, 2, 1] */
static void
hs_cos_3_15(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0000049f93edde27);

    /* V_3 = c + z V_4   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x0001a01a01a01a01), UWORD(0xa01a01a01a01a01a), h, l);

    /* V_2 = c + z V_3   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x005b05b05b05b05b), UWORD(0x05b05b05b05b05b0), h, l);

    /* V_1 = c + z V_2   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x0aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 3);

    /* cos = 1 - z V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, P, 3);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 3))
    {
        flint_mpn_zero(res, 3);
        res[3] = 1;
    }
    else
    {
        mpn_neg(res, T, 3);
        res[3] = 0;
    }
}

/* cos, n = 3, r = 16: 5 levels, widths [3, 2, 2, 1, 1] */
static void
hs_cos_3_16(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0000049f93edde27);

    /* V_3 = c + z V_4   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x0001a01a01a01a01) + v;

    /* V_2 = c + z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x005b05b05b05b05b), UWORD(0x05b05b05b05b05b0), h, l);

    /* V_1 = c + z V_2   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x0aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 3);

    /* cos = 1 - z V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, P, 3);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 3))
    {
        flint_mpn_zero(res, 3);
        res[3] = 1;
    }
    else
    {
        mpn_neg(res, T, 3);
        res[3] = 0;
    }
}

/* cos, n = 3, r = 17: 5 levels, widths [3, 2, 2, 1, 1] */
static void
hs_cos_3_17(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0000049f93edde27);

    /* V_3 = c + z V_4   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x0001a01a01a01a01) + v;

    /* V_2 = c + z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x005b05b05b05b05b), UWORD(0x05b05b05b05b05b0), h, l);

    /* V_1 = c + z V_2   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x0aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 3);

    /* cos = 1 - z V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, P, 3);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 3))
    {
        flint_mpn_zero(res, 3);
        res[3] = 1;
    }
    else
    {
        mpn_neg(res, T, 3);
        res[3] = 0;
    }
}

/* cos, n = 3, r = 18: 4 levels, widths [3, 2, 2, 1] */
static void
hs_cos_3_18(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0001a01a01a01a01);

    /* V_2 = c + z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x005b05b05b05b05b), UWORD(0x05b05b05b05b05b0), h, l);

    /* V_1 = c + z V_2   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x0aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 3);

    /* cos = 1 - z V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, P, 3);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 3))
    {
        flint_mpn_zero(res, 3);
        res[3] = 1;
    }
    else
    {
        mpn_neg(res, T, 3);
        res[3] = 0;
    }
}

/* cos, n = 3, r = 19: 4 levels, widths [3, 2, 2, 1] */
static void
hs_cos_3_19(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0001a01a01a01a01);

    /* V_2 = c + z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x005b05b05b05b05b), UWORD(0x05b05b05b05b05b0), h, l);

    /* V_1 = c + z V_2   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x0aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 3);

    /* cos = 1 - z V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, P, 3);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 3))
    {
        flint_mpn_zero(res, 3);
        res[3] = 1;
    }
    else
    {
        mpn_neg(res, T, 3);
        res[3] = 0;
    }
}

/* cos, n = 3, r = 20: 4 levels, widths [3, 2, 2, 1] */
static void
hs_cos_3_20(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0001a01a01a01a01);

    /* V_2 = c + z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x005b05b05b05b05b), UWORD(0x05b05b05b05b05b0), h, l);

    /* V_1 = c + z V_2   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x0aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 3);

    /* cos = 1 - z V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, P, 3);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 3))
    {
        flint_mpn_zero(res, 3);
        res[3] = 1;
    }
    else
    {
        mpn_neg(res, T, 3);
        res[3] = 0;
    }
}

/* cos, n = 3, r = 21: 4 levels, widths [3, 2, 1, 1] */
static void
hs_cos_3_21(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0001a01a01a01a01);

    /* V_2 = c + z V_3   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x005b05b05b05b05b) + v;

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x0aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 3);

    /* cos = 1 - z V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, P, 3);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 3))
    {
        flint_mpn_zero(res, 3);
        res[3] = 1;
    }
    else
    {
        mpn_neg(res, T, 3);
        res[3] = 0;
    }
}

/* cos, n = 3, r = 22: 4 levels, widths [3, 2, 1, 1] */
static void
hs_cos_3_22(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0001a01a01a01a01);

    /* V_2 = c + z V_3   (1 x 1 limbs) */
    v = n_mulhi(z[2], v);
    v = UWORD(0x005b05b05b05b05b) + v;

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x0aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 3);

    /* cos = 1 - z V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, P, 3);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 3))
    {
        flint_mpn_zero(res, 3);
        res[3] = 1;
    }
    else
    {
        mpn_neg(res, T, 3);
        res[3] = 0;
    }
}

/* cos, n = 3, r = 23: 3 levels, widths [3, 2, 1] */
static void
hs_cos_3_23(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x005b05b05b05b05b);

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x0aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 3);

    /* cos = 1 - z V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, P, 3);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 3))
    {
        flint_mpn_zero(res, 3);
        res[3] = 1;
    }
    else
    {
        mpn_neg(res, T, 3);
        res[3] = 0;
    }
}

/* cos, n = 3, r = 24: 3 levels, widths [3, 2, 1] */
static void
hs_cos_3_24(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x005b05b05b05b05b);

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x0aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 3);

    /* cos = 1 - z V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, P, 3);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 3))
    {
        flint_mpn_zero(res, 3);
        res[3] = 1;
    }
    else
    {
        mpn_neg(res, T, 3);
        res[3] = 0;
    }
}

/* cos, n = 3, r = 25: 3 levels, widths [3, 2, 1] */
static void
hs_cos_3_25(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x005b05b05b05b05b);

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x0aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 3);

    /* cos = 1 - z V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, P, 3);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 3))
    {
        flint_mpn_zero(res, 3);
        res[3] = 1;
    }
    else
    {
        mpn_neg(res, T, 3);
        res[3] = 0;
    }
}

/* cos, n = 3, r = 26: 3 levels, widths [3, 2, 1] */
static void
hs_cos_3_26(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x005b05b05b05b05b);

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x0aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 3);

    /* cos = 1 - z V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, P, 3);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 3))
    {
        flint_mpn_zero(res, 3);
        res[3] = 1;
    }
    else
    {
        mpn_neg(res, T, 3);
        res[3] = 0;
    }
}

/* cos, n = 3, r = 27: 3 levels, widths [3, 2, 1] */
static void
hs_cos_3_27(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x005b05b05b05b05b);

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x0aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 3);

    /* cos = 1 - z V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, P, 3);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 3))
    {
        flint_mpn_zero(res, 3);
        res[3] = 1;
    }
    else
    {
        mpn_neg(res, T, 3);
        res[3] = 0;
    }
}

/* cos, n = 3, r = 28: 3 levels, widths [3, 2, 1] */
static void
hs_cos_3_28(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x005b05b05b05b05b);

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x0aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 0, P, 3);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 3);

    /* cos = 1 - z V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z, P, 3);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 3))
    {
        flint_mpn_zero(res, 3);
        res[3] = 1;
    }
    else
    {
        mpn_neg(res, T, 3);
        res[3] = 0;
    }
}

/* cos, n = 3, r = 29: 3 levels, widths [2, 2, 1] */
static void
hs_cos_3_29(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x005b05b05b05b05b);

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x0aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* V_0 = c + z V_1   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x8000000000000000), UWORD(0x0000000000000000), h, l);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, P, 3);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 3))
    {
        flint_mpn_zero(res, 3);
        res[3] = 1;
    }
    else
    {
        mpn_neg(res, T, 3);
        res[3] = 0;
    }
}

/* cos, n = 3, r = 30: 3 levels, widths [2, 2, 1] */
static void
hs_cos_3_30(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x005b05b05b05b05b);

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x0aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* V_0 = c + z V_1   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[2], z[1], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x8000000000000000), UWORD(0x0000000000000000), h, l);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, P, 3);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 3))
    {
        flint_mpn_zero(res, 3);
        res[3] = 1;
    }
    else
    {
        mpn_neg(res, T, 3);
        res[3] = 0;
    }
}

/* cos, n = 3, r = 31: 2 levels, widths [2, 1] */
static void
hs_cos_3_31(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0aaaaaaaaaaaaaaa);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x8000000000000000), UWORD(0x0000000000000000), h, l);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, P, 3);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 3))
    {
        flint_mpn_zero(res, 3);
        res[3] = 1;
    }
    else
    {
        mpn_neg(res, T, 3);
        res[3] = 0;
    }
}

/* cos, n = 3, r = 32: 2 levels, widths [2, 1] */
static void
hs_cos_3_32(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0aaaaaaaaaaaaaaa);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x8000000000000000), UWORD(0x0000000000000000), h, l);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, P, 3);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 3))
    {
        flint_mpn_zero(res, 3);
        res[3] = 1;
    }
    else
    {
        mpn_neg(res, T, 3);
        res[3] = 0;
    }
}

/* cos, n = 3, r = 33: 2 levels, widths [2, 1] */
static void
hs_cos_3_33(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0aaaaaaaaaaaaaaa);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x8000000000000000), UWORD(0x0000000000000000), h, l);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, P, 3);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 3))
    {
        flint_mpn_zero(res, 3);
        res[3] = 1;
    }
    else
    {
        mpn_neg(res, T, 3);
        res[3] = 0;
    }
}

/* cos, n = 3, r = 34: 2 levels, widths [2, 1] */
static void
hs_cos_3_34(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0aaaaaaaaaaaaaaa);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x8000000000000000), UWORD(0x0000000000000000), h, l);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, P, 3);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 3))
    {
        flint_mpn_zero(res, 3);
        res[3] = 1;
    }
    else
    {
        mpn_neg(res, T, 3);
        res[3] = 0;
    }
}

/* cos, n = 3, r = 35: 2 levels, widths [2, 1] */
static void
hs_cos_3_35(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0aaaaaaaaaaaaaaa);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x8000000000000000), UWORD(0x0000000000000000), h, l);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, P, 3);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 3))
    {
        flint_mpn_zero(res, 3);
        res[3] = 1;
    }
    else
    {
        mpn_neg(res, T, 3);
        res[3] = 0;
    }
}

/* cos, n = 3, r = 36: 2 levels, widths [2, 1] */
static void
hs_cos_3_36(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0aaaaaaaaaaaaaaa);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x8000000000000000), UWORD(0x0000000000000000), h, l);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, P, 3);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 3))
    {
        flint_mpn_zero(res, 3);
        res[3] = 1;
    }
    else
    {
        mpn_neg(res, T, 3);
        res[3] = 0;
    }
}

/* cos, n = 3, r = 37: 2 levels, widths [2, 1] */
static void
hs_cos_3_37(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0aaaaaaaaaaaaaaa);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x8000000000000000), UWORD(0x0000000000000000), h, l);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, P, 3);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 3))
    {
        flint_mpn_zero(res, 3);
        res[3] = 1;
    }
    else
    {
        mpn_neg(res, T, 3);
        res[3] = 0;
    }
}

/* cos, n = 3, r = 38: 2 levels, widths [2, 1] */
static void
hs_cos_3_38(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0aaaaaaaaaaaaaaa);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x8000000000000000), UWORD(0x0000000000000000), h, l);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, P, 3);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 3))
    {
        flint_mpn_zero(res, 3);
        res[3] = 1;
    }
    else
    {
        mpn_neg(res, T, 3);
        res[3] = 0;
    }
}

/* cos, n = 3, r = 39: 2 levels, widths [2, 1] */
static void
hs_cos_3_39(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0aaaaaaaaaaaaaaa);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x8000000000000000), UWORD(0x0000000000000000), h, l);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, P, 3);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 3))
    {
        flint_mpn_zero(res, 3);
        res[3] = 1;
    }
    else
    {
        mpn_neg(res, T, 3);
        res[3] = 0;
    }
}

/* cos, n = 3, r = 40: 2 levels, widths [2, 1] */
static void
hs_cos_3_40(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0aaaaaaaaaaaaaaa);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x8000000000000000), UWORD(0x0000000000000000), h, l);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, P, 3);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 3))
    {
        flint_mpn_zero(res, 3);
        res[3] = 1;
    }
    else
    {
        mpn_neg(res, T, 3);
        res[3] = 0;
    }
}

/* cos, n = 3, r = 41: 2 levels, widths [2, 1] */
static void
hs_cos_3_41(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0aaaaaaaaaaaaaaa);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x8000000000000000), UWORD(0x0000000000000000), h, l);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, P, 3);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 3))
    {
        flint_mpn_zero(res, 3);
        res[3] = 1;
    }
    else
    {
        mpn_neg(res, T, 3);
        res[3] = 0;
    }
}

/* cos, n = 3, r = 42: 2 levels, widths [2, 1] */
static void
hs_cos_3_42(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0aaaaaaaaaaaaaaa);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x8000000000000000), UWORD(0x0000000000000000), h, l);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, P, 3);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 3))
    {
        flint_mpn_zero(res, 3);
        res[3] = 1;
    }
    else
    {
        mpn_neg(res, T, 3);
        res[3] = 0;
    }
}

/* cos, n = 3, r = 43: 2 levels, widths [2, 1] */
static void
hs_cos_3_43(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0aaaaaaaaaaaaaaa);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x8000000000000000), UWORD(0x0000000000000000), h, l);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, P, 3);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 3))
    {
        flint_mpn_zero(res, 3);
        res[3] = 1;
    }
    else
    {
        mpn_neg(res, T, 3);
        res[3] = 0;
    }
}

/* cos, n = 3, r = 44: 2 levels, widths [2, 1] */
static void
hs_cos_3_44(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0aaaaaaaaaaaaaaa);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x8000000000000000), UWORD(0x0000000000000000), h, l);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, P, 3);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 3))
    {
        flint_mpn_zero(res, 3);
        res[3] = 1;
    }
    else
    {
        mpn_neg(res, T, 3);
        res[3] = 0;
    }
}

/* cos, n = 3, r = 45: 2 levels, widths [2, 1] */
static void
hs_cos_3_45(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0aaaaaaaaaaaaaaa);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x8000000000000000), UWORD(0x0000000000000000), h, l);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, P, 3);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 3))
    {
        flint_mpn_zero(res, 3);
        res[3] = 1;
    }
    else
    {
        mpn_neg(res, T, 3);
        res[3] = 0;
    }
}

/* cos, n = 3, r = 46: 2 levels, widths [2, 1] */
static void
hs_cos_3_46(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    v = UWORD(0x0aaaaaaaaaaaaaaa);

    /* V_0 = c + z V_1   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[2], z[1], v);
    add_ssaaaa(w1, w0, UWORD(0x8000000000000000), UWORD(0x0000000000000000), h, l);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, P, 3);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 3))
    {
        flint_mpn_zero(res, 3);
        res[3] = 1;
    }
    else
    {
        mpn_neg(res, T, 3);
        res[3] = 0;
    }
}

/* cos, n = 3, r = 47: 1 levels, widths [2] */
static void
hs_cos_3_47(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    w1 = UWORD(0x8000000000000000); w0 = UWORD(0x0000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, P, 3);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 3))
    {
        flint_mpn_zero(res, 3);
        res[3] = 1;
    }
    else
    {
        mpn_neg(res, T, 3);
        res[3] = 0;
    }
}

/* cos, n = 3, r = 48: 1 levels, widths [2] */
static void
hs_cos_3_48(nn_ptr res, nn_srcptr x)
{
    ulong z[3], V[3], T[3], P[3];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 3);

    w1 = UWORD(0x8000000000000000); w0 = UWORD(0x0000000000000000);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z, P, 3);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 3))
    {
        flint_mpn_zero(res, 3);
        res[3] = 1;
    }
    else
    {
        mpn_neg(res, T, 3);
        res[3] = 0;
    }
}

/* cos, n = 4, r = 8: 11 levels, widths [4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2] */
static void
hs_cos_4_8(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    w1 = UWORD(0x0000000000000000); w0 = UWORD(0x04338e5b6dfe14a5);

    /* V_9 = c + z V_10   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x0000000000000007), UWORD(0x950ae900808941ea), h, l);

    /* V_8 = c + z V_9   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x0000000000000b41), UWORD(0x3c31dcbecbbdd802), h, l);

    /* V_7 = c + z V_8   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x00000000000d73f9), UWORD(0xf399dc0f88ec32b5), h, l);

    /* V_6 = c + z V_7   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xfd1f2754668c46d4);
    P[1] = UWORD(0x603e4e905d6f8a2e);
    P[2] = UWORD(0x000000000c9cba54);
    mpn_add_n(V, P, T, 3);

    /* V_5 = c + z V_6   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xf425f600e7ba5b3c);
    P[1] = UWORD(0x6c4bdaa26d4c3d67);
    P[2] = UWORD(0x00000008f76c77fc);
    mpn_add_n(V, P, T, 3);

    /* V_4 = c + z V_5   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xe392d8777c170b65);
    P[1] = UWORD(0xd71cbbc05b4fa999);
    P[2] = UWORD(0x0000049f93edde27);
    mpn_add_n(V, P, T, 3);

    /* V_3 = c + z V_4   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x01a01a01a01a01a0);
    P[1] = UWORD(0xa01a01a01a01a01a);
    P[2] = UWORD(0x0001a01a01a01a01);
    mpn_add_n(V, P, T, 3);

    /* V_2 = c + z V_3   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0xb05b05b05b05b05b);
    P[1] = UWORD(0x5b05b05b05b05b05);
    P[2] = UWORD(0x05b05b05b05b05b0);
    P[3] = UWORD(0x005b05b05b05b05b);
    mpn_add_n(V, P, T, 4);

    /* V_1 = c + z V_2   (4 x 4 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[3] = UWORD(0x0aaaaaaaaaaaaaaa);
    mpn_add_n(V, P, T, 4);

    /* V_0 = c + z V_1   (4 x 4 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x0000000000000000);
    P[3] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 4);

    /* cos = 1 - z V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, P, 4);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 4))
    {
        flint_mpn_zero(res, 4);
        res[4] = 1;
    }
    else
    {
        mpn_neg(res, T, 4);
        res[4] = 0;
    }
}

/* cos, n = 4, r = 9: 10 levels, widths [4, 4, 4, 3, 3, 3, 2, 2, 2, 2] */
static void
hs_cos_4_9(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    w1 = UWORD(0x0000000000000007); w0 = UWORD(0x950ae900808941ea);

    /* V_8 = c + z V_9   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x0000000000000b41), UWORD(0x3c31dcbecbbdd802), h, l);

    /* V_7 = c + z V_8   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x00000000000d73f9), UWORD(0xf399dc0f88ec32b5), h, l);

    /* V_6 = c + z V_7   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x000000000c9cba54), UWORD(0x603e4e905d6f8a2e), h, l);

    /* V_5 = c + z V_6   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xf425f600e7ba5b3c);
    P[1] = UWORD(0x6c4bdaa26d4c3d67);
    P[2] = UWORD(0x00000008f76c77fc);
    mpn_add_n(V, P, T, 3);

    /* V_4 = c + z V_5   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xe392d8777c170b65);
    P[1] = UWORD(0xd71cbbc05b4fa999);
    P[2] = UWORD(0x0000049f93edde27);
    mpn_add_n(V, P, T, 3);

    /* V_3 = c + z V_4   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x01a01a01a01a01a0);
    P[1] = UWORD(0xa01a01a01a01a01a);
    P[2] = UWORD(0x0001a01a01a01a01);
    mpn_add_n(V, P, T, 3);

    /* V_2 = c + z V_3   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0xb05b05b05b05b05b);
    P[1] = UWORD(0x5b05b05b05b05b05);
    P[2] = UWORD(0x05b05b05b05b05b0);
    P[3] = UWORD(0x005b05b05b05b05b);
    mpn_add_n(V, P, T, 4);

    /* V_1 = c + z V_2   (4 x 4 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[3] = UWORD(0x0aaaaaaaaaaaaaaa);
    mpn_add_n(V, P, T, 4);

    /* V_0 = c + z V_1   (4 x 4 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x0000000000000000);
    P[3] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 4);

    /* cos = 1 - z V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, P, 4);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 4))
    {
        flint_mpn_zero(res, 4);
        res[4] = 1;
    }
    else
    {
        mpn_neg(res, T, 4);
        res[4] = 0;
    }
}

/* cos, n = 4, r = 10: 9 levels, widths [4, 4, 3, 3, 3, 3, 2, 2, 2] */
static void
hs_cos_4_10(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    w1 = UWORD(0x0000000000000b41); w0 = UWORD(0x3c31dcbecbbdd802);

    /* V_7 = c + z V_8   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x00000000000d73f9), UWORD(0xf399dc0f88ec32b5), h, l);

    /* V_6 = c + z V_7   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x000000000c9cba54), UWORD(0x603e4e905d6f8a2e), h, l);

    /* V_5 = c + z V_6   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xf425f600e7ba5b3c);
    P[1] = UWORD(0x6c4bdaa26d4c3d67);
    P[2] = UWORD(0x00000008f76c77fc);
    mpn_add_n(V, P, T, 3);

    /* V_4 = c + z V_5   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xe392d8777c170b65);
    P[1] = UWORD(0xd71cbbc05b4fa999);
    P[2] = UWORD(0x0000049f93edde27);
    mpn_add_n(V, P, T, 3);

    /* V_3 = c + z V_4   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x01a01a01a01a01a0);
    P[1] = UWORD(0xa01a01a01a01a01a);
    P[2] = UWORD(0x0001a01a01a01a01);
    mpn_add_n(V, P, T, 3);

    /* V_2 = c + z V_3   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5b05b05b05b05b05);
    P[1] = UWORD(0x05b05b05b05b05b0);
    P[2] = UWORD(0x005b05b05b05b05b);
    mpn_add_n(V, P, T, 3);

    /* V_1 = c + z V_2   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[3] = UWORD(0x0aaaaaaaaaaaaaaa);
    mpn_add_n(V, P, T, 4);

    /* V_0 = c + z V_1   (4 x 4 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x0000000000000000);
    P[3] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 4);

    /* cos = 1 - z V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, P, 4);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 4))
    {
        flint_mpn_zero(res, 4);
        res[4] = 1;
    }
    else
    {
        mpn_neg(res, T, 4);
        res[4] = 0;
    }
}

/* cos, n = 4, r = 11: 9 levels, widths [4, 4, 3, 3, 3, 2, 2, 2, 1] */
static void
hs_cos_4_11(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0000000000000b41);

    /* V_7 = c + z V_8   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x00000000000d73f9), UWORD(0xf399dc0f88ec32b5), h, l);

    /* V_6 = c + z V_7   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x000000000c9cba54), UWORD(0x603e4e905d6f8a2e), h, l);

    /* V_5 = c + z V_6   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x00000008f76c77fc), UWORD(0x6c4bdaa26d4c3d67), h, l);

    /* V_4 = c + z V_5   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xe392d8777c170b65);
    P[1] = UWORD(0xd71cbbc05b4fa999);
    P[2] = UWORD(0x0000049f93edde27);
    mpn_add_n(V, P, T, 3);

    /* V_3 = c + z V_4   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x01a01a01a01a01a0);
    P[1] = UWORD(0xa01a01a01a01a01a);
    P[2] = UWORD(0x0001a01a01a01a01);
    mpn_add_n(V, P, T, 3);

    /* V_2 = c + z V_3   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5b05b05b05b05b05);
    P[1] = UWORD(0x05b05b05b05b05b0);
    P[2] = UWORD(0x005b05b05b05b05b);
    mpn_add_n(V, P, T, 3);

    /* V_1 = c + z V_2   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[3] = UWORD(0x0aaaaaaaaaaaaaaa);
    mpn_add_n(V, P, T, 4);

    /* V_0 = c + z V_1   (4 x 4 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x0000000000000000);
    P[3] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 4);

    /* cos = 1 - z V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, P, 4);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 4))
    {
        flint_mpn_zero(res, 4);
        res[4] = 1;
    }
    else
    {
        mpn_neg(res, T, 4);
        res[4] = 0;
    }
}

/* cos, n = 4, r = 12: 8 levels, widths [4, 4, 3, 3, 3, 2, 2, 1] */
static void
hs_cos_4_12(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x00000000000d73f9);

    /* V_6 = c + z V_7   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x000000000c9cba54), UWORD(0x603e4e905d6f8a2e), h, l);

    /* V_5 = c + z V_6   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x00000008f76c77fc), UWORD(0x6c4bdaa26d4c3d67), h, l);

    /* V_4 = c + z V_5   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xe392d8777c170b65);
    P[1] = UWORD(0xd71cbbc05b4fa999);
    P[2] = UWORD(0x0000049f93edde27);
    mpn_add_n(V, P, T, 3);

    /* V_3 = c + z V_4   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x01a01a01a01a01a0);
    P[1] = UWORD(0xa01a01a01a01a01a);
    P[2] = UWORD(0x0001a01a01a01a01);
    mpn_add_n(V, P, T, 3);

    /* V_2 = c + z V_3   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5b05b05b05b05b05);
    P[1] = UWORD(0x05b05b05b05b05b0);
    P[2] = UWORD(0x005b05b05b05b05b);
    mpn_add_n(V, P, T, 3);

    /* V_1 = c + z V_2   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[3] = UWORD(0x0aaaaaaaaaaaaaaa);
    mpn_add_n(V, P, T, 4);

    /* V_0 = c + z V_1   (4 x 4 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x0000000000000000);
    P[3] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 4);

    /* cos = 1 - z V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, P, 4);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 4))
    {
        flint_mpn_zero(res, 4);
        res[4] = 1;
    }
    else
    {
        mpn_neg(res, T, 4);
        res[4] = 0;
    }
}

/* cos, n = 4, r = 13: 8 levels, widths [4, 4, 3, 3, 2, 2, 2, 1] */
static void
hs_cos_4_13(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x00000000000d73f9);

    /* V_6 = c + z V_7   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x000000000c9cba54), UWORD(0x603e4e905d6f8a2e), h, l);

    /* V_5 = c + z V_6   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x00000008f76c77fc), UWORD(0x6c4bdaa26d4c3d67), h, l);

    /* V_4 = c + z V_5   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x0000049f93edde27), UWORD(0xd71cbbc05b4fa999), h, l);

    /* V_3 = c + z V_4   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x01a01a01a01a01a0);
    P[1] = UWORD(0xa01a01a01a01a01a);
    P[2] = UWORD(0x0001a01a01a01a01);
    mpn_add_n(V, P, T, 3);

    /* V_2 = c + z V_3   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5b05b05b05b05b05);
    P[1] = UWORD(0x05b05b05b05b05b0);
    P[2] = UWORD(0x005b05b05b05b05b);
    mpn_add_n(V, P, T, 3);

    /* V_1 = c + z V_2   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[3] = UWORD(0x0aaaaaaaaaaaaaaa);
    mpn_add_n(V, P, T, 4);

    /* V_0 = c + z V_1   (4 x 4 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x0000000000000000);
    P[3] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 4);

    /* cos = 1 - z V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, P, 4);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 4))
    {
        flint_mpn_zero(res, 4);
        res[4] = 1;
    }
    else
    {
        mpn_neg(res, T, 4);
        res[4] = 0;
    }
}

/* cos, n = 4, r = 14: 7 levels, widths [4, 4, 3, 3, 2, 2, 1] */
static void
hs_cos_4_14(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x000000000c9cba54);

    /* V_5 = c + z V_6   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x00000008f76c77fc), UWORD(0x6c4bdaa26d4c3d67), h, l);

    /* V_4 = c + z V_5   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x0000049f93edde27), UWORD(0xd71cbbc05b4fa999), h, l);

    /* V_3 = c + z V_4   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x01a01a01a01a01a0);
    P[1] = UWORD(0xa01a01a01a01a01a);
    P[2] = UWORD(0x0001a01a01a01a01);
    mpn_add_n(V, P, T, 3);

    /* V_2 = c + z V_3   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5b05b05b05b05b05);
    P[1] = UWORD(0x05b05b05b05b05b0);
    P[2] = UWORD(0x005b05b05b05b05b);
    mpn_add_n(V, P, T, 3);

    /* V_1 = c + z V_2   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[3] = UWORD(0x0aaaaaaaaaaaaaaa);
    mpn_add_n(V, P, T, 4);

    /* V_0 = c + z V_1   (4 x 4 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x0000000000000000);
    P[3] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 4);

    /* cos = 1 - z V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, P, 4);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 4))
    {
        flint_mpn_zero(res, 4);
        res[4] = 1;
    }
    else
    {
        mpn_neg(res, T, 4);
        res[4] = 0;
    }
}

/* cos, n = 4, r = 15: 7 levels, widths [4, 3, 3, 3, 2, 2, 1] */
static void
hs_cos_4_15(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x000000000c9cba54);

    /* V_5 = c + z V_6   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x00000008f76c77fc), UWORD(0x6c4bdaa26d4c3d67), h, l);

    /* V_4 = c + z V_5   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x0000049f93edde27), UWORD(0xd71cbbc05b4fa999), h, l);

    /* V_3 = c + z V_4   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x01a01a01a01a01a0);
    P[1] = UWORD(0xa01a01a01a01a01a);
    P[2] = UWORD(0x0001a01a01a01a01);
    mpn_add_n(V, P, T, 3);

    /* V_2 = c + z V_3   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5b05b05b05b05b05);
    P[1] = UWORD(0x05b05b05b05b05b0);
    P[2] = UWORD(0x005b05b05b05b05b);
    mpn_add_n(V, P, T, 3);

    /* V_1 = c + z V_2   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x0aaaaaaaaaaaaaaa);
    mpn_add_n(V, P, T, 3);

    /* V_0 = c + z V_1   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x0000000000000000);
    P[3] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 4);

    /* cos = 1 - z V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, P, 4);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 4))
    {
        flint_mpn_zero(res, 4);
        res[4] = 1;
    }
    else
    {
        mpn_neg(res, T, 4);
        res[4] = 0;
    }
}

/* cos, n = 4, r = 16: 6 levels, widths [4, 3, 3, 2, 2, 1] */
static void
hs_cos_4_16(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x00000008f76c77fc);

    /* V_4 = c + z V_5   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x0000049f93edde27), UWORD(0xd71cbbc05b4fa999), h, l);

    /* V_3 = c + z V_4   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x0001a01a01a01a01), UWORD(0xa01a01a01a01a01a), h, l);

    /* V_2 = c + z V_3   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5b05b05b05b05b05);
    P[1] = UWORD(0x05b05b05b05b05b0);
    P[2] = UWORD(0x005b05b05b05b05b);
    mpn_add_n(V, P, T, 3);

    /* V_1 = c + z V_2   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x0aaaaaaaaaaaaaaa);
    mpn_add_n(V, P, T, 3);

    /* V_0 = c + z V_1   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x0000000000000000);
    P[3] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 4);

    /* cos = 1 - z V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, P, 4);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 4))
    {
        flint_mpn_zero(res, 4);
        res[4] = 1;
    }
    else
    {
        mpn_neg(res, T, 4);
        res[4] = 0;
    }
}

/* cos, n = 4, r = 17: 6 levels, widths [4, 3, 3, 2, 2, 1] */
static void
hs_cos_4_17(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x00000008f76c77fc);

    /* V_4 = c + z V_5   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x0000049f93edde27), UWORD(0xd71cbbc05b4fa999), h, l);

    /* V_3 = c + z V_4   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x0001a01a01a01a01), UWORD(0xa01a01a01a01a01a), h, l);

    /* V_2 = c + z V_3   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5b05b05b05b05b05);
    P[1] = UWORD(0x05b05b05b05b05b0);
    P[2] = UWORD(0x005b05b05b05b05b);
    mpn_add_n(V, P, T, 3);

    /* V_1 = c + z V_2   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x0aaaaaaaaaaaaaaa);
    mpn_add_n(V, P, T, 3);

    /* V_0 = c + z V_1   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x0000000000000000);
    P[3] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 4);

    /* cos = 1 - z V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, P, 4);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 4))
    {
        flint_mpn_zero(res, 4);
        res[4] = 1;
    }
    else
    {
        mpn_neg(res, T, 4);
        res[4] = 0;
    }
}

/* cos, n = 4, r = 18: 6 levels, widths [4, 3, 3, 2, 2, 1] */
static void
hs_cos_4_18(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x00000008f76c77fc);

    /* V_4 = c + z V_5   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x0000049f93edde27), UWORD(0xd71cbbc05b4fa999), h, l);

    /* V_3 = c + z V_4   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x0001a01a01a01a01), UWORD(0xa01a01a01a01a01a), h, l);

    /* V_2 = c + z V_3   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5b05b05b05b05b05);
    P[1] = UWORD(0x05b05b05b05b05b0);
    P[2] = UWORD(0x005b05b05b05b05b);
    mpn_add_n(V, P, T, 3);

    /* V_1 = c + z V_2   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x0aaaaaaaaaaaaaaa);
    mpn_add_n(V, P, T, 3);

    /* V_0 = c + z V_1   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x0000000000000000);
    P[3] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 4);

    /* cos = 1 - z V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, P, 4);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 4))
    {
        flint_mpn_zero(res, 4);
        res[4] = 1;
    }
    else
    {
        mpn_neg(res, T, 4);
        res[4] = 0;
    }
}

/* cos, n = 4, r = 19: 5 levels, widths [4, 3, 3, 2, 1] */
static void
hs_cos_4_19(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0000049f93edde27);

    /* V_3 = c + z V_4   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x0001a01a01a01a01), UWORD(0xa01a01a01a01a01a), h, l);

    /* V_2 = c + z V_3   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5b05b05b05b05b05);
    P[1] = UWORD(0x05b05b05b05b05b0);
    P[2] = UWORD(0x005b05b05b05b05b);
    mpn_add_n(V, P, T, 3);

    /* V_1 = c + z V_2   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x0aaaaaaaaaaaaaaa);
    mpn_add_n(V, P, T, 3);

    /* V_0 = c + z V_1   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x0000000000000000);
    P[3] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 4);

    /* cos = 1 - z V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, P, 4);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 4))
    {
        flint_mpn_zero(res, 4);
        res[4] = 1;
    }
    else
    {
        mpn_neg(res, T, 4);
        res[4] = 0;
    }
}

/* cos, n = 4, r = 20: 5 levels, widths [4, 3, 3, 2, 1] */
static void
hs_cos_4_20(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0000049f93edde27);

    /* V_3 = c + z V_4   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x0001a01a01a01a01), UWORD(0xa01a01a01a01a01a), h, l);

    /* V_2 = c + z V_3   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x5b05b05b05b05b05);
    P[1] = UWORD(0x05b05b05b05b05b0);
    P[2] = UWORD(0x005b05b05b05b05b);
    mpn_add_n(V, P, T, 3);

    /* V_1 = c + z V_2   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x0aaaaaaaaaaaaaaa);
    mpn_add_n(V, P, T, 3);

    /* V_0 = c + z V_1   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x0000000000000000);
    P[3] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 4);

    /* cos = 1 - z V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, P, 4);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 4))
    {
        flint_mpn_zero(res, 4);
        res[4] = 1;
    }
    else
    {
        mpn_neg(res, T, 4);
        res[4] = 0;
    }
}

/* cos, n = 4, r = 21: 5 levels, widths [4, 3, 2, 2, 1] */
static void
hs_cos_4_21(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0000049f93edde27);

    /* V_3 = c + z V_4   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x0001a01a01a01a01), UWORD(0xa01a01a01a01a01a), h, l);

    /* V_2 = c + z V_3   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x005b05b05b05b05b), UWORD(0x05b05b05b05b05b0), h, l);

    /* V_1 = c + z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x0aaaaaaaaaaaaaaa);
    mpn_add_n(V, P, T, 3);

    /* V_0 = c + z V_1   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x0000000000000000);
    P[3] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 4);

    /* cos = 1 - z V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, P, 4);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 4))
    {
        flint_mpn_zero(res, 4);
        res[4] = 1;
    }
    else
    {
        mpn_neg(res, T, 4);
        res[4] = 0;
    }
}

/* cos, n = 4, r = 22: 5 levels, widths [4, 3, 2, 2, 1] */
static void
hs_cos_4_22(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0000049f93edde27);

    /* V_3 = c + z V_4   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x0001a01a01a01a01), UWORD(0xa01a01a01a01a01a), h, l);

    /* V_2 = c + z V_3   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x005b05b05b05b05b), UWORD(0x05b05b05b05b05b0), h, l);

    /* V_1 = c + z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x0aaaaaaaaaaaaaaa);
    mpn_add_n(V, P, T, 3);

    /* V_0 = c + z V_1   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x0000000000000000);
    P[3] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 4);

    /* cos = 1 - z V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, P, 4);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 4))
    {
        flint_mpn_zero(res, 4);
        res[4] = 1;
    }
    else
    {
        mpn_neg(res, T, 4);
        res[4] = 0;
    }
}

/* cos, n = 4, r = 23: 5 levels, widths [4, 3, 2, 2, 1] */
static void
hs_cos_4_23(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0000049f93edde27);

    /* V_3 = c + z V_4   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x0001a01a01a01a01), UWORD(0xa01a01a01a01a01a), h, l);

    /* V_2 = c + z V_3   (2 x 2 limbs) */
    _fixed_mulhi_2x2_sloppy(&h, &l, z[3], z[2], w1, w0);
    add_ssaaaa(w1, w0, UWORD(0x005b05b05b05b05b), UWORD(0x05b05b05b05b05b0), h, l);

    /* V_1 = c + z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x0aaaaaaaaaaaaaaa);
    mpn_add_n(V, P, T, 3);

    /* V_0 = c + z V_1   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x0000000000000000);
    P[3] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 4);

    /* cos = 1 - z V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, P, 4);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 4))
    {
        flint_mpn_zero(res, 4);
        res[4] = 1;
    }
    else
    {
        mpn_neg(res, T, 4);
        res[4] = 0;
    }
}

/* cos, n = 4, r = 24: 4 levels, widths [4, 3, 2, 1] */
static void
hs_cos_4_24(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0001a01a01a01a01);

    /* V_2 = c + z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x005b05b05b05b05b), UWORD(0x05b05b05b05b05b0), h, l);

    /* V_1 = c + z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x0aaaaaaaaaaaaaaa);
    mpn_add_n(V, P, T, 3);

    /* V_0 = c + z V_1   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x0000000000000000);
    P[3] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 4);

    /* cos = 1 - z V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, P, 4);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 4))
    {
        flint_mpn_zero(res, 4);
        res[4] = 1;
    }
    else
    {
        mpn_neg(res, T, 4);
        res[4] = 0;
    }
}

/* cos, n = 4, r = 25: 4 levels, widths [4, 3, 2, 1] */
static void
hs_cos_4_25(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0001a01a01a01a01);

    /* V_2 = c + z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x005b05b05b05b05b), UWORD(0x05b05b05b05b05b0), h, l);

    /* V_1 = c + z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x0aaaaaaaaaaaaaaa);
    mpn_add_n(V, P, T, 3);

    /* V_0 = c + z V_1   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x0000000000000000);
    P[3] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 4);

    /* cos = 1 - z V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, P, 4);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 4))
    {
        flint_mpn_zero(res, 4);
        res[4] = 1;
    }
    else
    {
        mpn_neg(res, T, 4);
        res[4] = 0;
    }
}

/* cos, n = 4, r = 26: 4 levels, widths [4, 3, 2, 1] */
static void
hs_cos_4_26(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0001a01a01a01a01);

    /* V_2 = c + z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x005b05b05b05b05b), UWORD(0x05b05b05b05b05b0), h, l);

    /* V_1 = c + z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x0aaaaaaaaaaaaaaa);
    mpn_add_n(V, P, T, 3);

    /* V_0 = c + z V_1   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x0000000000000000);
    P[3] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 4);

    /* cos = 1 - z V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, P, 4);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 4))
    {
        flint_mpn_zero(res, 4);
        res[4] = 1;
    }
    else
    {
        mpn_neg(res, T, 4);
        res[4] = 0;
    }
}

/* cos, n = 4, r = 27: 4 levels, widths [4, 3, 2, 1] */
static void
hs_cos_4_27(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0001a01a01a01a01);

    /* V_2 = c + z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x005b05b05b05b05b), UWORD(0x05b05b05b05b05b0), h, l);

    /* V_1 = c + z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x0aaaaaaaaaaaaaaa);
    mpn_add_n(V, P, T, 3);

    /* V_0 = c + z V_1   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x0000000000000000);
    P[3] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 4);

    /* cos = 1 - z V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, P, 4);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 4))
    {
        flint_mpn_zero(res, 4);
        res[4] = 1;
    }
    else
    {
        mpn_neg(res, T, 4);
        res[4] = 0;
    }
}

/* cos, n = 4, r = 28: 4 levels, widths [4, 3, 2, 1] */
static void
hs_cos_4_28(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0001a01a01a01a01);

    /* V_2 = c + z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x005b05b05b05b05b), UWORD(0x05b05b05b05b05b0), h, l);

    /* V_1 = c + z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x0aaaaaaaaaaaaaaa);
    mpn_add_n(V, P, T, 3);

    /* V_0 = c + z V_1   (4 x 3 limbs) */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z + 0, P, 4);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x0000000000000000);
    P[3] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 4);

    /* cos = 1 - z V_0 */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    P[3] = V[3];
    flint_mpn_mulhigh_n(T, z, P, 4);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 4))
    {
        flint_mpn_zero(res, 4);
        res[4] = 1;
    }
    else
    {
        mpn_neg(res, T, 4);
        res[4] = 0;
    }
}

/* cos, n = 4, r = 29: 4 levels, widths [3, 3, 2, 1] */
static void
hs_cos_4_29(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0001a01a01a01a01);

    /* V_2 = c + z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x005b05b05b05b05b), UWORD(0x05b05b05b05b05b0), h, l);

    /* V_1 = c + z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x0aaaaaaaaaaaaaaa);
    mpn_add_n(V, P, T, 3);

    /* V_0 = c + z V_1   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 3);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, P, 4);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 4))
    {
        flint_mpn_zero(res, 4);
        res[4] = 1;
    }
    else
    {
        mpn_neg(res, T, 4);
        res[4] = 0;
    }
}

/* cos, n = 4, r = 30: 4 levels, widths [3, 3, 2, 1] */
static void
hs_cos_4_30(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0001a01a01a01a01);

    /* V_2 = c + z V_3   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x005b05b05b05b05b), UWORD(0x05b05b05b05b05b0), h, l);

    /* V_1 = c + z V_2   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[1] = UWORD(0xaaaaaaaaaaaaaaaa);
    P[2] = UWORD(0x0aaaaaaaaaaaaaaa);
    mpn_add_n(V, P, T, 3);

    /* V_0 = c + z V_1   (3 x 3 limbs) */
    P[0] = V[0];
    P[1] = V[1];
    P[2] = V[2];
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 3);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, P, 4);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 4))
    {
        flint_mpn_zero(res, 4);
        res[4] = 1;
    }
    else
    {
        mpn_neg(res, T, 4);
        res[4] = 0;
    }
}

/* cos, n = 4, r = 31: 3 levels, widths [3, 2, 1] */
static void
hs_cos_4_31(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x005b05b05b05b05b);

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x0aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 3);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, P, 4);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 4))
    {
        flint_mpn_zero(res, 4);
        res[4] = 1;
    }
    else
    {
        mpn_neg(res, T, 4);
        res[4] = 0;
    }
}

/* cos, n = 4, r = 32: 3 levels, widths [3, 2, 1] */
static void
hs_cos_4_32(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x005b05b05b05b05b);

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x0aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 3);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, P, 4);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 4))
    {
        flint_mpn_zero(res, 4);
        res[4] = 1;
    }
    else
    {
        mpn_neg(res, T, 4);
        res[4] = 0;
    }
}

/* cos, n = 4, r = 33: 3 levels, widths [3, 2, 1] */
static void
hs_cos_4_33(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x005b05b05b05b05b);

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x0aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 3);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, P, 4);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 4))
    {
        flint_mpn_zero(res, 4);
        res[4] = 1;
    }
    else
    {
        mpn_neg(res, T, 4);
        res[4] = 0;
    }
}

/* cos, n = 4, r = 34: 3 levels, widths [3, 2, 1] */
static void
hs_cos_4_34(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x005b05b05b05b05b);

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x0aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 3);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, P, 4);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 4))
    {
        flint_mpn_zero(res, 4);
        res[4] = 1;
    }
    else
    {
        mpn_neg(res, T, 4);
        res[4] = 0;
    }
}

/* cos, n = 4, r = 35: 3 levels, widths [3, 2, 1] */
static void
hs_cos_4_35(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x005b05b05b05b05b);

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x0aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 3);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, P, 4);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 4))
    {
        flint_mpn_zero(res, 4);
        res[4] = 1;
    }
    else
    {
        mpn_neg(res, T, 4);
        res[4] = 0;
    }
}

/* cos, n = 4, r = 36: 3 levels, widths [3, 2, 1] */
static void
hs_cos_4_36(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x005b05b05b05b05b);

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x0aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 3);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, P, 4);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 4))
    {
        flint_mpn_zero(res, 4);
        res[4] = 1;
    }
    else
    {
        mpn_neg(res, T, 4);
        res[4] = 0;
    }
}

/* cos, n = 4, r = 37: 3 levels, widths [3, 2, 1] */
static void
hs_cos_4_37(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x005b05b05b05b05b);

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x0aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 3);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, P, 4);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 4))
    {
        flint_mpn_zero(res, 4);
        res[4] = 1;
    }
    else
    {
        mpn_neg(res, T, 4);
        res[4] = 0;
    }
}

/* cos, n = 4, r = 38: 3 levels, widths [3, 2, 1] */
static void
hs_cos_4_38(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x005b05b05b05b05b);

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x0aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 3);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, P, 4);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 4))
    {
        flint_mpn_zero(res, 4);
        res[4] = 1;
    }
    else
    {
        mpn_neg(res, T, 4);
        res[4] = 0;
    }
}

/* cos, n = 4, r = 39: 3 levels, widths [3, 2, 1] */
static void
hs_cos_4_39(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x005b05b05b05b05b);

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x0aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 3);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, P, 4);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 4))
    {
        flint_mpn_zero(res, 4);
        res[4] = 1;
    }
    else
    {
        mpn_neg(res, T, 4);
        res[4] = 0;
    }
}

/* cos, n = 4, r = 40: 3 levels, widths [3, 2, 1] */
static void
hs_cos_4_40(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x005b05b05b05b05b);

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x0aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 3);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, P, 4);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 4))
    {
        flint_mpn_zero(res, 4);
        res[4] = 1;
    }
    else
    {
        mpn_neg(res, T, 4);
        res[4] = 0;
    }
}

/* cos, n = 4, r = 41: 3 levels, widths [3, 2, 1] */
static void
hs_cos_4_41(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x005b05b05b05b05b);

    /* V_1 = c + z V_2   (2 x 1 limbs) */
    _fixed_mulhi_2x1(&h, &l, z[3], z[2], v);
    add_ssaaaa(w1, w0, UWORD(0x0aaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), h, l);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 3);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, P, 4);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 4))
    {
        flint_mpn_zero(res, 4);
        res[4] = 1;
    }
    else
    {
        mpn_neg(res, T, 4);
        res[4] = 0;
    }
}

/* cos, n = 4, r = 42: 2 levels, widths [3, 2] */
static void
hs_cos_4_42(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    w1 = UWORD(0x0aaaaaaaaaaaaaaa); w0 = UWORD(0xaaaaaaaaaaaaaaaa);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 3);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, P, 4);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 4))
    {
        flint_mpn_zero(res, 4);
        res[4] = 1;
    }
    else
    {
        mpn_neg(res, T, 4);
        res[4] = 0;
    }
}

/* cos, n = 4, r = 43: 2 levels, widths [3, 2] */
static void
hs_cos_4_43(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    w1 = UWORD(0x0aaaaaaaaaaaaaaa); w0 = UWORD(0xaaaaaaaaaaaaaaaa);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 3);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, P, 4);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 4))
    {
        flint_mpn_zero(res, 4);
        res[4] = 1;
    }
    else
    {
        mpn_neg(res, T, 4);
        res[4] = 0;
    }
}

/* cos, n = 4, r = 44: 2 levels, widths [3, 2] */
static void
hs_cos_4_44(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    w1 = UWORD(0x0aaaaaaaaaaaaaaa); w0 = UWORD(0xaaaaaaaaaaaaaaaa);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 3);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, P, 4);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 4))
    {
        flint_mpn_zero(res, 4);
        res[4] = 1;
    }
    else
    {
        mpn_neg(res, T, 4);
        res[4] = 0;
    }
}

/* cos, n = 4, r = 45: 2 levels, widths [3, 2] */
static void
hs_cos_4_45(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    w1 = UWORD(0x0aaaaaaaaaaaaaaa); w0 = UWORD(0xaaaaaaaaaaaaaaaa);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 3);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, P, 4);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 4))
    {
        flint_mpn_zero(res, 4);
        res[4] = 1;
    }
    else
    {
        mpn_neg(res, T, 4);
        res[4] = 0;
    }
}

/* cos, n = 4, r = 46: 2 levels, widths [3, 2] */
static void
hs_cos_4_46(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    w1 = UWORD(0x0aaaaaaaaaaaaaaa); w0 = UWORD(0xaaaaaaaaaaaaaaaa);

    /* V_0 = c + z V_1   (3 x 2 limbs) */
    P[0] = UWORD(0);
    P[1] = w0; P[2] = w1;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 3);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, P, 4);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 4))
    {
        flint_mpn_zero(res, 4);
        res[4] = 1;
    }
    else
    {
        mpn_neg(res, T, 4);
        res[4] = 0;
    }
}

/* cos, n = 4, r = 47: 2 levels, widths [3, 1] */
static void
hs_cos_4_47(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0aaaaaaaaaaaaaaa);

    /* V_0 = c + z V_1   (3 x 1 limbs) */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = v;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 3);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, P, 4);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 4))
    {
        flint_mpn_zero(res, 4);
        res[4] = 1;
    }
    else
    {
        mpn_neg(res, T, 4);
        res[4] = 0;
    }
}

/* cos, n = 4, r = 48: 2 levels, widths [3, 1] */
static void
hs_cos_4_48(nn_ptr res, nn_srcptr x)
{
    ulong z[4], V[4], T[4], P[4];
    ulong v, w1, w0, h, l;
    slong i;

    (void) v; (void) w1; (void) w0; (void) h; (void) l; (void) i;

    flint_mpn_sqrhigh(z, x, 4);

    v = UWORD(0x0aaaaaaaaaaaaaaa);

    /* V_0 = c + z V_1   (3 x 1 limbs) */
    P[0] = UWORD(0);
    P[1] = UWORD(0);
    P[2] = v;
    flint_mpn_mulhigh_n(T, z + 1, P, 3);
    P[0] = UWORD(0x0000000000000000);
    P[1] = UWORD(0x0000000000000000);
    P[2] = UWORD(0x8000000000000000);
    mpn_add_n(V, P, T, 3);

    /* cos = 1 - z V_0 */
    P[0] = UWORD(0);
    P[1] = V[0];
    P[2] = V[1];
    P[3] = V[2];
    flint_mpn_mulhigh_n(T, z, P, 4);
    /* 1 - T: below one the unit limb is zero */
    if (mpn_zero_p(T, 4))
    {
        flint_mpn_zero(res, 4);
        res[4] = 1;
    }
    else
    {
        mpn_neg(res, T, 4);
        res[4] = 0;
    }
}

typedef void (*hs_fn)(nn_ptr, nn_srcptr);
typedef struct { hs_fn f; const char * name; int n; int r; int N; } hs_entry;
static hs_entry hs_all[] = {
    { hs_atanh_1_8, "atanh", 1, 8, 3 },
    { hs_atanh_1_9, "atanh", 1, 9, 2 },
    { hs_atanh_1_10, "atanh", 1, 10, 2 },
    { hs_atanh_1_11, "atanh", 1, 11, 2 },
    { hs_atanh_1_12, "atanh", 1, 12, 2 },
    { hs_atanh_1_13, "atanh", 1, 13, 1 },
    { hs_atanh_1_14, "atanh", 1, 14, 1 },
    { hs_atanh_1_15, "atanh", 1, 15, 1 },
    { hs_atanh_1_16, "atanh", 1, 16, 1 },
    { hs_atanh_1_17, "atanh", 1, 17, 1 },
    { hs_atanh_1_18, "atanh", 1, 18, 1 },
    { hs_atanh_1_19, "atanh", 1, 19, 1 },
    { hs_atanh_1_20, "atanh", 1, 20, 1 },
    { hs_atanh_1_21, "atanh", 1, 21, 1 },
    { hs_atanh_1_22, "atanh", 1, 22, 1 },
    { hs_atanh_1_23, "atanh", 1, 23, 1 },
    { hs_atanh_1_24, "atanh", 1, 24, 1 },
    { hs_atanh_1_25, "atanh", 1, 25, 1 },
    { hs_atanh_1_26, "atanh", 1, 26, 1 },
    { hs_atanh_1_27, "atanh", 1, 27, 1 },
    { hs_atanh_1_28, "atanh", 1, 28, 1 },
    { hs_atanh_1_29, "atanh", 1, 29, 1 },
    { hs_atanh_1_30, "atanh", 1, 30, 1 },
    { hs_atanh_1_31, "atanh", 1, 31, 1 },
    { hs_atanh_1_32, "atanh", 1, 32, 1 },
    { hs_atanh_1_33, "atanh", 1, 33, 1 },
    { hs_atanh_1_34, "atanh", 1, 34, 1 },
    { hs_atanh_1_35, "atanh", 1, 35, 1 },
    { hs_atanh_1_36, "atanh", 1, 36, 1 },
    { hs_atanh_1_37, "atanh", 1, 37, 1 },
    { hs_atanh_1_38, "atanh", 1, 38, 1 },
    { hs_atanh_1_39, "atanh", 1, 39, 1 },
    { hs_atanh_1_40, "atanh", 1, 40, 1 },
    { hs_atanh_1_41, "atanh", 1, 41, 1 },
    { hs_atanh_1_42, "atanh", 1, 42, 1 },
    { hs_atanh_1_43, "atanh", 1, 43, 1 },
    { hs_atanh_1_44, "atanh", 1, 44, 1 },
    { hs_atanh_1_45, "atanh", 1, 45, 1 },
    { hs_atanh_1_46, "atanh", 1, 46, 1 },
    { hs_atanh_1_47, "atanh", 1, 47, 1 },
    { hs_atanh_1_48, "atanh", 1, 48, 1 },
    { hs_atanh_2_8, "atanh", 2, 8, 7 },
    { hs_atanh_2_9, "atanh", 2, 9, 6 },
    { hs_atanh_2_10, "atanh", 2, 10, 5 },
    { hs_atanh_2_11, "atanh", 2, 11, 5 },
    { hs_atanh_2_12, "atanh", 2, 12, 4 },
    { hs_atanh_2_13, "atanh", 2, 13, 4 },
    { hs_atanh_2_14, "atanh", 2, 14, 3 },
    { hs_atanh_2_15, "atanh", 2, 15, 3 },
    { hs_atanh_2_16, "atanh", 2, 16, 3 },
    { hs_atanh_2_17, "atanh", 2, 17, 3 },
    { hs_atanh_2_18, "atanh", 2, 18, 2 },
    { hs_atanh_2_19, "atanh", 2, 19, 2 },
    { hs_atanh_2_20, "atanh", 2, 20, 2 },
    { hs_atanh_2_21, "atanh", 2, 21, 2 },
    { hs_atanh_2_22, "atanh", 2, 22, 2 },
    { hs_atanh_2_23, "atanh", 2, 23, 2 },
    { hs_atanh_2_24, "atanh", 2, 24, 2 },
    { hs_atanh_2_25, "atanh", 2, 25, 2 },
    { hs_atanh_2_26, "atanh", 2, 26, 1 },
    { hs_atanh_2_27, "atanh", 2, 27, 1 },
    { hs_atanh_2_28, "atanh", 2, 28, 1 },
    { hs_atanh_2_29, "atanh", 2, 29, 1 },
    { hs_atanh_2_30, "atanh", 2, 30, 1 },
    { hs_atanh_2_31, "atanh", 2, 31, 1 },
    { hs_atanh_2_32, "atanh", 2, 32, 1 },
    { hs_atanh_2_33, "atanh", 2, 33, 1 },
    { hs_atanh_2_34, "atanh", 2, 34, 1 },
    { hs_atanh_2_35, "atanh", 2, 35, 1 },
    { hs_atanh_2_36, "atanh", 2, 36, 1 },
    { hs_atanh_2_37, "atanh", 2, 37, 1 },
    { hs_atanh_2_38, "atanh", 2, 38, 1 },
    { hs_atanh_2_39, "atanh", 2, 39, 1 },
    { hs_atanh_2_40, "atanh", 2, 40, 1 },
    { hs_atanh_2_41, "atanh", 2, 41, 1 },
    { hs_atanh_2_42, "atanh", 2, 42, 1 },
    { hs_atanh_2_43, "atanh", 2, 43, 1 },
    { hs_atanh_2_44, "atanh", 2, 44, 1 },
    { hs_atanh_2_45, "atanh", 2, 45, 1 },
    { hs_atanh_2_46, "atanh", 2, 46, 1 },
    { hs_atanh_2_47, "atanh", 2, 47, 1 },
    { hs_atanh_2_48, "atanh", 2, 48, 1 },
    { hs_atanh_3_8, "atanh", 3, 8, 11 },
    { hs_atanh_3_9, "atanh", 3, 9, 9 },
    { hs_atanh_3_10, "atanh", 3, 10, 8 },
    { hs_atanh_3_11, "atanh", 3, 11, 8 },
    { hs_atanh_3_12, "atanh", 3, 12, 7 },
    { hs_atanh_3_13, "atanh", 3, 13, 6 },
    { hs_atanh_3_14, "atanh", 3, 14, 6 },
    { hs_atanh_3_15, "atanh", 3, 15, 5 },
    { hs_atanh_3_16, "atanh", 3, 16, 5 },
    { hs_atanh_3_17, "atanh", 3, 17, 5 },
    { hs_atanh_3_18, "atanh", 3, 18, 4 },
    { hs_atanh_3_19, "atanh", 3, 19, 4 },
    { hs_atanh_3_20, "atanh", 3, 20, 4 },
    { hs_atanh_3_21, "atanh", 3, 21, 3 },
    { hs_atanh_3_22, "atanh", 3, 22, 3 },
    { hs_atanh_3_23, "atanh", 3, 23, 3 },
    { hs_atanh_3_24, "atanh", 3, 24, 3 },
    { hs_atanh_3_25, "atanh", 3, 25, 3 },
    { hs_atanh_3_26, "atanh", 3, 26, 3 },
    { hs_atanh_3_27, "atanh", 3, 27, 3 },
    { hs_atanh_3_28, "atanh", 3, 28, 2 },
    { hs_atanh_3_29, "atanh", 3, 29, 2 },
    { hs_atanh_3_30, "atanh", 3, 30, 2 },
    { hs_atanh_3_31, "atanh", 3, 31, 2 },
    { hs_atanh_3_32, "atanh", 3, 32, 2 },
    { hs_atanh_3_33, "atanh", 3, 33, 2 },
    { hs_atanh_3_34, "atanh", 3, 34, 2 },
    { hs_atanh_3_35, "atanh", 3, 35, 2 },
    { hs_atanh_3_36, "atanh", 3, 36, 2 },
    { hs_atanh_3_37, "atanh", 3, 37, 2 },
    { hs_atanh_3_38, "atanh", 3, 38, 1 },
    { hs_atanh_3_39, "atanh", 3, 39, 1 },
    { hs_atanh_3_40, "atanh", 3, 40, 1 },
    { hs_atanh_3_41, "atanh", 3, 41, 1 },
    { hs_atanh_3_42, "atanh", 3, 42, 1 },
    { hs_atanh_3_43, "atanh", 3, 43, 1 },
    { hs_atanh_3_44, "atanh", 3, 44, 1 },
    { hs_atanh_3_45, "atanh", 3, 45, 1 },
    { hs_atanh_3_46, "atanh", 3, 46, 1 },
    { hs_atanh_3_47, "atanh", 3, 47, 1 },
    { hs_atanh_3_48, "atanh", 3, 48, 1 },
    { hs_atanh_4_8, "atanh", 4, 8, 15 },
    { hs_atanh_4_9, "atanh", 4, 9, 13 },
    { hs_atanh_4_10, "atanh", 4, 10, 12 },
    { hs_atanh_4_11, "atanh", 4, 11, 10 },
    { hs_atanh_4_12, "atanh", 4, 12, 9 },
    { hs_atanh_4_13, "atanh", 4, 13, 9 },
    { hs_atanh_4_14, "atanh", 4, 14, 8 },
    { hs_atanh_4_15, "atanh", 4, 15, 7 },
    { hs_atanh_4_16, "atanh", 4, 16, 7 },
    { hs_atanh_4_17, "atanh", 4, 17, 6 },
    { hs_atanh_4_18, "atanh", 4, 18, 6 },
    { hs_atanh_4_19, "atanh", 4, 19, 6 },
    { hs_atanh_4_20, "atanh", 4, 20, 5 },
    { hs_atanh_4_21, "atanh", 4, 21, 5 },
    { hs_atanh_4_22, "atanh", 4, 22, 5 },
    { hs_atanh_4_23, "atanh", 4, 23, 4 },
    { hs_atanh_4_24, "atanh", 4, 24, 4 },
    { hs_atanh_4_25, "atanh", 4, 25, 4 },
    { hs_atanh_4_26, "atanh", 4, 26, 4 },
    { hs_atanh_4_27, "atanh", 4, 27, 4 },
    { hs_atanh_4_28, "atanh", 4, 28, 4 },
    { hs_atanh_4_29, "atanh", 4, 29, 3 },
    { hs_atanh_4_30, "atanh", 4, 30, 3 },
    { hs_atanh_4_31, "atanh", 4, 31, 3 },
    { hs_atanh_4_32, "atanh", 4, 32, 3 },
    { hs_atanh_4_33, "atanh", 4, 33, 3 },
    { hs_atanh_4_34, "atanh", 4, 34, 3 },
    { hs_atanh_4_35, "atanh", 4, 35, 3 },
    { hs_atanh_4_36, "atanh", 4, 36, 3 },
    { hs_atanh_4_37, "atanh", 4, 37, 2 },
    { hs_atanh_4_38, "atanh", 4, 38, 2 },
    { hs_atanh_4_39, "atanh", 4, 39, 2 },
    { hs_atanh_4_40, "atanh", 4, 40, 2 },
    { hs_atanh_4_41, "atanh", 4, 41, 2 },
    { hs_atanh_4_42, "atanh", 4, 42, 2 },
    { hs_atanh_4_43, "atanh", 4, 43, 2 },
    { hs_atanh_4_44, "atanh", 4, 44, 2 },
    { hs_atanh_4_45, "atanh", 4, 45, 2 },
    { hs_atanh_4_46, "atanh", 4, 46, 2 },
    { hs_atanh_4_47, "atanh", 4, 47, 2 },
    { hs_atanh_4_48, "atanh", 4, 48, 2 },
    { hs_atan_1_8, "atan", 1, 8, 3 },
    { hs_atan_1_9, "atan", 1, 9, 2 },
    { hs_atan_1_10, "atan", 1, 10, 2 },
    { hs_atan_1_11, "atan", 1, 11, 2 },
    { hs_atan_1_12, "atan", 1, 12, 2 },
    { hs_atan_1_13, "atan", 1, 13, 1 },
    { hs_atan_1_14, "atan", 1, 14, 1 },
    { hs_atan_1_15, "atan", 1, 15, 1 },
    { hs_atan_1_16, "atan", 1, 16, 1 },
    { hs_atan_1_17, "atan", 1, 17, 1 },
    { hs_atan_1_18, "atan", 1, 18, 1 },
    { hs_atan_1_19, "atan", 1, 19, 1 },
    { hs_atan_1_20, "atan", 1, 20, 1 },
    { hs_atan_1_21, "atan", 1, 21, 1 },
    { hs_atan_1_22, "atan", 1, 22, 1 },
    { hs_atan_1_23, "atan", 1, 23, 1 },
    { hs_atan_1_24, "atan", 1, 24, 1 },
    { hs_atan_1_25, "atan", 1, 25, 1 },
    { hs_atan_1_26, "atan", 1, 26, 1 },
    { hs_atan_1_27, "atan", 1, 27, 1 },
    { hs_atan_1_28, "atan", 1, 28, 1 },
    { hs_atan_1_29, "atan", 1, 29, 1 },
    { hs_atan_1_30, "atan", 1, 30, 1 },
    { hs_atan_1_31, "atan", 1, 31, 1 },
    { hs_atan_1_32, "atan", 1, 32, 1 },
    { hs_atan_1_33, "atan", 1, 33, 1 },
    { hs_atan_1_34, "atan", 1, 34, 1 },
    { hs_atan_1_35, "atan", 1, 35, 1 },
    { hs_atan_1_36, "atan", 1, 36, 1 },
    { hs_atan_1_37, "atan", 1, 37, 1 },
    { hs_atan_1_38, "atan", 1, 38, 1 },
    { hs_atan_1_39, "atan", 1, 39, 1 },
    { hs_atan_1_40, "atan", 1, 40, 1 },
    { hs_atan_1_41, "atan", 1, 41, 1 },
    { hs_atan_1_42, "atan", 1, 42, 1 },
    { hs_atan_1_43, "atan", 1, 43, 1 },
    { hs_atan_1_44, "atan", 1, 44, 1 },
    { hs_atan_1_45, "atan", 1, 45, 1 },
    { hs_atan_1_46, "atan", 1, 46, 1 },
    { hs_atan_1_47, "atan", 1, 47, 1 },
    { hs_atan_1_48, "atan", 1, 48, 1 },
    { hs_atan_2_8, "atan", 2, 8, 7 },
    { hs_atan_2_9, "atan", 2, 9, 6 },
    { hs_atan_2_10, "atan", 2, 10, 5 },
    { hs_atan_2_11, "atan", 2, 11, 5 },
    { hs_atan_2_12, "atan", 2, 12, 4 },
    { hs_atan_2_13, "atan", 2, 13, 4 },
    { hs_atan_2_14, "atan", 2, 14, 3 },
    { hs_atan_2_15, "atan", 2, 15, 3 },
    { hs_atan_2_16, "atan", 2, 16, 3 },
    { hs_atan_2_17, "atan", 2, 17, 3 },
    { hs_atan_2_18, "atan", 2, 18, 2 },
    { hs_atan_2_19, "atan", 2, 19, 2 },
    { hs_atan_2_20, "atan", 2, 20, 2 },
    { hs_atan_2_21, "atan", 2, 21, 2 },
    { hs_atan_2_22, "atan", 2, 22, 2 },
    { hs_atan_2_23, "atan", 2, 23, 2 },
    { hs_atan_2_24, "atan", 2, 24, 2 },
    { hs_atan_2_25, "atan", 2, 25, 2 },
    { hs_atan_2_26, "atan", 2, 26, 1 },
    { hs_atan_2_27, "atan", 2, 27, 1 },
    { hs_atan_2_28, "atan", 2, 28, 1 },
    { hs_atan_2_29, "atan", 2, 29, 1 },
    { hs_atan_2_30, "atan", 2, 30, 1 },
    { hs_atan_2_31, "atan", 2, 31, 1 },
    { hs_atan_2_32, "atan", 2, 32, 1 },
    { hs_atan_2_33, "atan", 2, 33, 1 },
    { hs_atan_2_34, "atan", 2, 34, 1 },
    { hs_atan_2_35, "atan", 2, 35, 1 },
    { hs_atan_2_36, "atan", 2, 36, 1 },
    { hs_atan_2_37, "atan", 2, 37, 1 },
    { hs_atan_2_38, "atan", 2, 38, 1 },
    { hs_atan_2_39, "atan", 2, 39, 1 },
    { hs_atan_2_40, "atan", 2, 40, 1 },
    { hs_atan_2_41, "atan", 2, 41, 1 },
    { hs_atan_2_42, "atan", 2, 42, 1 },
    { hs_atan_2_43, "atan", 2, 43, 1 },
    { hs_atan_2_44, "atan", 2, 44, 1 },
    { hs_atan_2_45, "atan", 2, 45, 1 },
    { hs_atan_2_46, "atan", 2, 46, 1 },
    { hs_atan_2_47, "atan", 2, 47, 1 },
    { hs_atan_2_48, "atan", 2, 48, 1 },
    { hs_atan_3_8, "atan", 3, 8, 11 },
    { hs_atan_3_9, "atan", 3, 9, 9 },
    { hs_atan_3_10, "atan", 3, 10, 8 },
    { hs_atan_3_11, "atan", 3, 11, 8 },
    { hs_atan_3_12, "atan", 3, 12, 7 },
    { hs_atan_3_13, "atan", 3, 13, 6 },
    { hs_atan_3_14, "atan", 3, 14, 6 },
    { hs_atan_3_15, "atan", 3, 15, 5 },
    { hs_atan_3_16, "atan", 3, 16, 5 },
    { hs_atan_3_17, "atan", 3, 17, 5 },
    { hs_atan_3_18, "atan", 3, 18, 4 },
    { hs_atan_3_19, "atan", 3, 19, 4 },
    { hs_atan_3_20, "atan", 3, 20, 4 },
    { hs_atan_3_21, "atan", 3, 21, 3 },
    { hs_atan_3_22, "atan", 3, 22, 3 },
    { hs_atan_3_23, "atan", 3, 23, 3 },
    { hs_atan_3_24, "atan", 3, 24, 3 },
    { hs_atan_3_25, "atan", 3, 25, 3 },
    { hs_atan_3_26, "atan", 3, 26, 3 },
    { hs_atan_3_27, "atan", 3, 27, 3 },
    { hs_atan_3_28, "atan", 3, 28, 2 },
    { hs_atan_3_29, "atan", 3, 29, 2 },
    { hs_atan_3_30, "atan", 3, 30, 2 },
    { hs_atan_3_31, "atan", 3, 31, 2 },
    { hs_atan_3_32, "atan", 3, 32, 2 },
    { hs_atan_3_33, "atan", 3, 33, 2 },
    { hs_atan_3_34, "atan", 3, 34, 2 },
    { hs_atan_3_35, "atan", 3, 35, 2 },
    { hs_atan_3_36, "atan", 3, 36, 2 },
    { hs_atan_3_37, "atan", 3, 37, 2 },
    { hs_atan_3_38, "atan", 3, 38, 1 },
    { hs_atan_3_39, "atan", 3, 39, 1 },
    { hs_atan_3_40, "atan", 3, 40, 1 },
    { hs_atan_3_41, "atan", 3, 41, 1 },
    { hs_atan_3_42, "atan", 3, 42, 1 },
    { hs_atan_3_43, "atan", 3, 43, 1 },
    { hs_atan_3_44, "atan", 3, 44, 1 },
    { hs_atan_3_45, "atan", 3, 45, 1 },
    { hs_atan_3_46, "atan", 3, 46, 1 },
    { hs_atan_3_47, "atan", 3, 47, 1 },
    { hs_atan_3_48, "atan", 3, 48, 1 },
    { hs_atan_4_8, "atan", 4, 8, 15 },
    { hs_atan_4_9, "atan", 4, 9, 13 },
    { hs_atan_4_10, "atan", 4, 10, 12 },
    { hs_atan_4_11, "atan", 4, 11, 10 },
    { hs_atan_4_12, "atan", 4, 12, 9 },
    { hs_atan_4_13, "atan", 4, 13, 9 },
    { hs_atan_4_14, "atan", 4, 14, 8 },
    { hs_atan_4_15, "atan", 4, 15, 7 },
    { hs_atan_4_16, "atan", 4, 16, 7 },
    { hs_atan_4_17, "atan", 4, 17, 6 },
    { hs_atan_4_18, "atan", 4, 18, 6 },
    { hs_atan_4_19, "atan", 4, 19, 6 },
    { hs_atan_4_20, "atan", 4, 20, 5 },
    { hs_atan_4_21, "atan", 4, 21, 5 },
    { hs_atan_4_22, "atan", 4, 22, 5 },
    { hs_atan_4_23, "atan", 4, 23, 4 },
    { hs_atan_4_24, "atan", 4, 24, 4 },
    { hs_atan_4_25, "atan", 4, 25, 4 },
    { hs_atan_4_26, "atan", 4, 26, 4 },
    { hs_atan_4_27, "atan", 4, 27, 4 },
    { hs_atan_4_28, "atan", 4, 28, 4 },
    { hs_atan_4_29, "atan", 4, 29, 3 },
    { hs_atan_4_30, "atan", 4, 30, 3 },
    { hs_atan_4_31, "atan", 4, 31, 3 },
    { hs_atan_4_32, "atan", 4, 32, 3 },
    { hs_atan_4_33, "atan", 4, 33, 3 },
    { hs_atan_4_34, "atan", 4, 34, 3 },
    { hs_atan_4_35, "atan", 4, 35, 3 },
    { hs_atan_4_36, "atan", 4, 36, 3 },
    { hs_atan_4_37, "atan", 4, 37, 2 },
    { hs_atan_4_38, "atan", 4, 38, 2 },
    { hs_atan_4_39, "atan", 4, 39, 2 },
    { hs_atan_4_40, "atan", 4, 40, 2 },
    { hs_atan_4_41, "atan", 4, 41, 2 },
    { hs_atan_4_42, "atan", 4, 42, 2 },
    { hs_atan_4_43, "atan", 4, 43, 2 },
    { hs_atan_4_44, "atan", 4, 44, 2 },
    { hs_atan_4_45, "atan", 4, 45, 2 },
    { hs_atan_4_46, "atan", 4, 46, 2 },
    { hs_atan_4_47, "atan", 4, 47, 2 },
    { hs_atan_4_48, "atan", 4, 48, 2 },
    { hs_sin_1_8, "sin", 1, 8, 2 },
    { hs_sin_1_9, "sin", 1, 9, 2 },
    { hs_sin_1_10, "sin", 1, 10, 2 },
    { hs_sin_1_11, "sin", 1, 11, 2 },
    { hs_sin_1_12, "sin", 1, 12, 1 },
    { hs_sin_1_13, "sin", 1, 13, 1 },
    { hs_sin_1_14, "sin", 1, 14, 1 },
    { hs_sin_1_15, "sin", 1, 15, 1 },
    { hs_sin_1_16, "sin", 1, 16, 1 },
    { hs_sin_1_17, "sin", 1, 17, 1 },
    { hs_sin_1_18, "sin", 1, 18, 1 },
    { hs_sin_1_19, "sin", 1, 19, 1 },
    { hs_sin_1_20, "sin", 1, 20, 1 },
    { hs_sin_1_21, "sin", 1, 21, 1 },
    { hs_sin_1_22, "sin", 1, 22, 1 },
    { hs_sin_1_23, "sin", 1, 23, 1 },
    { hs_sin_1_24, "sin", 1, 24, 1 },
    { hs_sin_1_25, "sin", 1, 25, 1 },
    { hs_sin_1_26, "sin", 1, 26, 1 },
    { hs_sin_1_27, "sin", 1, 27, 1 },
    { hs_sin_1_28, "sin", 1, 28, 1 },
    { hs_sin_1_29, "sin", 1, 29, 1 },
    { hs_sin_1_30, "sin", 1, 30, 1 },
    { hs_sin_1_31, "sin", 1, 31, 1 },
    { hs_sin_1_32, "sin", 1, 32, 1 },
    { hs_sin_1_33, "sin", 1, 33, 1 },
    { hs_sin_1_34, "sin", 1, 34, 1 },
    { hs_sin_1_35, "sin", 1, 35, 1 },
    { hs_sin_1_36, "sin", 1, 36, 1 },
    { hs_sin_1_37, "sin", 1, 37, 1 },
    { hs_sin_1_38, "sin", 1, 38, 1 },
    { hs_sin_1_39, "sin", 1, 39, 1 },
    { hs_sin_1_40, "sin", 1, 40, 1 },
    { hs_sin_1_41, "sin", 1, 41, 1 },
    { hs_sin_1_42, "sin", 1, 42, 1 },
    { hs_sin_1_43, "sin", 1, 43, 1 },
    { hs_sin_1_44, "sin", 1, 44, 1 },
    { hs_sin_1_45, "sin", 1, 45, 1 },
    { hs_sin_1_46, "sin", 1, 46, 1 },
    { hs_sin_1_47, "sin", 1, 47, 1 },
    { hs_sin_1_48, "sin", 1, 48, 1 },
    { hs_sin_2_8, "sin", 2, 8, 5 },
    { hs_sin_2_9, "sin", 2, 9, 5 },
    { hs_sin_2_10, "sin", 2, 10, 4 },
    { hs_sin_2_11, "sin", 2, 11, 4 },
    { hs_sin_2_12, "sin", 2, 12, 4 },
    { hs_sin_2_13, "sin", 2, 13, 3 },
    { hs_sin_2_14, "sin", 2, 14, 3 },
    { hs_sin_2_15, "sin", 2, 15, 3 },
    { hs_sin_2_16, "sin", 2, 16, 3 },
    { hs_sin_2_17, "sin", 2, 17, 2 },
    { hs_sin_2_18, "sin", 2, 18, 2 },
    { hs_sin_2_19, "sin", 2, 19, 2 },
    { hs_sin_2_20, "sin", 2, 20, 2 },
    { hs_sin_2_21, "sin", 2, 21, 2 },
    { hs_sin_2_22, "sin", 2, 22, 2 },
    { hs_sin_2_23, "sin", 2, 23, 2 },
    { hs_sin_2_24, "sin", 2, 24, 2 },
    { hs_sin_2_25, "sin", 2, 25, 1 },
    { hs_sin_2_26, "sin", 2, 26, 1 },
    { hs_sin_2_27, "sin", 2, 27, 1 },
    { hs_sin_2_28, "sin", 2, 28, 1 },
    { hs_sin_2_29, "sin", 2, 29, 1 },
    { hs_sin_2_30, "sin", 2, 30, 1 },
    { hs_sin_2_31, "sin", 2, 31, 1 },
    { hs_sin_2_32, "sin", 2, 32, 1 },
    { hs_sin_2_33, "sin", 2, 33, 1 },
    { hs_sin_2_34, "sin", 2, 34, 1 },
    { hs_sin_2_35, "sin", 2, 35, 1 },
    { hs_sin_2_36, "sin", 2, 36, 1 },
    { hs_sin_2_37, "sin", 2, 37, 1 },
    { hs_sin_2_38, "sin", 2, 38, 1 },
    { hs_sin_2_39, "sin", 2, 39, 1 },
    { hs_sin_2_40, "sin", 2, 40, 1 },
    { hs_sin_2_41, "sin", 2, 41, 1 },
    { hs_sin_2_42, "sin", 2, 42, 1 },
    { hs_sin_2_43, "sin", 2, 43, 1 },
    { hs_sin_2_44, "sin", 2, 44, 1 },
    { hs_sin_2_45, "sin", 2, 45, 1 },
    { hs_sin_2_46, "sin", 2, 46, 1 },
    { hs_sin_2_47, "sin", 2, 47, 1 },
    { hs_sin_2_48, "sin", 2, 48, 1 },
    { hs_sin_3_8, "sin", 3, 8, 8 },
    { hs_sin_3_9, "sin", 3, 9, 7 },
    { hs_sin_3_10, "sin", 3, 10, 7 },
    { hs_sin_3_11, "sin", 3, 11, 6 },
    { hs_sin_3_12, "sin", 3, 12, 6 },
    { hs_sin_3_13, "sin", 3, 13, 5 },
    { hs_sin_3_14, "sin", 3, 14, 5 },
    { hs_sin_3_15, "sin", 3, 15, 5 },
    { hs_sin_3_16, "sin", 3, 16, 4 },
    { hs_sin_3_17, "sin", 3, 17, 4 },
    { hs_sin_3_18, "sin", 3, 18, 4 },
    { hs_sin_3_19, "sin", 3, 19, 4 },
    { hs_sin_3_20, "sin", 3, 20, 3 },
    { hs_sin_3_21, "sin", 3, 21, 3 },
    { hs_sin_3_22, "sin", 3, 22, 3 },
    { hs_sin_3_23, "sin", 3, 23, 3 },
    { hs_sin_3_24, "sin", 3, 24, 3 },
    { hs_sin_3_25, "sin", 3, 25, 3 },
    { hs_sin_3_26, "sin", 3, 26, 2 },
    { hs_sin_3_27, "sin", 3, 27, 2 },
    { hs_sin_3_28, "sin", 3, 28, 2 },
    { hs_sin_3_29, "sin", 3, 29, 2 },
    { hs_sin_3_30, "sin", 3, 30, 2 },
    { hs_sin_3_31, "sin", 3, 31, 2 },
    { hs_sin_3_32, "sin", 3, 32, 2 },
    { hs_sin_3_33, "sin", 3, 33, 2 },
    { hs_sin_3_34, "sin", 3, 34, 2 },
    { hs_sin_3_35, "sin", 3, 35, 2 },
    { hs_sin_3_36, "sin", 3, 36, 2 },
    { hs_sin_3_37, "sin", 3, 37, 2 },
    { hs_sin_3_38, "sin", 3, 38, 1 },
    { hs_sin_3_39, "sin", 3, 39, 1 },
    { hs_sin_3_40, "sin", 3, 40, 1 },
    { hs_sin_3_41, "sin", 3, 41, 1 },
    { hs_sin_3_42, "sin", 3, 42, 1 },
    { hs_sin_3_43, "sin", 3, 43, 1 },
    { hs_sin_3_44, "sin", 3, 44, 1 },
    { hs_sin_3_45, "sin", 3, 45, 1 },
    { hs_sin_3_46, "sin", 3, 46, 1 },
    { hs_sin_3_47, "sin", 3, 47, 1 },
    { hs_sin_3_48, "sin", 3, 48, 1 },
    { hs_sin_4_8, "sin", 4, 8, 10 },
    { hs_sin_4_9, "sin", 4, 9, 10 },
    { hs_sin_4_10, "sin", 4, 10, 9 },
    { hs_sin_4_11, "sin", 4, 11, 8 },
    { hs_sin_4_12, "sin", 4, 12, 8 },
    { hs_sin_4_13, "sin", 4, 13, 7 },
    { hs_sin_4_14, "sin", 4, 14, 7 },
    { hs_sin_4_15, "sin", 4, 15, 6 },
    { hs_sin_4_16, "sin", 4, 16, 6 },
    { hs_sin_4_17, "sin", 4, 17, 6 },
    { hs_sin_4_18, "sin", 4, 18, 5 },
    { hs_sin_4_19, "sin", 4, 19, 5 },
    { hs_sin_4_20, "sin", 4, 20, 5 },
    { hs_sin_4_21, "sin", 4, 21, 4 },
    { hs_sin_4_22, "sin", 4, 22, 4 },
    { hs_sin_4_23, "sin", 4, 23, 4 },
    { hs_sin_4_24, "sin", 4, 24, 4 },
    { hs_sin_4_25, "sin", 4, 25, 4 },
    { hs_sin_4_26, "sin", 4, 26, 4 },
    { hs_sin_4_27, "sin", 4, 27, 3 },
    { hs_sin_4_28, "sin", 4, 28, 3 },
    { hs_sin_4_29, "sin", 4, 29, 3 },
    { hs_sin_4_30, "sin", 4, 30, 3 },
    { hs_sin_4_31, "sin", 4, 31, 3 },
    { hs_sin_4_32, "sin", 4, 32, 3 },
    { hs_sin_4_33, "sin", 4, 33, 3 },
    { hs_sin_4_34, "sin", 4, 34, 3 },
    { hs_sin_4_35, "sin", 4, 35, 2 },
    { hs_sin_4_36, "sin", 4, 36, 2 },
    { hs_sin_4_37, "sin", 4, 37, 2 },
    { hs_sin_4_38, "sin", 4, 38, 2 },
    { hs_sin_4_39, "sin", 4, 39, 2 },
    { hs_sin_4_40, "sin", 4, 40, 2 },
    { hs_sin_4_41, "sin", 4, 41, 2 },
    { hs_sin_4_42, "sin", 4, 42, 2 },
    { hs_sin_4_43, "sin", 4, 43, 2 },
    { hs_sin_4_44, "sin", 4, 44, 2 },
    { hs_sin_4_45, "sin", 4, 45, 2 },
    { hs_sin_4_46, "sin", 4, 46, 2 },
    { hs_sin_4_47, "sin", 4, 47, 2 },
    { hs_sin_4_48, "sin", 4, 48, 2 },
    { hs_cos_1_8, "cos", 1, 8, 3 },
    { hs_cos_1_9, "cos", 1, 9, 3 },
    { hs_cos_1_10, "cos", 1, 10, 2 },
    { hs_cos_1_11, "cos", 1, 11, 2 },
    { hs_cos_1_12, "cos", 1, 12, 2 },
    { hs_cos_1_13, "cos", 1, 13, 2 },
    { hs_cos_1_14, "cos", 1, 14, 2 },
    { hs_cos_1_15, "cos", 1, 15, 1 },
    { hs_cos_1_16, "cos", 1, 16, 1 },
    { hs_cos_1_17, "cos", 1, 17, 1 },
    { hs_cos_1_18, "cos", 1, 18, 1 },
    { hs_cos_1_19, "cos", 1, 19, 1 },
    { hs_cos_1_20, "cos", 1, 20, 1 },
    { hs_cos_1_21, "cos", 1, 21, 1 },
    { hs_cos_1_22, "cos", 1, 22, 1 },
    { hs_cos_1_23, "cos", 1, 23, 1 },
    { hs_cos_1_24, "cos", 1, 24, 1 },
    { hs_cos_1_25, "cos", 1, 25, 1 },
    { hs_cos_1_26, "cos", 1, 26, 1 },
    { hs_cos_1_27, "cos", 1, 27, 1 },
    { hs_cos_1_28, "cos", 1, 28, 1 },
    { hs_cos_1_29, "cos", 1, 29, 1 },
    { hs_cos_1_30, "cos", 1, 30, 1 },
    { hs_cos_1_31, "cos", 1, 31, 1 },
    { hs_cos_1_32, "cos", 1, 32, 1 },
    { hs_cos_1_33, "cos", 1, 33, 1 },
    { hs_cos_1_34, "cos", 1, 34, 1 },
    { hs_cos_1_35, "cos", 1, 35, 1 },
    { hs_cos_1_36, "cos", 1, 36, 1 },
    { hs_cos_1_37, "cos", 1, 37, 1 },
    { hs_cos_1_38, "cos", 1, 38, 1 },
    { hs_cos_1_39, "cos", 1, 39, 1 },
    { hs_cos_1_40, "cos", 1, 40, 1 },
    { hs_cos_1_41, "cos", 1, 41, 1 },
    { hs_cos_1_42, "cos", 1, 42, 1 },
    { hs_cos_1_43, "cos", 1, 43, 1 },
    { hs_cos_1_44, "cos", 1, 44, 1 },
    { hs_cos_1_45, "cos", 1, 45, 1 },
    { hs_cos_1_46, "cos", 1, 46, 1 },
    { hs_cos_1_47, "cos", 1, 47, 1 },
    { hs_cos_1_48, "cos", 1, 48, 1 },
    { hs_cos_2_8, "cos", 2, 8, 6 },
    { hs_cos_2_9, "cos", 2, 9, 5 },
    { hs_cos_2_10, "cos", 2, 10, 5 },
    { hs_cos_2_11, "cos", 2, 11, 4 },
    { hs_cos_2_12, "cos", 2, 12, 4 },
    { hs_cos_2_13, "cos", 2, 13, 4 },
    { hs_cos_2_14, "cos", 2, 14, 4 },
    { hs_cos_2_15, "cos", 2, 15, 3 },
    { hs_cos_2_16, "cos", 2, 16, 3 },
    { hs_cos_2_17, "cos", 2, 17, 3 },
    { hs_cos_2_18, "cos", 2, 18, 3 },
    { hs_cos_2_19, "cos", 2, 19, 3 },
    { hs_cos_2_20, "cos", 2, 20, 2 },
    { hs_cos_2_21, "cos", 2, 21, 2 },
    { hs_cos_2_22, "cos", 2, 22, 2 },
    { hs_cos_2_23, "cos", 2, 23, 2 },
    { hs_cos_2_24, "cos", 2, 24, 2 },
    { hs_cos_2_25, "cos", 2, 25, 2 },
    { hs_cos_2_26, "cos", 2, 26, 2 },
    { hs_cos_2_27, "cos", 2, 27, 2 },
    { hs_cos_2_28, "cos", 2, 28, 2 },
    { hs_cos_2_29, "cos", 2, 29, 2 },
    { hs_cos_2_30, "cos", 2, 30, 2 },
    { hs_cos_2_31, "cos", 2, 31, 1 },
    { hs_cos_2_32, "cos", 2, 32, 1 },
    { hs_cos_2_33, "cos", 2, 33, 1 },
    { hs_cos_2_34, "cos", 2, 34, 1 },
    { hs_cos_2_35, "cos", 2, 35, 1 },
    { hs_cos_2_36, "cos", 2, 36, 1 },
    { hs_cos_2_37, "cos", 2, 37, 1 },
    { hs_cos_2_38, "cos", 2, 38, 1 },
    { hs_cos_2_39, "cos", 2, 39, 1 },
    { hs_cos_2_40, "cos", 2, 40, 1 },
    { hs_cos_2_41, "cos", 2, 41, 1 },
    { hs_cos_2_42, "cos", 2, 42, 1 },
    { hs_cos_2_43, "cos", 2, 43, 1 },
    { hs_cos_2_44, "cos", 2, 44, 1 },
    { hs_cos_2_45, "cos", 2, 45, 1 },
    { hs_cos_2_46, "cos", 2, 46, 1 },
    { hs_cos_2_47, "cos", 2, 47, 1 },
    { hs_cos_2_48, "cos", 2, 48, 1 },
    { hs_cos_3_8, "cos", 3, 8, 8 },
    { hs_cos_3_9, "cos", 3, 9, 8 },
    { hs_cos_3_10, "cos", 3, 10, 7 },
    { hs_cos_3_11, "cos", 3, 11, 7 },
    { hs_cos_3_12, "cos", 3, 12, 6 },
    { hs_cos_3_13, "cos", 3, 13, 6 },
    { hs_cos_3_14, "cos", 3, 14, 5 },
    { hs_cos_3_15, "cos", 3, 15, 5 },
    { hs_cos_3_16, "cos", 3, 16, 5 },
    { hs_cos_3_17, "cos", 3, 17, 5 },
    { hs_cos_3_18, "cos", 3, 18, 4 },
    { hs_cos_3_19, "cos", 3, 19, 4 },
    { hs_cos_3_20, "cos", 3, 20, 4 },
    { hs_cos_3_21, "cos", 3, 21, 4 },
    { hs_cos_3_22, "cos", 3, 22, 4 },
    { hs_cos_3_23, "cos", 3, 23, 3 },
    { hs_cos_3_24, "cos", 3, 24, 3 },
    { hs_cos_3_25, "cos", 3, 25, 3 },
    { hs_cos_3_26, "cos", 3, 26, 3 },
    { hs_cos_3_27, "cos", 3, 27, 3 },
    { hs_cos_3_28, "cos", 3, 28, 3 },
    { hs_cos_3_29, "cos", 3, 29, 3 },
    { hs_cos_3_30, "cos", 3, 30, 3 },
    { hs_cos_3_31, "cos", 3, 31, 2 },
    { hs_cos_3_32, "cos", 3, 32, 2 },
    { hs_cos_3_33, "cos", 3, 33, 2 },
    { hs_cos_3_34, "cos", 3, 34, 2 },
    { hs_cos_3_35, "cos", 3, 35, 2 },
    { hs_cos_3_36, "cos", 3, 36, 2 },
    { hs_cos_3_37, "cos", 3, 37, 2 },
    { hs_cos_3_38, "cos", 3, 38, 2 },
    { hs_cos_3_39, "cos", 3, 39, 2 },
    { hs_cos_3_40, "cos", 3, 40, 2 },
    { hs_cos_3_41, "cos", 3, 41, 2 },
    { hs_cos_3_42, "cos", 3, 42, 2 },
    { hs_cos_3_43, "cos", 3, 43, 2 },
    { hs_cos_3_44, "cos", 3, 44, 2 },
    { hs_cos_3_45, "cos", 3, 45, 2 },
    { hs_cos_3_46, "cos", 3, 46, 2 },
    { hs_cos_3_47, "cos", 3, 47, 1 },
    { hs_cos_3_48, "cos", 3, 48, 1 },
    { hs_cos_4_8, "cos", 4, 8, 11 },
    { hs_cos_4_9, "cos", 4, 9, 10 },
    { hs_cos_4_10, "cos", 4, 10, 9 },
    { hs_cos_4_11, "cos", 4, 11, 9 },
    { hs_cos_4_12, "cos", 4, 12, 8 },
    { hs_cos_4_13, "cos", 4, 13, 8 },
    { hs_cos_4_14, "cos", 4, 14, 7 },
    { hs_cos_4_15, "cos", 4, 15, 7 },
    { hs_cos_4_16, "cos", 4, 16, 6 },
    { hs_cos_4_17, "cos", 4, 17, 6 },
    { hs_cos_4_18, "cos", 4, 18, 6 },
    { hs_cos_4_19, "cos", 4, 19, 5 },
    { hs_cos_4_20, "cos", 4, 20, 5 },
    { hs_cos_4_21, "cos", 4, 21, 5 },
    { hs_cos_4_22, "cos", 4, 22, 5 },
    { hs_cos_4_23, "cos", 4, 23, 5 },
    { hs_cos_4_24, "cos", 4, 24, 4 },
    { hs_cos_4_25, "cos", 4, 25, 4 },
    { hs_cos_4_26, "cos", 4, 26, 4 },
    { hs_cos_4_27, "cos", 4, 27, 4 },
    { hs_cos_4_28, "cos", 4, 28, 4 },
    { hs_cos_4_29, "cos", 4, 29, 4 },
    { hs_cos_4_30, "cos", 4, 30, 4 },
    { hs_cos_4_31, "cos", 4, 31, 3 },
    { hs_cos_4_32, "cos", 4, 32, 3 },
    { hs_cos_4_33, "cos", 4, 33, 3 },
    { hs_cos_4_34, "cos", 4, 34, 3 },
    { hs_cos_4_35, "cos", 4, 35, 3 },
    { hs_cos_4_36, "cos", 4, 36, 3 },
    { hs_cos_4_37, "cos", 4, 37, 3 },
    { hs_cos_4_38, "cos", 4, 38, 3 },
    { hs_cos_4_39, "cos", 4, 39, 3 },
    { hs_cos_4_40, "cos", 4, 40, 3 },
    { hs_cos_4_41, "cos", 4, 41, 3 },
    { hs_cos_4_42, "cos", 4, 42, 2 },
    { hs_cos_4_43, "cos", 4, 43, 2 },
    { hs_cos_4_44, "cos", 4, 44, 2 },
    { hs_cos_4_45, "cos", 4, 45, 2 },
    { hs_cos_4_46, "cos", 4, 46, 2 },
    { hs_cos_4_47, "cos", 4, 47, 2 },
    { hs_cos_4_48, "cos", 4, 48, 2 },
};
static int hs_count = 656;
