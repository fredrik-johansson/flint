#!/usr/bin/env python3
"""Emit mixed-precision hand-written series for the odd families used by
the bitwise log, atan and sin/cos, in the shape of exp_rs_opt_hand.inc.

All four are series in z = x^2, so an error at Horner level k is
attenuated by TWICE as many powers of x as in the exponential:

    atanh(x) = x + x^3 V_1,   V_1 = 1/3 + z(1/5 + z(1/7 + ...))
    atan(x)  = x - x^3 V_1,   V_1 = 1/3 - z(1/5 - z(1/7 - ...))
    sin(x)   = x - x^3 V_1,   V_1 = 1/3! - z(1/5! - z(1/7! - ...))
    cos(x)   = 1 - z  V_0,    V_0 = 1/2! - z(1/4! - z(1/6! - ...))

so that err(atanh) = x^3 z^(k-1) err(V_k), giving

    bits(V_k) = 64 n - 6 - r (2k + 1)          (odd families)
    bits(V_k) = 64 n - 6 - 2 r (k + 1)         (cos)

Every V_k is positive and decreasing (z < 2^-2r swamps the alternating
term), so the recurrences need no sign handling.  Levels of one and two
limbs stay in registers; wider ones use flint_mpn_mulhigh_n.
"""
import sys
from fractions import Fraction
from math import factorial, log2

FUNCS = ("atanh", "atan", "sin", "cos")


def coeff(func, k):
    """the coefficient of z^k in V (as a Fraction), all positive"""
    if func in ("atanh", "atan"):
        return Fraction(1, 2 * k + 3)          # V_1 starts at 1/3
    if func == "sin":
        return Fraction(1, factorial(2 * k + 3))
    return Fraction(1, factorial(2 * k + 2))   # cos: V_0 starts at 1/2!


def minN(func, n, r):
    """number of V levels needed"""
    k = 0
    while True:
        c = coeff(func, k)
        # the term contributes x^3 z^k (odd) or z^(k+1) (cos)
        if func == "cos":
            mag = -2 * r * (k + 1) + log2(float(c))
        else:
            mag = -r * (2 * k + 3) + log2(float(c))
        if mag < -64 * n:
            return max(1, k)
        k += 1


def width(func, n, r, k):
    if func == "cos":
        bits = 64 * n - 6 - 2 * r * (k + 1)
    else:
        bits = 64 * n - 6 - r * (2 * k + 3)
    if bits <= 0:
        return 1
    return max(1, -(-bits // 64))


def cbits(c, L):
    """Fraction c as an L-limb fixed value, low word first"""
    v = (c.numerator << (64 * L)) // c.denominator
    return [(v >> (64 * i)) & ((1 << 64) - 1) for i in range(L)]


def emit_sqrhigh(o, dst, src, L):
    """dst = high L limbs of src^2"""
    if L == 1:
        o("    %s[0] = n_mulhi(%s[0], %s[0]);" % (dst, src, src))
    elif L == 2:
        o("    _fixed_sqrhi_2x2(&%s[1], &%s[0], %s[1], %s[0]);"
          % (dst, dst, src, src))
    else:
        o("    flint_mpn_sqrhigh(%s, %s, %d);" % (dst, src, L))


def emit_mulhigh(o, dst, a, b, L):
    """dst = high L limbs of a*b (a, b are L-limb pointers)"""
    if L == 1:
        o("    %s[0] = n_mulhi(%s[0], %s[0]);" % (dst, a, b))
    elif L == 2:
        o("    _fixed_mulhi_2x2_sloppy(&%s[1], &%s[0], %s[1], %s[0], %s[1], %s[0]);"
          % (dst, dst, a, a, b, b))
    else:
        o("    flint_mpn_mulhigh_n(%s, %s, %s, %d);" % (dst, a, b, L))


def emit_addsub(o, sub, dst, u, v, L):
    """dst = u -+ v over L limbs"""
    if L == 1:
        o("    %s[0] = %s[0] %s %s[0];" % (dst, u, "-" if sub else "+", v))
    else:
        o("    NN_%s_%d(%s, %s, %s);"
          % ("SUB" if sub else "ADD", L, dst, u, v))


def emit(o, func, n, r, name):
    N = minN(func, n, r)
    L = [min(n, width(func, n, r, k)) for k in range(N)]
    for k in range(N - 1, 0, -1):              # inner levels never wider
        L[k - 1] = max(L[k - 1], L[k])

    sub = func in ("atan", "sin", "cos")       # alternating?
    o("/* %s, n = %d, r = %d: %d levels, widths %s */"
      % (func, n, r, N, L))
    o("static void")
    o("%s(nn_ptr res, nn_srcptr x)" % name)
    o("{")
    o("    ulong z[%d], V[%d], T[%d], P[%d], Z[2];" % (n, n, n, n))
    o("    ulong v, w1, w0, h, l;")
    o("    slong i;")
    o("")
    o("    (void) v; (void) w1; (void) w0; (void) h; (void) l;")
    o("    (void) i; (void) Z; (void) V; (void) T; (void) P;")
    o("")
    emit_sqrhigh(o, "z", "x", n)
    o("")

    # innermost level
    k = N - 1
    cur = L[k]
    c = cbits(coeff(func, k), cur)
    if cur == 1:
        o("    v = UWORD(0x%016x);" % c[0])
    elif cur == 2:
        o("    w1 = UWORD(0x%016x); w0 = UWORD(0x%016x);" % (c[1], c[0]))
    else:
        for i in range(cur):
            o("    V[%d] = UWORD(0x%016x);" % (i, c[i]))

    for k in range(N - 2, -1, -1):
        prev, cur = cur, L[k]
        c = cbits(coeff(func, k), cur)
        op = "-" if sub else "+"
        o("")
        o("    /* V_%d = c %s z V_%d   (%d x %d limbs) */" % (k, op, k + 1, cur, prev))
        if cur == 1:
            o("    v = n_mulhi(z[%d], v);" % (n - 1))
            o("    v = UWORD(0x%016x) %s v;" % (c[0], op))
        elif cur == 2 and prev == 1:
            o("    _fixed_mulhi_2x1(&h, &l, z[%d], z[%d], v);" % (n - 1, n - 2))
            if sub:
                o("    sub_ddmmss(w1, w0, UWORD(0x%016x), UWORD(0x%016x), h, l);"
                  % (c[1], c[0]))
            else:
                o("    add_ssaaaa(w1, w0, UWORD(0x%016x), UWORD(0x%016x), h, l);"
                  % (c[1], c[0]))
        elif cur == 2:
            o("    _fixed_mulhi_2x2_sloppy(&h, &l, z[%d], z[%d], w1, w0);"
              % (n - 1, n - 2))
            if sub:
                o("    sub_ddmmss(w1, w0, UWORD(0x%016x), UWORD(0x%016x), h, l);"
                  % (c[1], c[0]))
            else:
                o("    add_ssaaaa(w1, w0, UWORD(0x%016x), UWORD(0x%016x), h, l);"
                  % (c[1], c[0]))
        else:
            # promote the previous level into a cur-limb array (top aligned)
            if prev == 1:
                for i in range(cur - 1):
                    o("    P[%d] = UWORD(0);" % i)
                o("    P[%d] = v;" % (cur - 1))
            elif prev == 2:
                for i in range(cur - 2):
                    o("    P[%d] = UWORD(0);" % i)
                o("    P[%d] = w0; P[%d] = w1;" % (cur - 2, cur - 1))
            else:
                for i in range(cur - prev):
                    o("    P[%d] = UWORD(0);" % i)
                for i in range(prev):
                    o("    P[%d] = V[%d];" % (cur - prev + i, i))
            if cur >= 3:
                emit_mulhigh(o, "T", "z + %d" % (n - cur), "P", cur)
            else:
                for i in range(cur):
                    o("    Z[%d] = z[%d];" % (i, n - cur + i))
                emit_mulhigh(o, "T", "Z", "P", cur)
            for i in range(cur):
                o("    P[%d] = UWORD(0x%016x);" % (i, c[i]))
            emit_addsub(o, sub, "V", "P", "T", cur)

    # tail
    o("")
    if func == "cos":
        o("    /* cos = 1 - z V_0 */")
    else:
        o("    /* %s = x %s x^3 V_0 */" % (func, "-" if sub else "+"))
    if cur == 1:
        for i in range(n - 1):
            o("    P[%d] = UWORD(0);" % i)
        o("    P[%d] = v;" % (n - 1))
    elif cur == 2:
        for i in range(n - 2):
            o("    P[%d] = UWORD(0);" % i)
        o("    P[%d] = w0; P[%d] = w1;" % (n - 2, n - 1))
    else:
        for i in range(n - cur):
            o("    P[%d] = UWORD(0);" % i)
        for i in range(cur):
            o("    P[%d] = V[%d];" % (n - cur + i, i))

    if func == "cos":
        emit_mulhigh(o, "T", "z", "P", n)
        o("    /* 1 - T: below one the unit limb is zero */")
        o("    if (mpn_zero_p(T, %d))" % n)
        o("    {")
        o("        flint_mpn_zero(res, %d);" % n)
        o("        res[%d] = 1;" % n)
        o("    }")
        o("    else")
        o("    {")
        o("        mpn_neg(res, T, %d);" % n)
        o("        res[%d] = 0;" % n)
        o("    }")
    else:
        emit_mulhigh(o, "T", "z", "x", n)      # x^3
        emit_mulhigh(o, "V", "T", "P", n)
        emit_addsub(o, sub, "res", "x", "V", n)
    o("}")
    o("")
    return N


def emit_chain(o, func, n, r, zname):
    """emit the Horner chain for one family on an already-computed z,
    leaving the result in the register/array named by the return value"""
    N = minN(func, n, r)
    L = [min(n, width(func, n, r, k)) for k in range(N)]
    for k in range(N - 1, 0, -1):
        L[k - 1] = max(L[k - 1], L[k])
    sub = True                              # sin and cos both alternate
    pre = "s" if func == "sin" else "c"
    k = N - 1
    cur = L[k]
    c = cbits(coeff(func, k), cur)
    if cur == 1:
        o("    %sv = UWORD(0x%016x);" % (pre, c[0]))
    elif cur == 2:
        o("    %sw1 = UWORD(0x%016x); %sw0 = UWORD(0x%016x);"
          % (pre, c[1], pre, c[0]))
    else:
        for i in range(cur):
            o("    %sV[%d] = UWORD(0x%016x);" % (pre, i, c[i]))
    for k in range(N - 2, -1, -1):
        prev, cur = cur, L[k]
        c = cbits(coeff(func, k), cur)
        o("")
        o("    /* %s V_%d  (%d x %d) */" % (func, k, cur, prev))
        if cur == 1:
            o("    %sv = n_mulhi(%s[%d], %sv);" % (pre, zname, n - 1, pre))
            o("    %sv = UWORD(0x%016x) - %sv;" % (pre, c[0], pre))
        elif cur == 2 and prev == 1:
            o("    _fixed_mulhi_2x1(&h, &l, %s[%d], %s[%d], %sv);"
              % (zname, n - 1, zname, n - 2, pre))
            o("    sub_ddmmss(%sw1, %sw0, UWORD(0x%016x), UWORD(0x%016x), h, l);"
              % (pre, pre, c[1], c[0]))
        elif cur == 2:
            o("    _fixed_mulhi_2x2_sloppy(&h, &l, %s[%d], %s[%d], %sw1, %sw0);"
              % (zname, n - 1, zname, n - 2, pre, pre))
            o("    sub_ddmmss(%sw1, %sw0, UWORD(0x%016x), UWORD(0x%016x), h, l);"
              % (pre, pre, c[1], c[0]))
        else:
            if prev == 1:
                for i in range(cur - 1):
                    o("    P[%d] = UWORD(0);" % i)
                o("    P[%d] = %sv;" % (cur - 1, pre))
            elif prev == 2:
                for i in range(cur - 2):
                    o("    P[%d] = UWORD(0);" % i)
                o("    P[%d] = %sw0; P[%d] = %sw1;" % (cur - 2, pre, cur - 1, pre))
            else:
                for i in range(cur - prev):
                    o("    P[%d] = UWORD(0);" % i)
                for i in range(prev):
                    o("    P[%d] = %sV[%d];" % (cur - prev + i, pre, i))
            if cur >= 3:
                emit_mulhigh(o, "T", "%s + %d" % (zname, n - cur), "P", cur)
            else:
                for i in range(cur):
                    o("    Z[%d] = %s[%d];" % (i, zname, n - cur + i))
                emit_mulhigh(o, "T", "Z", "P", cur)
            for i in range(cur):
                o("    P[%d] = UWORD(0x%016x);" % (i, c[i]))
            emit_addsub(o, True, "%sV" % pre, "P", "T", cur)
    return cur


def emit_sincos(o, n, r, name):
    """sin and cos together: they are both series in z = x^2, so ONE
    squaring serves both"""
    o("/* sin_cos, n = %d, r = %d */" % (n, r))
    o("static void")
    o("%s(nn_ptr ysin, nn_ptr ycos, nn_srcptr x)" % name)
    o("{")
    o("    ulong z[%d], sV[%d], cV[%d], T[%d], P[%d], Z[2];" % (n, n, n, n, n))
    o("    ulong sv, sw1, sw0, cv, cw1, cw0, h, l;")
    o("")
    o("    (void) sv; (void) sw1; (void) sw0; (void) cv; (void) cw1;")
    o("    (void) cw0; (void) h; (void) l; (void) Z; (void) T; (void) P;")
    o("    (void) sV; (void) cV;")
    o("")
    o("    /* the one squaring shared by both series */")
    emit_sqrhigh(o, "z", "x", n)
    o("")
    cs = emit_chain(o, "sin", n, r, "z")
    o("")
    cc = emit_chain(o, "cos", n, r, "z")
    o("")
    o("    /* sin = x - x^3 V_0 */")
    if cs == 1:
        for i in range(n - 1):
            o("    P[%d] = UWORD(0);" % i)
        o("    P[%d] = sv;" % (n - 1))
    elif cs == 2:
        for i in range(n - 2):
            o("    P[%d] = UWORD(0);" % i)
        o("    P[%d] = sw0; P[%d] = sw1;" % (n - 2, n - 1))
    else:
        for i in range(n - cs):
            o("    P[%d] = UWORD(0);" % i)
        for i in range(cs):
            o("    P[%d] = sV[%d];" % (n - cs + i, i))
    emit_mulhigh(o, "T", "z", "x", n)
    emit_mulhigh(o, "sV", "T", "P", n)
    emit_addsub(o, True, "ysin", "x", "sV", n)
    o("    ysin[%d] = 0;" % n)
    o("")
    o("    /* cos = 1 - z V_0 */")
    if cc == 1:
        for i in range(n - 1):
            o("    P[%d] = UWORD(0);" % i)
        o("    P[%d] = cv;" % (n - 1))
    elif cc == 2:
        for i in range(n - 2):
            o("    P[%d] = UWORD(0);" % i)
        o("    P[%d] = cw0; P[%d] = cw1;" % (n - 2, n - 1))
    else:
        for i in range(n - cc):
            o("    P[%d] = UWORD(0);" % i)
        for i in range(cc):
            o("    P[%d] = cV[%d];" % (n - cc + i, i))
    emit_mulhigh(o, "T", "z", "P", n)
    o("    if (mpn_zero_p(T, %d))" % n)
    o("    {")
    o("        flint_mpn_zero(ycos, %d);" % n)
    o("        ycos[%d] = 1;" % n)
    o("    }")
    o("    else")
    o("    {")
    o("        mpn_neg(ycos, T, %d);" % n)
    o("        ycos[%d] = 0;" % n)
    o("    }")
    o("}")
    o("")

def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    out = []
    o = out.append
    ent = []
    for func in FUNCS:
        for n in range(1, 5):
            for r in range(4, 49):
                if r > 64 * n - 16:
                    continue
                nm = "hs_%s_%d_%d" % (func, n, r)
                N = emit(o, func, n, r, nm)
                ent.append((nm, func, n, r, N))
    o("typedef void (*hs_fn)(nn_ptr, nn_srcptr);")
    o("typedef struct { hs_fn f; const char * name; int n; int r; int N; } hs_entry;")
    o("static hs_entry hs_all[] = {")
    for nm, func, n, r, N in ent:
        o('    { %s, "%s", %d, %d, %d },' % (nm, func, n, r, N))
    o("};")
    o("static int hs_count = %d;" % len(ent))
    print("\n".join(out))


if __name__ == "__main__":
    main()
