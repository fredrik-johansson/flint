# HANDOFF — FLINT `fixed` module

Everything is in `/mnt/user-data/outputs/`:

- `0001..0015-*.patch` — apply **in order** to a FLINT tree; all verified clean
  (`git apply --check --whitespace=error`, then apply, in sequence)
- `fixed/` — a mirror of the module sources, `fixed.h`, `fixed.rst`, `test/`,
  `tune/`, and `dev/` (all 8 generators + 2 prototypes)
- `README.md` — the running project log; the per-patch entries have the
  reasoning, the measurements, and the negative results

Nothing is withheld. All 8 generators are committed inside the patches.

---

## 1. The generators (all in `dev/`, all committed)

| generator | emits | notes |
|---|---|---|
| `gen_fixed_exp_rs_hard.py` | `exp_rs_hard.inc` | RS series for exp. Extended by me: explicit `xn`, **arbitrary zbits**, `noshrink`, `WINNERS_OPT`. |
| `gen_fixed_trig_rs_hard.py` | `trig_rs_hard.inc` | RS series for atanh/atan/sin_cos/sinh_cosh. Extended: **alternating families at zbits=16** via `RS_*_ACC`. |
| `gen_fixed_exp_bitwise_small.py` | `exp_bitwise_rs_small.inc` | register reduction for exp, incl. constant-r `opt` variants |
| `gen_fixed_log1p_bitwise_small.py` | `log1p_bitwise_rs_small.inc` | ditto for log1p |
| `gen_fixed_atan_bitwise_small.py` | `atan_bitwise_rs_small.inc` | ditto for atan |
| **`gen_fixed_exp_hand_mixed.py`** | hand exp series, any `(n, r)` | **mixed precision**: level widths from `bits(V_k) = 64n − 6 − r·k` |
| **`gen_fixed_trig_hand_mixed.py`** | hand atanh/atan/sin/cos/**tan** series, any `(func, n, r)` | odd families: `bits(V_k) = 64n − 6 − r(2k+3)`. Also `emit_sincos()` (shared square) |
| **`gen_fixed_tan_rotate.py`** | `tan_rotate.inc` | the W rotation in **pure registers**, n ≤ 6 |

`gen_fixed_sin_cos_bitwise_small.py` was **deleted** in 0013 — the tangent
reconstruction made it dead.

---

## 2. The load-bearing ideas (carry these, they are not obvious)

**Error attenuation is what makes everything cheap.** Writing a series as
`f = ... + x^k·V_k`, an error in an inner Horner value reaches the result
multiplied by the powers of `x` above it. So

    exp:          bits(V_k) = 64n − 6 − r·k
    odd families: bits(V_k) = 64n − 6 − r(2k+3)      (series in z = x², so 2× faster)

Most levels collapse to 1–2 limbs. **Choose every operation by width**:
len 1 → `n_mulhi`; len 2 → the 2×2 forms in `hand_mulhi.inc`; len ≥3 →
`flint_mpn_mulhigh_n`/`sqrhigh`; adds → `NN_ADD_k`/`NN_SUB_k`. A length-1
`flint_mpn_sqrhigh` is *pure overhead* — inlining alone halved atan at n=1.

**A cheaper series wants a SHORTER reduction.** Every time a series got faster,
the optimal `r` moved down, sometimes a long way (atan: 16/32 → 4, 6, 22, 20).
**Always re-sweep r after touching a series.** The sweep harnesses are the
`/tmp/rsw_*.c` pattern in the log; rebuild them from the generators.

**Tan is scale-invariant — this is the whole trick.** `arg(W) = Σ atan(2^−i)`,
so `tan(Σ Aᵢ) = wy/wx`: **|W|'s growth cancels in the ratio.** No `|W|²`, no
normalization, no unimodularity. Everything else follows from
`t = tan(x/2) = (wy + wx·u)/(wx − wy·u)`.

**Bit-exactness harnesses catch what tests miss.** Every register path was
diffed against its generic twin over 400k cases. That found the atan ordering
bug (X updated from the *new* Y) in minutes.

---

## 3. Traps that bit me more than once

- **Unnormalized divisors.** `mpn_tdiv_qr` requires the divisor's top limb
  nonzero. Every denominator in this module sits *near 1* (`wx ∈ (0.72,1.17)`,
  `1±t²`), so the top limb is often zero. **Trim the divisor length.** This bug
  appeared three separate times (stress test in 0002, the tan prototype, and
  the divisor rework below).
- **`2t` exceeds 1** (`t < tan(½) = 0.5464`) — outputs need the unit limb.
- **The greedy reduce's window model requires `tabᵢ < 2^−i`.** `log(1+2^−i)` and
  `atan(2^−i)` satisfy it; the *doubled* angles `2·atan(2^−i) ≈ 2^{1−i}` do not.
  Hence the undoubled table + reducing `x/2`.
- **Truncation creep**: floored table entries under-subtract, so the boundary
  and final steps must *repeat* the subtraction (`FIXED_BITWISE_REDUCE_USED_ALLOC`).
- Container noise is ±30–50% between passes. **Only compare within one process.**

---

## 4. Open work, in the order I would do it

1. **Re-tune `r` for sin/cos at n=3,4.** Cheap, and they are the only sizes below
   arb (0.86×, 0.67×) while n=5,6 sit at 1.03×. Their `r` was picked for the
   *old* conjugate-ratio tail and never re-swept. Harness already exists.

2. **Divisor normalization — ~20%, and I measured it working.**
   `1 + t² ∈ [1, 1.2984)` is the worst possible divisor: n+1 limbs to hold *one
   bit* above the point. Halve it:

       D = (1+t²)/2 ∈ [0.5, 0.6492)   →  n limbs, TOP BIT SET

   then `sin = t/D`, `cos = ((1−t²)/2)/D`. **One limb shorter *and* normalized.**
   For tan, `1−t² ∈ (0.7015, 1]` already has its top bit set (t=0 excepted).

   With this in, **sin/cos passed every test and measured ~20% faster**
   (148.8→118.3, 255.6→201.2, 412.0→324.2, 465.4→385.3 for n=1..4). **Only the
   `ytan` branch failed**, and I reverted rather than ship it broken.
   Suspects I did *not* rule out: the quotient is **n+2 limbs** (tan(1)<1.56
   needs n+1, but GMP writes n+2 for a (2n+1)/(n) division); and `mpn_tdiv_qr`
   overlap rules. A separate quotient buffer did *not* fix it, so look at the
   `2t` left-shift and the numerator length first.

3. **Benchmark 2×`mpn_tdiv_qr` vs 1 reciprocal + 2 mulhighs** for sin/cos.

4. **Specialize the sin/cos-from-tan reconstruction for n=1..4** (everything but
   the division in registers), same shape as `gen_fixed_tan_rotate.py`.

5. `mpn_tdiv_qr` itself can be beaten (Fredrik's note) — not urgent.

6. Cosmetic: patch 0002's commit message and the `fixed.rst` sin/cos section
   still describe the conjugate-ratio construction, which 0013 deleted. If you
   want clean history, 0002 could be rewritten to introduce the tangent form
   directly — a rebase, so I left it.

---

## 5. Verification recipe (what I ran every time)

**Regenerated 2026-07: patch 0011 originally shipped a stray scratch file,
`src/fixed/trighand.c` — a bare dump of generated series with no includes.
FLINT's build globs every `.c` in the module directory, so it broke the build.
It is removed and the patches are regenerated. The lesson is in the recipe
below: I had been compiling explicit file lists, which never exercised the
glob. ALWAYS build the module the way the build system does.**

    # every module .c must compile STANDALONE -- the build globs them.
    # (An explicit file list will NOT catch a stray scratch .c in the dir.)
    for f in src/fixed/*.c; do gcc -O2 -Wall -I src -I . -I src/fixed -c $f -o /tmp/o.o || echo "FAILS: $f"; done

    # 64-bit
    gcc -O2 -march=native -Wall -I src -I . -I src/fixed -o t \
        src/fixed/test/main.c src/fixed/*.c -L. -Wl,-rpath,$(pwd) -lflint -lmpfr -lgmp -lm
    FLINT_TEST_MULTIPLIER=10 ./t          # expect 8 PASS

    # 32-bit (separate tree, ABI=32; the hand paths must compile out)
    gcc -m32 -O2 -DFLINT_WANT_ASSERT -Wall ...   # expect 8 PASS

    # patches
    cd <clean tree>; for p in 0*.patch; do git apply --check --whitespace=error $p && git apply $p; done

Coverage: `-O0 --coverage`, `FLINT_TEST_MULTIPLIER=25`, then `gcov`.
**Rebuild *all* objects** — a stale `.o` once made me think coverage had dropped.
