# HANDOFF — FLINT `fixed` module

The branch itself is the source of truth: everything below is committed.
`/mnt/user-data/outputs/` carries `0016..0018-*.patch` (this session's
commits, `git format-patch` off the squashed base, verified to apply
clean and to reproduce the working tree exactly), a mirror of
`src/fixed/`, `src/fixed.h`, `doc/source/fixed.rst`, and — new —
`dev/harnesses/` with every harness this session used ready to rebuild:
`rsw_tan.c` (the r sweep), `canary.c` (out-of-contract writes),
`tailbench.c` (division-strategy comparison), `errcheck.c`/`prodcheck.c`
(MPFR ulp validation), `spec1.c` (the n = 1 register prototype),
`bench_arb.c` (arb ratios), and `gen_tan_sweep.py`/`gen_prod_tan.py`
(drive `gen_fixed_trig_hand_mixed.py` for the sweep series and for
regenerating `trig_rs_opt_hand.inc`'s tan section — the tuned (n, r)
PICK list lives in the latter).  The per-entry log with all
measurements and negative results is `src/fixed/README.md`.

---

## 1. What earlier sessions shipped

1. **Divisor normalization — DONE, and the old ytan bug explained.**
   It was never the normalization: the old tan division fed a
   (2n+1)-limb numerator to a divisor trimmed to n limbs, and
   `mpn_tdiv_qr` writes exactly `nn − dn + 1 = n + 2` quotient limbs —
   ONE PAST the documented `(res, n+1)` contract.  The tests allocate
   `res[66]` and never saw it; a canary at `res[n+1]` trips on the old
   code.  The fix and the speedup are one change: **halve numerator and
   denominator together** so every divisor is a normalized n-limb value
   and every quotient fits n+1 limbs by construction
   (`sin = t/D`, `cos = ((1−t²)/2)/D`, `D = (1+t²)/2`; `tan` as twice
   `t/(1−t²)`).  The same conditional halving now also normalizes the
   FIRST division `t = (wy + wx·u)/(wx − wy·u)` (numerator < 0.64 never
   uses its unit limb), and the large-n `tan = sin/cos` fallback got the
   same quotient-length fix.

2. **r re-swept for ALL n = 1..6** (was: just 3, 4).  New
   `_fixed_tan_opt_r = {., 4, 6, 9, 14, 15, 16}` (was 5, 6, 15, 27, 27,
   30); series regenerated in `trig_rs_opt_hand.inc`.  As always, the
   cheaper tail moved every optimum down.  n = 6 tied 15/16; 16 taken
   for cleaner measured error.

3. **One reciprocal + two mulhighs replaces the two sin/cos divisions**
   (open item 3, measured all three ways): `R = 1/(1+t²) = 1/(2D)` by
   one division of `2^(128n−1) − 1` (NOT `2^(128n−1)`: when t² is one
   ulp, `D = 2^(64n−1)` exactly and the quotient would hit `2^(64n)`),
   then `sin = 2tR`, `cos = (1−t²)R`.  A third to a half off the tail
   at every n, ≤ 2 ulp from the divisions.  `flint_mpn_preinvn` +
   `divrem_preinvn` loses even to two plain divisions here — Newton
   setup only amortizes over many divisions per divisor.  Used only
   when BOTH outputs are wanted; single outputs keep one direct
   division; the t² = 0 pair case is now division-free (sin = 2t,
   cos = 1).

4. **n = 1 reconstruction in pure registers** (open item 4, first
   size): after the shared reduction, everything is umul_ppmm chains
   and one udiv_qrnnd per division, each divisor rounded to 64
   normalized bits by the same conditional halving.  99 → 42 ns; the
   mpn dispatch and buffer traffic WAS two thirds of the n = 1 call.
   Dispatched only for `(n == 1, r == 0)`, so the explicit-r test paths
   still exercise the generic route.

Current standing vs arb (one process, ratio = arb/fixed):

    n            1     2     3     4     5     6     7     8
    sin_cos   4.91  2.01  1.56  1.64  1.70  1.52  0.91  1.16
    tan       6.91  3.03  2.19  2.27  2.22  1.99  1.17  1.38

Production errors at n = 1..6, 20k adversarial-biased points per size:
sin ≤ 27, cos ≤ 18, tan ≤ 72 ulp (budgets 6r+128, 8r+256 with the
per-size r).

---

## 2. The load-bearing ideas (carried over; still not obvious)

**Error attenuation is what makes everything cheap.** Writing a series
as `f = ... + x^k·V_k`, an inner error reaches the result multiplied by
the powers of `x` above it: `bits(V_k) = 64n − 6 − r·k` (exp),
`64n − 6 − r(2k+3)` (odd families).  Choose every operation by width:
len 1 → `n_mulhi`; len 2 → `hand_mulhi.inc`; len ≥3 → `flint_mpn_
mulhigh_n`/`sqrhigh`; adds → `NN_ADD_k`/`NN_SUB_k`.

**A cheaper tail wants a shorter reduction.** Confirmed AGAIN this
session (r fell at every n).  Always re-sweep r after touching a series
OR the reconstruction.  Harness pattern: `/tmp/rsw_tan.c` in the log —
a parametrized copy of `_fixed_tan_halfangle` + generated `hs_tan_n_r`
series, best-of-5 timing, two independent processes to confirm ranking.

**Tan is scale-invariant — this is the whole trick.** `arg(W) = Σ Aᵢ`,
so `tan(Σ Aᵢ) = wy/wx`: |W|'s growth cancels in the ratio.

**Every divisor in this module can be made n limbs with its TOP BIT
SET** by conditionally halving numerator and denominator together —
the numerators are all safely below their denominators' scale.  This
is worth more than it looks: no trimming, no internal normalization
shifts, and (the actual bug class) quotient lengths become predictable.

**Quotient length is part of the contract.** `mpn_tdiv_qr` writes
exactly `nn − dn + 1` limbs.  Sentinel-canary tests at `res[n+1]`
catch what the value tests miss; keep `/tmp/canary.c` (in the log)
alive.

**Container noise is ±30–50% between passes** — and a background build
eats the timings entirely (bitten this session).  Only compare within
one process, on an idle machine; the arb ratio is the durable metric.

---

## 3. Traps that bit (old and new)

- **Unnormalized divisors / n+2-limb quotients** — see above; struck a
  fourth time this session before being killed at the root.
- **`2t` exceeds 1** (`t < tan(½) = 0.5464`) — outputs carry a unit limb.
- **Reciprocal numerators must be all-ones** (`2^(128n−1) − 1`), or the
  `D = 2^(64n−1)` corner overflows R into limb n.
- **Zero the WHOLE numerator buffer**: a partial `flint_mpn_zero(N, n)`
  before setting the top limb left stale middle limbs and produced
  errors that looked like broken algebra; exact arithmetic (a 10-line
  Python check) pinned it on the buffer in minutes.  Cheap and decisive
  — reach for it before staring at identities.
- **The greedy reduce's window model requires `tabᵢ < 2^−i`** — hence
  the undoubled table and reducing `x/2`.
- **Truncation creep**: boundary and final reduce steps may repeat an
  index (`FIXED_BITWISE_REDUCE_USED_ALLOC`).
- **32-bit environment** (not module) potholes: Ubuntu's fat i386
  libgmp exports `__gmpn_add_nc` etc. only with CPU suffixes but
  configure misdetects them ("none required") — comment the
  `FLINT_HAVE_NATIVE_mpn_*` lines out of **`src/config.h`** (NOT
  `src/flint-config.h`, which doesn't carry these).  The merged-object
  step needs `make LD="ld -m elf_i386"`.

---

## 4. Open work, in the order I would do it

(Old items 1 and 2 shipped this session: register reconstruction runs
through n = 4 -- `dev/gen_fixed_tan_halfangle_small.py` ->
`tan_halfangle_small.inc`, dispatched at r = 0 -- and n = 7 has its
own tangent series plus register rotation, with n = 8 getting a series
too.  The tail was also restructured so that NO tangent is ever
divided out: A = TT DE, B = DE^2, C = TT^2, sin = 2A/(B+C),
cos = 1 - 2C/(B+C), tan = 2A/(B-C), one reciprocal per output, and the
n > 8 fallback forms TT and DE from sin/cos of the residual WITHOUT
dividing them -- cos t' cancels.  See the last three README entries.)

1. **n = 8..12 register rotation** needs 9-wide and wider add/sub
   chains that `longlong.h` stops short of (NN_ADD_8/NN_SUB_8 are the
   widest).  Either add the macros, or emit the unit limb's borrow by
   a comparison; the n = 7 -> 8 sin_cos step (+38% for one limb, when
   9 -> 10 costs ~18%) says the generic rotation loop is worth ~150 ns
   at n = 8.  Tangent series now run through n = 12 (r = 16, 19, 23,
   25 for 9..12); n = 13 is the fallback boundary (+63%, ratio 1.13
   against arb, recovering to 1.28 by n = 16).  Extending the series
   further is mechanical if ever wanted -- the _TAN table holds 65
   coefficients and the generator's add/sub emitter handles any width
   -- but past 12 the returns thin as arb's own overhead amortizes.

2. **`mpn_tdiv_qr` itself can be beaten** (Fredrik's note) -- after
   this session exactly ONE tdiv_qr per requested reciprocal remains,
   so the ceiling is small; a 2n/n schoolbook with a cached top-limb
   preinverse might still shave the call overhead at n = 2..4, where
   the division is now the majority of the tail.

3. Cosmetic: patch 0002's commit message still describes the
   conjugate-ratio construction (a rebase, so left).

---


## 5. Verification recipe (run every time)

    # every module .c must compile STANDALONE -- the build GLOBS them
    for f in src/fixed/*.c; do gcc -O2 -Wall -I src -I . -I src/fixed -c $f -o /tmp/o.o || echo "FAILS: $f"; done

    # 64-bit
    gcc -O2 -march=native -Wall -I src -I . -I src/fixed -o t \
        src/fixed/test/main.c -L. -Wl,-rpath,$(pwd) -lflint -lmpfr -lgmp -lm
    FLINT_TEST_MULTIPLIER=10 ./t          # expect 8 PASS

    # canary for out-of-contract writes (source in the log): expect clean
    # 32-bit tree: ABI=32, see the environment potholes above; expect 8 PASS
    # patches: git apply --check --whitespace=error, in order

    # benchmarks: idle machine, one process, best-of-5; ratio vs arb
