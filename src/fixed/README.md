# Optimized Taylor series evaluation for FLINT/arb

Drop-in replacements for `_arb_exp_taylor_rs`, `_arb_atan_taylor_rs`,
`_arb_sin_cos_taylor_rs` (FLINT 3.3.x):

```c
void _arb_exp_taylor_rs2(nn_ptr y, ulong * error, nn_srcptr x, slong xn, ulong N);
void _arb_atan_taylor_rs2(nn_ptr y, ulong * error, nn_srcptr x, slong xn, ulong N, int alternating);
void _arb_sin_cos_taylor_rs2(nn_ptr ysin, nn_ptr ycos, ulong * error, nn_srcptr x, slong xn, ulong N, int sinonly, int alternating);
```

Same calling conventions, preconditions (`0 <= x <= 1/16`, same N limits)
and output semantics as the classical functions; a rigorous ulp error
bound is returned in `error` (typically 8–16 ulp instead of the classical
2 ulp — callers absorb this in the radius as usual, costing at most a few
guard bits). On 32-bit systems, and on small inputs where the classical
code is faster, the functions transparently dispatch to the originals.

## Architecture

Each entry point dispatches between four implementations of the
rectangular-splitting sum `S = sum_k (+-1)^k c_k z^k` (`z = x` for exp,
`z = x^2` otherwise; `zb` = leading zero bits of `z`):

1. **hardcoded** small cases (exp `N <= 4`, atan `N <= 2`, sin/cos
   `N <= 1`), using sqrhigh/mulhigh;
2. **mid**: the classical fixed-precision structure with cheaper
   primitives, for small operands (`xn < 12` and `N*zb < 560`);
3. **L**: a whole-limb specialization for `zb >= 64` and `xn >= 8`
   (`L = zb >> 6` limbs of precision change per term);
4. **general**: bit-level variable-precision windows, for everything
   else.

All variants share the inverse tables, the preinv division helper and
the boundary-multiplication helper.

## The general variant

1. **Variable working precision.** All quantities live on a common
   fixed-point grid of `xn` fractional limbs. With `zb` a lower bound on
   the leading zero bits of `z` (`2*zbx` is used for `z = x^2`), a partial
   sum that is still going to be multiplied by `z^k` is kept in a *window*
   of the top `wn(k) = xn - max(0, ((k*zb) >> 6) - 2)` fractional limbs
   plus the units limb. The schedule is maintained by addition-only
   recurrences on `k*zb`; the only division in the control flow is one
   upfront `Nmax = 64*(xn+1)/zb + 2` used to drop terms lying entirely
   below the final ulp. The two guard limbs make every truncation
   committed at a reduced window worth `< 2^-1` final ulp *after* the
   pending multiplication by `z^k` (see error analysis in `taylor_rs2.c`).

2. **mulhigh everywhere.** `flint_mpn_mulhigh_n` / `flint_mpn_sqrhigh`
   (precise, one-sided, `< 2` ulp of the bottom returned limb) replace all
   full products: the power table `z^1..z^m`, the Horner boundary
   multiplications by `z^m`, and divisions by precomputed inverses.
   Powers are stored *top-aligned* in zero-initialized slots of `xn + 2`
   limbs with `floor(i*zb/64)` implied leading-zero limbs and one guard
   limb below the grid; “the top `h` limbs of `z^i`, zero-extended” is
   then just `slot + pslot - h`, giving padding-free unequal-length
   mulhighs (generalizing the zero-padding tricks of the
   proof-of-concept).

3. **Divisions by one-limb constants.** Windows of at most
   `DIV_PREINV_CUTOFF (= 10) + 1` limbs are divided by a single mulhigh
   with a precomputed `DINV_LIMBS`-limb inverse of the block denominator
   (`floor(2^(64*16)/den)`); the window is padded with one zero limb below
   so that truncating the inverse costs `< 2` window-bottom ulp even when
   the value has a full units limb. Longer windows use `mpn_divrem_1`.
   The inverse tables (60 KB) are computed lazily on first use; the
   initialization is idempotent so the race is benign, but production
   code should emit them as constant tables or use one-time init.

4. **Hardcoded small cases** for exp `N <= 4`, atan `N <= 2`,
   sin/cos `N <= 1`, using `sqrhigh`/`mulhigh`.

5. **Micro-optimized inner loop**: running power pointer (no multiply per
   term), direct carry into the units limb when the power has no implied
   leading-zero limbs (the common case), branch-free window bookkeeping.

## The mid variant (small operands)

For `xn < 12` with `N*zb < 560` the window bookkeeping costs more than it
saves, so the mid variant keeps the classical layout exactly (xn-limb
power table, full (xn+1)-limb accumulator, same loop) and only swaps
primitives: sqrhigh/mulhigh powers, one mulhigh against a zero-padded
`z^m` per Horner boundary, preinv divisions, and the upfront term clamp
(applied only when `N*zb >= 64*(xn+1)`, so the small-input path never
pays the division).  Error bounds: 11 (exp), 18 (atan), 12 (sin/cos)
ulp, derived in the source (the power-chain error converges to <= 2.28
ulp because each halving-and-squaring step is damped by `z <= 1/16`).

Measured honestly: the gain over the classical code here is only
~1.0-1.1x, occasionally more (e.g. exp at `zb = 64, xn = 4`: ~1.5-2x from
the clamp + preinv final division).  The reason is that FLINT's hardcoded
assembly `mul_n`/`sqr` small cases are already nearly as fast as mulhigh
at `n <= 3-4`, and at small N the first denominator group covers the
whole sum so there are almost no divisions to replace.  The mid variant
is kept because it is never slower beyond a few ns of dispatch overhead
and it wins whenever any of its levers engages.

## The L variant (z < 2^-64)

When `zb >= 64` we may take `lz_i = L*i` with `L = zb >> 6` (a valid
lower bound; the resulting windows are at least as wide as the general
schedule, so the general error analysis applies verbatim).  Alignment
becomes exact-limb, which yields the simplifications of the
proof-of-concept generalized to arbitrary `xn`, `N`, block structure:

- the shift deltas between power slots vanish, so `z^2..z^m` are
  sqrhigh/mulhigh'd **directly into their top-aligned slots** — no
  scratch pass, no copies;
- the inner loop maintains the addmul length `un = wn - L*power` with a
  single addition per term (no lz loads);
- terms whose windowed length is `<= 0` are skipped outright.

This is a consistent 3-18% improvement over the general path in its
regime (largest around `xn = 8-16`), and it takes over the
`zb >= 64, xn >= 8` corner where the general path only broke even.

## Error bound

The reported bound is assembled from a constant part (5 for exp and
sin/cos, 6 for atan) covering the clamped tail, the power-table guard
truncations and all reduced-window contributions, plus dynamically
counted full-precision events (denominator changes, boundary
multiplications touching the final grid, the final division, and the
final multiplication by `x` for atan/sin). The full accounting is in the
header comment of `taylor_rs2.c`. All contributions are one-sided
(downward) except through alternating subtraction; the bound is valid
two-sided.

The test program verifies `|computed - S| <= error * 2^(-64*xn)`
*rigorously* (evaluating `S` in ball arithmetic at `64*xn + 192` bits and
including the reference radius in the comparison) on randomized inputs
covering `xn` up to 60, the full range of `zb`, all `N`, and all flag
combinations. Across 25,000 instances (exercising all four dispatch tiers) the maximum
observed error was 1.7 ulp against reported bounds of 8–18, and the max
observed/reported ratio was 0.5.

## Hardcoded generated evaluators (gen_exp_series.py)

`hardcoded/gen_exp_series.py` generates fully unrolled evaluators

    void mpn_exp_series_{N}(mp_ptr res, mp_srcptr x)

for `sum_{k<N} x^k/k!` with `xn = N` fractional limbs and
`0 <= x < 2^-64` (top limb never read), generalizing the hand-written
`mpn_exp_series_10_rs4` proof of concept, which was first verified by
full alignment tracing and 60,000 rigorous arb comparisons (max
observed 3.11 ulp, one-sided, derived bound 6 ulp).

The generated functions add 1, x, x^2/2 exactly and evaluate k = 3..N-1
by rectangular splitting with block length m over the common denominator
20!: coefficients 20!/k! all fit one limb (hence N <= 21) and there is a
single final division, emitted as a mulhigh by a truncated inverse of
20! or as mpn_divrem_1.  Because x < 2^-64, every alignment is
limb-exact (windows shrink one limb per term), each addmul is exactly
one limb longer than the previous so carry limbs chain with no
propagation, and the division by 20! damps all splitting-internal
truncation errors by 2^-61, making guard-free windows sound.  The RS
sum is a pure fraction < c_3 x^3, so its top limb is zero, the division
length is always N-2 and reading the top N-2 limbs of the inverse is
valid independent of N.  One-sided error <= 8 ulp (11 for m = 2),
verified for all 132 (N, m) x {mulhigh, divrem} variants against arb
(global max observed 3.66 ulp).

Profiling every m per N (three interleaved passes, per-candidate
minimum) gives: m = 3 optimal for N = 5-8, m = 4 for N >= 10 (m = 5-6
tie within noise near N = 9, 12); the final division is fastest as
mulhigh for N <= 11 (length <= 9, inside FLINT's hardcoded mulhigh
assembly table) and as divrem_1 for N >= 12.  `winners.json` +
`gen_exp_series.py best` emit `exp_series_hard.c` with one winner per N
plus dispatch tables `mpn_exp_series_tab[]` / `mpn_exp_series_err_tab[]`.

Measured against `_arb_exp_taylor_rs` at N = xn, zb = 64: 2.7-4.4x
standalone.  Building `taylor_rs2.c` with `-DEXP_SERIES_HARD` and
linking `exp_series_hard.c` adds these as a dispatch tier in
`_arb_exp_taylor_rs2` (taken when `zb >= 64`, `N >= xn`,
`5 <= xn <= 21`; extra terms beyond a larger requested N differ by
< x^xn/xn! << 1 ulp), lifting those cells from 1.0-2.1x to 2.6-3.8x.

## Hardcoded evaluators for x < 2^-32 (gen_exp_series.py, zbits = 32)

The same generator emits `mpn_exp_series32_{N}` for N = 1..21 with
xn = ceil(N/2) fractional limbs and 0 <= x < 2^-32, targeting the
few-limb precision range of aggressive table-based argument reduction.
All bookkeeping stays limb-granular via lz(j) = floor(j*zbits/64):
power products pick factor pairs with lz(a) + lz(b) = lz(j) (odd j
from (j//2, j//2+1); j = 0 mod 4 by squaring; j = 2 mod 4 from
(j/2-1, j/2+1); x^2 is a full sqrhigh of x with one always-zero top
limb), block bases stay aligned only for even m, and since the window
grows a limb every other term, an addmul whose carry limb is already
occupied adds into it (at most two carries < 2*20! share a limb, so no
overflow).  The RS sum occupies xn limbs, division length xn.  Same
error bounds (8 ulp, 11 for m = 2, 3 for N <= 3); all 151 generated
variants verified against arb (max observed 2.43 ulp, one-sided).

Both families cover N = 1..21 (N <= 3 without rectangular splitting;
the 64-bit N = 4 winner is m = 2 with mulhigh division).  The 64-bit
dispatch tier in `_arb_exp_taylor_rs2` now takes any 1 <= xn <= 21.

For zbits = 32, N = 3, 4, 5 the generator emits hand-optimized
functions using inline limb operations (umul_ppmm/add_ssaaaa) instead
of mpn calls.  Tricks: for xn = 2 the top limb x[1] < 2^32, so x[1]^2
is a plain single multiply with no high word and x^2*2^128 fits one
limb; floor(v/6) or floor(v/24) of a two-limb v with small high limb
v1 is computed branch-free and one-sided as
(v1 + ch)*floor(2^64/D) + floor(cl/D) with ch:cl = v0 + (2^64 mod D)*v1,
short by at most 1; for N = 5, x^3/6 + x^4/24 = (4x^3 + x^4)/24 costs
one division; and the x^3 assembly folds partial products with carried
128-bit adds because floor((x2*b0 + x1*b1)/2^64) can exceed 2^64 (a
single-limb intermediate silently wrapped at maximal x -- caught by the
boundary stress test).  Measured: N=3: 1.8 ns, N=4: 2.7 ns (4.5x and
6.3x over their generated counterparts; 7.4x and 22.9x over
`_arb_exp_taylor_rs`), N=5: 11.5 ns (2.5x over generated, 7.7x over
classical).  Bounds 3/4/5 ulp, verified one-sided against arb
including x -> 2^-32.

Profiled winners for the generated remainder: m = 2 for N = 6, m = 4-6
above (middle band within noise), divrem_1 for N >= 17 (division
length xn >= 9).  Speedups vs `_arb_exp_taylor_rs` at xn = ceil(N/2):
2.4-3.2x for N >= 6.

## General limb-windowed evaluator (exp_series_gen.c)

`mpn_exp_series_gen(y, error, x, xn, N)` is the production path for
larger precision with x < 2^(-64*zL), zL >= 1 counted as whole leading
zero limbs of x (the limb-granular counterpart of _arb_exp_taylor_rs2's
exp path).  Structure: direct-write powers into top-aligned slots (all
shift deltas vanish), rectangular splitting with even m ~ sqrt(N/2)
(power computations cost full-length mulhighs while boundary
multiplications shrink with the window, so the optimum sits below
sqrt(N)), windows of wn = xn - (zL*base - GUARD_LIMBS) limbs with
additive per-term read lengths, one mulhigh per block boundary landing
exactly at the new window bottom, and mpn_divrem_1 (in place) for all
divisions.

The window guard is dynamic: guard = 1 when N <= 21 and guard = 2
otherwise.  With limb-exact zL, every truncation event committed at a
shrunk window (term reads, boundary mulhighs, divisions) is later
scaled by the pending x^base <= 2^(-64*zL*base) while the window
bottom sits zL*base - guard limbs down, so each event is suppressed by
2^(-64*guard) independently of the block divisors; N <= 21 means a
single coefficient block (divisor (N-1)!) and only the final division,
giving the guard-1 budget: full-window term drops sum 1/k! <= 1.72,
last boundary <= 3/(N-1)! <= 0.13, final division <= 1, reported error
4 ulp (verified over 4000 dedicated instances, max observed 1.67).
Guard 1 saves 4-7% at zL <= 2 (noise-level at zL >= 3).  Note the
boundary reads extend up to guard limbs below x^m's stored region, so
exactly guard pad limbs are zeroed.

Coefficients are generated on the fly, removing the FACTORIAL_TAB_SIZE
limit: descending k, the running numerator c = c_k and block divisor
dsf = max(k,1)*(k+1)*...*k_top are maintained in ulongs; when
dsf*max(k-1,1) would overflow, the block is closed by one window
division by dsf and restarted with c = 1.  For N <= 21 this reproduces
a single block with divisor (N-1)!.  c <= dsf < 2^64 always, so every
addmul coefficient fits a limb.  GUARD_LIMBS = 2 is kept (unlike the
hardcoded family, whose uniform 20! denominator damps window
truncations by 2^-61, the on-the-fly final block divisor can be small,
so guard-free windows would be unsound); error accounting matches
_arb_exp_taylor_rs2 (base 5, +2 per full-window division, +3 per
full-window boundary, +2 final).

Verified against arb for xn up to 339 limbs, zL in 1..3, N beyond 288
(dozens of on-the-fly blocks), max observed 1.68 ulp, one-sided.
Measured with exactly zL leading zero limbs (benchmark inputs pin the
top significant bit; n_randtest otherwise sprinkles zero limbs that
contaminate min-of-passes timings): 2.1-3.0x over
`_arb_exp_taylor_rs` (which aborts beyond N = 288) and 1.04-1.12x
faster than `_arb_exp_taylor_rs2` in nearly all configurations while
being far simpler.  At xn = N = 20, zL = 1 it remains ~1.4x slower
than `mpn_exp_series_20` even with guard 1 (the hardcoded family's
remaining edge: guard-free windows via the uniform 20! denominator,
exact 1/x/x^2-half terms, no loop or setup overhead), so the hardcoded
functions are still worth keeping.  Intended dispatch: x < 2^-32
few-limb -> exp_series32 tab; x < 2^-64, 5 <= xn <= 21, N >= xn ->
exp_series tab; otherwise -> mpn_exp_series_gen (which now also serves
x < 2^-128 with N <= 21 at its guard-1 speed).

## Dispatch heuristic

```
N small                          -> hardcoded
xn < 12 && N*zb < 560            -> mid   (checked both before and after
                                           the term clamp)
zb >= 64 && xn >= 8              -> L
otherwise                        -> general
```

plus fallbacks to the classical functions for 32-bit limbs, `zb < 4`
(outside the documented `x <= 1/16` domain) and degenerate corners.
`DISPATCH_XN`, `DISPATCH_NZB`, `DISPATCH_L_XN`, `DIV_PREINV_CUTOFF` and
`GUARD_LIMBS` are `#define`s meant to be tuned per architecture; building
with `-DNO_L_VARIANT` disables the L tier for A/B comparison.

## Benchmarks

Single-core container, gcc 13.3 `-O2 -march=native`, FLINT 3.3.1;
numbers indicative (noisy environment), each cell uses the term count N a
real caller would use at precision `64*xn` for an argument with `zb`
leading zero bits. Speedup = classical time / new time:

```
exp                     atan (zb = zbx)          sin+cos
zb\xn  4    8    16   32    4    8    16   32    4    8    16   32
 4    1.0x 1.0x 1.3x 1.3x  1.0x 1.0x 1.5x 1.9x  1.0x 1.0x 1.3x 1.4x
 8    1.0x 0.9x 1.3x 1.5x  1.0x 1.0x 1.5x 1.8x  1.4x 1.0x 1.3x 1.5x
16    1.0x 1.0x 1.5x 1.7x  1.0x 1.2x 1.3x 1.8x  1.0x 1.0x 1.5x 1.8x
32    1.0x 1.0x 1.6x 1.9x  0.9x 1.2x 1.8x 2.1x  1.0x 1.0x 1.7x 1.9x
64    2.2x 1.0x 1.8x 2.3x  0.9x 1.3x 1.8x 2.2x  0.9x 1.1x 1.9x 2.0x
128   1.0x 1.6x 2.0x 2.4x  1.0x 1.2x 2.6x 2.5x  1.0x 1.1x 1.8x 2.5x
```

Cells at small `xn` are the mid variant (~1.0-1.1x, occasionally more
where the clamp/preinv engage); the sub-1.0 entries are 20-90 ns calls
within the noise band of a few ns of dispatch cost.  Gains grow with both
precision (mulhigh advantage, quadratic work in shrunken windows) and
argument smallness (window shrinkage): 1.2–2x typical at `xn >= 16`, up
to ~2.4–2.6x for very small arguments at high precision.  Small
hardcoded cases (e.g. exp `N = 4`) gain 1.5–2x.

## Files

- `taylor_rs2.h`, `taylor_rs2.c` — implementation (self-contained except
  for linking against libflint for the coefficient tables and fallbacks)
- `test_taylor_rs2.c` — randomized rigorous correctness tests
  (`./test_taylor_rs2 [iters]`)
- `bench_taylor_rs2.c` — old-vs-new grid benchmark
- `Makefile` — set `FLINT_DIR` to a built FLINT 3.3.x source tree

Building `taylor_rs2.c` with `-DPHASE_TIMING` exposes cycle counters
(`taylor_rs2_phase[]`) for the setup/main-loop/final phases of the exp
evaluator, used during tuning.

## Notes and possible follow-ups

- The block length is the classical `m = smallest even m with m^2 >= N`;
  since window shrinkage makes late blocks cheaper than early ones, a
  slightly larger m (fewer boundary multiplications, which now dominate
  at high precision) may pay off; not explored.
- The inverse tables could be emitted as static const data by a
  generator, removing the lazy init.
- `_arb_exp_taylor_bound`-style callers could additionally exploit that
  the reported error is usually much smaller than the worst case by
  taking the returned dynamic bound instead of a static constant.


## atan/atanh and sin/cos/sinh/cosh series evaluators

Companions of the exp families for the odd/even series after argument
reduction, generated by gen_atan_sincos.py: `mpn_atan_series[32]_{xn}`
and `mpn_atanh_series[32]_{xn}` compute sum (-+1)^k x^(2k+1)/(2k+1)
(xn = 1..35 resp. 1..17); `mpn_sin_cos_series[32]_{xn}` and
`mpn_sinh_cosh_series[32]_{xn}` compute the odd/even factorial sums
(xn = 1..22 resp. 1..10), with either output pointer NULL-able: the
series summations are grouped after the shared z-power block and each
is skipped entirely when its output is not requested.

All four substitute z = x^2 (power z^k carries exactly 2*zbits/64
leading zero limbs) and run descending plain mpn_addmul_1 /
mpn_submul_1 rectangular splitting with even m.  For the alternating
functions the accumulator may wrap negative between the dominant
units terms; the fresh gap and carry limbs then use a running
sign-extension limb updated branch-free from each op's borrow/carry
(RS_SUBMUL / RS_ADDMUL / RS_MUL_NEG macros; the sign cannot be
tracked statically since a subtracted term can be exactly zero).
All-positive stretches and the non-alternating functions emit pure
exp-style chains with zero overhead.  Term counts are per series (odd:
2k+1 < xn; even: 2k < xn at zbits = 64), the even k <= 1 head is
exact (ycos = 1 - z/2 by shift+negate, ycosh = 1 + z/2, splitting
over k >= 2 only), and the non-alternating series exclude k = 0 from
the splitting (their sums exceed 1), finishing with an unconditional
res = x + mulhigh(q, x).  mpn_sin_cos_series_3 is exactly a copy, a
1x1 sqrhigh, a shift and a negation.

The single division by the shared denominator (3*5*...*33 ~ 2^62.5
for atan/atanh, 20! for the factorial families, whose damping keeps
guard-free windows sound) is either mpn_divrem_1 or a mulhigh by a
truncated inverse, swept per size: the top n limbs of
floor(2^(64*NMAX)/D) equal floor(2^(64*n)/D) exactly, so one constant
array serves every length with a one-sided error <= 3 units.  The
truncated inverse strictly underestimates, so the quotient never
reaches 2^(64*xn); the divrem variant subtracts 1 from the
accumulator first (error < 2^-61 quotient units) for the same
guarantee -- all output assembly is unconditional and branch-free.
The inverse wins at small sizes (also curing the previous small-size
regressions vs classical: atan32 xn=2 went from 0.96x to 1.50x),
divrem at large.

The final assemblies are minimal: the quotient's zero limbs are
never multiplied (the coefficient damping makes its top buffer limb
provably zero), and the closing multiply by x runs over only the
product's significant limbs -- n' = occ + xsig - xn limbs, i.e. xn-3
at zbits = 64 for the positive odd series -- with x's upper limbs
folded in by a single unconditional mpn_add_1 whose carry-out becomes
the top output limb (x*(sum) may cross 2^-64 when x is within an ulp
of it, a case exercised by the all-ones + garbage-top test patterns).
The alternating finals shorten by one limb since x*S <= x < 2^-64.
As a result no zbits = 64 function ever reads x's top limb, and the
contract only requires the value itself to be below 2^-64.

Accuracy: max observed 3.03 ulp (64-bit) and 2.10 ulp (32-bit) over
all 884 generated functions against arb (bounds 12/10 two-sided for
the alternating, one-sided for the positive families), plus NULL-call
consistency checks and a 2M-input fuzz of the sparse-x corner.
Speed at equal precision: atan 2.5-3.2x over `_arb_atan_taylor_rs`
(xn=18: 420 vs 1327 ns; xn=27: 1025 vs 3299 ns), atanh 2.8-3.2x
(xn=10: 114 vs 368 ns), sin/cos 2.1-3.3x and sinh/cosh up to 3.5x
over `_arb_sin_cos_taylor_rs` (xn=17: 490 vs 1734 ns), 2.1-3.6x for
the 32-bit families, and 12-23x at the tiny sizes that reduce to
copies.

`mpn_atan_series_gen` and `mpn_sin_cos_series_gen` (atan_sincos_gen.c)
are the general limb-granular versions, mirroring the classical
signatures (alternating/sinonly flags) with divrem_1 divisions, dynamic
guard (1 when a single denominator block covers the run, else 2), and
on-the-fly coefficients: sin/cos numerators chain integrally walking k
downward (c_{k-1} = dsf_k, dsf_{k-1} = dsf_k * f(k-1)(f(k-1)-1),
closing a block by dividing the window on overflow), while odd
reciprocals do not chain downward, so atan runs a cheap ascending
prepass collecting blocks of consecutive odds into a temporary array
and evaluates descending with one exact ulong division D/(2k+1) per
term plus the classical multiply/divide rescale at block changes.  The
sign bookkeeping follows the L-variants: even m keeps every boundary
multiplication positive, and block divisions in the negative phase
offset the units limb by the old denominator before dividing and
correct by the new denominator (resp. 1) after.  Measured max error
~1 ulp against exact partial sums (bounds 9-15); speed 1.4-3.6x over
the classical evaluators for xn = 8..256, zL = 1..4.


## fixed module (FLINT integration)

`fixed/` packages the evaluators behind the FLINT-facing interface in
fixed.h: `fixed_exp_rs_pre32/64` (res, n+1), the sin/cos families
`fixed_{sin,cos,sin_cos,sinh,cosh,sinh_cosh}_rs_pre32/64` (outputs
n+1 with a units limb; the combined versions accept NULL for either
output), and `fixed_{atan,atanh}_rs_pre32/64` (res, n).  Errors are
bounded by the public FIXED_*_RS_MAX_ERR(n) macros -- flat constants
(10 for exp, 15 for the others), justified because the windowed
general routines report O(1) error for these series: their
full-window increments never fire, since the first coefficient block
closes (exp: k = 20; atan: k = 16; factorial families: k ~ 10) after
the windows have already shrunk.

All hardcoded routines are private statics in exp_rs_hard.inc /
trig_rs_hard.inc (renamed _fixed_*_rs[32]_*, ulong/nn_ptr types,
unused error tables removed), included on 64-bit machines only; the
windowed general routines are inlined as statics into exp_rs.c and
trig_rs.c, with the sin/cos general routine adapted to n+1-limb
outputs carrying a real units limb (replacing the classical 1-ulp
cos clamp, which would be wrong for cosh).  Dispatch: pre32 uses
hardcoded routines for n <= 11 (exp; the 32-bit exp table is indexed
by term count, N = 2n, with N = 21 serving n = 11), n <= 17
(atan/atanh), n <= 10 (sin/cos families), else the fallbacks; pre64
routes x < 2^-128 or out-of-table n to the general routines, else
hardcoded (exp n <= 21, atan n <= 35, sin/cos n <= 22).  Cos-only
calls on the pre64 general path compute a discarded sin (the general
routine has no cos-only mode; possible future optimization).

The fallbacks _fixed_exp_rs_fallback, _fixed_sin_cos_rs_fallback and
_fixed_atan_rs_fallback run at constant full precision with
coefficients generated on the fly (odd reciprocals via an ascending
block prepass), determine the term count from the actual leading zero
bits of x, and are FLINT_BITS-generic: 32-bit machines route
everything through them.  Tests (FLINT test-harness style,
test/t-exp_rs.c, test/t-trig_rs.c, test/main.c) stress every public
function against its fallback at one extra limb of precision, check
single-output vs combined consistency, and include an absolute arb
comparison for small n (with a comment noting it becomes circular
once arb is rebuilt on this module, to be replaced by mpfr); pass at
10x the default iteration counts with n up to 250.


## fixed_exp_bitwise_rs

exp((x, n)) for any 0 <= x < 1 via the bitwise argument reduction of
the paper: subtract in turn each L_i = log(1 + 2^-i), i = 0..r, for
which L_i <= x (the classical single-pass if-reduction, sound on
[0, 1) since L_{i-1} < 2 L_i keeps the remainder below L_i after
step i, base case 1 < 2 log 2, each factor used at most once),
evaluate the Taylor series below 2^-r via fixed_exp_rs_pre32/64, and
multiply back the used (1 + 2^-i) factors by shift-and-add.  The
stored logarithms are truncations, letting the remainder creep at
most one guard ulp per step past the exact invariant; one extra
conditional subtraction of L_r after the loop restores a remainder
strictly below L_r < 2^-r.  Negative arguments are deliberately not
handled here (a second -log(1 - 2^-i) table would double the cache
footprint); they will later be folded into the large-argument
reduction by adding a multiple of log 2.  The reconstruction
processes the increasing used indices in chunks sharing the same limb
offset, with the i < FLINT_BITS bulk running branch-free full-length
shift-adds (no carry can leave the units limb since the value stays
below e).  One guard limb runs through reduction, series and
reconstruction: FIXED_EXP_BITWISE_RS_MAX_ERR(n) = 3 ulp, with r
clamped to FLINT_BITS n.  The single log table is arb-generated on
demand and cached per thread at the largest precision/index range
requested (top-limb reads serve smaller requests exactly); a cleanup
function is registered with FLINT.  Timings vs arb_exp at equal
precision (best r): n=8: 1.9x, n=16: 2.3x, n=32: 3.3x, n=64: 2.5x;
at n <= 4 arb's lookup-table exp remains faster.


### sinh reconstruction for long series

When the reduced-argument series gets long, fixed_exp_bitwise_rs
evaluates sinh(t) instead (the odd series has half the terms, over
z = t^2 whose powers also carry twice the leading zeros) and
reconstructs exp(t) = sinh(t) + sqrt(1 + sinh(t)^2), costing one
truncated squaring and one mpn_sqrtrem on the (2 wn + 1)-limb integer
(1 + sinh^2) 2^(2 FLINT_BITS wn).  All the extra error (sinh's 15 ulp
at the guarded precision, the squaring truncation, the square-root
floor) sits below one guard ulp, so FIXED_EXP_BITWISE_RS_MAX_ERR
stays 3.  Measured crossovers: for r < 64 (series through the pre32
path) sinh wins from ~45 terms (n ~ 24 at r = 32, saving 15-20% by
n = 64+); for r >= 64 the windowed pre64 series is more efficient per
term and the square root sets a higher floor, moving the crossover to
~128 terms (n ~ 2r).  Encoded as EXP_USE_SINH(wn, r) with those two
constants; margins near the crossovers are within a few percent.
Future improvement: for t < 2^-r the square root argument satisfies
sinh^2 < 2^-2r, so sqrt(1 + w) could be evaluated by a short binomial
series (a quarter of the exp terms) instead of mpn_sqrtrem, which
would pull the r >= 64 crossover substantially earlier.

End-to-end at the best r per size, vs arb_exp at equal precision:
n=16: 1404 vs 3925 ns (2.8x), n=32: 3359 vs 12537 (3.7x), n=64: 13.6
vs 31.6 us (2.3x), n=128: 48.4 vs 98.6 us (2.0x), n=256: 186 vs 300
us (1.6x).


### register small path (n <= 5)

The argument reduction and reconstruction loops of
fixed_exp_bitwise_rs run entirely on register operands for n <= 5,
via per-size functions generated by gen_exp_bitwise_small.py into
exp_bitwise_rs_small.inc (named scalar limbs, the fixed-size
add/sub chains from longlong.h, MPN_RIGHT_SHIFT_LOW funnel shifts
from mpn_extras.h, and an advancing table pointer).  The reduction
compares and subtracts branch-free: a (wn+1)-limb subtraction with
zero top limbs yields the borrow as a full mask selecting between
difference and old value and advancing the used count; the
reconstruction factors with i < FLINT_BITS (all of them at the
typical r) are funnel shifts plus one carry chain, with rare larger
i spilling to the generic chunk loop.  Two failed designs are worth
recording: a single top-aligned fixed-width (6/7 limb) body was
slower than the mpn loops at n <= 3 (3x the raw width at n = 1, and
the branchless step cannot exploit mpn_cmp's top-limb early exit),
and a WN-templated inline body was sabotaged by the compiler
(auto-vectorization of the select loop forcing GPR<->XMM traffic,
and incomplete scalar replacement leaving the difference limbs in
stack slots).  The generated named-scalar form is deterministic:
1.9x/1.5x/1.3x/1.4x/1.0x over the mpn loops at n = 1..5 (r = 32) and
1.1-2.0x at r = 64, verified bit-identical to the generic path
(FIXED_EXP_BITWISE_NO_SMALL disables the dispatch for A/B testing).
With it, fixed_exp_bitwise_rs beats arb_exp at every precision:
n=1: 108 vs 148 ns, n=2: 150 vs 222, n=3: 197 vs 287, n=4: 257 vs
378, n=5: 306 vs 540, on through the large-n numbers above.


### windowed single-limb reduction

The reduction loop decides its subtractions one FLINT_BITS window at
a time on a single-limb model h of the remainder.  h is maintained as
a lower bound on the remainder's current window limb (the model
subtracts lt + 1, since each table limb is a truncation), with upper
slack e growing by one per model subtraction; then h > lt proves
L_i <= t and h <= lt - e - 1 proves t < L_i, so decisions outside the
ambiguity band [lt - e, lt] are exact and branch-free at one limb per
step.  Steps inside the band -- a ~2^-57 fraction of random steps,
plus automatically the last few steps of each window where lt itself
is small, plus the first step of each window (the remainder may still
carry one bit above it) -- fall back to the exact full-width
compare-subtract with a model resync.  Decided subtractions are
recorded and applied to the full-precision remainder in unconditional
batches (branch-free register carry chains in the small path).  The
decisions agree exactly with the plain greedy, verified bit-identical
over 40000 random (n, r).

Results match the prediction that motivated the scheme: the generic
path is neutral (mpn_cmp's top-limb early exit was already cheap;
kept for the ~7% gain at r = 128), the register path loses ~10% at
n = 1 where the windowing overhead cannot amortize over a 3-limb
chain, and wins 13-22% at n = 3..5 (n=5, r=64: 398 vs 490 ns; r=320:
2256 vs 2706).  The generator therefore emits the direct per-step
reduction for n <= 2 and the windowed one for n >= 3, making the
hybrid at least as fast as the previous code at every (n, r).


### accuracy vs r, and removal of the internal guard limb

Measured accuracy of the padded implementation (max ulp over 400
samples per cell, arb reference): 1.00 flat for every n in 1..64 and
every r in 32..320 -- the r-dependent error sources (table and shift
truncations, ~3r of them) all sit at ~3r 2^-FLINT_BITS output ulp
behind the guard limb, so guard limbs would only need to grow around
r ~ 2^60.  Since callers of this function will typically have padded
the precision already, the internal guard limb was double padding:
fixed_exp_bitwise_rs now works entirely at the output precision of n
limbs, with the r-linear bound
FIXED_EXP_BITWISE_RS_MAX_ERR(n, r) = 9 r + 100 (derived: at most
r + 2 subtracted table entries each short by < 2 ulp and r + 2
reconstruction shifts each truncating < 1 ulp, amplified by the
remaining product < e, plus series/sinh terms; measured maxima
0.5-0.7 r, a 12-20x margin).  r is now clamped to FLINT_BITS n - 16,
which keeps the reduction slop below L_r so the single extra
subtraction still suffices.  The narrower working width speeds up
every size, most visibly the register path: n=1: 73 vs 102 ns
(r=32) and 117 vs 216 (r=64), n=2: 109 vs 165, with ~2-8% at large n.

The binomial-series square root idea is withdrawn: sqrt(1 + w) with
w = sinh^2 < 2^-2r needs terms until 2rk > FLINT_BITS n, i.e. ~N/2 of
them -- the same count as the cosh series it secretly is -- so it is
equivalent to evaluating a second series as long as the first and can
never beat a 2-3 M(n) mpn_sqrtrem in the regime where the sinh path
runs at all.


### unified public interface

The pre32/pre64 split is gone from the public interface: the series
evaluation functions are now fixed_exp_rs, fixed_{sin,cos,sin_cos,
sinh,cosh,sinh_cosh}_rs and fixed_{atan,atanh}_rs, documented as
requiring x < 2^-32 and checking it with
FLINT_ASSERT((x[n-1] >> (FLINT_BITS - 32)) == 0) (which degenerates
to x[n-1] == 0 on 32-bit machines).  Each dispatches on the top limb:
nonzero selects the 32-bit internal versions, zero the 64-bit
windowed/hardcoded ones -- so a caller who happens to pass an
argument below 2^-64 through the "32-bit" contract now automatically
gets the faster path.  A pleasant consequence: fixed_exp_bitwise_rs
and its generated small functions no longer branch on r to pick a
series variant; they just call fixed_exp_rs / fixed_sinh_rs and the
top limb of the reduced argument makes the same choice implicitly.
Tests reshaped accordingly (the zbits parameter now only shapes the
random inputs so both internal paths are exercised); everything
passes at 10x iterations with timings unchanged.


## FLINT patch

`0001-Add-fixed-module-efficient-low-level-fixed-point-rea.patch` is
a complete git patch on top of the FLINT tree adding the module: 19
files, ~15200 lines.  Contents: src/fixed.h; src/fixed/{exp_rs.c,
trig_rs.c, exp_bitwise_rs.c} with the generated exp_rs_hard.inc,
trig_rs_hard.inc, exp_bitwise_rs_small.inc and  src/fixed/test/{main.c, t-exp_rs.c,
t-trig_rs.c, t-exp_bitwise_rs.c}; src/fixed/profile/
p-exp_bitwise_rs.c (the n-vs-r comparison table against arb_exp,
using profiler.h timing); doc/source/fixed.rst (module conventions,
error macros, all public functions, r tuning guidance) listed in
doc/source/index.rst; and the module registered in Makefile.in
(HEADER_DIRS) and CMakeLists.txt (_BUILD_DIRS).  Validated in-tree:
make compiles the three module objects warning-free under FLINT's
flags (-std=c11 -pedantic -Wall -Werror=implicit-function-declaration
-O3), make -n confirms the test and profile programs are wired into
the tests and profile targets, the test suite passes at 10x
iterations built against the tree, and tests and profile also compile
pedantic-clean.


### code generators in dev/

The patch includes self-contained regeneration scripts:
dev/gen_fixed_exp_rs_hard.py, dev/gen_fixed_trig_rs_hard.py and
dev/gen_fixed_exp_bitwise_small.py.  The first two embed the original
rectangular-splitting generators plus the benchmark-selected
(m, use_divrem) winner tables (formerly external JSON files) and the
mpn->fixed transform (renames to _fixed_*_rs[32]_*, static
qualifiers, ulong/nn_ptr types, error tables dropped), so that
`python3 dev/gen_fixed_X.py > src/fixed/X.inc` reproduces each
shipped file byte-for-byte -- verified by diff for all three.


### patch rebased on FLINT main context

The patch now applies cleanly to FLINT main: the pristine base for
Makefile.in, CMakeLists.txt, doc/source/index.rst and
doc/source/references.rst was replaced by the current main versions
(fetched from GitHub), with `fixed` registered in main's tab-column
module lists; verified with git apply --check --whitespace=error on a
simulated main tree, plus an actual apply.  Citations corrected
throughout (commit message, fixed.rst, source comments): the
algorithms follow "Fast multiple precision exp(x) with
precomputations" (van der Hoeven and Johansson, ARITH 2024) and
"Efficient implementation of elementary functions in the
medium-precision range" (Johansson, ARITH 2015); references.rst gains
the [HJ2024] entry and [Joh2014c] is updated from preprint-only to
its ARITH 22 publication data.  The module intro in fixed.rst now
reads: low-level, low-overhead functions for fixed-point real
arithmetic, intended as kernels for arbitrary-precision
floating-point arithmetic and ball arithmetic.  Trailing-blank-line
whitespace errors were removed and the dev/ generators updated so all
three .inc files remain byte-for-byte reproducible.


### n = 6, 7 cliff investigation

Measured: at n = 6 the full function and a generic-only build run at
the same speed (375 vs 361 ns, r = 32), so the doubling from n = 5 was
entirely the end of the hardcoded register path, not compiler
inlining -- beyond it, the mpn loops are dominated by per-call
overhead at these widths, which is also why n = 6 and n = 7 were so
close.  Since longlong.h provides carry chains up to 8 limbs, the
generated small path now extends to n <= 7: n = 6 improves from 564
to 360 ns and n = 7 from 719 to 453 (r = 64; 1.25-1.3x at r = 32),
with n <= 5 unchanged and outputs bit-identical over 30000 random
(n, r).  n = 8 is the natural remaining boundary (no 9-limb chains,
and register pressure would force spills anyway).

The generator now also emits the near-identical per-window reduction
blocks through a single per-function WINDOW(cc, hreg) macro (the
model limb's source register is a macro parameter, hence still
compile-time), invoked once per window: the expansion is unchanged,
so the machine code and timings are identical, while the generated
file stays at ~1000 lines despite growing from five functions to
seven.  The stale sinh-path error comment (written for the earlier
internally-padded version, where the 15 ulp sinh error lived at the
guarded (n+1)-limb grid and hence below one output ulp) was rewritten
for the no-padding code: those ~20 ulp now land directly in the
output and form part of the +100 constant in
FIXED_EXP_BITWISE_RS_MAX_ERR.


## log via the dual (L-mode BKM) reduction: worked out

log_proto.c (prototype, not yet in the patch) validates the design:
greedily multiply P (from 1) by factors 1 + 2^-i, i = 1..r, while
P (1 + 2^-i) <= X = 1 + x.  This greedy IS the exp reduction applied
to log(X/P), so everything ports verbatim: each factor used at most
once (L_{i-1} < 2 L_i), one extra conditional step at i = r absorbs
truncation creep, the SAME cached L table serves (no new
precomputation, no cache growth), and the windowed single-limb model
applies (decisions compare the deficit D = X - P's window limb
against (P >> i)'s window limb, computable in O(1) from P's top two
limbs since P is in [1, 2)).  Each accepted factor is one
shift-and-add on P, one exact subtraction from D, and one table add
into the L-accumulator.

The fusion: log(X/P) = 2 atanh((X - P)/(X + P)).  The numerator is
the exact deficit D, so ONE division t = D/S with S = X + P in
[2, 4) simultaneously normalizes by P and produces the atanh
argument; D < P 2^-r and S >= 2 give t < 2^-r, meeting the
fixed_atanh_rs contract out of the box.  Two structural wins fall
out: no per-factor division and no separate division by P at all;
and the odd atanh series (t ~ u/2) needs HALF the terms of
log(1 + u) -- log's analogue of exp's sinh trick, with no square
root.  Errors are O(r) exactly as in exp_bitwise: table reads < 2
ulp per entry, P truncations < 1 ulp per accept entering as
log(Pi/P) < r ulp, division floor, doubled atanh error.

Prototype measurements: max error ~0.2 r ulp across n in 1..64 and
r in 32..256 (a (4r + 64)-style macro would hold with 5x margin);
timings vs arb_log at equal precision: n=16: 2408 vs 4068 ns
(1.7x), n=32: 7391 vs 15818 (2.1x), n=64: 20722 vs 50967 (2.5x),
n=128: 69624 vs 149470 (2.1x), losing only below n ~ 8 where the
prototype's naive per-step mpn calls dominate (a rejected step still
pays zero + rshift; an accept pays three mpn ops; the division is a
generic mpn_tdiv_qr).  The remedies are exactly the exp playbook:
register small path, windowed single-limb decisions, short division
by precomputed inverse.  Domain [1, 2) (log1p on [0, 1)) covers the
mantissa case; the (1 - 2^-i) direction is deliberately omitted (it
would double the table cache, as for exp's dropped xneg), with X < 1
handled by the caller through the exponent.


## fixed_log1p_bitwise_rs (production)

log1p_bitwise_rs.c lands the dual reduction in the module, mirroring
exp: same shared L-table cache (now exposed internally through
fixed.h), same r-linear error philosophy
(FIXED_LOG1P_BITWISE_RS_MAX_ERR(n, r) = 3r + 64, measured maxima
~0.2r), FLINT-style test (t-log1p_bitwise_rs.c) and profile program
(p-log1p_bitwise_rs.c), docs in fixed.rst, commit message updated.
Beyond the prototype: (1) the reduction decisions run on a
single-limb model of the deficit D = X - P, with
lt = MPN_RIGHT_SHIFT_LOW(P[n], P[n-1], i mod 64) computed from P's
top two limbs (window-independent since P < 2); because accepts
update D immediately, the model resyncs exactly after every accept
and the ambiguity band degenerates to the single value h == lt, so
rejected steps cost a handful of register operations; (2) one
correctness subtlety absent from exp: the deficit obeys only
D < P 2^-i < 2^(1-i), so a bit above the current window can persist
-- but then P >> i <= D is a certain accept, handled by one extra
limb check per step; (3) P grows gradually as suggested: each
accepted factor 1 + 2^-i moves P's exact lowest bit down by i, so
the significant length starts at one limb and all accept-side
operations (shift, D-subtract, P-add, table add) run on significant
slices with carry/borrow ripple, P being truncated only once its
exact length would exceed n limbs.  The production version matches
the prototype bit-exactly over 3000 random (n, r) and runs 1.4-2.8x
faster than it; vs arb_log at equal precision: n=8: 737 vs 1121 ns
(1.5x), n=16: 1475 vs 4025 (2.7x), n=32: 4293 vs 12705 (3.0x),
n=64: 14800 vs 49653 (3.4x), n=128: 49113 vs 136899 (2.8x), losing
only below n ~ 6 (arb's lookup-table log; a register small path as
for exp is the known remedy).  Interestingly log's cost now roughly
matches exp_bitwise at equal (n, r): the single division costs about
what exp's sinh square root does.


### r = 16 unlock and n <= 2 specializations for log

fixed_log1p_bitwise_rs now accepts r >= 16.  Values below 32 take
effect in new register implementations for n = 1 and n = 2 (the
general path clamps r up to 32, keeping the fixed_atanh_rs
contract): their inline atanh evaluations only need t < 2^-16, since
at these precisions the series collapses -- t + t^3/3 at n = 1 (t^5
is below the 64-bit ulp already for r = 16), four terms at n = 2.
The n = 1 reduction is a fully branch-free ~8-operation register
loop (funnel shift of P's top limbs, compare-select mask, masked
updates of deficit/product/accumulator, all carry-free by the
invariants P + v <= X < 2 and acc < log 2), with the final division
a hand-written normalized 3-by-2 limb step; n = 2 uses two-limb
register windows with sub_dddmmmsss borrow masks (one boundary bug
caught by the tests: the extra step at exactly r = 64 must use
v = (1 : p1), not a zero-shift funnel).  Measured: n = 1: 78 ns at
r = 16 (vs 115 at r = 32, 303 for the previous generic path, 213
for arb_log -- 2.7x over arb, 3.9x over generic); n = 2: 182 ns at
r = 16 (vs 232/510/306 -- 1.7x over arb).  r = 16 is decisively
optimal at both sizes, confirming that the reduction dominated
there.  Accuracy at r = 16: ~9-11 ulp, far under the 3r + 64 = 112
bound.  n = 3..4 remain on the generic path (~540/630 ns): extending
the register treatment to n <= 7 as for exp is the known next step.


### generated atanh16 series (zbits = 16 in the generator)

The hand-inlined n = 2 atanh series in log's small path (three
divisions, divrem_1 instead of inverse multiplication) is replaced by
generator output: dev/gen_fixed_trig_rs_hard.py now supports
zbits = 16, producing _fixed_atanh_rs16_{1..4} (in trig_rs_hard.inc)
with the full optimization set -- one denominator division by
truncated-inverse multiplication, squaring chains, rectangular
splitting with swept (m, dr) winners (dr = 1 everywhere, m = 2 except
m = 4 at xn = 4).  Getting there surfaced a real structural limit:
the generator assumes each power z^j has a whole number of leading
zero limbs, and at zbits = 16 the 32 j zero bits sit at half-limb
offsets -- most strikingly, computing z^2 from z's top limb alone
returns exactly zero (that limb is below 2^32), silently deleting the
t^5 term.  The shipped zbits = 16 mode therefore disables the
per-power shrink (full-width windows and powers), which in turn
required handling carry-slot collisions absent at 32/64: with g == w,
the units coefficient and subsequent addmul carries must add into,
not overwrite, the occupied top slot (two additive emission paths,
provably inert for 32/64 -- the regenerated 32/64 sections stay
byte-identical).  The alternating (RS_* macro) paths are excluded
from 16 pending the same treatment, so atan16 is deferred.  A proper
sub-limb extension (shift-normalized z) is noted as the way to
recover the remaining 1-limb-operation savings.  Validated: all 14
candidate configurations max 1.00 ulp against arb over 3000 samples
each; internal dispatch _fixed_atanh_rs16 (hardcoded n <= 4, fallback
beyond, fallback contract relaxed to x < 2^-16); tested directly in
t-trig_rs.c.  Effect on log: n = 2 improves from 182 to 158 ns at
r = 16 (n = 1 keeps its two-term inline: a call would cost more than
the four operations it replaces).


### inline emission for the atanh16 series and cleanups

Following review: the generated _fixed_atanh_rs16_* functions now
emit inline umul_ppmm/add_ssaaaa sequences instead of mpn calls for
all 1- and 2-limb operations (base and power squarings with the
2 x1 x0 cross term, mul_1/addmul_1 coefficient steps, the
truncated-inverse multiply at length 2, the final multiply by x),
plus a dedicated 5-multiplication sequence for the zero-padded z^m
boundary (effectively a 3x2 product; the dense 3x3 inverse multiply
remains a mulhigh_n kernel call).  The length-0 mpn_add_1 in the
final combine (undefined behavior per GMP's n >= 1 requirement,
though its carry was provably zero) is eliminated: when q x covers
all limbs the combine is an inline carry-free add, and the then-dead
cy and unused scratch declarations are filtered out.  The n = 2 code
of fixed_log1p_bitwise_rs likewise replaces its length-2 mpn_add_n
(and a compound-literal temporary) with one add_sssaaaaaa.  Emission
mishaps along the way reinforced two working rules: always assert
(or verify per-edit) that automated source replacements matched, and
never trust cross-pass benchmark comparisons on this container
(whole-run speed varies ~1.5x).  Within one pass: atanh16_2 improves
about 24% (31.5 to 24.0 ns), atanh16_3 about 15%, xn = 1 and 4
roughly neutral; correctness re-verified against arb for all 14
configurations (max 2.92 ulp, inline mulhigh error slightly above
the kernels', well under the 12 ulp budget); module warning- and
pedantic-clean, full suite at 10x.


### log register small path extended to n <= 7

dev/gen_fixed_log1p_bitwise_small.py generates
_fixed_log1p_bitwise_rs_{3..7} (log1p_bitwise_rs_small.inc, ~2500
lines), mirroring the exp small-path generator.  The design
consolidates what the hand-written n <= 2 versions taught: deficit,
product fraction (units implicitly 1: accepted factors keep
P + v <= X < 2, so the fraction never carries out) and accumulator
all live in named scalars; every step -- window boundary, model tie,
certain accept, and the extra step at i = r -- is ONE masked
operation, a borrow-capturing subtraction chain whose borrow doubles
as the exact comparison and, as an all-ones/zero mask, selects the
new deficit and gates the additions into P and the accumulator; the
fast path rejects on h < lt (one funnel shift of P's top limbs) plus
the stray-bit check for windows past the first.  The runtime extra
step is a switch on r / 64 with per-window straight-line bodies,
including the b == 0 boundary form (the r = 64 bug class from the
n = 2 work, now handled by construction).  One generator bug caught
in code review before testing: the P update initially used a wrong
index mapping and span (v occupies value limbs 0..top, mirroring the
deficit chain, not a slice starting at limb c).  r >= 16 now takes
effect through n = 4 (atanh16 coverage); n >= 5 clamps to r >= 32.

Validation: bit-exact against the generic path over 60000 random
(n, r) at r >= 32 -- the masked-step decisions provably reproduce
the generic greedy -- plus arb accuracy at r in [16, 32) (8-13 ulp
vs the 3r + 64 bound).  Timings (one pass, same regime): n = 3:
167 ns at r = 16 vs 468 generic and 300 arb_log; n = 4: 229 vs
460/364; n = 5: 313 vs 504/474; n = 6: 325 vs 472/598; n = 7: 473
vs 654/692.  fixed_log1p_bitwise_rs now beats arb_log at every
precision; r = 16 is optimal through n = 4, r = 32 for n = 5..7
(extending the atanh16 family past xn = 4 would unlock smaller r
there -- untested whether it pays).


### r < 32 at n = 5-7: measured, does not pay

The atanh16 family was extended to xn <= 7 (winners m6_d1 at 5 and 6,
m6_d1 at 7; ~123/147/220 ns; all max <= 12 ulp vs arb by the sweep's
correctness gate) and the log small path temporarily unlocked to
r >= 16 for all n.  Measurement (one pass):

  n |  r=16   r=20   r=24   r=28   r=32   r=48
  5 |   295    305    336    352    285    354
  6 |   320    360    373    417    325    408
  7 |   437    439    458    492    387    476

r = 16 loses to r = 32 at n = 5 (295 vs 285) and n = 7 (437 vs 387)
and ties at n = 6 -- and the whole [16, 32) band is dominated,
rising with r because the 2^-16 series has a fixed term count while
the reduction grows, until the cheap 32-range series takes over at
r = 32.  An r-adaptive series (an atanh24 family, or zero-count
dispatch inside atanh16) is what it would take to change this.  The
unlock is reverted (n >= 5 clamps r to 32, as documented); the
atanh16_{5..7} functions stay in the tree -- correct, tested through
_fixed_atanh_rs16 (now dispatching hardcoded up to n = 7), and
available to future adaptive-r or atan-reduction work.


### 32-bit correctness audit

Method: no 32-bit toolchain is available in the container, so the
audit combined (a) an inventory of every FLINT_BITS == 64 guard and
every include of generated code, (b) a semantic trace of the paths a
32-bit build takes, and (c) a compile-level shim check per module
source (#include flint.h, redefine FLINT_BITS to 32, #include the
.c with -Wall -fsyntax-only) that catches references to gated-out
symbols and unused statics -- all four now pass warning-free.

Findings and fixes: (1) neither small-path include was guarded --
exp_bitwise_rs_small.inc and log1p_bitwise_rs_small.inc (whose
generated window bounds use literal 64s) now sit under
FLINT_BITS == 64, along with log's hand-written n = 1, 2 functions
and all the dispatches; (2) at 32-bit, n = 1 clamps r to
FLINT_BITS n - 16 = 16, below the 2^-32 contract of the exp series:
fixed_exp_bitwise_rs now asserts n >= 2 on 32-bit (documented),
while fixed_log1p_bitwise_rs keeps n = 1 working through the generic
path by branching to _fixed_atanh_rs16 (whose fallback is
FLINT_BITS-generic with the zb >= 16 contract) when r < 32;
(3) the audit caught an unrelated regression: exp's small-path
dispatch had silently reverted to n <= 5 while the generated inc
serves n <= 7 (n = 6, 7 were falling through to the generic path) --
fixed and spot-checked (258/332 ns vs the n = 8 cliff at 457).
The generic paths themselves audited parametric: exp_rs/trig_rs
dispatch to fallbacks at 32-bit (zb asserts consistent: exp and
sin_cos at 32, atan/atanh at 16), the exp_bitwise and log1p generic
bodies use FLINT_BITS throughout, and EXP_USE_SINH's constants are
tuning-only.  Remaining honest caveat (documented in fixed.rst):
the 32-bit paths are compile-checked and semantically audited but
have not run on real 32-bit hardware.


### automated r tuning and r = 0 defaults

src/fixed/tune/tune-bitwise-r.c (picked up by FLINT's `make tune`
directory rules) finds the optimal-r crossovers for both bitwise
functions: candidate parameters form a ladder (16 for log, then 32,
64 and every further multiple of 64 up to a user-specified rmax),
and each crossover is located by binary search over n (up to nmax,
defaults 160/320) under the stated assumption that the optimal r is
nondecreasing in n.  Two robustness lessons are baked in: (1) the
input must be uniform random limbs -- flint_mpn_rrandom's structured
bit runs degenerate the greedy accept pattern and skewed every
crossover on the first attempt; (2) the comparator must be
deterministic at the degenerate ends, or noise cascades: where r_old
is internally clamped up (log's r = 16 beyond n = 4) the two calls
are identical and the tie must count as a win for r_new by fiat
(before this fix, the 16->32 'crossover' landed at n = 56 on pure
jitter, and monotone chaining pushed every later crossover after
it), and r_new counts as unavailable while r > FLINT_BITS n - 16.
A 1% tolerance absorbs the flat margins.  Output is a pasteable pair
of tables per function.

Both functions now accept r = 0, selecting the tuned default via the
built-in tables (largest tabulated value for larger n); measured
x86-64 defaults shipped: exp r = 32/64/128/192/256/320 from
n = 1/12/39/115/137/154, log r = 16/32/64/128 from n = 1/4/18/64.
The error macros treat r = 0 as a conservative 512.  Tests exercise
r = 0 every fourth iteration; docs updated; tuner pedantic-clean and
32-bit shim-clean (r_floor made FLINT_BITS-aware).


### tuner hardening

The tuning program's hardcoded 400-limb stack buffers (input in
tune_func, output in time_call) crashed for nmax beyond them (e.g.
nmax = 10000, rmax = 512), and the crossover tables were fixed at 32
entries.  All sizing is now dynamic (input nmax limbs, output
nmax + 1, tables rmax/FLINT_BITS + 3 entries, flint_malloc/free),
the output buffer is allocated once and passed down, and the
arguments are validated (2 <= nmax <= 32767 so crossovers fit the
short output tables, rmax >= 32) with a usage message.  Verified
with an AddressSanitizer build run at the formerly crashing
parameters -- probes at n = 10000 (640k-bit calls) complete with no
sanitizer reports -- plus a small-parameter regression run and
pedantic-clean compilation.


### test coverage

Method: the module and test sources compiled with gcc --coverage at
-O0, the suite run at FLINT_TEST_MULTIPLIER = 25, gcov -b per object
(the .inc files report as included sources).  First results exposed
two real issues: (1) exp_rs_hard.inc sat at 75.5% because the odd-N
_fixed_exp_rs32_* functions are DEAD CODE -- the dispatch computes
N = min(2n, 21), so odd N below the top entry is unreachable; the
generator now skips them (10 functions, 369 lines removed, NULL tab
entries); (2) log1p_bitwise_rs_small.inc sat at 83.8% because the
tests' r range capped at 321, never reaching the n = 6, 7 high
windows, and never hitting exact multiples of FLINT_BITS (the
extra-step b == 0 boundary forms -- including the historically buggy
r = 64 case).  The tests now start with a deterministic (n, r)
sweep -- all window counts, both boundary forms, smallest and
largest legal r per n, r = 0 -- and sample random r up to 448.

Final coverage: exp_rs_hard.inc, trig_rs_hard.inc,
log1p_bitwise_rs_small.inc 100.0%; exp_bitwise_rs_small.inc 99.7%
(one chunked-reconstruction call at a rare used-factor shape); the
.c files 97.2-98.4% with every uncovered line classified: sixteen
top-limb carry corrections in the series windows (each needs a
~2^-64 coincidence), the exp extra-step accept (fires only from
accumulated truncation creep -- required for worst-case soundness,
probabilistically unreachable under random inputs), two rare
division-correction breaks in log's hand-written small sizes, and
log's generic atanh16 branch (reachable only on 32-bit limbs).

### 32-bit execution: all tests pass

A 32-bit toolchain was installed (gcc-multilib plus i386 gmp/mpfr
from the Ubuntu archive) and FLINT 3.3.1 built at ABI=32 in a
separate tree.  The build needed several environment fixes, none of
them fixed-module issues: the tree's Makefile.in (main-derived from
the patch rebase) lists modules the 3.3.1 tarball lacks and omits
ones it needs (mpfr_vec/mpfr_mat, required by fmpz_lll); the ld -r
partial links need LD='ld -m elf_i386'; fft_small requires 64-bit
AVX intrinsics and must be disabled at ABI=32 (configure cache
overrides); and configure misdetects five GMP internal symbols
(__gmpn_add_nc, sub_nc, rsh1add_n, rsh1sub_n, modexact_1_odd) as
native although Ubuntu's i386 libgmp does not export them -- flipped
in the generated config headers.  One test fix belongs to the patch:
t-exp_bitwise_rs now draws n >= 2 on 32-bit limbs, matching the
documented contract.

Result: all four fixed tests PASS as genuine ELF 32-bit executables
at FLINT_TEST_MULTIPLIER = 10, and again with the module compiled
under FLINT_WANT_ASSERT (all contract asserts hold at runtime).
The 32-bit paths are no longer merely audit-verified.

(Process note: a recurring tool-call interruption traced to
pkill -f 'make -j2' matching the invoking shell's own command line
and killing it mid-script; pkill -x make is the safe form.)


### adversarial stress test

t-bitwise_rs_stress.c feeds the bitwise reductions inputs built from
their own internal structure: for exp, sparse sums of the tabulated
L_i = log(1 + 2^-i) read from the shared cache (with independent
random indices, so repeats occur, and every fifth iteration
explicitly duplicating the entry at the reduction limit); for log,
sparse products of the factors 1 + 2^-i assembled by exact
shift-and-adds under the X < 2 constraint (reverting any factor that
would cross 2 -- itself a mirror of the greedy).  Both are then
perturbed by up to +-3 ulp to probe each side of every decision
boundary.  These inputs force what random data essentially never
produces: exact greedy matches, single-limb-model ties (h == lt),
deficits driven through zero, and residuals landing at the extra
correction step.  Verified effective by coverage:
exp_bitwise_rs.c reaches 100.0% -- the extra-step accept, previously
classified as probabilistically unreachable (it fires only from
accumulated truncation creep), is now exercised by the duplicated-L_r
construction.  Remaining uncovered across the module: sixteen
~2^-64-probability carry corrections in the series windows, two
division-correction breaks, one rare reconstruction shape, and the
32-bit-only generic atanh16 branch.  The stress test passes at
multiplier 10 on 64-bit and multiplier 5 on the real 32-bit build.


## design sketch: sin/cos and atan mirroring exp/log

Shared table: ONE cached vector alpha_i = 2 atan(2^-i)
(thread-local like _fixed_exp_logs, arb-generated floors; the
top-limb-exactness argument needs restating for atan).  The scale is
dictated by the rotation direction, where the table drives the
reduction comparisons against theta through the windowed model; the
vectoring direction's decisions never touch the table (pure
Y vs X >> i geometry), so it merely ACCUMULATES the doubled entries
-- into an (n+1)-limb accumulator, since the sum reaches pi/2 > 1,
or batched from a used[] list -- and halves the sum once at the end
(one exact shift, <= 1 ulp).  One table, one exactness proof.

atan (the log mirror, vectoring mode): state (X, Y) = (1, t),
t in [0, 1), output atan(t) in [0, pi/4).  Greedy for i = 1..r:
apply (1 - i 2^-i), i.e. (X, Y) <- (X + Y 2^-i, Y - X 2^-i) and add
A_i, whenever Y >= X 2^-i (accumulating the doubled table entry, see
above).  This is a greedy on the angle with steps
A_i = alpha_i / 2; concavity of atan gives A_{i-1} < 2 A_i, so single use per
factor, one extra step at i = r, and a final angle below A_r ~ 2^-r
-- the exp/log convergence analysis verbatim.  The decision
Y >= X >> i is structurally log's D >= P >> i (Y = deficit,
X = product): the single-limb windowed model ports directly (h from
Y's window limb, lt one funnel of X's top two limbs; X in [1, K),
K < 1.21, unit limb stays 1), with exact resync after immediate
accepts and a single-valued tie band.  Fused residual: ONE division
t' = Y/X < 2^-r feeds fixed_atan_rs directly (r >= 32; r >= 16
needs the deferred atan16/RS-macro work).  Differences from log:
accepts cost 2 shifts + 2 add/subs (cross terms, old values), about
1.5x log's; and Y is not an exact deficit -- truncations perturb the
angle by <= ulp each, so the bound is two-sided but keeps the
C r + D shape.  Growing significance applies to both components;
register small paths follow the log small generator (X, Y, acc
scalars ~ 3n live, cap around n <= 5).

sin_cos (the exp mirror, rotation mode): theta in [0, 1), include
i = 0 (sum A_i ~ 1.68 covers the domain).  Phase 1: reduce theta by
greedy A_i subtraction -- literally _fixed_exp_reduce with the table
pointer swapped (theta is a known real accumulator like exp's t);
refactor the exp reduce to take the table as a parameter and share.
Phase 2: series on theta_res < 2^-r: sine is natively odd, so exp's
sinh half-term trick is intrinsic; cos(theta_res) = sqrt(1 - s^2)
via mpn_sqrtrem mirrors exp's reconstruction with a sign flip, same
threshold logic.  Phase 3 (per the identity from
arXiv:2207.02501: with alpha = 2 atan(b/a),
exp(i x) = exp(i (x - k alpha)) ((a + b i)/(a - b i))^k, using
a = 2^i, b = 1 so the table is alpha_i = 2 atan(2^-i)): each
rotation unit is a ratio of conjugates, hence EXACTLY unimodular --
no normalizing square root exists in this formulation at all.
Deferred over the used factors, with W = prod (1 + i 2^-i)^{k_j}
(the real 2^i denominators cancel between numerator and conjugate,
so W is the same two-shift-adds-per-factor fixed-point
reconstruction), the answer is (c + i s) W^2 / |W|^2.  Tail: two
squarings yield x^2 - y^2 and x^2 + y^2 simultaneously (Re W^2 and
|W|^2 from the same pair), one cross product for 2xy, a complex
multiply by (c + i s), and ONE real reciprocal of |W|^2 -- no
square root, and no scale tracking during the loop (one fewer
shift-add per accepted factor than the earlier draft).  The deep
win is structural: the ratio of the SAME truncated W~ against its
own conjugate is exactly unimodular up to the final rounding, so
the ~r reconstruction truncations perturb only the ANGLE of the
correction (<= r ulp, absorbed by the linear budget) and cannot
perturb its modulus -- the scale error cancels by construction,
where a separately tracked M would have needed its own accounting.
The extra step at i = r simply allows k_r = 2 (W takes the factor
twice).  The sin + sqrt(1 - s^2) half-term trick survives only as
the conditional large-N optimization mirroring exp's sinh path;
fixed_sin_cos_rs on shared z-powers is the simple default.

Co-billed variation for sin_cos (less compelling now that the
restoring tail is nearly CORDIC-simple): nonrestoring CORDIC -- every
i <= r with a sign choice (the decision is the sign bit of the
theta accumulator: branchless, no model), making the scale a
CONSTANT K_r = prod_{i<=r} sqrt(1 + 4^-i): tabulate 1/K_r in the
cache, no M tracking, no invsqrt, final 2-mulhigh scale.  Trade:
2-shift-add vector work at every step instead of ~half, ~2.7x the
vector arithmetic -- restoring should win at large n, nonrestoring
plausibly at small n; benchmark decides the dispatch, as
sinh-vs-direct did for exp.

Everything else carries mechanically: C r + D error macros with the
caller-pads-one-limb advice, r = 0 tuning by extending
tune-bitwise-r.c with two more ladders, register small paths from
the log small generator, the 32-bit gating policy, and the
adversarial stress family (sparse sums of A_i for sin_cos; sparse
rotation products under |Z| constraints for atan).  Cost estimates:
atan ~ 1.2-1.5x log_bitwise; sin_cos ~ 1.5-2x exp_bitwise.  Open
items: the alpha-table exactness proof, atan16 for r < 32 at tiny
n, and the restoring/nonrestoring crossover (shifted toward
restoring by the conjugate-ratio tail).

## Patch 0002 progress (trig)

Implemented and passing at 10x (64-bit): atan_bitwise_rs.c (vectoring,
shared alpha table, aliasing-safe series output, r=0 tab {32,64,128}/
{1,18,64}), sin_cos_bitwise_rs.c (rotation via [Joh2022] conjugate-ratio
identity, alpha cache, SINCOS_USE_SQRT threshold, sqrt cos path via
mpn_neg+mpn_sqrtrem, exact-bookkeeping tail, r=0 tab {32,64,128,192}/
{1,12,40,115}). Shared _fixed_bitwise_reduce (exp refactored onto it,
istart/tab params). Validated vs arb 3000 cases: sin 85.8 cos 90.3
atan 67.9 ulp, all in budget (6r+128, 4r+64). Tests t-sin_cos/t-atan
registered; stress extended with alpha-sum and rotation-product
adversaries. All 7 tests pass.

Key fixes during bringup: product lengths in the tail (2wn+1 not 2wn+2,
read a garbage limb); reduce to r+1 for sin_cos because DOUBLED angles
halve the effective reduction (residual < alpha_r ~ 2^{1-r}, violated
the 2^-r series contract at r+0 -> 2^-64 relative errors).

Remaining: coverage on new files; tune ladders for sin_cos+atan (needs
2-output wrapper); real r=0 tables; bench vs arb; docs (fixed.rst +
[Joh2022] in references.rst); pedantic+32-bit shim; run on flint32;
commit 2, format-patch --start-number 2, verify sequential apply.

## Patch 0002 COMPLETE

Committed (bd26fff) and emitted as 0002-fixed-sin-cos-and-atan-*.patch;
0001 then 0002 apply clean in sequence to a pristine tree.  All 7 tests
pass at 10x on 64-bit and at 5x on 32-bit with assertions on.  Coverage:
sin_cos 100%, exp_bitwise 100%, atan 98.4% (sole gap = the 32-bit-only
series fallback).  Tuned r tables: sin_cos {32,64,128,192}/{1,13,83,160},
atan {32,64,128}/{1,21,104}.

THE DESIGN CHANGE THAT MATTERED.  The table stores the UNDOUBLED angles
A_i = atan(2^-i) and the rotation reduces x/2, rather than tabulating
2 atan(2^-i) and reducing x.  This is not just the economy of one shared
table: the windowed single-limb decision model REQUIRES tab_i < 2^-i.
A_i satisfies it (as log(1+2^-i) does); the doubled angles ~2^(1-i) do
not, so after a window-boundary step the residual keeps one bit above
the window limb, the model desynchronizes, and the greedy silently skips
indices.  Caught on 32-bit at n=2 (residual escaped to ~2^-31, tripping
the series contract assert); invisible on 64-bit, where the windows are
twice as wide.

THE LATENT BUG IT EXPOSED.  The table entries are floors, so each accept
under-subtracts by up to an ulp; this creep can push the remainder above
the next entry whenever that entry's deficit from 2^-i falls below an
ulp (roughly 3i > FLINT_BITS*wn for atan, 2i > FLINT_BITS*wn for log).
A single boundary subtraction then leaves the model's precondition
false.  Fixed by repeating the boundary and final subtractions (an index
used twice = its factor applied twice, which every reconstruction
absorbs; used[] sized by FIXED_BITWISE_REDUCE_USED_ALLOC).  exp was
checked to be NOT reachable-buggy: 3M targeted adversarial cases at the
window boundaries, identical worst error (117.6 ulp) with and without
the fix -- its log deficit 2^(-2i-1) is far larger than atan's
2^(-3i)/3.  So patch 0001 stands as frozen; the hardening rides along in
0002, where the reduce becomes shared, and is what makes sharing safe
for any table.

BENCHMARK (r = 0, same-pass ratios vs arb).  atan wins from n ~ 8
(1.1x) rising to 3.3-4.6x at n = 64-128; sin_cos wins from n ~ 24
(1.2x) to 1.6x at n = 64.  BELOW those sizes both LOSE to arb (sin_cos
0.3-0.7x, atan 0.4-0.5x): unlike exp and log, the trig functions have no
generated register small paths yet, so even n = 1 runs the full generic
machinery.  That is the obvious next piece of work.

## Patch 0002, second increment: r = 16 and the atan small path

r = 16 FOR BOTH TRIG FUNCTIONS.  This needed the hardcoded series for
x < 2^-16, which the generator refused to emit for the ALTERNATING
families (atan, sin/cos) -- the deferred "carry-slot collision" item.
Cause: at zbits = 16 the terms shrink by half a limb, so the
leading-zero shrink is disabled (SGen.zl) and every window is full;
the fresh top slot that the RS_* macros ASSIGN to (a[n] = fl +- cy)
is then already holding the previous term's carry limb.  Fix:
accumulating RS_SUBMUL_ACC / RS_ADDMUL_ACC / RS_SUBMUL_Z_ACC fold the
borrow into that slot instead of overwriting it; and the units
coefficient is simply ADDED into the wrapped accumulator -- correct
even with a live sign extension, because the discarded carry out of
limb w is exactly what cancels fl_s = -1 (the block value is
nonnegative, so the sign extension dies there).

The reach is set by the SHARED DENOMINATORS, not the emission: the
series clamp to SC_NMAX / ATAN_NMAX terms while x < 2^-16 wants ~2n of
them, so the first dropped term bounds the precision -- sin/cos
x^21/21! < 2^-401 covers 64n <= 384 (n <= 6), atan x^35/35 goes well
past n = 7.  n = 7 sin/cos was wrong by ~2^46 ulp, matching the
predicted 2^(448-401.5); observed 1.016e14.

Measured (tuner agrees with head-to-head): r = 16 pays to n = 4 for
sin_cos, to n = 6 for atan.  Tables: sin_cos {16,32,64,128,192,256} /
{1,5,11,38,87,128}; atan {16,32,64,128} / {1,7,18,92}.  The tuner's
r_floor is now a per-function parameter (log clamps its generic path
up to 32, so r=16 lives only in its small paths; the trig functions
route r<32 to the wider series in the generic path too).

ATAN SMALL PATH (dev/gen_fixed_atan_bitwise_small.py, n = 1..7).
Mirrors the log1p generator, but a vectoring step needs TWO
funnel-shifted vectors read from the OLD components (v = X>>i to
compare/subtract from Y, w = Y>>i to add into X) where log needs one.
Bug found by the bit-exactness harness: the masked step overwrote Y
before using it to build w, so X got the NEW Y>>i -- failures began
exactly at r = 64 (first c>=1 boundary, where wt names the y
registers rather than temporaries) and n=1 (r <= 48) never failed.
Fixed by updating X and the accumulator BEFORE the Y select.  Now
bit-exact vs the generic over 400k cases, all (n, r); 100% coverage
on all 1953 generated lines.

Speedup vs generic (same pass): n=1 4.4x, n=2 1.45x, n=3 1.25x,
n=4 1.8x, n=5 1.1x, n=6 1.35x, n>=7 ~1.0x.  atan n=1 now 60 ns
(was 580 before r=16 + small path).  vs arb it is jagged, since arb's
own atan is jagged in precision (fast at n=2,5,7, slow at n=4,6).

REMAINING: the sin_cos small path.  Same masked-step machinery, but
its rotation state differs from atan's -- wx DECREASES from 1 (to
> 0.72), so unlike atan's X it needs an explicit units limb, while wy
stays below 0.92 and fits n limbs.  The W factors can be applied on
the fly (no used[] list), and the table-driven boundary step must be
emitted twice to absorb the creep.  The tail (two squarings, complex
multiply, two divisions) stays as mpn calls.

## Patch 0002 (final) + patch 0003

SIN_COS SMALL PATH (dev/gen_fixed_sin_cos_bitwise_small.py, n=1..7).
Combines the two earlier ones: table-driven DECISIONS (as in exp) with
a rotating-vector STATE (as in atan), the W factors applied on the fly
so no used[] list is needed.  W's real part needs an explicit units
limb -- unlike atan's X, which only grows from 1, wx DECREASES (to
~0.72), so wu falls 1 -> 0.  With the UNDOUBLED table the residual is
strictly below the window after the boundary step, so no stray-bit
check is needed (atan's geometric residual, bounded only by
X 2^-i < 2^(1-i), does need one).  Boundary and terminal steps are
emitted twice to absorb the truncation creep.  Bit-exact vs generic
over 400k cases FIRST TRY (the X/acc-before-Y-select ordering
discipline learned from atan carried over).  100% coverage of all 2723
generated lines -- after adding r = FLINT_BITS*k - 1 to the
deterministic sweep: the reduction runs to R = r+1, so the boundary
form of the terminal step needs r one BELOW a multiple of FLINT_BITS,
which the old sweep (multiples and multiples+1) never produced.

Speedup only 1.0-1.7x, vs atan's 4.4x, and sin_cos still loses to arb
at small n.  MEASURED WHY: the reconstruction, not the reduction, is
the bulk of the work -- the two divisions by |W|^2 alone cost 89 ns at
n=1 and 366 ns at n=6, i.e. a third to a half of the whole call, while
the series is only 20-133 ns.  So the tail is where the remaining win
is.  Two options: (a) one reciprocal + two multiplies instead of two
divisions (the shared _fixed_sin_cos_recon helper now makes this a
one-place change benefiting both paths); (b) NONRESTORING rotation,
where every factor is applied with a sign, making the scale factor a
CONSTANT K_r that can be tabulated -- eliminating |W|^2, both
squarings and both divisions, at the price of ~2x the vector work in
a reduction that is now nearly free.  The measurement makes (b) look
much more attractive than it did when the reduction dominated.

PATCH 0003: r = 16 FOR EXP.  The exp series generator extends to
zbits = 16 almost trivially: exp's series runs in x, not x^2, so terms
lose 16 bits apiece and lz(j) = j*zbits/64 stays limb-aligned (4 terms
per limb) -- none of the fractional-limb trouble the trig families
hit.  Only constraint: m must be a multiple of 64/zbits (block bases
need lz(b)-lz(m) = lz(b-m)), which leaves N=4 with no legal split
below m=N, so a single block (Horner) is admitted.  Reach set by the
shared denominator: N = 4n terms, 20! supports N <= 21 -> n <= 5.  The
64/32 families regenerate byte-identically.

VERDICT (measured, not assumed): like-for-like in the generic code
r=16 beats r=32 at n=1,2,3,5.  BUT exp's register small path at r=32
beats generic-r=16 nearly everywhere -- so the small path had to be
taught the wider series too.  Only then does r=16 win: ~23% at n=1,
16% at n=3, 8% at n=4; it LOSES slightly at n=2 and n=5 (the series
cost is jagged in N).  Tuner picks r=16 up to n=4.  Default table
prepends r=16 and keeps 0001's tuned entries above unchanged.
Accuracy 12-15 ulp (budget 244).  32-bit clamps back to 32.

## Patch 0004: fully specialized exp for small n

The idea (spend the 1/N! instead of banking it) is confirmed and shipped.
The generator tied precision to xn = ceil(N*zbits/64), i.e. it demanded
N*r >= 64n and treated log2(N!) as slack.  Gen now takes xn explicitly
and accepts ANY zbits, so N is the least integer with
N*r + log2(N!) >= 64n.  At r=16 that is N = 11, 14, 17 for n = 3, 4, 5
where 4n = 12, 16, 20 terms were used.  n=2, r=20 -> N=6 exactly as
predicted; the landmarks shift too (N=4 at n=2 needs only r=31).

Selected by measurement over a bit-granular sweep (267 series, all
validated vs arb):
  n=1: r=12, N=5 (single block)   n=2: r=16, N=8
  n=3: r=16, N=11                 n=5: r=16, N=17
~5-10% over the previous defaults; accuracy 12-15 ulp.  Dispatch is on
r == 0 (the caller leaving r to us); explicit r keeps the old paths.

TWO FINDINGS THAT SHAPED THE RESULT.

(1) Away from the divisors of 64, rectangular splitting is ILLEGAL: the
block bases need lz(b) - lz(m) = lz(b-m), which holds only when
64 % zbits == 0.  So a bit-granular series must be a single block, and
that costs about what the shorter reduction saves -- which is exactly
why the landmarks sit at 16/32/64.  Disabling the leading-zero shrink to
legalize splitting everywhere is a WORSE trade (lost up to 49%: the
shrink is worth more than the splitting).  Hence the flat basin and the
modest gains.

(2) A compile-time-constant r pays only while r is SMALL.  Specializing
at r=32 was markedly SLOWER than passing 32 at run time (194 vs 141 ns
at n=4): gcc fully unrolls the 32-iteration window into more code than
the folding saves.  All specialized sizes therefore sit at r <= 16.

n=4 is deliberately NOT specialized: there r=32 measured fastest, and at
zbits=32 the count N=2n is already factorial-minimal, so there is
nothing to spend -- its default simply moves to 32.

## Patch 0005: hand-written 1- and 2-limb exp series

The n=1 and n=2 specialized series are now hand-written
(exp_rs_opt_hand.inc); the generator skips them.  At those sizes the RS
framework (sqrhigh/mulhigh_n/addmul_1 calls, shift-split machinery,
fac20_inv rescale) dwarfs the actual arithmetic.

Plain Horner over c_k = 1/k! (hardcoded n-limb fixed constants), high
products keeping only the top n limbs.  Structure per your suggestion:
1 + x + x^2*(C2 + x*(C3 + ...)).  N cut to the true minimum (4 for n=1,
7 for n=2) since the first dropped term is already 2^-67 / 2^-143.

n=2 uses a SLOPPY 2x2 mulhigh dropping the a0*b0 partial entirely
(MULH2): 3.0 ulp vs 1.5 for the exact fold, both trivially inside the
~4708 budget, and 18 vs 22 ns -- so sloppy wins.  Your instinct was
right.

RESULTS (same-process): n=1 full call 45 -> 26 ns (4.2x arb, was 3.2x);
n=2 series 33 -> 18 ns.  Robustness 12/15 ulp over 5M adversarial cases.
100% coverage of the hand code; 32-bit correctly excludes it (7/7).

n=3 was PROTOTYPED and REJECTED: a 3x3 sloppy mulhigh is 5-6 partials
per multiply x 11 Horner steps = genuine arithmetic, and the hand
version tied the generated one (53 vs 54 ns).  That is the boundary --
once the series does real work the framework overhead is noise, so
hand-writing only pays at n <= 2.

All 5 patches apply clean in sequence; 7/7 tests at 10x (64-bit) and 5x
(32-bit, asserts on).

## Patch 0006: mixed-precision hand series

The n=1/n=2 hand evaluators did all work at full width; they need not.
Error attenuation: exp(x) = 1 + x + x^2*P, P = C2 + x*Q, Q = C3 + x*R...
err(P) enters scaled by x^2 < 2^-2r, each level below by another x.
At n=2 (ulp 2^-128, x < 2^-16), holding series error ~2^-122:
    P 90 bits, Q 74, R 58, S 42, U 26
So R, S, U -- THREE of the five Horner steps -- fit in ONE limb, computed
on x's top limb alone.  Plus x^2 is a squaring (cross products coincide,
2 umuls not 3) and x*R multiplies by a single limb (a column vanishes).

RESULT: n=2 series 21 -> 13 word multiplies, 16.2 -> 8.5 ns (1.9x);
full n=2 call 82 -> 73 ns.  n=1 restructured as 1 + x + x^2(C2 + x*C3 +
x^2*C4) (reuses x^2, independent terms instead of a serial chain):
26 -> 23 ns, now 4.4x arb_exp.  Accuracy UNCHANGED (12/15 ulp over 5M
cases, budget 4708) -- the series rounding is invisible beside the
reduction's own error.

TRIED AND REJECTED:
- plain non-widening muls on pre-shifted operands for the small terms
  (they need only 20 and 6 bits): no faster on x86-64, imul == mulx.
  Would matter on a target lacking a fast widening multiply.
- __uint128_t high word instead of umul_ppmm: identical codegen.
NOTE: n_mulhi and ulong_extras/ll_is_prime.c are NOT in the 3.3.1
tarball (they are newer main-branch additions), so they could not be
used or compared here.

Also fixed: the MULH1 macro used a GNU braced-group expression; now in
ISO C statement form (pedantic-clean).

## Patch 0007: n_mulhi + sloppy 2x2 (main-branch helpers)

Adopted both helpers you supplied.  mulhi_2x2_sloppy is SLOPPIER than my
MULH2: it discards the low words of BOTH cross products (not just
a0*b0), giving 3 muls + 2 adds vs my 3 muls + 4 adds.  Series error
5.1 -> 6.1 ulp (budget 4708), n=2 series 8.5 -> 7.2 ns, full n=2 call
73 -> 71 ns.

n_mulhi's __uint128_t form generates the SAME code as umul_ppmm (measured
identical), so it is adopted purely as the pre-existing, clearly-named
function -- exactly the reason you gave.

KEPT two specializations over the generic 2x2, both measured better:
  - multiply by a ONE-limb value (a whole column vanishes): 2 muls
  - squaring (the two cross products coincide): 2 muls
Routing the squaring through the generic sloppy 2x2 was SLOWER
(7.4 vs 7.2 ns), so the dedicated form stays.

Compat: n_mulhi postdates 3.3, so a definition guarded on
__FLINT_RELEASE < 3.4.0 is included (inert on main).  mulhi_2x2_sloppy
is file-local to ll_is_prime.c on main, so it is duplicated here --
promoting it to a shared header would suit both callers.

CUMULATIVE small-n exp: n=1 45 -> 23 ns (4.5x arb), n=2 82 -> 71 ns
(3.5x arb).  100% coverage of the hand code; 7/7 tests both word sizes.

## Patch 0008: hand-written 3-limb series

n=3 now hand-written too; the generator emits only n=5.

WHY THE EARLIER ATTEMPT FAILED: it did all the work at 3 limbs and tied
the generated code.  The same error attenuation as n=2 splits the
11-term Horner:
    V_2 154 bits, V_3 138, V_4 122, V_5 106, V_6 90, V_7 74,
    V_8 58, V_9 42, V_10 26
-> only V_2,V_3 at 3 limbs; V_4..V_7 at 2; V_8..V_10 at 1.
Just FOUR genuine 3-limb multiplies remain (x*V_4, x*V_3, x^2, x^2*V_2)
-- flint_mpn_mulhigh_n / flint_mpn_sqrhigh for those, registers above.
Plus a dedicated 3x2 (V_4 embeds with a zero low limb: 3 full + 2
high-only products vs 6 for a general 3x3), worth ~4%.

RESULT: series 57 -> 36 ns; full n=3 call ~104 -> 84 ns (3.3x arb, was
2.7x).  Error 15 ulp / 5M cases (budget 4708).

WEIGHED AND REJECTED:
- rectangular splitting: with only two 3-limb Horner steps left there is
  nothing to split.
- C2 = 1/2 as a free shift (frac = x + x^2/2 + x^2*(x*V_3)): needs an
  EXTRA 3-limb multiply, measured slower (42 vs 36 ns).

CUMULATIVE (vs the generated specialized series):
  n=1  45 -> 23 ns   n=2  82 -> 71 ns   n=3  ~104 -> 84 ns
100% coverage of the hand code; 7/7 both word sizes; pedantic-clean.

## Patch 0009: hand 4-limb series -- and r MOVES

Your hypothesis confirmed, and it matters: a cheaper series can afford
more terms, so it wants a SHORTER reduction.

n=4: with the GENERATED series r=32 was fastest (which is why 0004 left
n=4 unspecialized).  With a mixed-precision HAND series the optimum
moves down to a broad basin r ~ 15-21, and the full call drops ~13%.
Now specialized at r=16 like its neighbours; it no longer consults the
default table.  Level widths (N=14): V2:4 V3:4 V4-V7:3 V8-V11:2
V12-V13:1 -- only the two outermost steps at full width.

n=5: swept the same way, NO GAIN (206 vs 202 ns) -> KEPT GENERATED.
With 17 terms the chain is long enough that the generator's tuned
RECTANGULAR SPLITTING beats a serial Horner at any precision.  That is
the real boundary: mixed-precision Horner wins while the chain is
short; splitting wins once it is long.  Combining them (splitting with
per-block precision) is the obvious next step and is NOT attempted.

NEW TOOL: dev/gen_fixed_exp_hand_mixed.py emits hand-style mixed
evaluators for ANY (n,r) from bits(V_k) = 64n - 6 - r*k.  This let me
sweep r with hand-quality series at every candidate instead of guessing
-- which is what made the n=4 result trustworthy.

CUMULATIVE small-n exp (same process, vs arb_exp):
  n=1 21.5 ns (4.7x)   n=2 66 ns (3.5x)   n=3 74 ns (3.3x)
  n=4 106 ns (1.6x)    n=5 198 ns (2.7x)
Errors 12-16 ulp (budget 4708); 100% coverage; 7/7 both word sizes.

## Rectangular splitting with per-block precision (n=5): TRIED, NOT SHIPPED

Implemented and validated; it works, but the gain is ~2.5% on the series
and ~1.4% on the full call, which does not justify a hand-written 5-limb
evaluator plus its own error analysis.  The tree stays at 9 patches.

CONSTRUCTION.  The two ideas do combine: RS gets its speed from
SINGLE-LIMB coefficients (c_k = 20!/k!, so each term is an addmul_1),
which is orthogonal to precision.  So keep the addmul_1 blocks and let
each block's accumulator be only as wide as it needs to be.  With
S_b = sum_{k>=b} c_k x^(k-b) and S_b = acc_b + x^m * S_{b+m}, the
binding constraint is truncating the powers to the accumulator width:
    64 W_b >= 314 - 16b - log2((b+1)!)
For n=5, r=16, N=17, m=4 that gives accumulator widths 5,4,3,2,1 for
b = 0,4,8,12,16, against 5,5,5,5,5 today.

KEY TRICK.  The first attempt was SLOWER (139 vs 114 ns) because the
boundary step z^4 * S was a full mpn_mul.  But x^4 < 2^-64 exactly when
r*m >= 64, so z^4's top limb is zero: shift it up one limb and the
boundary becomes an EXACT flint_mpn_mulhigh_n at the width of S.  The
boundaries then run at 2,3,4,5 limbs instead of 5,5,5,5 -- that is where
essentially all of the gain comes from.

MEASURED (same process): generated 110.6/112.5 ns; blocked m=4
108.7/108.7; blocked m=5 116.3/120.2 (more powers, fewer boundaries --
a bad trade); mixed Horner 134.5.  Error 2.6 vs 3.0 ulp.
m >= 4 is forced: the shift trick needs r*m >= 64.

WHY SO SMALL.  Per-block precision only saves work in the DEEP blocks,
and those were already the cheap ones.  The dominant costs are shared
and unchanged: three 5-limb power multiplies (z2, z3, z4), the b=0 block
at full width, and the 6-limb final rescale by 1/20!.  The savings are
bounded by the fraction of work sitting in the deep blocks, which is
small.

This also confirms the boundary found in 0009 is robust: mixed-precision
Horner for n <= 4, rectangular splitting for n >= 5 -- and the splitting
is already close to optimal at n=5, with little left on the table.

## Patch 0010: hand atan series + re-optimized r  (log/sin_cos: NEXT)
## [amended] every op now width-appropriate; NN_ADD/NN_SUB macros

The exp treatment carried to atan.  These are series in z = x^2, so the
attenuation is TWICE as strong per level as exp's:
    atan(x) = x - x^3 V_0,  V_0 = 1/3 - z(1/5 - z(1/7 - ...))
    bits(V_k) = 64n - 6 - r(2k+3)
Level widths collapse fast (n=4: 4,3,3,2,2,1,1).  All V_k stay positive
(z < 2^-2r swamps the alternating term) -> no sign handling.

r MOVES A LOT.  With the cheap series the optimum falls to
    r = 4, 6, 21, 18   for n = 1..4   (generic code wants 16 or 32)
Gains on the WHOLE call: +46%, +30%, +30%, +15%.  Error 2-13 ulp
(budget 2112).  Now 3.0x / 2.5x / 1.1x / 1.8x arb_atan.

INFRASTRUCTURE (dev/gen_fixed_trig_hand_mixed.py) emits mixed-precision
series for ANY (function, n, r) over atanh, atan, sin, cos.  All 656
generated variants VALIDATED against arb (a real bug it caught: cos is
alternating too -- my first version emitted V = c + z*V).
Helpers factored into src/fixed/hand_mulhi.inc, shared by exp and trig.

REMAINING (infrastructure is ready, series already generated+validated):
  - log1p: the atanh series is emitted and validated; needs the same
    constant-r opt path in gen_fixed_log1p_bitwise_small.py + r sweep.
  - sin_cos: sin and cos series emitted and validated; needs a two-output
    opt path (they SHARE z = x^2, so one sqrhigh serves both) + r sweep.
Expect similar gains: the attenuation argument is identical, and both
have the same cheap-series-wants-smaller-r dynamic.


## 0010 amendment: inline everything narrow

Fixed the gap you spotted: the emitter still bottomed out in
flint_mpn_sqrhigh / flint_mpn_mulhigh_n at LENGTH 1 AND 2, and used
mpn_add_n / mpn_sub_n throughout.  Now every op is chosen by width:
  len 1  -> n_mulhi                         (and n_mulhi(x,x) for sqr)
  len 2  -> _fixed_mulhi_2x2_sloppy / _fixed_sqrhi_2x2
  len 3+ -> flint_mpn_mulhigh_n / flint_mpn_sqrhigh
  add/sub-> NN_ADD_k / NN_SUB_k  (mpn_extras.h)
At n=1 there is now NO mpn call at all.

WORTH A LOT: inlining alone HALVES n=1 (21.7 -> 10.2 ns), on top of what
the mixed precision already gave; n=2 +15%, n=4 +9%.  That is the same
lesson as exp's n=1: when the arithmetic is a handful of word ops, the
call/framework overhead IS the cost.

Re-swept r afterwards (the cheaper series shifts the optimum again):
r = 4, 6, 22, 20 for n = 1..4.  All 720 emitted series re-validated.
atan now 3.4x / 2.9x / 1.2x / 1.9x arb_atan; error 2-19 ulp (budget 2112).

STILL TO DO (unchanged): log1p (atanh) and sin_cos.  Series already
emitted and validated by the same generator; each needs the constant-r
opt path + r sweep.  sin_cos additionally wants a two-output path --
sin and cos SHARE z = x^2, so one squaring serves both.

## Patch 0011: log1p + sin/cos hand series (with the shared square)

Completes the module.  Both are series in z = x^2, same collapse of level
widths, same width-appropriate ops (inline <=2 limbs, mpn >=3,
NN_ADD/NN_SUB).

SHARED SQUARE: sin and cos are emitted TOGETHER from ONE sqrhigh -- a
saving neither separate generated series can take (each computes its own
x^2).  emit_sincos() in dev/gen_fixed_trig_hand_mixed.py.

r MOVES AGAIN (cheap series -> shorter reduction):
  log1p    r = 10, 26        (n=3,4)      was 16, 32
  sin_cos  r = 5, 8, 12, 19  (n=1..4)     was 16, 32
GAINS: log1p +22% (n=3), +4% (n=4); sin_cos +14..21% (all n).
Errors 5-46 ulp vs budgets 91-242.  360 more series validated vs arb.

log1p n=1,2 LEFT ALONE: their register code already evaluates atanh
INLINE (patch 0001 hand work), so there is no framework overhead to
remove -- only r could move, which alone does not justify replacing
hand-tuned code.

FINAL vs arb (n=1..4):
  log1p   4.1x  2.0x  2.2x  2.5x
  atan    4.2x  1.0x  1.3x  1.9x
  sin_cos 0.9x  0.5x  0.5x  0.7x   <- still LOSES; its reconstruction
    (two divisions by |W|^2) dominates, as measured back in patch 0002.
    The series work cannot fix that; the nonrestoring rotation (constant
    scale factor K_r, no |W|^2, no divisions) is the open lead.

11 patches, all apply clean; 7/7 tests both word sizes.

## tan half-angle for sin/cos: PROTOTYPE WINS BIG (not yet wired)

Your idea works, and the key is that TAN IS SCALE-INVARIANT.  Since
arg(W) = sum atan(2^-i) with W = prod (1 + i 2^-i), we have
tan(arg W) = wy/wx -- W's modulus growth CANCELS IN THE RATIO.  So the
whole reason for the conjugate-ratio trick (avoiding normalization)
evaporates: no |W|^2, no squarings, no unimodularity argument.

    x/2 = sum A_i + t'          (the SAME reduction, no doubling needed)
    tan(x/2) = (wy + wx u) / (wx - wy u),      u = tan(t')
    sin = 2t/(1+t^2),  cos = (1-t^2)/(1+t^2),  tan(x) = 2t/(1-t^2)

Well-conditioned throughout: t < tan(0.5) = 0.5464, denominator
wx - wy*u ~ wx in (0.72,1.17) -- no cancellation.  Tan's coefficients
(tangent numbers, no hypergeometric recurrence) are just tabulated;
at n<=4 only a handful are needed, exactly as you predicted.

MEASURED (prototype still uses GENERIC mpn reduction + W rotation, while
"shipped" uses the register small paths -- so this is a LOWER bound):
  n | sin+cos | shipped | gain | sin only | tan only | err
  1 |  136.4  |  232.7  | +41% |  140.2   |  106.8   |  6
  2 |  240.3  |  373.8  | +36% |  248.0   |  190.7   | 14
  3 |  354.8  |  534.1  | +34% |  278.5   |  195.5   | 96
  4 |  453.9  |  671.4  | +32% |  354.8   |  255.6   | 39
Budgets are 150-250 ulp, so accuracy is fine (n=3's 96 wants a slightly
larger r or one more term -- to be tuned when wired).

AND standalone functions fall out: tan alone is ~45% cheaper than
sin+cos, sin alone ~25% cheaper -- no more computing two functions and
throwing one away.

BUGS FOUND EN ROUTE (both classic): 2t can EXCEED 1 (t up to 0.546) so
the outputs need n+1 limbs; and DE ~ wx in (0.72,1.17) often has a ZERO
TOP LIMB, so mpn_tdiv_qr needs its divisor length trimmed -- the same
unnormalized-divisor trap as the stress-test bug in 0002.

NOT YET WIRED: needs the recon in sin_cos_bitwise_rs replaced, the small
paths regenerated for the new tail, r re-tuned, and public fixed_tan_*
added.  180 tan series already emitted and validated.
Prototype: fixed/dev/prototype_tan_halfangle_sincos.c

## Patch 0012: tangent half-angle -- sin, cos AND tan

Shipped.  The conjugate-ratio tail is gone from the specialized sizes.

WHY IT WORKS: arg(W) = sum A_i, so tan(sum A_i) = wy/wx -- and TAN IS A
RATIO, so |W|'s growth CANCELS.  No |W|^2, no normalization, no
unimodularity.  The reduction also loses its doubling (t' wanted as is).
    t = tan(x/2) = (wy + wx u)/(wx - wy u),   u = tan(t')     [1 division]
    sin = 2t/(1+t^2), cos = (1-t^2)/(1+t^2), tan = 2t/(1-t^2) [share t^2]

WINS AT EVERY n = 1..8 (+33..+52%).  No crossover in range, so the
direct tangent is used for n <= 6 (series tabulated); beyond that tan =
sin/cos with one division.  (Fredrik's fallback idea -- get u as
sin(t')/cos(t') with one extra division -- is the route if the tan
series ever gets unattractive at larger n; not needed here.)

RESULTS vs arb:
  sin_cos  1.09 0.96 0.69 0.64 0.76 0.80   (was 0.89 0.49 0.51 0.67)
  tan      2.13 2.31 1.28 1.35 1.21 1.18   (NEW public function)
Errors 7-107 ulp (budgets 3200/4352).  8/8 tests both word sizes.

TRAPS: 2t EXCEEDS 1 (t < tan(1/2) = 0.5464) -> outputs need a unit limb;
denominators sit near 1 so their top limb is often ZERO -> mpn_tdiv_qr
needs its divisor length TRIMMED (the same unnormalized-divisor trap as
in 0002 and the prototype).

12 patches, all apply clean.

## Patch 0013: tangent half-angle EVERYWHERE; old sin/cos tail DELETED

Answer to "is it simply better at every n": YES.
  n     8    16    32    64   128   192
  gain +32%  +25%  +13%  +13%  +18%  +10%
(plus +33..52% below n=8).  No crossover anywhere, so it becomes the ONLY
reconstruction.

Large n uses Fredrik's own fallback: tan(t') = sin(t')/cos(t') with ONE
extra division, reusing the existing sin/cos series -- no need to
tabulate tangent numbers indefinitely (they have no tidy recurrence).
That same route serves 32-bit, where the hand tan series don't exist.

DELETED:
  - |W|^2 normalization: 2 squarings + 2 WIDE divisions (3n+3 over 2n+1)
  - the doubled-angle reduction (the tangent wants t' as it stands)
  - sin_cos_bitwise_rs_small.inc + its generator
  - SINCOS_USE_SQRT and the sinh-style sqrt path
  - sin_cos_bitwise_rs.c: 380+ -> 162 lines; fixed_sin_cos_bitwise_rs is
    now a dozen lines calling the reconstruction.

What remains is: the reduction, the rotation, and the observation that
arg(W) = sum A_i so tan(sum A_i) = wy/wx -- and since tan is a RATIO,
|W|'s growth CANCELS.  No normalization anywhere.

13 patches, all apply clean; 8/8 tests both word sizes, warning-free.

## Patch 0014: rotation in registers; avoidable copying removed

Fredrik spotted that the W rotation was carrying dead work: each factor
materialized wx>>i and wy>>i into SCRATCH -- zeroing them, filling with
mpn_rshift or copyi, then mpn_sub + mpn_add -- when the components fit
in registers, the shift is two funnel shifts, and the updates are single
sub_ddmmss / add_ssaaaa chains.

dev/gen_fixed_tan_rotate.py emits it in PURE REGISTERS for n<=6 (no
scratch arrays at all): wx keeps the unit limb (it falls from 1 toward
0.72), wy needs none (stays < 0.92); the runtime index picks one of a
handful of emitted cases by limb offset.  The generic loop keeps buffers
but now zeroes ONLY the limbs the shift leaves stale.

RESULT vs arb:
  sin_cos  1.26 1.15 0.86 0.67 1.03 1.03   (was 1.09 0.96 0.69 0.64 0.76 0.80)
  tan      2.47 2.82 1.70 1.34 1.60 1.61
~15-25% off the call.  14 patches, all apply clean; 8/8 both word sizes.

## Patch 0015 + OPEN WORK (honest status)

SHIPPED: flint_mpn_mul in place of mpn_mul (your catch).  15 patches,
all apply clean, 8/8 both word sizes.

NOT SHIPPED -- divisor normalization.  The idea is RIGHT and measured
FASTER, but I hit a bug in the tan branch and ran out of budget to chase
it safely, so I REVERTED rather than ship a broken tan.

  The insight: 1 + t^2 in [1, 1.2984) is the worst possible divisor --
  n+1 limbs to hold ONE BIT above the point.  Halve it:
      D = (1+t^2)/2 in [0.5, 0.6492)  -> n limbs, TOP BIT SET
  then  sin = t/D,  cos = ((1-t^2)/2)/D.  One limb shorter AND
  normalized, which is what mpn_tdiv_qr wants.  For tan, 1-t^2 in
  (0.7015, 1] already has its top bit set (t=0 is the exception).

  MEASURED with this in (sin/cos passing, tan failing):
    sin_cos  148.8 -> 118.3 (n=1), 255.6 -> 201.2 (n=2),
             412.0 -> 324.2 (n=3), 465.4 -> 385.3 (n=4)
  i.e. ~20% off the call.  Worth finishing.

  The bug is in the ytan branch only (sin/cos passed all tests).  Likely
  suspects I did not get to rule out: the n+2-limb quotient (tan(1)<1.56
  needs n+1, GMP writes n+2), and mpn_tdiv_qr overlap rules between the
  quotient buffer and the numerator.

ALSO OPEN (your list):
  - benchmark 2x mpn_tdiv_qr vs 1 reciprocal + 2 mulhighs for sin/cos
  - specialize the sin/cos-from-tan reconstruction for n=1..4 (everything
    but the division in registers)
  - re-tune r for n=3,4 (still the only sizes below arb: 0.86x, 0.67x)

## FIX: stray scratch file in patch 0011

`src/fixed/trighand.c` -- a bare dump of generated series, no includes, only
ever meant to be #included into a test harness -- was accidentally copied into
src/fixed/ and swept up by `git add -A src/fixed` in patch 0011.  FLINT's build
GLOBS every .c in the module dir, so it broke the build.

Removed from history (rebase of 0011..0015); patches regenerated.  All 15 still
apply clean; 8/8 tests both word sizes.

WHY I MISSED IT: I was compiling EXPLICIT FILE LISTS in every check, which never
exercises the glob.  Now verified with:
    for f in src/fixed/*.c; do gcc ... -c $f || echo FAILS; done
All 7 module .c files compile standalone.  Audited every file the 15 patches
touch: no other strays.

## Divisor normalization SHIPPED (the ytan bug was a quotient-length bug)

The reverted ~20% idea from 0015 is now in, and the ytan failure is
explained: it was never about the normalization.  The old tan division
fed a (2n+1)-limb numerator (2t) to a divisor TRIMMED to n limbs
(1 - t^2 < 1 has a zero top limb), and mpn_tdiv_qr writes exactly
nn - dn + 1 = n + 2 quotient limbs -- ONE PAST the documented
(res, n + 1) contract.  The tests never saw it (they allocate res[66]);
a canary at res[n + 1] trips on the old code and is clean on the new.
The same latent overwrite sat in the large-n fallback (sin/cos then
divide).

The fix and the speedup are the same change: HALVE numerator and
denominator so that every divisor is a normalized n-limb value and
every quotient fits n + 1 limbs by construction:

    sin = t/D,  cos = ((1 - t^2)/2)/D,   D = (1 + t^2)/2 in [0.5, 0.6492)
    tan = 2 (t/(1 - t^2))                1 - t^2 in (0.7015, 1)

- D has its TOP BIT SET: no trimming, no normalization shifts inside
  mpn_tdiv_qr, one limb shorter than 1 + t^2.
- numerators become n-limb (t < 0.5464, (1 - t^2)/2 <= 1/2), so nn = 2n
  and nn - dn + 1 = n + 1 exactly.
- tan divides t/(1 - t^2) = tan/2 < 0.78 and doubles by a final 1-bit
  shift of the (n+1)-limb quotient (costs <= 1 ulp extra, budget 8r+256).
- t^2 = 0 at working precision: cos numerator is 1/2 exactly (cos comes
  out exactly 1, no special-cased division), tan short-circuits to 2t.
- fallback tan = sin/cos: numerator n limbs (sin < 0.842 never uses its
  unit limb); when cos rounds to 1 the divisor keeps its unit limb, the
  quotient is n limbs, and res[n] is pre-zeroed.

8/8 tests, 25x multiplier on the trig ones; 200k-case canary clean.

## r RE-SWEPT for the tangent half-angle (all n = 1..6)

With the register rotation (0014) and the normalized tail both cheaper
than what the old parameters were tuned for, the optima moved DOWN, as
they always have ("a cheaper tail wants a shorter reduction"):

    n            1    2    3    4    5    6
    r (old)      5    6   15   27   27   30
    r (new)      4    6    9   14   15   16

Sweep: /tmp/rsw_tan.c pattern (parametrized copy of _fixed_tan_halfangle
+ generated hs_tan_n_r series for n <= 6, r = 4..48; best-of-5 timing,
two independent processes agreeing on the ranking).  n = 6 tied between
r = 15 and 16; 16 taken for its cleaner measured error (9/9 vs 65/43
ulp).  Errors validated at 12000 adversarial-biased points per size:
sin <= 26, cos <= 18, tan <= 56 ulp, budgets 6r+128 / 8r+256.

MEASURED (one process, best-of-5, ratio = arb time / fixed time):
    n            1     2     3     4     5     6     7     8
    sin_cos   2.06  1.65  1.29  1.33  1.42  1.40  0.85  1.03
    tan       3.39  2.98  2.06  2.04  2.12  2.02  1.13  1.32
(previous session, different machine: sin_cos 1.26 1.15 0.86 0.67 1.03
1.03 -- cross-machine ratios are only indicative, but n = 3, 4 were the
only sizes below arb and are now comfortably above, and every specialized
size beats arb for both functions.)  n = 7 is the first fallback size
(no tan series, generic rotation) and is the remaining soft spot.

## One reciprocal + two mulhighs beats two divisions for sin/cos

Open item 3 measured and shipped.  Three tails compared in isolation
(from t and t^2 to sin and cos, /tmp/tailbench.c pattern):

    n          1     2     3     4     5     6     8    10    12
    2 x tdiv  95   141   231   283   377   448   621   833  1068
    reciprocal 34    99   169   174   231   317   386   622   848
    preinvn   220   291   336   436   496   556   675  1148  1468

The reciprocal route computes R = 1/(1 + t^2) = 1/(2D) with ONE
division, then sin = 2 t R and cos = (1 - t^2) R are each a single
mulhigh: a third to a half off the tail at every size, at most 2 ulp
from the division results.  flint_mpn_preinvn + divrem_preinvn loses
even to two plain divisions here -- the Newton setup only amortizes
over many divisions by one divisor.

Two traps found by the harness itself:
  - the reciprocal numerator must be 2^(128n-1) - 1, not 2^(128n-1):
    when t^2 is one ulp, D = 2^(64n-1) EXACTLY and the quotient hits
    2^(64n), overflowing the n limbs the mulhighs read.  (All-ones
    numerator caps R at 2^(64n) - 1.)
  - a first version zeroed only n limbs of the numerator buffer and
    left stale limbs in the middle -- the errors looked like a broken
    identity until exact arithmetic confirmed the algebra was fine and
    pointed back at the buffer.

The reciprocal is used only when BOTH outputs are wanted (that is the
common call); a single output keeps its one direct division.  The
t^2 = 0 pair case now short-circuits to sin = 2t, cos = 1 with no
division at all.  R lives in Q + n (t occupies Q[0..n-1], the buffer
has 2n + 2).

MEASURED, whole call, one process (ratio = arb / fixed):
    n            1     2     3     4     5     6     7     8
    sin_cos   2.25  1.93  1.45  1.56  1.65  1.53  0.89  1.15
    tan       3.36  2.99  2.14  2.18  2.14  2.01  1.15  1.37

## 32-bit: 8/8 with asserts (plus two environment notes)

The full suite passes on an ABI=32 tree with FLINT_WANT_ASSERT.  Two
issues in the ENVIRONMENT, not the module, for whoever repeats this:
Ubuntu's fat (cpu-dispatched) i386 libgmp exports __gmpn_add_nc etc.
only with CPU suffixes, but configure's AC_SEARCH_LIBS misdetects them
as present ("none required") -- comment the corresponding
FLINT_HAVE_NATIVE_mpn_* lines out of src/config.h (NOT
src/flint-config.h, which does not carry these).  And the merged-object
ld -r step needs make LD="ld -m elf_i386" or it emits elf64 stubs.

## First division normalized too; n = 1 reconstruction in REGISTERS

Two more shipped on top of the reciprocal tail.

FIRST DIVISION HALVED (all n).  t = (wy + wx u)/(wx - wy u) had the
same disease as the old tail: denominator in (0.72, 1.17), so its top
limb is 0 or 1 and mpn_tdiv_qr either gets trimmed or normalizes
internally.  The numerator is t < tan(1/2) = 0.5464 times the
denominator, hence < 0.64 and never uses ITS unit limb -- so when the
denominator reaches one, halve BOTH: the divisor becomes an n-limb
top-bit-set value in every case, at <= 2 ulp of extra floor.
_fixed_divide now serves only the u = sin/cos fallback division.

n = 1 IN REGISTERS (open item 4, first size).  After the reduction the
whole call is straight-line: the rotation and series already were, and
the reconstruction becomes two umul_ppmm chains (products against wx
and wy, unit limbs 0/1 handled by conditional adds) and one udiv_qrnnd
per division, each divisor rounded to 64 normalized bits by the same
conditional halving.  Errors 7/5/16 ulp (sin/cos/tan) at 200k points
against budgets 152/288.

MEASURED, whole call (one process, ratio = arb / fixed):
    n            1     2     3     4     5     6     7     8
    sin_cos ns  42.4 144.4 253.9 329.4 418.2 539.9
    sin_cos   4.91  2.01  1.56  1.64  1.70  1.52  0.91  1.16
    tan       6.91  3.03  2.19  2.27  2.22  1.99  1.17  1.38
n = 1 fell from 99 to 42 ns -- the mpn dispatch and buffer traffic WAS
two thirds of the call.  Extending the register reconstruction to
n = 2..4 (via a generator in the shape of gen_fixed_tan_rotate.py,
using the 2x2/3x3 forms of hand_mulhi.inc and preinv 3/2 division
steps) is the natural next step; the n = 1 result bounds what it is
worth.  8/8 both word sizes; canary clean; production errors at
n = 1..6 measured 5..72 ulp at 20k points per size.

### one reciprocal per output: the tangent never materializes

The reconstruction used to compute t = TT/DE, square it, and divide
again per output.  With A = TT DE, B = DE^2, C = TT^2, S = B + C the
outputs are

    sin = 2A/S,   cos = 1 - 2C/S,   tan = 2A/(B - C)

-- each ONE reciprocal division (the divisor normalizes to n limbs
with its top bit set by a shift of at most two bits, whose exponent
folds into the final doubling shift) and one mulhigh.  Every special
case (t = 0, t^2 = 0) dissolves; cos comes out exactly 1 through the
generic path.  The doubling shifts amplify the last mulhigh's
truncation by at most 8 (sin/cos) resp. 16 (tan) ulp -- noise against
the 6r + 128 and 8r + 256 budgets.

The fallback's u-division is gone too.  For n beyond the tangent
series, u = tan(t') = ss/cc was one more division; writing
cos t' = 1 - g with g < 2^(-2r-1),

    TT = wy + wx ss - wy g,    DE = wx - wx g - wy ss,

four mulhighs against the tiny ss and g -- cc cancels in every ratio.
No difference can go negative (wy g < wy; DE's two subtrahends are
below wx).  The old chain spent FOUR divisions per sin_cos call at
n = 7 (u, t, one per output); now one, and n = 7 sin_cos went
1020 -> 811 ns on this alone.  _fixed_divide and the dead
sin/cos-then-divide fallback in fixed_tan_bitwise_rs are deleted: the
half-angle path handles every n >= 1 on both word sizes.

### r re-swept (third time), tan series to n = 8, rotation to n = 7

Every cheaper tail moves the optimum: final
r = {4, 5, 9, 14, 15, 18, 16, 16} for n = 1..8.  The sweep's fastest
candidates at n = 5 (r = 17) and n = 7 (r = 15) measured fine at 300
points and BLOW THE TAN BUDGET at 8000 (431 > 392, 440 > 376): the
sweep's error column is an estimate, not a certificate -- always
big-sample validate before locking.  n = 8's register rotation would
need 9-wide add/sub chains that longlong.h doesn't have; it keeps the
generic loop, and the n = 7 -> 8 step shows it (691 -> 937 ns).

### register reconstruction for n = 2..4

dev/gen_fixed_tan_halfangle_small.py -> tan_halfangle_small.inc: the
whole tail in registers -- hand 2x2 sloppy forms at n = 2,
assembly-backed mulhigh dispatch at n = 3, 4, NN_ADD/NN_SUB chains and
unrolled funnel shifts -- except the one mpn_tdiv_qr per reciprocal.
The n = 1 scalar path was rewritten to the same A/B/C form (one
udiv_qrnnd per reciprocal).  n = 2: 157 -> 116 ns.

MEASURED, whole call (one process, ratio = arb / fixed):
    n            1     2     3     4     5     6     7     8
    sin_cos ns  51.2 116.3 240.5 313.9 448.7 558.4 701.4 954.8
    sin_cos     4.97  3.02  1.89  2.00  1.88  1.81  1.70  1.71
    tan         7.52  4.90  2.82  2.90  2.53  2.39  2.23  2.16

(The session opened with n = 7 at 0.85-0.91x.)  Errors at 3000
adversarial points: sin <= 33, cos <= 40, tan <= 70 ulp against
budgets 416/640 at r = 0.  Both word sizes 8/8; canary clean over
200k cases.

### tan series to n = 12: the fallback cliff pushed out

With the reconstruction settled, the drop past n = 8 had two causes:
the tangent series ended there, and fixed_sin_cos_rs's hand routines
end at n = 10, so the fallback paid its generic path from n = 11 --
one process measured sin_cos ratios against arb of 1.43 / 1.37 /
1.10 / 1.15 at n = 9..12, with a +48% step at 10 -> 11.  Emitting
tangent series through n = 12 removes BOTH (the half-angle path stops
touching fixed_sin_cos_rs below n = 13):

    n            8     9    10    11    12    13    14    15    16
    sin_cos ns 875  1028  1210  1388  1621  2650  2774  3015  3345
    sin_cos   1.64  1.62  1.58  1.61  1.66  1.13  1.20  1.22  1.28
    tan       2.02  1.92  1.86  1.88  1.92  1.29  1.36  1.37  1.44

Per-limb steps of 15-20% through n = 12, no cliff; n = 12 -> 13 is
the new fallback boundary (+63%), past which the ratio recovers on
its own as arb's overhead amortizes.  The remaining outsized step,
n = 7 -> 8 (+38%), is the register rotation ending at 7 (9-wide
longlong chains; deliberately left).

r for the new sizes: {16, 19, 23, 25}, big-sample validated as
always -- and the validation again vetoed the sweep's fastest picks:
n = 9 r = 18 reaches 90% of the tan budget at 6000 points, n = 12
r = 22 blows it outright (469 > 432).  The emitted widths stay
honest: 14-16 levels each, trig_rs_opt_hand.inc grows to ~2500
lines.  gen_fixed_trig_hand_mixed.py's add/sub emitter now falls
back to one mpn call for widths past NN_ADD_8/NN_SUB_8, which is
what unlocked n > 8 (their functions are dozens of mulhigh rows;
one dispatched add is noise).

Errors at 2000 adversarial points, n = 9..12: sin <= 42, cos <= 17,
tan <= 59 ulp against budgets 416/640.  Both word sizes 8/8; canary
clean.

### fully specialized atan to n = 7, log1p to n = 6

The register survey after the tan work: every bitwise function already
runs its reduction in registers to n = 7 (the _small variants), but
the FULLY specialized tier -- compile-time r plus the hand series
built for exactly that r, the tier the tan half-angle paths live in --
stopped at exp <= 5, atan <= 4, log1p <= 4.  Measuring the tier delta
per size (r = 0 dispatch against the same call with r passed
explicitly, which forces the _small path):

    n           1    2    3    4    5    6    7
    exp        30%  26%  15%   7%   3%   1%   1%
    atan       34%  42%  23%  23%   -    -    -
    log1p       2%  31%  22%  13%   -    -    -

exp's curve has decayed to nothing by its last specialized size --
extending it buys ~1-3% and is NOT done.  atan's had NOT (23% at
n = 4), and log1p's not quite (13%).  Extended: atan OPT_R += {5: 18,
6: 20, 7: 19}, log1p OPT_R += {5: 31, 6: 30}, with the matching
_fixed_atan_rs_opt_5..7 / _fixed_atanh_rs_opt_5..6 hand series.
log1p at n = 7 measured a 0% tier delta in the sweep proxy and is
left alone; likewise the sweep vetoed nothing this time but the
big-sample rule was applied anyway (all picks <= 30 ulp at 8000
points against budgets of 136-157).

Measured (one process, ratio = arb / fixed):

    n              5      6      7
    atan   345 (1.78) 428 (1.72) 574 (1.61)   was 1.49 / 1.44 / 1.41
    log1p  342 (1.83) 402 (1.82)              was 1.69 / 1.61

The n = 7 -> 8 cliff (+64-73%) remains across exp, log1p and atan for
the same reason as sin/cos/tan: the register _small paths end where
the 8-wide longlong chains do.

Harness note: the sweep drives the SAME generator bodies as
production (gen_one with the series call swapped for a function
pointer), so the sweep's r optimum transfers directly; the residual
divisor in the atan/log1p bodies still carries its unit limb into
mpn_tdiv_qr (the pre-normalization trick from the tan tail has not
been applied there -- a possible next nibble).

### restructure: one file per specialized size, tuned by one script

The dispatch had grown byzantine: each function threaded r = 0 calls
through up to three tiers (fully specialized opt variants, register
_small variants taking a runtime r, and the generic mpn path), the
generated .inc files mixed live and dead series families, and the
2^-16 window series existed for explicit reduction parameters nobody
sensible passes.  Reorganized:

- **The public contract is r = 0 or r >= 32.**  Everything the r < 32
  window family served is covered better by the specialized tier;
  the whole rs16 machinery died (about 3300 lines across
  trig_rs_hard.inc, exp_rs_hard.inc, their wrappers and tables), as
  did the runtime-r register _small variants (three .inc files,
  another 12600 lines) -- explicit r now takes the generic mpn path,
  which is a tuning/testing interface, not a fast path.

- **One source file per (function, size)**: exp_opt_1..7.c,
  log1p_opt_1..7.c, atan_opt_1..7.c, trig_opt_1..12.c (the last
  defines _fixed_trig_opt_n plus the fixed_sin_cos_opt_n and
  fixed_tan_opt_n wrappers -- sin/cos and tan share the half-angle
  reconstruction, so they share the file).  Each file is complete:
  its hand series (static, built for exactly its compile-time r),
  its reduction and reconstruction.  Shared helpers stay shared:
  _fixed_bitwise_reduce, the angle/log tables, _fixed_exp_recon,
  _fixed_tan_halfangle_mid (the generic half-angle body with the
  residual series behind a function pointer; series = NULL takes
  sin/cos of the residual and never divides them), hand_mulhi.inc
  (now include-guarded), tan_rotate.inc.

- **dev/tune_fixed.py generates, builds, runs, selects and emits.**
  For FUNC in {exp, log1p, atan, trig} and a size n it writes an
  out-of-tree source with one fully specialized candidate per r,
  each the SAME body production will use; builds it against the
  in-tree libflint; validates every candidate against MPFR over
  --points adversarial inputs and times it; selects the fastest
  whose error stays within --margin (default 0.5) of the budget,
  breaking near-ties (1.5%) toward the CLEANEST error -- the sweep
  blowout lesson, now institutional; and with --emit writes
  src/fixed/FUNC_opt_n.c.  --pin R skips the sweep; --pin at the
  shipped r reproduces the shipped file byte for byte (verified),
  which is the reproducibility check.  Hand-written pieces that
  predate the series generator (exp n <= 2 series, the log1p n <= 2
  bodies, the trig n = 1 scalar body) are archived in dev/ and
  carried over verbatim by the tuner.

- The old per-function generators remain as body-emitter libraries
  for the tuner; their inc-emitting mains are gone.
  gen_fixed_exp_rs_hard.py and gen_fixed_trig_rs_hard.py now emit
  only the surviving 2^-64/2^-32 families and reproduce the pruned
  .inc files.

Net: 32085 -> 25976 lines in src/fixed (-19%) while ADDING 33 files.
Performance parity verified across every function and size (log1p
n = 7 gained 8% incidentally: its new dedicated r = 25 beats the old
runtime-32 path).  8/8 tests both word sizes at multiplier 25;
canary clean; error validation at 3000+ points per trig size.

### specialization pruned at generation time, not left to the compiler

The per-size opt files carried the full runtime-r bodies with a
const int r on top -- correct, and the compiler eliminated the dead
windows, but the SOURCE still contained them: log1p_opt_7.c had six
guarded windows and a seven-case trailing switch of which only
window 0 and case 0 were reachable at r = 25.  The body emitters now
prune when the reduction parameter is a compile-time constant: only
windows c <= r/64 are emitted (guard-free), the trailing boundary
step collapses to the single live case with its shift count b = r
mod 64 a literal, and window-0 loop bounds resolve numerically.  The
hand-carried log1p n <= 2 bodies in dev/log1p_hand.inc were
specialized the same way (their r >= 64 branches pruned, with the
header noting they are the r = 16 specialization of the originals).

17 emitted files shrank by 3900 lines (atan_opt_7.c: 1000 -> 383);
src/fixed is now 22099 lines, from 32085 before the restructure
began.  Object code is unchanged by construction -- the pruned
branches were compiler-dead -- and the full verification (tests at
multiplier 25 on both word sizes, canary, MPFR validation, benchmark
parity) agrees.

### large-n tuning to r = 768; the profile directory rebuilt

The tuning process now has two documented tiers.  Small sizes:
dev/tune_fixed.py, per (function, size), arbitrary r, emitting the
*_opt_<n>.c files.  Large sizes: src/fixed/tune/tune-bitwise-r,
choosing from the ladder of 32 and the multiples of 64 where the
general series evaluation exists, by binary-searching each
crossover.  The tuner lost its dead r = 16 machinery, gained an
untimed warm-up call before every measurement (the angle/log tables
used to leak into the first pass), and now also covers the
trigonometric functions: sin_cos and tan share the half-angle path,
so they share one new table (_fixed_trig_bitwise_rs_r_tab, tuned on
sin_cos with a tangent cross-check run whose crossovers agree at the
low rungs) -- the generic half-angle path previously HARDWIRED
r = 32 at every size, and the tuned ladder lifts the n = 13+
fallback region visibly (n = 16 sin_cos: 1.32x -> 1.40x over arb).

Freshly generated tables to r = 768 (nmax 640) are installed in all
four dispatch files; exp tops out at r = 448 from n = 343, log1p
uses the full ladder to 768 by n = 578, atan reaches 320 by n = 322,
trig 448 by n = 395.  The r = 0 selection is queryable through
fixed_<func>_bitwise_rs_default_r(n), which the dispatch itself, the
tuner and the profiler all share (the small-n constants are mirrored
there from the *_opt_<n>.c pins, with a keep-in-sync comment).

src/fixed/profile: the two stale per-function programs are replaced
by one p-fixed taking the function as argument, running n = 1..12
and then geometric ~4/3 steps until nmax or until arb wins twice in
a row, printing n, bits, digits, the selected r, both times and the
ratio.  One warm-up call per size keeps precomputation out of the
measurement, and the timing loop cycles over 64 random inputs so the
branchy reductions pay their real misprediction costs.  A full exp
run holds 2.0-4.0x over arb through n = 477 (30528 bits, 9190
digits) without arb ever winning.

### authoritative large-n tables (full-depth run)

The container-tuned tables were replaced by a run on real hardware
with nmax deep enough that every function reaches the r = 768 end of
the ladder (exp crosses at n = 567, log1p at 1004, atan at 1483,
trig at 1450).  The container run had truncated several ladders --
its nmax of 640 hid the upper rungs -- and its noisy timings placed
some crossovers early (exp r = 256: 163 vs 100 on quiet hardware).
Selection spot-checked across all boundaries: the opt-file constants
serve n <= 7 (<= 12 for trig), the ladder takes over exactly at the
tabulated thresholds, and everything saturates at 768.

### two-tier table precomputation

The log(1 + 2^-i) and atan(2^-i) tables were built with one arb call
per entry; per the prototype (newfastexp2.c, prolog_1) they are now
built in two tiers, with the split at about sqrt(precision) bits
(everything through binary splitting above 65536 bits):

- **Small i: binary splitting.**  For the logarithms, i <= 6 factors
  over the primes up to 17, so seven bsplit logs
  (arb_log_primes_vec_bsplit) serve all of them, and beyond that
  log(1 + 2^-i) = 2 atanh(1/(2^(i+1) + 1)) through
  arb_atan_frac_bsplit.  For the arctangents, i <= 3 falls out of
  the first four Gaussian-prime angles {atan 1, atan 2, atan(3/2),
  atan 4} -- atan(1/2) = 2 theta_0 - theta_1, atan(1/4) = 2 theta_0
  - theta_3, atan(1/8) = theta_1 - theta_2 (from 8 + i =
  -i(1+2i)(3+2i)) -- and the direct series bsplit takes over.  TODO
  (noted in both files): reimplement the binary splitting natively
  with mpn arithmetic, and restore the prototype's direct
  log(1 + p/q) series form for very high precision; arb is fine for
  now and this tier is not the bottleneck.

- **Large i: fixed-point mpn multi-summation, directly in the target
  storage.**  Every table entry now carries one guard limb below its
  value limbs (the stride variables _fixed_exp_logs_n /
  _fixed_atans_n count it; consumers, which read the top n limbs of
  each entry, needed NO changes).  The series are summed at the full
  storage width: grouping the log series by the odd part of the
  index makes every term a shift of floor(2^wp / k) for odd k -- and
  that ONE division serves every i at once, so the whole tier costs
  about wp / (2 iter_start) single-limb divisions plus shifted
  additions.  The atan series is purely alternating over odd k and
  is simpler still.

  One-sidedness bug found and fixed on the way: for i past wp/2
  every correction term falls below the guard limb and the raw sum
  lands at exactly 2^-i, violating the strict tab_i < 2^-i invariant
  the reductions rely on (t-log1p caught it at n = 3, r = 176 via
  entry i = 128).  Subtracted terms now round UP (one extra guard
  ulp for the dropped bits) and the leading term carries a blanket
  two-guard-ulp bias covering the sub-guard tail, so every entry is
  the exact floor or one ulp below -- the historical contract.

Fresh-build timings (this container, us): logs/atans at (nv, rc):
(4,32): 687/118 -> 320/52; (64,192): 8642/8804 -> 1313/1424;
(256,448): 106995/119893 -> 15503/12825; (1024,768):
2096931/2470287 -> 211453/156494.  That is 2-16x, growing with
precision.  Verified: table entries against arb (floor or floor-1,
one-sided) to rc = 176; 8/8 tests at multiplier 25 on both word
sizes; canary; MPFR product check; explicit-r path timing unchanged
(852 vs 858 ns, atan n = 3 r = 48, old vs new library).

### prototype bsplit helpers carried over verbatim

arb_atan1_frac_bsplit and arb_log1_frac_bsplit from the newfastexp2
prototype now live in src/fixed/frac_bsplit.inc as static helpers
(verbatim: the only edits are static linkage and a renamed LOG2
macro).  Compared with arb_atan_frac_bsplit they fold pairs of terms
at the leaves and strip two-powers from the denominators as they go.
The logarithm table's bsplit tier uses them exactly as prolog_1 did:
2 atanh(1/(2^(i+1)+1)) through atan1, switching to the direct
log(1 + 2^-i) series of log1 once the precision passes 6000 bits and
i passes 30.  At (nv, rc) = (1024, 768) -- 65536-bit entries, the
all-bsplit regime -- the logs build drops from ~211 ms to ~140 ms
(three-run average).

One reading note: atan1_bsplit ignores its alternate parameter and
its leaves assume p = 1, so the helper is currently valid only for
the hyperbolic p = 1 case it serves; the ALTERNATING atan table
entries therefore stay on the stock arb_atan_frac_bsplit (a comment
in sin_cos_bitwise_rs.c records this).  Extending the leaves with
the (-1)^a signs would let the atan table share the optimization.

Verified: both tables against arb (exact floor or one ulp below,
one-sided) at (3,176) and at (128,100) -- the latter exercises the
log1 path from i = 30 -- plus 8/8 tests at multiplier 25 on both
word sizes and the canary.

### alternating flag fixed in atan1_bsplit; fixed.rst restructured

atan1_bsplit ignored its alternate parameter (fine for the
hyperbolic logarithm use, wrong for atan).  The signs now enter
ABSOLUTELY at the leaves -- the single leaf carries (-1)^a and the
paired leaf becomes (2a+3) qpow2 - (2a+1) with the (-1)^a sign --
so the merges need no sign logic at all; everything else stays
verbatim.  Verified against arb_atan/arb_atanh over i = 1..12 at
200..20000 bits, both modes, power-of-two and odd denominators.
The atan table's bsplit tier now uses the helper in alternating
mode; at 65536-bit entries the build time is unchanged within noise
(for q = 2^i the two-power stripping is nearly free in fmpz, so the
log1 series was where the gain lived), but the tables share one
code path.

doc/source/fixed.rst was restructured: after the introduction, one
section "Elementary functions" with subsections "Series evaluation"
(absorbing the fallbacks) and "Bitwise argument reduction"
(absorbing the five per-function sections and the tuning/profiling
text).  The tangent documentation lost its references to past
iterations of the code ("as before", the comparison against the
conjugate-ratio reconstruction) and now describes only what is
implemented.

### CI fixes: three real bugs, and the reason they were missed

All three CI failures were genuine, and all three were invisible in
this container because THE LOCAL BUILD HAD ASSERTIONS OFF.  FLINT
compiles the library with -DBUILDING_FLINT, which reads src/config.h
(not src/flint-config.h); FLINT_WANT_ASSERT is undefined there, so
every FLINT_ASSERT in the module compiled to nothing all session.
The fix for the process: verification now includes an
assertion-enabled binary, built by compiling the module sources into
the test program,

  gcc -O2 -DFLINT_WANT_ASSERT=1 -I src -I . -I src/fixed \
      -o ttestA src/fixed/*.c src/fixed/test/main.c -lflint ...

which reproduces all three CI failures on the unfixed tree, message
for message.  (A whole-library assertion build also works but takes
~25 minutes here and leaves the shared object inconsistent unless
fully rebuilt.)

1. MinGW link failure ("undefined reference to _fixed_exp_logs").
   The stress test read the cached tables through their thread-local
   pointers, and Windows DLLs cannot export thread-local data.  The
   globals are no longer declared in fixed.h; they live in the new
   library-internal src/fixed/impl.h (the name FLINT uses for
   module-private headers, e.g. src/acb/impl.h; only src/*.h are
   installed), which every module TU includes, and external code goes through two new accessors,
   _fixed_exp_logs_entry(i, n) / _fixed_atans_entry(i, n) (plus
   _..._max_index()), which return the top n limbs of an entry.  The
   emitters emit the new include, so re-emission still reproduces
   the shipped files byte for byte.

2. 64-bit assertion "r == 0 || r >= 32" in fixed_atan_bitwise_rs.
   t-bitwise_rs_stress still drew the trigonometric r from
   16 + n_randint(...) -- the floor from before the r >= 32 contract
   -- so it passed r in 16..31 to sin_cos and atan.  Now 32 + ....
   (t-log1p_bitwise_rs likewise still probed r = 16 in its
   deterministic sweep.)

3. 32-bit assertion in _fixed_at_rs (residual not below 2^-32).  The
   clamp r <= FLINT_BITS n - 16 yields r = 16 at n = 1 on 32-bit,
   below the residual series contract -- exp had carried an
   n >= 2 assertion for exactly this since it was written, but
   log1p, atan, sin_cos and tan never got one, and their tests
   happily ran n = 1.  All four now assert it, the four tests use
   the nmin = (FLINT_BITS == 64) ? 1 : 2 pattern the exp test
   already had, and fixed.rst states the rule once for all the
   bitwise functions instead of only for exp.

Also fixed while here: dev/tune_fixed.py could not actually emit
log1p_opt_1.c / log1p_opt_2.c (--pin 16 --emit raised KeyError; the
hand bodies had been installed by hand and the tuner never learned
to carry them, contradicting the documented contract).  Log1pFunc
now carries them from dev/log1p_hand.inc, which is parametrized with
@FN@/@SERIES@ markers like the trig one, and --pin 16 --emit
reproduces both files.
