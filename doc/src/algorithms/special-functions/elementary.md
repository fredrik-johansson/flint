---
title: Elementary functions
summary: How the elementary functions of arf, nfloat and arb share the fixed-point kernels of mp_real - the floating-point entry points of the kernels, the register fast paths at one and two limbs, correct rounding by certification of a fixed-point frame (arf), balls from a frame with propagation through derivative bounds (arb) - and how the derived functions are composed.
complexity: that of the `mp_real` kernels, $O(M(n) \log^2 n)$ or better at $n$ bits, plus a rounding; strict rounding modes of `arf` retry with more guard bits with probability about $2^{-20}$ per retry.
order: 10
---

The real elementary functions of `arf`, `nfloat` and `arb` -- the
exponential, the logarithm, sine, cosine and tangent, the arctangent,
and the functions derived from them -- are all evaluated by the
fixed-point kernels of `mp_real`, described in
[](algorithms/floating-point/elementary-functions). This page describes
the shared layer between those kernels and the three number types: how a
floating-point argument enters a kernel, how a fixed-point result with an
error count becomes a correctly rounded number or a ball, and the short
cuts for arguments that fit in one or two words.

## The kernels and their floating-point entry points {#kernels}

A kernel of `mp_real` evaluates a function of a fixed-point argument in
$[0, 1)$, $n$ limbs long, and returns a fixed-point *frame*: an array
of limbs $Y$ and an error count $\varepsilon$ such that the true value
lies within $\varepsilon$ units of $Y 2^{e}$, where $2^e$ is the unit of
the frame (a power of $2^{-\bits}$, possibly shifted by the reduction).
For a floating-point argument $x = \pm d \cdot 2^{e_x}$ with an
msb-aligned mantissa $d$ of $n_x$ limbs, the entry points
`_mp_real_exp_fp`, `_mp_real_log_fp`, `_mp_real_sin_cos_fp`,
`_mp_real_atan_fp`, `_mp_real_atan2_fp`, `_mp_real_sinh_cosh_fp`,
`_mp_real_tanh_fp` and `_mp_real_pow_fp` (`mp_real.h`) perform the
argument reduction of that floating-point number and call the kernels
at a requested number of limbs $n_k$:

- The exponential reduces $x = q \log 2 + t$ with an exact integer $q$,
  returning the frame of $e^t$ and the binary shift $q$; a quotient
  beyond the word range is returned as a multiprecision shift by
  `_mp_real_exp_fp_huge`, so that $e^x$ is computed as long as the
  caller can represent the exponent.
- The logarithm writes $x = 2^{e} (1 + t)$ and returns the frame of
  $\log(1 + t) + e \log 2$ (for $\log(1 + x)$ of a small $x$, the frame
  is relative to $x$ so that no accuracy is lost to cancellation).
- Sine and cosine reduce $x$ modulo $\pi/4$ with the cached floor of
  $\pi/4$ (and of $2/\pi$) to the precision needed by the exponent of
  $x$, with the octant selecting which of $\pm \sin$, $\pm \cos$ is
  which; the tangent is the quotient of the two frames.
- The arctangent of $|x| > 1$ is $\pi/2 - \operatorname{atan}(1/|x|)$
  with a reciprocal to the working precision; `_mp_real_atan2_fp`
  forms the quotient of two mantissas first and handles the octants of
  $\operatorname{atan2}$ (a tiny quotient is computed exactly, with its
  remainder, so that the bounds are one-sided).
- The hyperbolic functions come from $e^{|x|}$ and its reciprocal, the
  tangent as a quotient; $x^y$ is $\exp(y \log x)$ with the logarithm at
  the precision increased by the size of $y \log x$, the exact product
  of the frame with the mantissa of $y$, and the exponential of the
  result (with the codes for $x^y$ rounding to 1, to a huge or tiny
  number).

The number of limbs is chosen by the caller from the target precision
and a few guard bits: $n_k = \lceil (p + g) / \bits \rceil$, with $g =
12$ (`MP_REAL_ELEM_GUARD_BITS`) for `nfloat` and `arb`, and $g = 20$ or
$24$ for `arf` ([](#rounding)). Each entry point returns its error
count with the frame, so that every caller performs its own rounding
step with a rigorous bound; the error counts are small (tens of units)
and are dominated by the reduction's use of the tabulated constants.

Arguments near 0 bypass the kernels: $e^x$, $\log(1 + x)$, $\sin x$,
$\operatorname{atan} x$ and so on are evaluated by a Taylor polynomial of
degree one, two or three in register arithmetic with a known sign of
the remainder (`arf/elem.h`), which gives a one-sided error bound. This
preserves the relative accuracy of results near 0 and, for `arf`, lets
the rounding of a value such as $x (1 - \delta)$ with $0 < \delta <
x^2/3$ be decided without a loop (`_arf_set_tiny_rel_err`).

## Register fast paths {#fast}

On 64-bit machines, arguments with a mantissa of one or two limbs at a
target precision of one or two limbs are evaluated by the functions
`_arf_small_exp_1`, `_arf_small_exp_2`, `_arf_small_log_1`,
`_arf_small_sin_cos_1`, `_arf_small_atan_1`, `_arf_small_sinh_cosh_1`
and their siblings in `arf/elem.h`, shared by `arf`, `nfloat` and
`arb`. They work on three- or four-word fixed-point numbers with the
`double`-based reductions and the dense tables of `mp_real` (see
[](algorithms/floating-point/elementary-functions#small)), return the
frame in registers with an error count of a few units, and fall back to
the kernels (return value nonzero) for arguments outside their range.
These paths give the elementary functions of `nfloat` and of exact
`arb` balls their speed at 64 and 128 bits, and the same code serves
`arf` with its rounding layer on top.

## Correct rounding in arf {#rounding}

`arf_exp`, `arf_log`, `arf_sin_cos`, `arf_atan`, ... return the function
value rounded to `prec` bits in the direction `rnd` (see
[](arf#elementary)). In the strict modes (`ARF_RND_DOWN`, `ARF_RND_UP`,
`ARF_RND_FLOOR`, `ARF_RND_CEIL`, `ARF_RND_NEAR`) the result is correctly
rounded, which requires more than a bound on the error: it requires
knowing on which side of a rounding boundary the true value lies.

The certification step `_arf_set_mpn_err` (`arf/elem.c`) receives a
frame $Y 2^e$ with lower and upper error counts $\varepsilon_-$ and
$\varepsilon_+$ (the value lies in $[(Y - \varepsilon_-) 2^e, (Y +
\varepsilon_+) 2^e]$; one of them may be zero for a one-sided bound)
and attempts to round. Let $D$ be the bits of $Y$ between the rounding
point for `prec` bits and the position of the error bound. The rounding
is decided when $D$ is neither all zeros nor all ones (for the directed
modes; for rounding to nearest, the same test on the bits past the
rounding bit): then both ends of the enclosure round to the same number.
If the window is ambiguous, the function returns `ARF_ELEM_UNCERTAIN`
and the caller retries with more guard bits.

::: algorithm Rounding by retries (Ziv's loop) {#alg-ziv}
**Input:** $x$, `prec`, a strict rounding mode.

**Output:** $f(x)$ correctly rounded.

1. Set $\mathit{extra} = 0$ and $n_k = \lceil (\mathit{prec} + 24 +
   \mathit{extra}) / \bits \rceil$.
2. Evaluate the entry point of the kernel at $n_k$ limbs and certify the
   frame. If the rounding is decided, return.
3. Set $\mathit{extra} = 32$ if it was 0, else double it, and go to 2.
:::

With $24$ initial guard bits the first attempt fails with probability
about $2^{-20}$ for a value whose bits are random, and each retry
squares the chance of failure. Values whose binary expansion is not
random near the rounding point -- exact results such as $e^0$, $\log 2^k$
or $2^{1/2}$ raised to 2, and values $x (1 + \delta)$ with a tiny
$\delta$ -- would make the loop run for ever with two-sided bounds, so
the entry points and the drivers detect them: exact cases are returned
directly (`_arf_pow_exact` for powers), and the tiny-argument,
near-one and huge-argument expansions use one-sided bounds, which decide
a rounding whenever the value is not exactly a representable number.

In the relaxed modes the first attempt is final: `ARF_RND_ACCURATE`
rounds the frame to nearest, with an error below $0.51$ ulp from the
$20$ guard bits, and `ARF_RND_FAST` truncates the frame plus its upper
error count, within one ulp. In these modes the inexact flag is set
whenever a kernel was evaluated; the exact cases (such as $e^0$, $\log
1$ and the powers of two of `arf_exp2` at integers) are detected before
any evaluation and reported as exact in every mode.

The quotient variants `_arf_set_mpn_err_quot` (for the tangents and the
reciprocal functions, where the result is a frame divided by an exact
mantissa or its reciprocal, at the working precision of the current
attempt) and `_arf_set_tiny_recip_rel_err` (for $(1/x)(1 + \delta)$,
which decides by one exact division with remainder) round composite
values without forming a floating-point intermediate.

Results with a bignum exponent (the exponential of a large argument,
powers) are rounded at exponent 0 and shifted by the exact integer
afterwards; the elementary functions of `arf` support the same range of
exponents as those of `arb`, with the same magnitude limits ($|x| >
2^{\max(128, 2 \,\mathit{prec})}$ for the exponential and $|x| >
2^{\max(65536, 4 \,\mathit{prec})}$ for the trigonometric functions give
NaN rather than a result whose cost would not be polynomial in `prec`).

## Balls in arb {#ball}

A point evaluation such as `arb_exp_arf`, `arb_log_arf`, `arb_atan_arf`,
`_arb_sin_cos_arf` or `_arb_sinh_cosh_arf` converts a frame into a ball
with `_arb_set_mpn_err` (`arb/elem.c`): the midpoint is $Y 2^e$
rounded down to `prec` bits and the radius is $\max(\varepsilon_-,
\varepsilon_+) 2^e$ plus one ulp of the midpoint. With the $12$ guard
bits of the kernel the radius stays within a few ulps, so the result
has close to `prec` bits of accuracy.

A ball function such as `arb_exp` receives $[m \pm r]$ and proceeds as
follows.

1. Special values (zero, infinities, NaN, an infinite radius) are handled
   directly, and so are arguments beyond the magnitude limits: `arb_exp`
   returns $[0 \pm \infty]$ or a tiny ball when $|m| > 2^{\max(128,
   2p)}$, and `arb_sin_cos` returns $[0 \pm 1]$ when $|m| > 2^{\max(65536,
   4p)}$ (see [](guide/types/ball-arithmetic#precision)).
2. If $r > 0$, the precision is lowered to about the relative accuracy of
   the input plus 20 to 30 bits; more bits would be hidden by the
   propagated error. If the input has less than about 20 bits of
   accuracy, the function is evaluated on the interval instead, using
   its monotonicity (for $\sin$ and $\cos$, the endpoints and the extrema
   in between, in `double` arithmetic: `_arb_sin_cos_wide`; for the
   hyperbolic functions `arb_sinh_cosh_wide`).
3. Otherwise the function of the exact midpoint is computed by the point
   evaluation, and the propagated error is added.

The propagated error uses a bound on the derivative over the ball
([](algorithms/ball-arithmetic/error-propagation)), which the point
evaluation provides nearly for free: `_arb_sin_cos_arf` returns upper
bounds of $|\sin m|$ and $|\cos m|$ read from the top bits of the two
frames (`_mag_set_frame`), so that the error of $\sin$ is bounded by
$r |\cos m| + r^2/2$ and that of $\cos$ by $r |\sin m| + r^2/2$; the
hyperbolic functions use $r \cosh m \,(1 + 2r)$ for $r \le 1/2$ from the
frame of $\cosh m$, and the exponential $e^m (e^r - 1)$ with
`mag_expm1`. For the functions of $\pi x$ the bounds are multiplied by
$\pi$. The logarithm bounds the error of $\log [m \pm r]$, $r < m$, by
$\log(1 + r/(m - r))$ and the arctangent by $r / (1 + (|m| - r)^2)$.

## Constants {#constants}

The constants $\pi$ and $\log 2$ at up to about 4000 bits are read from
the static tables of `mp_real` (the floors of $\pi/4$ and $\log 2$ at
64 limbs, `_mp_real_const_ptr`), rounded to nearest with half an ulp
plus the table's error as the radius; beyond, `arb_const_pi` and
`arb_const_log2` use the series of `mp_real` with a per-thread cache
([](algorithms/special-functions/constants)). The same floors serve the
argument reductions of the kernels, `nfloat_pi`, and `arf_const_pi`
(correctly rounded by certification of the floor at increasing
precision).

## Derived functions {#derived}

- `arb_expm1` and `arb_log1p` use the entry points with the subtraction
  of 1 built in, at a precision increased by the bits that would be lost
  to cancellation for small arguments.
- `arb_sinh_cosh`, `arb_tanh` and `arb_coth` are computed from the
  hyperbolic entry points; `arb_asinh`, `arb_acosh` and `arb_atanh` from
  `arb_log`, `arb_log1p` and `arb_sqrt`.
- `arb_asin` and `arb_acos` use $\operatorname{atan}(x / \sqrt{1 -
  x^2})$, and `arb_atan2` combines `arb_atan` with $\pi$.
- `arb_sin_pi` and `arb_cos_pi` reduce $x$ modulo $1/2$ exactly (by
  subtracting the integer nearest to $2x$) before multiplying by $\pi$, so
  they are exact at integers and half-integers and accurate near the
  zeros. `arb_sin_cos_pi_fmpq` evaluates small denominators as algebraic
  numbers, by refining a root of the minimal polynomial of $\cos(2\pi/q)$
  with Newton iteration.
- `arb_pow` with exact arguments and a positive base (or an integer
  exponent) uses `_mp_real_pow_fp` in one evaluation; otherwise it is
  $\exp(y \log x)$ in ball arithmetic, with the special cases of
  [](arb#powers).
- The complex functions of `acb` (`acb_exp`, `acb_log`, `acb_sin_cos`,
  `acb_atan`, ...) are built from the real ones: for example $e^{a + bi}
  = e^a (\cos b + i \sin b)$ and $\log z = \log|z| + i \arg z$ with
  `arb_log_hypot` and `arb_atan2`.
- In `arf`, the compositions (`arf_asin`, `arf_acosh`, `arf_log10`,
  `arf_sinc_pi`, ...) evaluate the inner functions with the relaxed
  rounding on temporaries of a few more limbs and bound the error of the
  composition explicitly, then certify the final rounding (`arf/elem_misc.c`);
  the one-sided expansions near the singular points (such as
  $\operatorname{acos} x$ near $x = 1$, or $\operatorname{acoth} x$ for
  huge $x$) avoid the stalls of two-sided bounds.
- Logarithms of smooth integers at high precision are assembled from
  cached logarithms of primes (`arb_log_ui_from_prev`,
  `arb_log_primes_vec_bsplit`).

## Use in FLINT {#use}

::: functions arf_exp arf_log arf_sin_cos arf_atan arf_pow arb_exp arb_log arb_sin_cos arb_atan arb_exp_arf arb_log_arf arb_atan_arf _mp_real_exp_fp _mp_real_log_fp _mp_real_sin_cos_fp _mp_real_atan_fp _mp_real_pow_fp :::
