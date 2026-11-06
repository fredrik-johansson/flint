---
title: acf.h
header: acf.h
description: Complex floating-point numbers as pairs of arf_t, for approximate computation
category: Real and complex numbers
order: 5
---

An `acf_t` is a complex floating-point number $a + bi$ whose real and
imaginary parts are `arf_t` numbers, without error bounds. It is the
complex counterpart of `arf_t` and is used where complex arithmetic is
needed but error bounds are not: in approximate (non-rigorous)
algorithms such as the floating-point stage of eigenvalue or root
computations that are certified afterwards, and as the element type of the
generic ring of complex floating-point numbers
(`gr_ctx_init_complex_float_acf`).

Addition, subtraction and multiplication are
correctly rounded in each part, like the `arf` operations; division,
square root and dot products are only approximate. All functions allow
aliasing of inputs and outputs.

## Types {#types}

::: type
acf_struct
acf_t
acf_ptr
acf_srcptr

An `acf_struct` has two `arf_struct` fields, the real and imaginary
parts. `acf_t` is an array of one `acf_struct`, passed by reference;
`acf_ptr` and `acf_srcptr` point to arrays.
:::

::: macro
acf_realref(x)
acf_imagref(x)

Pointers to the real and imaginary parts of `x` as `arf_t`.
:::

::: function
arf_ptr acf_real_ptr(acf_t z)
arf_ptr acf_imag_ptr(acf_t z)

Function versions of `acf_realref` and `acf_imagref`, for language
bindings.
:::

## Memory management {#memory}

::: function
void acf_init(acf_t x)
void acf_clear(acf_t x)

Initialise `x` to zero, and free the memory used by `x`.
:::

::: function
acf_ptr _acf_vec_init(slong n)
void _acf_vec_clear(acf_ptr v, slong n)

Allocate a vector of `n` numbers initialised to zero, and free such a
vector.
:::

::: function
void acf_swap(acf_t z, acf_t x)

Swaps `z` and `x` without copying.
:::

::: function
slong acf_allocated_bytes(const acf_t x)

Returns the number of bytes allocated on the heap by `x`, excluding
`sizeof(acf_struct)`.
:::

## Basic operations {#basic}

::: function
void acf_set(acf_t z, const acf_t x)
void acf_neg(acf_t z, const acf_t x)

Set `z` to `x`, respectively $-x$, exactly.
:::

::: function
int acf_set_round(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_neg_round(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)

Set `res` to `x`, respectively $-x$, with both parts rounded to `prec`
bits in the direction `rnd`. Bit 0 of the return value is the inexact
flag of the real part and bit 1 that of the imaginary part.
:::

::: function
int acf_equal(const acf_t x, const acf_t y)

Returns whether `x` and `y` are equal (both parts equal in the sense of
`arf_equal`).
:::

::: function
slong acf_bits(const acf_t x)

Returns the larger of `arf_bits` of the two parts.
:::

::: function
void acf_get_mag(mag_t res, const acf_t x)

Sets `res` to an upper bound for $|x|$.
:::

::: function
void acf_printd(const acf_t x, slong n)

Prints `x` as `a + b*I` with `n` significant decimal digits for each part
(see `arf_printd`).
:::

::: function
void acf_randtest(acf_t x, flint_rand_t state, slong bits, slong mag_bits)

Sets both parts of `x` to random numbers as `arf_randtest` does.
:::

## Arithmetic {#arithmetic}

::: function
int acf_add(acf_t res, const acf_t x, const acf_t y, slong prec, arf_rnd_t rnd)
int acf_sub(acf_t res, const acf_t x, const acf_t y, slong prec, arf_rnd_t rnd)
int acf_mul(acf_t res, const acf_t x, const acf_t y, slong prec, arf_rnd_t rnd)

Set `res` to $x + y$, $x - y$, respectively $xy$, with each part correctly
rounded to `prec` bits in the direction `rnd`. Bit 0 of the return value
is the inexact flag of the real part and bit 1 that of the imaginary part.

::: details
`acf_mul` uses `arf_complex_mul`, or `arf_complex_sqr` when `x` and `y`
are the same object; each part of the result is the exact $ac - bd$ or
$ad + bc$ rounded once.
:::
:::

## Approximate arithmetic {#approx}

The following functions are not correctly rounded: `rnd` is the direction
of the final rounding, but the intermediate roundings are not specified
and no error bound is given.

::: function
void acf_approx_inv(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
void acf_approx_div(acf_t res, const acf_t x, const acf_t y, slong prec, arf_rnd_t rnd)

Set `res` to an approximation of $1/x$, respectively $x/y$. A real or
purely imaginary divisor is handled by real divisions of the parts;
otherwise, above two limbs of precision, the complex kernels of
[](arf#mpn) compute the parts within about an ulp (relative to the
magnitude of the result for the quotient, whose numerator parts may
cancel) at $\lceil \mathit{prec} / B \rceil$ limbs, rounded to `prec`
bits in the direction `rnd`, for operands with moderate exponents;
otherwise $1/(a + bi) = (a - bi)/(a^2 + b^2)$ with $a^2 + b^2$ computed
by `arf_sosq`.
:::

::: function
void acf_approx_sqrt(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)

Sets `res` to an approximation of the principal square root of `x`,
through the complex kernel of [](arf#mpn) for operands with
nonzero parts and moderate exponents.
:::

::: function
void acf_approx_rsqrt(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)

Sets `res` to an approximation of $1 / \sqrt{x}$ (the principal root),
through the complex kernel of [](arf#mpn) for operands with nonzero
parts and moderate exponents, otherwise as the reciprocal of the square
root; `rsqrt(0)` is $+\infty$.
:::

::: function
void acf_approx_dot(acf_t res, const acf_t initial, int subtract, acf_srcptr x, slong xstep, acf_srcptr y, slong ystep, slong len, slong prec, arf_rnd_t rnd)

Sets `res` to an approximation of $s + (-1)^{\mathit{subtract}}
\sum_{i=0}^{\mathit{len}-1} x_i y_i$, where $s$ is `initial` (0 if
`NULL`), with the conventions of `arb_dot`, using a fixed-point
accumulation like that of `acb_approx_dot`.
:::

## Rounded arithmetic and elementary functions {#elem}

The following functions round each part of the result to `prec` bits.
In the strict modes (`ARF_RND_DOWN`, `ARF_RND_UP`, `ARF_RND_FLOOR`,
`ARF_RND_CEIL`, `ARF_RND_NEAR`) each part is correctly rounded, as the
functions of [](arf#elementary) round their results, whatever the size
of the part against the other one. The relaxed modes guarantee the
relative accuracy of the complex number instead: in `ARF_RND_FAST` each
part is within 1 ulp and in `ARF_RND_ACCURATE` within 0.51 ulp of the
larger part (at `prec` bits) of the exact value, so that a part much
smaller than the other may have a worse relative accuracy or be zero.
They are not slowed down by a part that is exponentially small against
the other one, such as the real part of $\cot(2^k i - 1)$, which the
strict modes must resolve. The return value carries the inexact flags
of the parts in the strict modes (bit 0 the real part, bit 1 the
imaginary part; bits 2 and 3 those of a second result); in the relaxed
modes it is not a reliable inexact flag.

::: details
The functions compose the real functions of [](arf#elementary) ($\exp(a + bi)
= e^a (\cos b + i \sin b)$, $\log(a + bi) = \log(a^2 + b^2)/2 + i
\operatorname{atan2}(b, a)$ with $a^2 + b^2$ formed exactly near the
unit circle, $\sin(a + bi) = \sin a \cosh b + i \cos a \sinh b$, and so
on), each part being a product of two accurately known real numbers or
a real function of an exactly known one, so that each part comes out
with its own relative accuracy. In the strict modes the compositions are
evaluated at an increasing working precision until the bits between
the error bound and the rounding point decide the rounding (a Ziv
loop); the relaxed modes take the first evaluation whose error bound is
small enough. At one and two limbs (64-bit, precision up to 128 bits)
the parts come from register compositions of the kernels of the real
functions, shared with `nfloat_complex` and `acb` (`acf/elem.h`).

Exact results are a concern for the strict modes, since no window of
bits ever decides the rounding of a dyadic value: the compositions
detect the exact zeros of the sine and cosine of $\pi x$ at the
half-integers (the imaginary part of $(-2)^3$, say), and `acf_pow`
computes integer powers (and small dyadic powers of exact roots) exactly
when the result has at most $2^{24}$ bits, and the powers of $2^j(\pm 1
\pm i)$ for any integer exponent. A dyadic value the tests miss (as
$(2i)^{1/2} = 1 + i$) would make the retries fail to converge: beyond a
working precision of $4 \mathit{prec} + 65536$ bits, the approximation
is rounded as it stands, without the guarantee of correct rounding.
:::

::: function
int acf_exp(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_exp_pi_i(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_log(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_sin(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_cos(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_sin_cos(acf_t res1, acf_t res2, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_sin_pi(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_cos_pi(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_sin_cos_pi(acf_t res1, acf_t res2, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_sinh(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_cosh(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_sinh_cosh(acf_t res1, acf_t res2, const acf_t x, slong prec, arf_rnd_t rnd)

Set `res` (or the pair `res1`, `res2`) to the exponential, $e^{\pi i x}$,
the principal logarithm, the sine and cosine (of $x$ and of $\pi x$) and
the hyperbolic sine and cosine of `x`. A nan or infinite part gives nan
parts, as do the arguments beyond the ranges of the real functions (the
exponential beyond $2^{\max(128, 2 \mathit{prec})}$ in the real part,
the trigonometric functions beyond $2^{\max(65536, 4 \mathit{prec})}$);
`log(0)` is $-\infty$. Real and purely imaginary arguments are computed
by the real functions (with $\log(-x) = \log x + \pi i$, $\log(ix) =
\log|x| \pm \pi i/2$).
:::

::: function
int acf_pow(acf_t res, const acf_t x, const acf_t y, slong prec, arf_rnd_t rnd)

Sets `res` to $x^y = \exp(y \log x)$ with the principal logarithm.
$x^0 = 1$ (including $0^0$); $0^y = 0$ for $\operatorname{Re}(y) > 0$
and nan otherwise. Real powers of positive reals are computed by
`arf_pow`; the powers of real and purely imaginary numbers with real
exponents (or with $|x| = 1$) from $|x|^{\operatorname{Re} y}$ and the
sine and cosine of $\pi k \operatorname{Re}(y) / 2$, so that their
exact values (such as $(-2)^3$ and $i^{1000}$) are recognised; see
the details above for the other exact cases.
:::

::: function
int acf_tan(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_cot(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_sec(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_csc(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_tan_pi(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_cot_pi(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_sec_pi(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_csc_pi(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_tanh(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_coth(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_sech(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_csch(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)

Set `res` to the tangent, cotangent, secant and cosecant of `x`
(respectively of $\pi x$) and to the hyperbolic tangent, cotangent,
secant and cosecant of `x`. Real arguments are computed by the real
functions (whose results at the poles are kept, with a zero imaginary
part); a nan or infinite part gives nan parts.

::: details
With $s = \sin a$, $c = \cos a$, $\mathit{sh} = \sinh b$ and
$\mathit{ch} = \cosh b$ (of $\pi a$ and $\pi b$ for the $\pi$ versions),
$\tan(a + bi) = (s c + \mathit{sh}\,\mathit{ch}\, i) / (c^2 +
\mathit{sh}^2)$, $\cot(a + bi) = (s c - \mathit{sh}\,\mathit{ch}\, i) /
(s^2 + \mathit{sh}^2)$, $\sec(a + bi) = (c\, \mathit{ch} + s\,
\mathit{sh}\, i) / (c^2 + \mathit{sh}^2)$ and $\csc(a + bi) = (s\,
\mathit{ch} - c\, \mathit{sh}\, i) / (s^2 + \mathit{sh}^2)$: the
denominators $|\cos z|^2$ and $|\sin z|^2$ are sums of squares, so that
each part is a quotient of accurately known numbers, also when it is
exponentially small against the other (the real part of $\tan(a + bi)$
is about $4 s c e^{-2|b|}$ for large $|b|$). The imaginary part of $\tan$
and $\cot$ is then $\pm(1 + \delta)$ with $\delta = -(\cos 2a +
e^{-2|b|}) / (2(c^2 + \mathit{sh}^2))$ resp. $(\cos 2a - e^{-2|b|}) /
(2(s^2 + \mathit{sh}^2))$, from which the directed modes are rounded
when $|\delta|$ is far below an ulp; for example, $\cot(2^k i - 1)$
takes about a microsecond at 53 bits in every mode, whatever $k$. The
hyperbolic functions are computed as $\tanh z = -i \tan(iz)$, $\coth z
= -i \cot(-iz)$, $\operatorname{sech} z = \sec(iz)$ and
$\operatorname{csch} z = -i \csc(-iz)$.
:::
:::

::: function
int acf_expm1(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_log1p(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_sinc(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_sinc_pi(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_asin(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_acos(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_atan(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_acot(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_asec(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_acsc(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_asin_pi(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_acos_pi(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_atan_pi(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_asinh(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_acosh(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_atanh(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_acoth(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_asech(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_acsch(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)

Set `res` to $e^x - 1$, $\log(1 + x)$, $\sin(x)/x$, $\sin(\pi x)/(\pi
x)$, and the principal inverse trigonometric and hyperbolic functions
of `x` (and the inverse trigonometric functions divided by $\pi$), with
the branch cuts of [](acb#inverse-trig) and [](acb#inverse-hyperbolic); the reciprocal functions are
$\operatorname{acot}(z) = \operatorname{atan}(1/z)$ and so on. At the
poles and the singular points (such as $\operatorname{atanh}(1)$ and
$\operatorname{acot}(0)$) the parts are nan.

::: details
These functions round the balls of the corresponding functions of
`acb`, evaluated at increasing working precisions: in the strict modes
until both endpoints of the ball of each part round to the same number,
in the relaxed modes until the radii lie below $2^{-\mathit{prec} - 4}$
resp. $2^{-\mathit{prec} - 8}$ times the larger part. Off the axes the
parts are transcendental at dyadic arguments, so that the strict loop
terminates. Real arguments in the real domain of a function are
computed by the real function of `arf`; on the axes outside it, the
parts that are exactly zero or a multiple of $\pi/2$ (of $1/2$ for the
$\pi$ versions) are set directly (for example $\operatorname{acos}(x)$
has the real part 0 for $x > 1$ and $\pi$ for $x < -1$). In the strict
modes an argument very close to an axis ($z = z_0 + \eta$ with
$|\eta| < 2^{-\mathit{prec} - 64} |z_0|$, $z_0$ on the axis and away from
the branch cuts) is evaluated from $f(z_0) + f'(z_0) \eta$ with the
remainder bounded by the Cauchy estimates, so that the small part of
the value comes out with its own relative accuracy at the usual working
precision: $\operatorname{asin}(1/2 + 2^{-100000} i)$ takes about two
microseconds at 53 bits.
:::
:::

::: function
int acf_inv(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_div(acf_t res, const acf_t x, const acf_t y, slong prec, arf_rnd_t rnd)

Set `res` to $1/x$, respectively $x/y$, with each part correctly rounded
(the relaxed modes as described above): $1/(a + bi) = (a -
bi)/(a^2 + b^2)$ and $(a + bi)/(c + di) = ((ac + bd) + (bc - ad) i) /
(c^2 + d^2)$ with the numerators and the denominator formed exactly, so
that each part is one correctly rounded quotient of exact numbers; when
the terms of a sum lie more than $2 \mathit{prec} + 128$ bits apart (so
that they cannot cancel) the sum is rounded to a working precision
instead and the quotient certified as for the elementary functions.
A real or purely imaginary divisor is handled by real divisions of the
parts; a zero divisor gives the results of `arf_div` (infinities and
nan).
:::

::: function
int acf_sqrt(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_rsqrt(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)

Set `res` to the principal square root of `x`, respectively its
reciprocal, with each part correctly rounded (the relaxed modes as
described above). The dyadic roots (and their reciprocals, when
$|x|$ is a power of two) are computed exactly; otherwise the complex
kernels of [](arf#mpn) at an increasing working precision, as for the
elementary functions (the compositions of real operations for extreme
exponents). Real and purely imaginary arguments are handled by the real
roots; `rsqrt(0)` is $+\infty$.
:::
