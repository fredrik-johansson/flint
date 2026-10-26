---
title: Arithmetic and precision
functions: arb_add arb_sub arb_mul arb_div arb_inv arb_addmul arb_fma arb_fmma arb_sosq arb_set_round_mpn_2exp_si
---

The arithmetic functions are the building blocks of every computation
with balls. Each computes the result of the operation on the midpoints
with one `arf` operation rounded to `prec` bits, and sets the radius to a
bound for the propagated error plus the rounding error.

## What is computed {#formulas}

For $x = [m_x \pm r_x]$ and $y = [m_y \pm r_y]$:

| Operation | Midpoint | Radius before rounding |
|---|---|---|
| `arb_add`, `arb_sub` | $m_x \pm m_y$ | $r_x + r_y$ |
| `arb_mul` | $m_x m_y$ | $\lvert m_x\rvert r_y + \lvert m_y\rvert r_x + r_x r_y$ |
| `arb_div` | $m_x / m_y$ | $\dfrac{\lvert m_x\rvert r_y + \lvert m_y\rvert r_x}{\lvert m_y\rvert (\lvert m_y\rvert - r_y)}$ |
| `arb_addmul` | $m_z + m_x m_y$ | $r_z + \lvert m_x\rvert r_y + \lvert m_y\rvert r_x + r_x r_y$ |
| `arb_fmma` | $m_a m_b + m_c m_d$ | the two product terms added |
| `arb_sosq` | $m_a^2 + m_b^2$ | likewise |

For `arb_add`, `arb_sub`, `arb_mul`, `arb_addmul`, `arb_submul`,
`arb_fma`, `arb_fmma`, `arb_sosq`, `arb_div`, `arb_inv`, `arb_sqrt`,
`arb_rsqrt`, `arb_set_round`, `arb_neg_round` and `arb_set_str` (and the
variants taking `arf`, integer or `fmpz` operands), the midpoint is
rounded to `prec` bits to nearest and the actual rounding error (bounded
to within $2^{-\mathrm{MAG\_BITS}}$ of an ulp) is added to the radius, by
`arb_set_round_mpn_2exp_si`; when the propagated error exceeds the
rounding error by more than `MAG_BITS` bits, the midpoint is rounded at
a correspondingly lower precision, since further bits would not be
meaningful. The quotients and roots are computed with a guard limb
(`flint_mpn_divhigh`, `flint_mpn_sqrthigh`, `flint_mpn_rsqrthigh`, see
[](mpn_extras#fphigh)), one more unit of $2^{-\mathrm{MAG\_BITS}}$ ulp
covering the error of the approximation; an exact quotient or root
(detected when the guard limb is near a limb boundary) is found exactly.
The remaining operations round the midpoint towards zero (`ARB_RND`) and
add one unit in the last place to the radius when the rounding was
inexact. All radius computations round upward.
The bounds are proved in
[](algorithms/ball-arithmetic/error-propagation); for addition and
multiplication of real balls they are the best possible up to the
rounding of the radius.

## Precision and exactness {#precision}

The precision bounds the length of the midpoint of the *result*. It is
not a property of the variables: an `arb_t` holds whatever midpoint the
last operation produced, and the inputs of an operation can be longer or
shorter than `prec`. Three consequences:

- If the exact result has at most `prec` bits, it is returned exactly,
  with no error added. Integers, dyadic fractions and sums and products of
  them stay exact as long as they fit, which is what makes exact input
  cheap at any precision.
- A rounded result has an error of at most one unit in the last place at
  `prec` bits, so its relative accuracy is about `prec` bits, *unless* the
  input radii are larger. The radius of the result is never smaller than
  the propagated input error, whatever the precision.
- Addition and multiplication do not lower the precision for inexact
  input: `arb_mul` of two balls with 20 correct bits at `prec` = 1000
  computes a 1000-bit midpoint. `arb_div`, the elementary functions and
  `arb_dot` do lower it, to the accuracy of the input plus about 30 bits.
  To save time and memory in long computations with inaccurate data, use
  a lower precision or `arb_trim`.

`prec` may be `ARF_PREC_EXACT` for addition, subtraction and
multiplication to get the exact result, if it is known to fit in memory.

```c example name=arb_arithmetic_exactness
#include <flint/flint.h>
#include <flint/arb.h>

int main(void)
{
    arb_t x, y, z;
    arb_init(x); arb_init(y); arb_init(z);

    /* exact input, exact result: no error */
    arb_set_ui(x, 3);
    arb_set_ui(y, 4);
    arb_mul(z, x, y, 2);
    flint_printf("3 * 4 at 2 bits:  "); arb_print(z); flint_printf("\n");

    /* exact input, inexact result: one rounding error */
    arb_add(z, x, y, 2);
    flint_printf("3 + 4 at 2 bits:  "); arb_print(z); flint_printf("\n");
    arb_div(z, x, y, 53);
    flint_printf("3 / 4 at 53 bits: "); arb_print(z); flint_printf("\n");
    arb_set_ui(y, 7);
    arb_div(z, x, y, 53);
    flint_printf("3 / 7 at 53 bits: "); arb_printn(z, 20, 0); flint_printf("\n");

    /* inexact input: the propagated error dominates */
    arb_set_str(x, "[3 +/- 1e-10]", 53);
    arb_mul(z, x, x, 53);
    flint_printf("x * x:            "); arb_printn(z, 20, 0); flint_printf("\n");

    /* division by a ball containing zero */
    arb_set_str(y, "[0.1 +/- 0.2]", 53);
    arb_div(z, x, y, 53);
    flint_printf("x / y:            "); arb_printn(z, 20, 0); flint_printf("\n");

    arb_clear(x); arb_clear(y); arb_clear(z);
    flint_cleanup();
    return 0;
}
```
```output
3 * 4 at 2 bits:  (3 * 2^2) +/- (0)
3 + 4 at 2 bits:  (1 * 2^3) +/- (536870912 * 2^-29)
3 / 4 at 53 bits: (3 * 2^-2) +/- (0)
3 / 7 at 53 bits: [0.4285714285714285 +/- 7.15e-17]
x * x:            [9.000000000 +/- 6.01e-10]
x / y:            nan
```

$12 = 3 \cdot 2^2$ fits in 2 bits and is exact, while 7 does not and
becomes $[8 \pm 2]$ (rounded to nearest), a correct but coarse
enclosure. $3/4$ is exact and
$3/7$ is rounded. The square of $[3 \pm 10^{-10}]$ has radius $6 \cdot
10^{-10} + 10^{-20}$, whatever the precision.

## Division and special values {#division}

If the divisor contains zero, the quotient is indeterminate
($[\operatorname{NaN} \pm \infty]$, printed as `nan`), even when the
dividend is 0 and even if only an endpoint of the divisor is 0: the set
of possible quotients is unbounded or undefined. Code that divides by a
quantity that may be close to zero should test `arb_is_nonzero` first or
be prepared for this result. For a divisor whose relative accuracy is
below about 20 bits, `arb_div` computes the quotient from the endpoints,
which is tighter than the first-order bound.

Infinite midpoints are handled as limits: $[+\infty] + [1 \pm 1] =
[+\infty \pm 1]$, which represents $+\infty$, and $[+\infty] \cdot [2
\pm 1] = [+\infty]$, while $[+\infty] \cdot [1 \pm 1]$ is the whole
extended real line $[+\infty \pm \infty]$. $\infty - \infty$ and $0 \cdot
\infty$ give a NaN midpoint.

For exact input the radius is zero exactly when the result is exact: an
exact quotient or the square root of a perfect square gets radius 0,
and every other result a radius of about half an ulp. Above a few
thousand bits the kernels switch to the Newton iterations of
[mp_real](mp_real) internally.

## Aliasing and fused operations {#aliasing}

Any output may be the same variable as any input. `arb_addmul(z, x, y,
prec)` computes $z + xy$ with a single rounding of the midpoint, which is
both faster and more accurate than a multiplication followed by an
addition; `arb_fma` is the same with a separate output. For sums of many
products, `arb_dot` is better still: it rounds only once in total and
lowers the working precision to what the input justifies.
`arb_sub(z, x, x, prec)` and `arb_mul(z, x, x, prec)` treat the two
arguments as independent balls: $x - x$ is $[0 \pm 2r]$, not 0 (see
[](guide/types/ball-arithmetic#conventions)).
