---
title: Roots and polynomials in castle arithmetic
summary: How castle arithmetic computes square roots and powers of elements, gcds of polynomials (subresultants), factorizations over towers (Trager's algorithm with shifted norms) and roots of polynomials — formulas, recognized small roots, exact factorization over the generators the coefficients involve, no-roots certificates at places, adjunction — and how splitting towers emerge from Cauchy moduli when conjugate roots meet.
complexity: Roots of a polynomial of degree $n$ over a tower of degree $D$ cost a factorization of a norm of degree $nD$ within the limits (64, or 24 with transcendental generators), and otherwise one expression search per root; the $k$-th Cauchy modulus of a splitting tower of a polynomial of degree $n$ has $\binom{n+1}{k}$ terms.
order: 37
---

This page describes how the castle fields ([](gr_ca)) compute roots of
elements and of polynomials, and the polynomial algorithms behind them.
The generators they create (radicals, roots of unity) follow the
policies of [](algorithms/exact/field-towers#structured); the
algorithms on towers they use (expression searches, Trager's method,
dynamic evaluation) are described in [](algorithms/exact/field-towers)
and [](algorithms/exact/dynamic-evaluation).

## Square roots and radicals of elements {#radicals}

Roots of rational numbers are canonical products of roots of primes and
roots of unity of prime power orders
([](algorithms/exact/field-towers#radicals)), so that, for instance,
$\sqrt{-5}\sqrt{-7} = -\sqrt{35}$ is found by merging the towers
$\mathbb{Q}(\sqrt{-5})$ and $\mathbb{Q}(\sqrt{-7})$, whose compositum
contains $\sqrt{35}$. Each root of prime power order is the generator
of a tower of its own registered in the context
(`gr_tower_adjoin_root_of_unity`), $i$ being $\exp(2\pi i/4)$.

Square roots of elements whose flat numerator and denominator are
squares (or minus squares) as polynomials are read off
($\sqrt{(1-\pi)^2} = \pi - 1$, $\sqrt{-(\pi+1)^2} = i(\pi + 1)$);
otherwise square roots are found exactly through quadratic steps
(`gr_tower_sqrt`, [](algorithms/exact/field-towers#express)) before any
lattice search, repeated roots of the same element are identified, and
roots of likely powers are found by factoring (see
`GR_TOWER_OPT_POWER_CHECK_DEGREE_LIMIT`). A positive rational content is
taken out of an $n$-th root ($\sqrt{4\pi} = 2\sqrt\pi$), and the root of
a monomial in roots of unity and roots of integers is a monomial in such
roots of higher orders ($\sqrt{4\sqrt 2} = 2 \cdot 2^{1/4}$). A root
which is not found in the tower of its argument is adjoined to it as a
new step, proven when Capelli's theorem can be certified at a place
([](algorithms/exact/dynamic-evaluation#capelli)) and dynamic otherwise.

## Powers and logarithms {#powers}

`gr_pow` handles rational exponents through roots and other exponents
as $\exp(y \log x)$. Logarithms of rational and Gaussian rational
numbers are decomposed over the logarithms of primes and of Gaussian
primes ([](algorithms/exact/field-towers#logarithms)), and $\log(c
\exp(u)) = \log c + u$ is taken from the definition when
$\operatorname{Im} u \in (-\pi, \pi]$. `gr_set_str` parses expressions
with `^`: `examples/huge_expr.c -castle` evaluates two algebraic numbers
given by huge expressions in nested radicals (degree 64) and proves them
equal in about 0.06 s, against 3 s for the Calcium version of the same
example and 4 s for the `qqbar` version.

## Greatest common divisors {#gcd}

The generic Euclidean algorithm is exponential over fields such as
$\mathbb{Q}(t)$: the coefficients of the remainders are rational
functions whose degrees double at each step. The castle fields and the
tower field contexts therefore register the subresultant polynomial
remainder sequence [@Col1967; @Duc2000] for `gr_poly_gcd`,
`gr_poly_xgcd` and `gr_poly_resultant`; a gcd of two degree-six
polynomials over $\mathbb{Q}(\pi, e)(\sqrt 2, \sqrt 3)$ with a common
quadratic factor went from not finishing to a few milliseconds. For
polynomials over the quotient rings themselves (the nested
representation), a subresultant method was tried and dropped: it was
3.6 times faster for random polynomials over quadratic and quartic
fields, but twice as slow on the small gcds that dominate the
radical-formula benchmarks, so the Euclidean algorithm remains the
default there. Modular gcd algorithms over number fields [@LM1989;
@Enc1995; @vHM2002] are not implemented
([](algorithms/exact/castle-design#open)).

*Example (Euclid against subresultants).* Over $\mathbb{Q}(\pi, e)(\sqrt
2, \sqrt 3)$, let $q = x^2 + \sqrt 2 x + \pi$ and

$$
f = q \, \bigl(x^4 + (\pi + \sqrt 3)x^3 + ex^2 + \sqrt 2 x + 1\bigr), \qquad
g = q \, \bigl(x^4 + (e - \sqrt 2)x^3 + \pi x + \sqrt 3\bigr).
$$

In the complex castle field, $\gcd(f, g)$ by the subresultant sequence
returns $q$ in about 0.03 s. The Euclidean remainders of degrees 5, 4, 3,
2 have largest coefficients of 76, 2963, 26002 and 609723 characters in
printed form, and the last division, which has to recognize that the
next remainder is zero, takes about a minute, as does
`gr_poly_gcd_euclidean`; in the monic gcd it returns, the coefficient of
$x$, which is $\sqrt 2$, is a fraction of $2.4 \cdot 10^5$ characters
whose denominator involves $\sqrt 2$ and $\sqrt 3$: an inverse beyond the
size limit of rationalization
([](algorithms/exact/castle-representations#canonical)), valid but not
canonical, so that its equality with $\sqrt 2$ is decided by the zero
test.

## Factorization {#factorization}

Polynomials over $F_k$ are factored by Trager's algorithm [@Tra1976],
for towers over $\mathbb{Q}$ and over $\mathbb{Q}(t_1, \ldots, t_r)$:

1. the steps up to $k$ are proven (modular places, then Trager's method
   on the steps themselves, refining reducible steps;
   [](algorithms/exact/dynamic-evaluation#certificates));
2. the polynomial is made squarefree by Yun's algorithm [@Yun1976] with
   subresultant gcds (`gr_poly_factor_squarefree`, which had to be
   taught that quotient rings pretending to be fields are unique
   factorization domains);
3. for each squarefree part $h$, the norm of $h(x - \theta)$ is computed
   by iterated resultants (`gr_poly_quotient_norm_poly`), with $\theta =
   \sum_j s_j a_j$ taken from a fixed sequence of twelve shifts until
   all norm factors are simple; the norm is factored by
   `fmpz_poly_factor` over $\mathbb{Q}$ or, after clearing denominators,
   as a polynomial in $t_1, \ldots, t_r, x$ by `fmpz_mpoly_factor`;
4. the factors are $\gcd(h(x - \theta), N_i)(x + \theta)$, taken by
   increasing degree of $N_i$, each divided out of $h(x - \theta)$ before
   the next, the last one being the cofactor (3.4 times faster than
   independent gcds on $\Phi_{24}$).

The work is limited by a bound on $\deg(h) [F_k : F_0]$
(`GR_TOWER_OPT_FACTOR_DEGREE_LIMIT`). During development, $\Phi_{24}$
over $\mathbb{Q}(\sqrt 2, \sqrt 3, i)$ split into eight linear factors
in 28 ms, $x^{12} - 4$ over $\mathbb{Q}(2^{1/6})$ into factors of degrees
1, 1, 2, 2, 2, 4 in 4 ms, and $x^4 + \pi^4$ over
$\mathbb{Q}(\pi)(\sqrt 2)$ into two quadratics in 1 ms. Factorization is
exposed through the `gr` polynomial methods, so that generic code can
factor over towers.

*Example (shifts for $x^3 - 2$ over $\mathbb{Q}(2^{1/3})$).* The norms
of [](algorithms/exact/dynamic-evaluation#trager) show the mechanism on
$x^4 + 1$ over $\mathbb{Q}(\sqrt 2)$, where the second shift succeeds.
For $h = x^3 - 2$ over $\mathbb{Q}(a)$, $a = \sqrt[3]{2}$, the roots of
$h(x - \theta)$ for $\theta = sa$ are $a(\zeta_3^i + s\zeta_3^j)$, and the
norm is squarefree only if these nine numbers are distinct. For one
generator the sequence of shifts begins $\theta = 0, -a, -a, -a, a,
-2a$ (it repeats $-a$), with norms

$$
(x^3 - 2)^3, \quad x^3(x^6 + 108), \quad (x^3 + 2)^2(x^3 - 16), \quad (x^3 + 2)(x^6 + 40x^3 + 1372)
$$

for $\theta = 0, -a, a, -2a$; the first three are rejected by the
squarefree test before any factorization. With $\theta = -2a$ the norm
$x^9 + 42x^6 + 1452x^3 + 2744$ is factored over $\mathbb{Q}$;
$\gcd(h(x + 2a), x^3 + 2) = x + a$, which shifted back is $x - a$, and
the cofactor gives $x^2 + ax + a^2$. Over $\mathbb{Q}(i)$ the same
algorithm gives $x^4 + 1 = (x^2 - i)(x^2 + i)$, and over
$\mathbb{Q}(\zeta_8)$ four linear factors $x \mp \zeta_8$, $x \mp
\zeta_8^3$.

## Roots {#strategies}

Root finding (`gr_poly_roots`) combines several strategies, ordered by
cost.

1. *Rational and structured roots.* Linear factors over $\mathbb{Q}$,
   quadratics by the formula with structured square roots of the
   discriminant, binomials $ax^n + b$ by principal radicals times roots
   of unity (so that the roots of $x^3 - 2$ are $\sqrt[3]{2}\zeta_3^j$),
   and factors over $\mathbb{Q}[t_1, \ldots, t_r]$ when the coefficients
   involve transcendental generators only (the characteristic polynomial
   of $\log M$ for a matrix $M$ has the roots $\log a$, $\log b$ and
   $-\log a - \log b$, found by multivariate factoring).
2. *Guesses.* Small algebraic roots are recognized by `qqbar_guess` and
   verified exactly (a quadratic one as a square root). A root of a
   degree-$n$ polynomial over a field of degree $D$ is generically of
   degree $nD$ over $\mathbb{Q}$; guesses are limited to degrees below
   $nD$ (and at most six), since the exact factorization finds generic
   roots more cheaply.
3. *Exact factorization over the involved generators.* The coefficients
   are mapped to their common tower, the generators they involve are
   marked (with the closure of their dependencies), and when that is not
   the whole tower the roots are sought over a new subset tower. Without
   this, $x^5 - \sqrt 2 x - 1$ solved after three other polynomials over
   $\mathbb{Q}(\sqrt 2)$ built a splitting tower of degree 5760 (40 s);
   with it, a tower of degree 240 in 15 ms. The factorization is used
   when the steps of the tower are proven and the norm degree is at most
   64, or 24 with transcendental generators
   (`GR_TOWER_OPT_ROOTS_FACTOR_DEGREE_LIMIT`).
4. *A no-roots certificate.* Before the lattice searches that try to
   express a root in the existing tower
   ([](algorithms/exact/field-towers#express)), a modular place is sought
   at which the reduction of the polynomial has no root in the residue
   field; by [](algorithms/exact/dynamic-evaluation#prop-place) the
   polynomial then has no root in the tower. Transcendental generators are
   specialized to random residues at the place.
5. *Adjunction.* An irreducible factor of degree $\ge 3$ gets a new
   algebraic step (proven, when the factorization was exact), and the
   cofactor is factored again over the extension. Beyond the limits of
   the exact methods, each root is expressed in the tower when it lies
   there ($\pi$ and $e$ for $(x - \pi)(x - e)$) and adjoined as a dynamic
   step otherwise, the linear factor being divided out before the next
   root is sought and a remaining quadratic factor being solved by the
   formula.

Multiplicities are returned as usual. `gr_factor` of a polynomial gives
the linear factors $x - r$ for the roots, and over the real fields the
quadratic factors $(x - r)(x - \overline r)$ for the pairs of nonreal
roots ([](algorithms/exact/castle-real-fields#subfields)).

*Example (roots of a product).* In a fresh complex castle field, the
nine roots of $(2x - 1)(x^2 - 2)(x^3 - 2)(x^3 - x - 1)$ come out as

$$
\tfrac12, \quad \pm\sqrt 2, \quad a_2, a_3, a_4, \quad a_5, \; a_5 a_6, \; -a_5a_6 - a_5 \quad (a_5 = \sqrt[3]{2}, \; a_6 = \zeta_3),
$$

where $a_2, a_3, a_4$ are the roots of $x^3 - x - 1$ near $1.32472$ and
$-0.662359 \pm 0.562280i$, each adjoined in a tower of its own with the
modulus $x^3 - x - 1$, and $-a_5a_6 - a_5 = \sqrt[3]{2}\zeta_3^2$. Over a
fixed tower the roots are the linear factors of the factorization: over
$\mathbb{Q}(\sqrt 2, \sqrt 3)$, $x^4 - 10x^2 + 1$ has the roots $\pm\sqrt
2 \pm \sqrt 3$.

## Splitting towers {#splitting}

The roots of an irreducible factor of degree $d \ge 3$ over $\mathbb{Q}$
are individual algebraic numbers, each adjoined with its *origin*
polynomial. When conjugate roots meet in an operation, the tower grows
by the *Cauchy moduli*: the $k$-th root is a root of $p(x)/\prod_{i<k}(x
- \alpha_i)$, computed by synthetic division (`gr_poly_div_root`)
without certification of irreducibility, since an enclosure isolating a
root of $p$ isolates a root of the quotient; the last root is obtained
from Vieta's formula. Symmetric functions of the roots then reduce to
the coefficients. The moduli have polynomially many terms in $d$ for
fixed $k$ ($\binom{d+1}{k}$), and they are converted to the flat
representation only when a reduction needs them
([](algorithms/exact/field-towers#representations)). The product $a_1
\cdots a_d$ cascades through every modulus, because $a_d = c -
\sum_{i<d} a_i$ produces $a_{d-1}^2$, whose reduction brings in
$h_2(a_1, \ldots, a_{d-2})$, and so on down.

With this, `examples/hilbert_matrix_ca.c -castle` verifies the trace
and the determinant of the Hilbert matrix against its eigenvalues for
$n = 16$ in 0.4 s, $n = 20$ in 10 s and $n = 22$ in 45 s, where Calcium
needs its Vieta option to do $n = 10$ ([](tutorials/exact/linear-algebra)).
Beyond that, the size of the normal forms modulo the Cauchy tower is the
limit, in either representation. (The reductions in the nested
representation use ordinary division by the moduli; with the Newton
inverse that the generic quotient rings would otherwise choose, $n = 20$
took 138 s.)

*Example (splitting towers of $x^3 - x - 1$ and of a Hilbert matrix).*
Let $a_1, a_2, a_3$ be the roots of $p = x^3 - x - 1$ near $1.32472$ and
$-0.662359 \pm 0.562280i$, as above. The sum $a_1 + a_2$ lives in the
tower

$$
a_1 : x^3 - x - 1 \;(\text{proven}), \qquad a_2 : p(x)/(x - a_1) = x^2 + a_1 x + a_1^2 - 1 \;(\text{dynamic}),
$$

of degree 6, in which $a_3 = -a_1 - a_2$. Then $a_1 + a_2 + a_3 = 0$ and
$a_1a_2a_3 = 1$, and

$$
(a_1 - a_2)(a_1 - a_3)(a_2 - a_3) = 6a_1^2a_2 + 2a_1 - 2a_2 + 3,
$$

which castle arithmetic proves equal to $\sqrt{-23}$ (its square reduces
to $-23$, the discriminant of $p$). The Galois group of $p$ is $S_3$, so
the quadratic modulus is in fact irreducible. For the $4 \times 4$
Hilbert matrix $H_4$, the characteristic polynomial $x^4 -
\frac{176}{105}x^3 + \frac{3341}{12600}x^2 - \frac{41}{23625}x +
\frac{1}{6048000}$ is irreducible with Galois group $S_4$; its roots
$\lambda_1, \ldots, \lambda_4$ generate a Cauchy tower of moduli of
degrees 4, 3, 2 (degree 24, the degree of the splitting field), and
$\lambda_1 + \cdots + \lambda_4$ and $\lambda_1 \cdots \lambda_4$ reduce
to the rationals $\operatorname{tr} H_4 = 176/105$ and $\det H_4 =
1/6048000$.

## Conjugates {#conjugates}

The complex conjugate of an algebraic generator is adjoined as a
conjugate root of the generator's origin polynomial, or, for a generator
that is not an absolute algebraic number (a root of $x^5 - \pi x - 1$),
as the root of the modulus with conjugated coefficients whose enclosure
meets the conjugate enclosure; the operation is restarted when this
changes the tower. Realness, ordering, real and imaginary parts follow
([](algorithms/exact/castle-real-fields#conj)).

*Example (conjugates).* For the roots of $x^3 - x - 1$ above,
$\overline{a_2}$ is the existing generator $a_3$, a root of the origin
polynomial, and $\operatorname{Re} a_2 = (a_2 + a_3)/2$ is computed in
the tower $\mathbb{Q}(a_2)(a_3)$ with the Cauchy modulus $x^2 + a_2x +
a_2^2 - 1$; it compares equal to $-a_1/2$. The roots of $x^5 - \pi x - 1$
over $\mathbb{Q}(\pi)$ form a splitting tower of degree 120: $a_1 \approx
-0.319367$ and $a_2 \approx -1.23581$ with proven moduli of degrees 5 and
4, $a_3 \approx 0.0769670 + 1.34249i$ with a dynamic cubic modulus, and
$a_4 = \sqrt{-3a_1^2 - 2a_1a_2 - 3a_2^2 - 2a_1a_3 - 2a_2a_3 - 3a_3^2}$
from the last quadratic, solved by the formula. The fourth and fifth
roots are $(\pm a_4 - a_1 - a_2 - a_3)/2$, and $\overline{a_3} = (-a_4 -
a_1 - a_2 - a_3)/2$. The real roots, sorted, are $a_2 < a_1 < (a_4 - a_1
- a_2 - a_3)/2 \approx 1.40124$; the last is real although it is written
with the nonreal generators $a_3$, $a_4$. The real view computes and
sorts them in a few hundredths of a second.

## Printing and reading roots {#printing}

A root of a polynomial which is not a structured generator prints as
`root(f, approx)`, with $f$ written out and the decimal approximation
nearest to the root, with as many digits as are needed to separate it
from the other roots. For $f = x^5 - 2(10^6x - 1)^2$, which has two
roots $10^{-6}(1 \pm 7.07\ldots \cdot 10^{-16})$, these print with 17
significant digits (`1.0000000000000007e-6` and
`9.9999999999999929e-7`), while the other three roots print with 6; all
five round trip through strings, and convert structurally to another
context as five distinct elements.

## Use in FLINT {#use}

- `gr_sqrt`, `gr_rsqrt`, `gr_pow` and `gr_root_of_unity` (through the
  `gr` methods), `gr_poly_roots`, `gr_poly_factor_squarefree` and
  `gr_factor` of the castle fields implement these methods
  (`src/gr_ca/radical.c`, `src/gr_ca/roots.c`).
- The tower layer provides `gr_tower_sqrt`, `gr_tower_poly_factor`,
  `gr_tower_poly_roots` and `gr_tower_adjoin_root_of_unity`;
  `gr_poly_div_root` and `gr_poly_quotient_norm_poly` are the generic
  primitives behind the Cauchy moduli and the norms.
- The options `GR_TOWER_OPT_POWER_CHECK_DEGREE_LIMIT`,
  `GR_TOWER_OPT_FACTOR_DEGREE_LIMIT`,
  `GR_TOWER_OPT_ROOTS_FACTOR_DEGREE_LIMIT`,
  `GR_TOWER_OPT_SAME_ORIGIN_DEGREE_LIMIT` and the expression limits bound
  the effort ([](gr_tower#tuning-options)).
