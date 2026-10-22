---
title: Real and algebraic fields in castle arithmetic
summary: How the real and algebraic castle fields restrict the complex one — four fields as views of one shared state, membership checks where values enter, sign tests; real forms of elements involving roots of unity and the tangent normal form of real trigonometric constants; real trigonometric functions through tangents and arctangents with a real angular relation search; real-first merging; and conjugation, real and imaginary parts.
complexity: A realness check costs one enclosure, and an exact test $\overline{x} = x$ when the imaginary part cannot be separated from zero; the tangent normal form at level $M$ costs arithmetic modulo the minimal polynomial of $\tan(\pi/M)$, of degree $\varphi(c)/2$.
order: 38
---

Castle arithmetic has four fields: the complex field of all the numbers
the towers can represent, its real subfield, and the complex and real
algebraic numbers ([](gr_ca#contexts)). This page describes how the
subfields are enforced and how their numbers are represented in real
terms.

## Four fields, one implementation {#subfields}

`gr_ctx_init_ca` creates the complex field; flags restrict it to the
real field, the field of algebraic numbers $\overline{\mathbb{Q}}$, or
the real algebraic numbers. These are not separate implementations but
*views* (`gr_ctx_init_ca_view`) of one shared state: the registry of
towers, the caches and the lock are shared and reference counted, and a
view adds a membership condition. Converting an element between views
of one state is a copy followed by a membership check: exact realness
for the real fields, and, for the algebraic fields, that all generators
the element involves (with their definitions) are algebraic.

In a real context, an operation whose result would not be real fails
with `GR_DOMAIN`: $i$, roots and logarithms of negative numbers, inverse
trigonometric functions outside their real domains, conversions of
nonreal numbers. Nonreal polynomial roots are dropped by `gr_poly_roots`
(for rational polynomials, only the real roots are isolated: 0.1 s
rather than 1.4 s for $x^{100} - 2(101x - 1)^2$, which has 4), and
`gr_factor` gives quadratic factors for the pairs of nonreal roots;
since all elements are real, the tests are sign tests. In an algebraic
context, $\pi$ is not available and the elementary functions only at
the arguments where their values are algebraic ($\exp(0)$, $\log(1)$,
$\cos(0)$, $1^y$, $x^{p/q}$, and so on — decided by the
Lindemann–Weierstrass and Gelfond–Schneider theorems), the others failing
with `GR_DOMAIN`. The restrictions apply to the *outermost* operation
only (a depth counter in the locking layer), not to the intermediate
values of the implementations: $\operatorname{acos}(1/2)$ in the real
field may pass through complex logarithms and still return $\pi/3$.

Following the `gr` convention that operands are valid elements of their
context, a real view does not check that its inputs are real: `gr_floor`
of $2 + i$ evaluated in the real view returns 2, which is the value FLINT
defines for the floor of a complex number (it acts on the real part).
Checking realness on every input would make every field operation pay
for an exact test.

*Example (one state, four views).* With $C$ a complex castle field and
$R$, $A$, $A_{\mathbb{R}}$ its real, algebraic and real algebraic views,
the same calls give (– for `GR_DOMAIN`):

| | $C$ | $R$ | $A$ | $A_{\mathbb{R}}$ |
|---|---|---|---|---|
| $\sqrt{-1}$ | $i$ | – | $i$ | – |
| $(-8)^{1/3}$ | $2\zeta_3 + 2$ | – | $2\zeta_3 + 2$ | – |
| $\operatorname{acos}(1/2)$ | $\pi/3$ | $\pi/3$ | – | – |
| $\log(-1)$ | $\pi i$ | – | – | – |
| $e^0$, $1^{\sqrt 2}$ | 1 | 1 | 1 | 1 |
| $\sqrt 2^{\sqrt 2}$ | $t_2 = \exp(\sqrt 2\,t_1)$, $t_1 = \log\sqrt 2$ | same | – | – |

The principal cube root of $-8$ is $1 + \sqrt{-3} = 2\zeta_3 + 2$, which
is not real. In $A$, $\operatorname{acos}(1/2)$ fails because its value
$\pi/3$ is transcendental. Converting $2 + i$ from $C$ to $R$ fails the
membership check, but `gr_floor` called with the context $R$ on the
element $2 + i$ of the shared state returns 2.

## Views {#views}

A view shares the towers, registry, generator names and lock of its
parent (and hence its elements), restricted further by its flags
(`GR_CA_REAL`, `GR_CA_ALGEBRAIC`, added to those of the parent). The
shared state is reference counted: the parent and its views may be
cleared in any order. `gr_set_other` between contexts of the same state
(`gr_ca_ctx_same_state`) copies the element and checks that it lies in
the target field; between different states it transfers the element
structurally (the generators it involves are recreated in the target
state from their definitions, with the roots matched by their
enclosures, so that even roots closer than their printed approximations
keep their identity). The real field is thus the real subfield of a
complex field whose computations it shares: its values are computed
through the complex numbers where that is the natural way, and only the
values entering it are checked.

Membership of the subfield is enforced where values enter it from
outside the operations: conversions, `gr_set_str` (whose definitions may
introduce auxiliary nonreal generators, as long as the value is real)
and `gr_gens`, which lists only the definitions lying in the subfield (in
a real field, not $i$ or $\exp(i)$). The realness test is numerical when
the imaginary part is separated from zero, and otherwise the exact test
$\overline{x} = x$ ([](#conj)).

## Real forms {#real-forms}

A real element whose representation involves roots of unity $\zeta_j =
e^{2\pi i/n_j}$ and otherwise only real generators is converted to a
*real form*: the roots of unity are cleared from the denominator
(through the nested representation if needed), and $x = \operatorname{Re}
x = \sum A \cos(2\pi \sum_j e_j/n_j)$ is computed over the terms $A
\prod_j \zeta_j^{e_j}$ of the numerator, the cosines being written in the
tangent normal form ([](#tangent-form)). This is applied to elementary
function results, real roots of polynomials, factor coefficients and
values entering a real field: $\cos(\pi/5) = (\sqrt 5 + 1)/4$,
$\sin(\pi/12) = \sqrt{2 - \sqrt 3}/2$, and the real factors of $x^4 + 1$
are $x^2 \mp \sqrt 2 x + 1$.

*Example (real forms).* In the complex field, $\zeta_8 + \zeta_8^{-1}$ is
the element $\zeta_8 - \zeta_8^3$; converted to the real view it becomes
$2\cos(\pi/4) = \sqrt 2$, a quadratic level in square root form.
Likewise $\zeta_9 + \zeta_9^{-1} = \zeta_9 - \zeta_9^2 - \zeta_9^5$ becomes

$$
2\cos\tfrac{2\pi}{9} = \frac{\tau^4 - 34\tau^2 + 29}{16}, \qquad \tau = \tan\tfrac{\pi}{9}.
$$

With a nonreal denominator: $q = (\zeta_8 + \sqrt 3)/(\zeta_8 - \sqrt 3)$
is rationalized in the complex field to $-(\sqrt 3\zeta_8^3 + 3\zeta_8^2 +
3\sqrt 3\zeta_8 + 4)/5$, and $q + \overline q = 2\sqrt 3(\zeta_8^3 -
\zeta_8)/5 - 8/5$ becomes $-(2\sqrt 2\sqrt 3 + 8)/5$ in the real view.

## The tangent normal form {#tangent-form}

For angles in $(\pi/D)\mathbb{Z}$ the values of $\sin$, $\cos$ and $\tan$
lie in the maximal real subfield of $\mathbb{Q}(\zeta_c)$, with $c = 4D$
for odd $D$ and $c = 2D$ for even $D$. Let $M$ be the smallest positive
integer with $\operatorname{lcm}(M, 4) = c$ ($M = D$ for odd $D$, $D/2$
for $D \equiv 2 \bmod 4$, $2D$ for $4 \mid D$), and $\tau = \tan(\pi/M)$.
Since

$$
\frac{(1 + i\tau)^2}{1 + \tau^2} = \frac{1 + i\tau}{1 - i\tau} = e^{2\pi i/M},
$$

the field $\mathbb{Q}(i, \tau)$ contains $e^{2\pi i/M}$ and $i$, hence all
roots of unity of order dividing $c$, and $\tau \in \mathbb{Q}(\zeta_c)$.
Hence $\mathbb{Q}(i, \tau) = \mathbb{Q}(\zeta_c)$, and since $\tau$ is
real and $i \notin \mathbb{Q}(\tau)$, $\mathbb{Q}(\tau) =
\mathbb{Q}(\zeta_c) \cap \mathbb{R}$. A real element written as a
polynomial in roots of unity is thus a polynomial in $\tau$, computed
modulo the minimal polynomial of $\tau$ with Gaussian arithmetic (binary
powering of $1 + ix$, sparse Horner evaluation over the nonzero terms).
The generator $\tan(\pi/M)$ is a definition kind of its own
(`GR_TOWER_TAN_PI`, printed `tan(pi/M)`); quadratic $\tau$ use their square
root forms ($\tan(\pi/3) = \sqrt 3$, $\tan(\pi/8) = \sqrt 2 - 1$,
$\tan(\pi/12) = 2 - \sqrt 3$), and other values in a quadratic subfield
also go through $\tau$ ($\cos(\pi/5) = (7 - \tau^2)/8$), so that each
level has exactly one generator. Tangents are computed at the level of
the doubled angle, $\tan(r\pi) = \sin(2r\pi)/(1 + \cos(2r\pi))$, whose
field is smaller ($\tan(\pi/12) = 2 - \sqrt 3$, not through
$\tan(\pi/24)$). This is the algebraic counterpart of the half-angle
substitution used for the trigonometric functions at other real
arguments ([](#trig)), and the real counterpart of the principle "one
number field, one generator"
([](algorithms/exact/field-towers#radicals)). The form is used by the
real trigonometric functions at rational multiples of $\pi$
(denominators up to 240), by the conversion of real elements to real
forms ([](#real-forms)), and by the sines, cosines and cotangents that
occur in the special function formulas: $\Gamma(3/5) = \pi(15\tau -
\tau^3)/(10\,\Gamma(2/5))$, $\psi'(4/5) = (15 - \tau^2)\pi^2/5 -
\psi'(1/5)$ with $\tau = \tan(\pi/5)$, and $\psi(1/7)$ through
$\tan(\pi/7)$ and the logarithms of its polynomials
([](algorithms/exact/castle-special-functions#gamma)).

*Example (levels 5, 7 and 12).* For $D = 5$, $c = 20$ and $M = 5$: $\tau
= \tan(\pi/5)$ has the minimal polynomial $x^4 - 10x^2 + 5$, and the real
field gives

$$
\cos\frac{\pi}{5} = \frac{7 - \tau^2}{8}, \qquad \sin\frac{\pi}{5} = \frac{7\tau - \tau^3}{8}, \qquad \sin\frac{2\pi}{5} = \frac{11\tau - \tau^3}{8},
$$

which are $(1 + \sqrt 5)/4$, $\sqrt{(5 - \sqrt 5)/8}$ and $\sqrt{(5 + \sqrt
5)/8}$, since $\tau^2 = 5 - 2\sqrt 5$. For $D = 7$, $c = 28$, $M = 7$, and
$\tau = \tan(\pi/7)$ has degree $6 = \varphi(28)/2$ (minimal polynomial
$x^6 - 21x^4 + 35x^2 - 7$); then $\cos(2\pi/7) = (\tau^4 - 22\tau^2 +
25)/32$, $\sin(\pi/7) = (\tau^5 - 20\tau^3 + 19\tau)/16$, and the product
$\sin\frac\pi7 \sin\frac{2\pi}7 \sin\frac{3\pi}7$ is computed in the tower
of the single generator $\tau$ as $(\tau^5 - 22\tau^3 + 49\tau)/64$, which
compares equal to $\sqrt 7/8$. For $D = 12$, $c = 24$ and $M = 24$:
$\sin(\pi/12) = (\tau^2 + 8\tau + 1)/8$ with $\tau = \tan(\pi/24)$,
whereas $\tan(\pi/12) = 2 - \sqrt 3$ is computed at the level of the
doubled angle.

## Real trigonometric functions {#trig}

Values such as $\cos 1$ or $\operatorname{atan} 2$ are not expressible
through roots of unity; in the complex field they are written through
$e^i$ and $\log(1 \pm 2i)$. Real views use real transcendental
generators instead: $\tan(u)$ for real $u$ (with $\sin x = 2q/(1 + q^2)$
and $\cos x = (1 - q^2)/(1 + q^2)$ for $q = \tan(x/2)$) and
$\operatorname{atan}(u)$; $\arcsin y = \arctan(y/\sqrt{1 - y^2})$ and
$\arccos = \pi/2 - \arcsin$. Odd functions are normalized to positive
arguments, and arguments affine in $\pi$ with rational coefficient are
split by the addition formulas with the algebraic values at $r\pi$, so
that no generator has $\pi$ in its argument (otherwise $\tan((\pi/3 -
a)/2)$ produced a relation with coefficient 3 and a degree-six dynamic
step, and the identity $\sin 3a = 4\sin a\sin(\pi/3 - a)\sin(\pi/3 + a)$
was undecidable). $\arctan(y)$ is recognized as $r\pi$ for a rational
$r$ with a small denominator when $\tan(r\pi) = y$ (verified exactly).

*Example ($\cos 1$ and $\operatorname{atan} 2$ in both fields).* In the
real view, $\sin 1 = 2q/(1 + q^2)$ and $\cos 1 = (1 - q^2)/(1 + q^2)$ with
the single generator $q = \tan\frac12$, and $\operatorname{atan} 2$ is a
generator. In the complex field, $\cos 1 = (t^2 + 1)/(2t)$ with $t =
e^i$, and $\operatorname{atan} 2 = \frac\pi2 + \frac i2\bigl(2\log(2 + i) -
\log 5\bigr)$ in the Gaussian logarithm basis of
[](algorithms/exact/field-towers#logarithms).

### The real angular relation search {#angular}

Without $i$ in the tower, Richardson's search
([](algorithms/exact/zero-testing#richardson)) runs on the *angles* $A =
2u$ of $\tan u$ and $A = 2\operatorname{atan} u$ of $\operatorname{atan}
u$, together with $\pi$ (included numerically when absent, and adjoined
in front when a candidate needs it; this is what makes
$\operatorname{atan}(\tan 2)$ print as $2 - \pi$ at once). Verification:
with a tangent on top, an additive zero test below (the angles live
there); with an arctangent on top, the product of the points $((1 -
q^2)/(1 + q^2), 2q/(1 + q^2))$ is checked to be $(1, 0)$ by pair
arithmetic, together with the branch condition $|\sum_j c_jA_j| < 1$.
Eliminations:

- an arctangent on top becomes linear: $\theta = -\sum_{j \ne g}
  c_jA_j/(2c_g)$;
- a tangent with coefficient $\pm 1$ becomes $S/(1 + C)$ for the point
  $(C, S)$ of the other angles;
- a tangent with $|c_g| > 1$ is swapped with another tangent of unit
  coefficient when the order allows; otherwise a basis of the lattice of
  angles is computed (Hermite normal form of the scaled coordinate
  matrix, with the unimodular transform inverted), new generators
  $\tan(B_r/2)$ are inserted before the tangents involved, and these
  become rational functions of them, so that the degree of the tower does
  not grow. The last resort is the degree-$2|c|$ real polynomial $C
  \operatorname{Im}((1 + ix)^{2c}) - S\operatorname{Re}((1 + ix)^{2c})$ as a
  dynamic step.

The basis step terminates because each insertion is paid for by a
decrease in the number of transcendental generators, unless an
arctangent is involved, in which case new tangents of old arguments are
refused. Without the basis step, random rational arguments (whose
tangents are all algebraic over each other) produced degree-12 dynamic
steps and multi-megabyte norms.

*Example (angular relations).* For $\operatorname{atan}\frac12 +
\operatorname{atan}\frac13$ the angles are $A_1 = 2\operatorname{atan}
\frac12$ and $A_2 = 2\operatorname{atan}\frac13$, with the relation
$2A_1 + 2A_2 - \pi = 0$. Its verification multiplies the points
$(\frac35, \frac45)$ and $(\frac45, \frac35)$ of $A_1$ and $A_2$ as
complex numbers: $(3 + 4i)^2(4 + 3i)^2/625 = -1$, the point of $\pi$. The
arctangent on top becomes the linear step $\operatorname{atan}\frac13 =
\pi/4 - \operatorname{atan}\frac12$, and the sum is $\pi/4$; in the same
way Machin's formula makes $\operatorname{atan}\frac{1}{239} =
4\operatorname{atan}\frac15 - \pi/4$ a linear step. For tangents:

- $\tan\frac12$ then $\tan 1$ (relation $A_2 - 2A_1 = 0$, unit
  coefficient on top): $\tan 1 = 2q/(1 - q^2)$ with $q = \tan\frac12$.
  Created in the opposite order, $\tan\frac12$ is on top with coefficient
  2; the two tangents are swapped, with the same result.
- $\tan\frac23$ and $\tan\frac12$: the angles $\frac43$ and 1 satisfy
  $3A_1 - 4A_2 = 0$, and neither coefficient is $\pm 1$. The lattice
  $\frac43\mathbb{Z} + \mathbb{Z}$ has the basis $\frac13$, the generator
  $\tau = \tan\frac16$ is inserted, and

  $$
  \tan\tfrac23 = \frac{4\tau(1 - \tau^2)}{1 - 6\tau^2 + \tau^4}, \qquad \tan\tfrac12 = \frac{3\tau - \tau^3}{1 - 3\tau^2},
  $$

  so that the tower keeps one transcendental generator and degree 1.
- $\operatorname{atan}(\tan 2)$: the angles 4 of $\tan 2$ and
  $2\operatorname{atan}(\tan 2)$, with $\pi$, give
  $\operatorname{atan}(\tan 2) = 2 - \pi$.

A complex field computing with the same state keeps its exponential
forms; relations between the two forms are found when both occur in a
tower (which then contains $i$).

::: remark Real forms do not mix with complex forms {#rem-forms}
Real forms must not be mixed into the complex field. Richardson's search
relates complex exponentials of the same angle only when the forms
agree; with tangent generators in the complex field, relations with
$e^{ix}$ were missed during development. The complex field uses complex
forms throughout, and real views use real forms; an option
(`GR_TOWER_OPT_TRIG_FORM`, after Calcium's option of the same name)
selects the tangent forms for real arguments in the complex field.
:::

## Real-first merging {#real-first}

In a state that has real views, a real algebraic generator that a merge
could only express through nonreal generators is adjoined and moved
before the first of them, provided its minimal polynomial lies below
that position (`GR_TOWER_MERGE_REAL_FIRST`); the steps after it become
dynamic and are proven or refined by Trager's method
([](algorithms/exact/dynamic-evaluation#trager)), and a root of unity
placed in front of a tower by its definition is placed after the tangent
generators already there. Thus the real factors of $x^5 - 3$ are $x^2 +
a(1 \mp \sqrt 5)x/2 + a^2$ with $a = \sqrt[5]{3}$, rather than
expressions in $\zeta_5$.

*Example ($\zeta_5$ and $\sqrt 5$).* Let $C$ be a complex castle field
with a real view $R$, $z = \zeta_5$ computed in $C$ and $s = \sqrt 5$ in
$R$. In $z + s$, the merge finds $\sqrt 5$ only through $\zeta_5$;
instead of the Gauss sum of
[](algorithms/exact/field-towers#gauss-sums), $\sqrt 5$ is moved in
front, and Trager's method refines the modulus of $\zeta_5$:

$$
\sqrt 5 : x^2 - 5, \qquad \zeta_5 : x^2 + \tfrac{1 - \sqrt 5}{2}x + 1,
$$

both proven, so that $z + s$ prints as $\zeta_5 + \sqrt 5$. (The result
depends on the direction of the merge: in a fresh state, $s + z$ gives
the Gauss sum form, with $\sqrt 5$ a linear step after $\zeta_5$.) For
$x^5 - 3$ the real factors are returned in tangent normal form,

$$
x^2 + \frac{(\tau^2 - 3)a}{4}x + a^2, \qquad x^2 + \frac{(7 - \tau^2)a}{4}x + a^2, \qquad a = \sqrt[5]{3}, \; \tau = \tan\tfrac\pi5,
$$

which are the factors above since $(\tau^2 - 3)/4 = (1 - \sqrt 5)/2$;
their coefficients live in a tower that still contains $\zeta_5$, with
$\tan(\pi/5)$ as a dynamic step after it, and involve only the real
generators $a$ and $\tau$.

## Conjugation and parts {#conj}

`gr_conj` fixes real generators, sends $i$ to $-i$ and roots of unity to
their inverses, conjugates $\exp$, $\log$ and root generators through
their definitions (away from the branch cut of $\log$, or with the cut
handled when the argument is a negative real number), and conjugates
any other algebraic generator as an algebraic number (the conjugate, a
root of the same polynomial, is adjoined next to the generator, building
a splitting tower, [](algorithms/exact/castle-roots#conjugates); the
operation is restarted when this changes the tower). Special function
generators are conjugated through their arguments, $F(u) \mapsto
F(\overline u)$ ($W_k(u) \mapsto W_{-k}(\overline u)$) off the branch
cuts, and the values at real points are real generators, which the real
fields accept. `gr_re`, `gr_im`, `gr_arg`, `gr_abs` (as $\sqrt{x\overline
x}$) and `gr_sgn` follow. The trigonometric and hyperbolic functions and
their inverses are expressed through $\exp$, $\log$ and $\sqrt{\cdot}$ in
the complex field (the inverse functions being checked against the
principal branch numerically, with `GR_UNABLE` on the branch cuts), and
`gr_floor`, `gr_ceil`, `gr_trunc`, `gr_nint` and `gr_cmp` are decided
for real elements by combining enclosures with exact zero tests.
Relations found in one tower are shared: when towers are merged,
generators are identified by their definitions (arguments compared
exactly), so that $\exp(ia)$ computed twice is one generator.

## Use in FLINT {#use}

- `gr_ctx_init_ca` with `GR_CA_REAL` and `GR_CA_ALGEBRAIC`, and
  `gr_ctx_init_ca_view`, create the subfields ([](gr_ca#contexts)).
- The generators $\tan(\pi/M)$ have the kind `GR_TOWER_TAN_PI`; merges
  in a state with a real view use `GR_TOWER_MERGE_REAL_FIRST`.
- `src/gr_ca/real.c` and `src/gr_ca/conj.c` implement the real forms and
  conjugation.
