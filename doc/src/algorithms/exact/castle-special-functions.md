---
title: Special functions in castle arithmetic
summary: How castle arithmetic represents values of the gamma, polygamma, zeta, Dirichlet L, Lerch, polylogarithm, error, elliptic integral and Lambert W functions exactly: canonical arguments chosen by functional equations, normal forms at rational arguments with respect to the distribution relations, and the relations between generators found by the zero test.
complexity: The normal forms at rational arguments of level $q$ cost a reduced row echelon form of a relation matrix of size about $q \times q$, cached per thread; the relations found by the zero test cost one such form per level and weight involved, plus a Hermite normal form for the gamma function.
order: 50
---

The values of special functions other than $\exp$ and $\log$ enter a
castle field ([](gr_ca)) as transcendental generators of their own kinds
(see [](gr_ca#special)). This page describes how they are made
canonical and which relations between them castle arithmetic knows.

## Principles {#principles}

Values of special functions (the gamma function, the error function, the
Lambert W function, polygamma functions, polylogarithms, the zeta
function, complete elliptic integrals, and constants such as Euler's
$\gamma$) are transcendental generators of their own kinds. Unlike
$\exp$ and $\log$, most of them are not covered by Schanuel's conjecture:
their algebraic relations come from functional equations (shifts,
reflection, multiplication theorems), special values and transformations.
The principle is the same as for Richardson's algorithm: every relation
used is an identity, applied exactly, and numerical evaluation only
proves nonvanishing and chooses representatives. Generators of these
kinds have the status `GR_TOWER_STATUS_CONJECTURAL`: distinct
generators are assumed algebraically independent unless a relation is
known, so that the zero test is sound (an element which is numerically
indistinguishable from zero without an explaining relation gives
`T_UNKNOWN`) but complete only under the corresponding conjectures
(for instance the Rohrlich–Lang conjecture for the gamma function at
rational arguments).

Castle arithmetic applies the identities before creating a generator, so
that values related by them are expressed through the same generator
(*canonical arguments*). The canonical arguments are intrinsic — they
depend on the value of the argument, not on its representation — so
that different towers create the same generators, which merges then
identify; the zero test then finds the remaining relations between the generators present
([](algorithms/exact/zero-testing#special-relations)).

The relations used are of two kinds. Identities applied *before* a
generator is created (shifts, reflection, multiplication theorems,
special values) express the value through existing generators or
choose its canonical argument. Relations *between* generators already
present, which may come from different towers or computations, are
found by the zero test when a number is not separated from zero, and
eliminate a generator, as for exponentials and logarithms
([](algorithms/exact/zero-testing#elimination)).

The guiding principle of castle arithmetic becomes concrete here: a
relation between special values enters a tower only if it is an instance
of a *known identity* (a functional equation, a special value, a
distribution or period relation), checked exactly. Numerics prove
nonvanishing, choose branches and representatives, and propose candidate
relations, which a family-specific procedure must then derive from the
identities. Each family has a *completeness conjecture* saying that the
known identities generate all algebraic relations: Schanuel's conjecture
for $\exp$, $\log$ and $W$; the Rohrlich–Lang conjecture for $\Gamma$ at
rationals [@Lan1978; @Wal2006]; a conjecture of Chowla–Milnor type for
Hurwitz zeta values [@Mil1983; @GMR2011]; conjectures on polylogarithms
[@Zag2007]; and Grothendieck's period conjecture for elliptic integrals
[@KZ2001]. The status of a generator records which one is relied upon
(`GR_TOWER_STATUS_SCHANUEL`, `GR_TOWER_STATUS_CONJECTURAL`), but never
decides an answer.

Most identities used in practice relate $F(z)$ to $F(z')$ for $z'$ in a
finite orbit of $z$ under a small group. The generator is created only
for a canonical representative, and the requested value is returned as
a cofactor (multiplicative or additive) times the canonical generator.
The canonical representative must be *intrinsic*, a function of the
value rather than of its representation; it is decided from enclosures,
with exact sign and floor tests only where enclosures touch. A first
version of the normalizations split off the rational constant term of
the flat representation of the argument. That term depends on the
generators chosen ($\Gamma(\sqrt 2)$ versus $\Gamma(1 + a)$ with $a =
\sqrt 2 - 1$), so different towers created unrelated generators for the
same value. Calcium simplifies $\Gamma$ at integers and half-integers
and $\operatorname{erf}$ at special points; the rules below subsume that.

## Gamma and polygamma functions {#gamma}

### Canonical arguments {#gamma-canonical}

$\Gamma(z)$, $\psi^{(m)}(z)$: the argument is shifted by an integer
into the strip $0 \le \operatorname{Re}(z) < 1$ (Pochhammer products,
respectively sums of $(z + k)^{-m-1}$) and reflected ($z \to 1 - z$,
with $\pi / \sin(\pi z)$, respectively derivatives of
$\pi \cot(\pi z)$) into $0 < \operatorname{Re}(z) < 1/2$, or
$\operatorname{Re}(z) \in \{0, 1/2\}$ with $\operatorname{Im}(z) > 0$.
Integers and half-integers give factorials and $\sqrt{\pi}$, the
digamma function at rationals is given by Gauss's digamma theorem
(through $\gamma$, logarithms and $\pi$), and $\psi^{(m)}$ at $1$ and
$1/2$ by zeta values. Gauss's digamma theorem [@NIST2012, 5.4.19],

$$
\psi(p/q) = -\gamma - \log(2q) - \tfrac\pi2\cot(\pi p/q) + 2\sum_{k=1}^{\lfloor (q-1)/2 \rfloor} \cos(2\pi kp/q) \log\sin(\pi k/q),
$$

reduces $\psi$ at rationals to $\gamma$ and logarithms of algebraic
numbers, written with the tangent normal form in real fields
([](algorithms/exact/castle-real-fields#tangent-form)).

*Example (canonical arguments of $\Gamma$).* In a fresh context,
$\Gamma(5/2)$ prints as `3*a1/4 {a1 = sqrt(pi)}`, and $\Gamma(2/3) =
2\pi/(\sqrt 3\,\Gamma(1/3))$ by reflection, with the generator
$\Gamma(1/3)$. The canonical point of both $\sqrt 2$ and $1 + \sqrt 2$ is
$\sqrt 2 - 1$, so that

$$
\Gamma(\sqrt 2) = (\sqrt 2 - 1)\Gamma(\sqrt 2 - 1), \qquad \Gamma(1 + \sqrt 2) = (2 - \sqrt 2)\Gamma(\sqrt 2 - 1)
$$

share the generator `t1 = gamma(a1-1)`. For $\sqrt 2 + \frac12 \approx
1.914$ the shift gives $0.914$ and the reflection $(3 - 2\sqrt 2)/2
\approx 0.086$:

$$
\Gamma\bigl(\sqrt 2 + \tfrac12\bigr) = \frac{(\sqrt 2 - \frac12)\pi}{\sin\bigl(\pi(\sqrt 2 - \frac12)\bigr)\Gamma\bigl(\frac{3 - 2\sqrt 2}{2}\bigr)},
$$

with the sine written through a generator $\exp(\pi i(\sqrt 2 -
\frac12))$. Intrinsic canonical points matter when one value reaches the
implementation in two representations. Let $a \approx 1.3247$ be the real
root of $x^3 - x - 1$ and $b \approx 0.3247$ the real root of $x^3 + 3x^2
+ 2x - 1$, created separately (so that $b = a - 1$ is not known). Then
$\Gamma(a) = (a - 1)\Gamma(a - 1)$ with the generator $\Gamma(a - 1)$,
and $\Gamma(1 + b) = b\Gamma(b)$ with the generator $\Gamma(b)$. In the
merged tower of $\Gamma(a)\Gamma(1 + b)$, $b$ maps to $a - 1$, the
argument of the second generator maps to that of the first, and the two
generators are identified: the product is $(a - 1)^2\Gamma(a - 1)^2$, in
a tower of degree 3 over $\mathbb{Q}(\Gamma(a - 1))$.

### Gamma at rational arguments {#gamma-rational}

$\Gamma(p/q)$ at rationals with denominator $q \le 36$: a normal form
with respect to all the relations which follow from the reflection
formula and Gauss's multiplication theorem (the distribution relations,
which by the Rohrlich–Lang conjecture are all the algebraic relations
between these values and $\pi$). A reduced row echelon form of the
relations at level $q$ (in the logarithms of the values, of $\pi$, of
integers and of the numbers $\sin(\pi k/q)$), with the values ordered
by preference (smallest denominator, then smallest numerator), selects
$\varphi(q)/2$ basis values at level $q$, which contain the basis values
at every level dividing $q$; every other value is a monomial in the
basis values times rational powers of $\pi$, of integers and of products
of sines. The reflection relations $x_k + x_{q-k} = \log\pi -
\log\sin(\pi k/q)$ and Gauss's multiplication relations

$$
\Gamma(nz) = (2\pi)^{(1-n)/2} n^{nz - 1/2} \prod_{j=0}^{n-1} \Gamma\Bigl(z + \frac jn\Bigr), \qquad n \mid q,
$$

are linear in the $x_k = \log\Gamma(k/q)$ with right-hand sides in
$\log\pi$, $\log 2$, $\log n$ and $\log\sin(\pi k/q)$; their reduced row
echelon form leaves $\varphi(q)/2$ free columns for $q \ge 3$, in
agreement with Kubert's theorem on the universal ordinary distribution
[@Kub1979] taken modulo the odd relations. Because smaller denominators
are preferred, the basis at level $q$ contains the bases at all divisors
of $q$, so values computed at different levels are expressed in one
basis (level 6 gives $\{\Gamma(1/3)\}$, level 12 gives $\{\Gamma(1/3),
\Gamma(1/4)\}$). Chudnovsky proved the algebraic independence of $\pi$
with $\Gamma(1/3)$ and with $\Gamma(1/4)$ [@Chu1984]. For example
$\Gamma(1/6) = 2^{-1/3} \sqrt{3/\pi}\, \Gamma(1/3)^2$
and $\Gamma(2/3) = 2\pi / (\sqrt{3}\, \Gamma(1/3))$. The products of sines
are written in terms of the root of unity $e^{\pi i/q}$ (square roots of
rationals for $q = 3, 4, 6$), so that values at one level share a
cyclotomic field rather than introducing independent algebraic
generators. Values with larger denominators (or when the normal form
cannot be built) use the reflection formula only and become
generators, possibly dependent on each other: the zero test then
computes the distribution relations between the gamma generators at
rational arguments of each level $N \le 240$ exactly (no numerical
relation search is needed, since they are theorems) and eliminates
the latest generator involved, which becomes algebraic of degree one
over the others, $\pi$, roots of integers (of prime power orders, as
for radicals, or Gauss sums for square roots) and of $\pi$, and
roots of products of sines (written through $e^{\pi i/N}$). The
relation is made an integral combination of the distribution relations
through a Hermite normal form at the level, cached per thread.

*Example (levels 5, 6 and 12).* $\Gamma(1/6)$ prints as
`a3*a2^2*a1*t1^2/(2*pi)` with $t_1 = \Gamma(1/3)$, $a_1 = \sqrt\pi$, $a_2 =
2^{1/3}$, $a_3 = \sqrt 3$, that is, $2^{-1/3}(3/\pi)^{1/2}\Gamma(1/3)^2$;
$\Gamma(5/6) = 2^{4/3}\pi^{3/2}/(\sqrt 3\,\Gamma(1/3)^2)$, and the product
of the two simplifies to $2\pi$. At level 12 the basis is $\{\Gamma(1/3),
\Gamma(1/4)\}$, $\Gamma(3/4) = \pi\sqrt 2/\Gamma(1/4)$, and

$$
\Gamma\bigl(\tfrac{1}{12}\bigr) = \frac{3^{1/8}\alpha\,\Gamma(\frac13)\Gamma(\frac14)}{\sqrt\pi}, \qquad
\alpha = \sqrt{-\zeta_3\zeta_8^3 - \zeta_3\zeta_8 - 2\zeta_8^3 + \zeta_8} = 2\Bigl(\sin\frac\pi3\sin\frac{5\pi}{12}\Bigr)^{1/2},
$$

the radical of a sine product being a square root over
$\mathbb{Q}(\zeta_8)(\zeta_3)$; this is the classical
$2^{-1/4}3^{3/8}(1 + \sqrt 3)^{1/2}\Gamma(\frac13)\Gamma(\frac14)/\sqrt\pi$. At
level 5 the basis has $\varphi(5)/2 = 2$ elements: $\Gamma(1/5)$ and
$\Gamma(2/5)$ are generators, and $\Gamma(3/5) = \pi/(\sin(2\pi/5)
\Gamma(2/5))$ with $1/\sin(2\pi/5)$ written in $\mathbb{Q}(i, \zeta_5)$.

Two lessons shaped the relation round of the zero test. First, a
relation is only useful in a good form: row reduction over $\mathbb{Q}$
expresses the constants in a basis of sines with fractional exponents,
so an integral combination of the original relations is computed by a
Hermite normal form of $[A \mid I]$, which gives the sines integer
exponents (one example went from 78 s to 1 s). Second, eliminations must
avoid algebraic denominators: the latest generator is written as

$$
g = \prod_j \Gamma_j^{a_j} \cdot \pi^{b/M} \cdot \prod_p p^{e_p/M} \cdot (\text{sine product})^{1/d},
$$

with the radicals as root generators (adjoined one prime power at a time
and combined with existing roots of the same radicand), inverse sines as
polynomials in the root of unity ($1/(w - 1) = \frac1r\sum_j jw^j$ for
$w^r = 1$), and $\pi$ and roots of unity adjoined and moved to the front
as needed. When the normal form failed in an early version, creating
$\Gamma(p/q)$ as a generator made the zero test return `T_UNKNOWN` on
true identities, because the generator depended on the basis; the first
fix refused to create it. The right fix was the opposite: the design
tolerates dependent generators, and what was missing was a relation
search for them. Once the zero test computed $\Gamma$ relations, the
fallback generator was harmless. The same pattern (allow, then discover)
is the dynamic evaluation principle applied to every kind of generator.

*Example (a duplication formula at level 38).* The duplication formula
at $z = 1/38$ reads $\Gamma(\frac{1}{38})\Gamma(\frac{10}{19}) =
2^{18/19}\sqrt\pi\,\Gamma(\frac{1}{19})$. Level 38 exceeds the limit 36 of
the normal form, so $\Gamma(1/38)$, $\Gamma(9/19)$ and $\Gamma(1/19)$ are
generators, and $\Gamma(10/19) = \pi/(\sin(9\pi/19)\Gamma(9/19))$ with
$1/\sin(9\pi/19)$ written in $\mathbb{Q}(i, \zeta_{19})$. The zero test
computes the relations among the three generators at level 38 and proves
the identity. With the option `GR_TOWER_OPT_SPECIAL_RELATION_LEVEL_LIMIT`
set to 37 the same comparison is undecided; setting the option to 240 on
the same context makes the same two elements compare equal. The
triplication formula $\Gamma(\frac{1}{40}) = (2\pi)^{-1}3^{-19/40}
\Gamma(\frac{1}{120})\Gamma(\frac{41}{120})\Gamma(\frac{81}{120})$ is
proved the same way at level 120.

### Polygamma and Hurwitz zeta at rational arguments {#hurwitz-rational}

$\psi^{(m)}(p/q)$ ($m \ge 1$) and $\zeta(s, p/q) = (-1)^s \psi^{(s-1)}(p/q) / (s-1)!$
at rationals with denominator $q \le 240$: the additive analogue of
the gamma normal form. The values $x_k = \zeta(s, k/q)$,
$0 < k \le q$ ($x_q = \zeta(s)$), satisfy the reflection relations
$x_k + (-1)^s x_{q-k} = (-1)^{s-1} \pi^s P_{s-1}(\cot(\pi k/q)) / (s-1)!$
(with $P_j(\cot x) = (d/dx)^j \cot x$) and the distribution relations
$\sum_{j<m} x_{b + jq/m} = m^s x_{bm}$ ($m \mid q$), which by a
conjecture of the Chowla–Milnor type are all the linear relations
between these values over the algebraic numbers. A reduced row echelon
form with the same preference order gives a basis of values at each
level (with $\zeta(s)$ for odd $s$, and Catalan's constant in place of
$\psi'(1/4) = \pi^2 + 8G$); every other value is a rational combination
of the basis values plus $\pi^s$ times an element of a cyclotomic field
(the cotangents written through roots of unity). For example
$\psi'(3/4) = \pi^2 - 8G$, $\psi''(1/3) = -4\pi^3/(3\sqrt 3) - 26\zeta(3)$,
$\psi'(1/6) = 5\psi'(1/3) - 4\pi^2/3$. The cotangent polynomials satisfy
$P_0(c) = c$, $P_{j+1} = -(1 + c^2)P_j'$, and the elementary parts
$\sum_k c_kP(\cot\pi k/N)$ are computed modulo $x^m - 1$, reduced modulo
$\Phi_m$ once, and written as monomials in the prime-power roots of unity
(an early version that built each cotangent as a flat element ran out of
memory at level 105). The linear independence of Hurwitz zeta values at
rationals is the subject of the Chowla–Milnor conjecture [@Mil1983;
@GMR2011]; for $\psi$ itself, the transcendence results of Murty and
Saradha [@MS2007] and Okada's theorem on the linear independence of
cotangent and digamma values [@Oka1981] are related. The row echelon forms are
cached per thread by weight and level. Values with larger denominators
are generators (after reflection), and the zero test computes the
relations between the Hurwitz generators (polygamma values at
rationals, $\zeta(s)$ for odd $s$, Catalan's constant) at each weight
and level $N \le 240$ exactly, through the normal form of the cached
level; every relation whose pivot is a
generator makes that generator a rational combination of earlier ones
plus $\pi^s$ times an element of $\mathbb{Q}(\zeta_{\operatorname{lcm}(N,4)})$
(computed modulo the cyclotomic polynomial, then written through the
roots of unity of the tower, with no algebraic denominator). On a
common rational line $a w + b$ with irrational $w$ (as for
$\Gamma$), the zero test uses the shift, reflection and multiplication
relations along the line: they are exact (no $2\pi i$ ambiguity, no
exponential factors), with constants $(a w + c)^{-s}$ and the
reflection terms, whose cotangents are written through
$\exp(\pi i w)$ and roots of unity as
$i (Y + 1)(1 + Y + \cdots + Y^{r-1}) / (Y^r - 1)$ with
$Y = e^{2\pi i (a w + c)}$, so that
$\psi'(\sqrt 2) + \psi'(\sqrt 2 + 1/2) = 4 \psi'(2 \sqrt 2)$ is decided.
The digamma function takes part with $s = 1$ ($-\psi$ in place of
$\zeta(1, z)$, $P_0(x) = x$), its multiplication theorem
$\sum_{j<n} \psi(z + j/n) = n \psi(nz) - n \log n$ contributing the
logarithms of the primes: $\psi(2\sqrt 2) = (\psi(\sqrt 2) +
\psi(\sqrt 2 + 1/2))/2 + \log 2$.

*Example (values at levels 3, 4 and 6).* Gauss's theorem gives
$\psi(1/4) = -\gamma - \pi/2 - 3\log 2$, printed as
`(-pi-6*t2-2*t1)/2 {t1 = euler; t2 = log(2)}`. At level 4, $\zeta(2,
1/4) = \pi^2 + 8G$ and $\psi'(3/4) = \pi^2 - 8G$ with Catalan's constant
$G$ as the generator, and $\zeta(3, 1/4) = \pi^3 + 28\zeta(3)$; at level
6, $\zeta(3, 1/6) = 2\sqrt 3\pi^3 + 91\zeta(3)$. At level 3 and $s = 2$,
$\psi'(1/3)$ itself is the generator:

$$
\zeta\bigl(2, \tfrac23\bigr) = \tfrac43\pi^2 - \psi'\bigl(\tfrac13\bigr), \qquad
\operatorname{Li}_2(e^{2\pi i/3}) = -\frac{\pi^2}{18} + \frac{i\sqrt 3}{27}\bigl(3\psi'\bigl(\tfrac13\bigr) - 2\pi^2\bigr).
$$

### The logarithmic gamma function {#lgamma}

$\log \Gamma(z)$ (`lgamma`, the analytic continuation from the
positive reals with the cut on the negative reals, as
`acb_lgamma`) is $\log(\Gamma(z)) + 2\pi i k$, the integer $k$
determined numerically: the relations of $\log \Gamma$ follow from
those of $\Gamma$ and the structure of the logarithms
($\log\Gamma(\pi + 1) = \log\Gamma(\pi) + \log \pi$).

### Rational lines {#lines}

$\Gamma$ at irrational arguments on a common rational line
$a w + b$ ($a, b$ rational: $\Gamma(\sqrt 2)$, $\Gamma(\sqrt 2 + 1/2)$,
$\Gamma(2 \sqrt 2)$): the zero test finds the line (an integer relation
between the arguments, verified exactly) and the shift, reflection and
multiplication relations along it, and eliminates a generator in terms
of the others, $\pi$, $\exp(\pi i w)$, $\exp(w \log p)$ (from
$n^{n z}$), roots and linear factors $a w + c$. The result is checked
numerically, since the logarithmic relations hold modulo $2 \pi i$.
Since $\Gamma(a w + b)$ is, by these relations, a constant times the
product of the $\Gamma(w + r)$ for $r$ in $\{(b + i)/a \bmod 1\}$
(likewise for the Hurwitz zeta function, with a sum), values each of
whose sets of residues has an element of its own (as two values from
unrelated computations which happen to lie on a line) have no
relation, and the linear algebra is skipped.

## Zeta and L-functions {#zeta}

### The Riemann zeta function {#riemann}

$\zeta(n)$ at integers: Bernoulli numbers (rational multiples of
$\pi^n$ for even $n$); $\zeta(3), \zeta(5), \ldots$ are generators.
At other $s$, the functional equation
$\zeta(s) = \pi^{s - 1/2} \Gamma((1-s)/2) / \Gamma(s/2)\, \zeta(1-s)$
gives the canonical half $\operatorname{Re}(s) > 1/2$ (or
$\operatorname{Re}(s) = 1/2$, $\operatorname{Im}(s) \ge 0$): $\zeta(1/3)$
is expressed through $\zeta(2/3)$.

### Dirichlet L-functions {#dirichlet}

Dirichlet $L$-functions (`gr_dirichlet_l`): $L(s, \chi) =
L(s, \chi^*) \prod_{p \mid q} (1 - \chi^*(p) p^{-s})$ for the primitive
character $\chi^*$ of conductor $f$ inducing $\chi$ mod $q$. At
integers, $L(s, \chi^*) = f^{-s} \sum_{a \le f} \chi^*(a) \zeta(s, a/f)$ in
the Hurwitz normal form ($-(1/f) \sum \chi^*(a) \psi(a/f)$ at $s = 1$,
Bernoulli polynomials at $s \le 0$), so that $L(1, \chi_{-4}) = \pi/4$
and $L(2, \chi_{-4}) = G$. Otherwise the functional equation
$L(s, \chi) = \varepsilon(\chi) (f/\pi)^{1/2 - s}
\frac{\Gamma((1 - s + a)/2)}{\Gamma((s + a)/2)} L(1 - s, \bar\chi)$
($\varepsilon(\chi) = \tau(\chi)/(i^a \sqrt f)$ with the Gauss sum
$\tau(\chi)$ written through roots of unity, $a$ the parity; conductors
$f \le 240$, otherwise `GR_UNABLE` at the other half) gives the
canonical half as for $\zeta$ (at $s = 1/2$, the smaller Conrey label
of $\chi$, $\bar\chi$), where $L(s, \chi^*)$ is a generator of kind
`GR_TOWER_DIRICHLET_L` (printed as `dirichlet_l(chi(f, k), s)`);
its conjugate is $L(\bar s, \bar\chi)$.

### The Hurwitz zeta function {#hurwitz}

$\zeta(s, a)$ for $s$ not an integer $\ge 2$ (those are polygamma
values): $-B_{1-s}(a)/(1-s)$ at integers $s \le 0$ (the continuation
in $a$, also at the nonpositive integers $a$, as
`acb_hurwitz_zeta` at $a = 0$; `GR_DOMAIN` there for the other
$s$, a term $(a + k)^{-s}$ being infinite); at rational $a$
(shifted into $(0, 1]$, denominators $q \le 240$) the character sum
$\zeta(s, p/q) = \frac{q^s}{\varphi(q)} \sum_{\chi \bmod q} \bar\chi(p) L(s, \chi)$,
so that the distribution relations
$\sum_{j<m} \zeta(s, (a + j)/m) = m^s \zeta(s, a)$ and the relations
with $L$-functions and $\zeta(s)$ hold by construction; at irrational
$a$ a generator of kind `GR_TOWER_HURWITZ_ZETA` (two arguments) at
the shift with $0 < \operatorname{Re}(a) \le 1$ (with the terms
$(a + k)^{-s}$, principal powers). The multiplication relations on
rational lines $a w + b$ are not used at non-integer $s$.

### The Lerch transcendent {#lerch}

The Lerch transcendent $\Phi(z, s, a) = \sum_n z^n (n + a)^{-s}$
(`gr_lerch_phi`): $\zeta(s, a)$ at $z = 1$,
$v^{-s} \sum_{k<v} z^k \zeta(s, (a + k)/v)$ at the roots of unity $z$
of order $v \le 240$, $\operatorname{Li}_s(z)/z$ at $a = 1$ with an
integer $s$; and $\operatorname{Li}_s(z) = z \Phi(z, s, 1)$ at the roots
of unity for non-integer $s$.

## Polylogarithms {#polylog}

$\operatorname{Li}_s(z)$: rational functions for $s \le 0$,
$-\log(1 - z)$ for $s = 1$, zeta values at $z = \pm 1$,
$\operatorname{Li}_2(1/2) = \pi^2/12 - \log^2(2)/2$; at a root of unity
$w = e^{2\pi i a/N}$ ($3 \le N \le 240$, recognized numerically and
verified exactly),
$\operatorname{Li}_s(w) = N^{-s} \sum_{k=1}^{N} w^k \zeta(s, k/N)$ in terms
of Hurwitz zeta values (in their normal form), so that
$\operatorname{Li}_2(i) = -\pi^2/48 + iG$ and
$\operatorname{Li}_3(i) = -3\zeta(3)/32 + i\pi^3/32$, and the
distribution relations of polylogarithms at roots of unity follow.
$\operatorname{Li}_2$ at other arguments is reduced under the
anharmonic group $\{z, 1-z, 1/z, 1/(1-z), 1-1/z, z/(z-1)\}$ by the
reflection and inversion formulas (with logarithms and $\pi^2$, and
$-i\pi\log z$ on the cut $z > 1$, as for `acb_polylog`) to a
canonical argument: in $[0, 1/2]$ for real $z$, otherwise the orbit
point in the region $\operatorname{Re}(w) \le 1/2$, $|w| \le 1$,
$|w - 1| \le 1$ (with $\operatorname{Im}(w) > 0$ on its boundary).
The zero test finds the distribution relations
$\operatorname{Li}_2(y^n) = n \sum_{k<n} \operatorname{Li}_2(\zeta_n^k y)$
($|y| \le 1$, $n = 2, 3, 4$) between the generators present (matching
the arguments through the anharmonic orbits, with the special values
at $0, \pm 1, 1/2, 2$), and eliminates the latest generator occurring
linearly in the relation after the eliminated ones are substituted. A
rational value without a generator is adjoined at its canonical point
when its height is smaller than that of the generator to eliminate
(which bounds the process: the descent is well founded, so the rounds
terminate). The relation is built as a flat element, reduced modulo the
eliminated generators, and solved for the latest generator occurring
linearly. The values $\operatorname{Li}_2((3 - \sqrt 5)/2) = \pi^2/15 -
\log^2((\sqrt 5 - 1)/2)$ and $\operatorname{Li}_2((\sqrt 5 - 1)/2) = \pi^2/10
- \log^2((\sqrt 5 - 1)/2)$ follow from the relation with $n = 2$ and $y =
(\sqrt 5 - 1)/2$, which connects three points of one orbit
[@Lew1981; @Zag2007]. Abel's five-term relation
$\operatorname{Li}_2(x) + \operatorname{Li}_2(y) - \operatorname{Li}_2(xy)
- \operatorname{Li}_2(\frac{x(1-y)}{1-xy}) - \operatorname{Li}_2(\frac{y(1-x)}{1-xy})
= \log\frac{1-x}{1-xy} \log\frac{1-y}{1-xy}$ is used likewise for
real $x, y \in (0, 1)$ among the orbit points of the generators'
arguments, when the five values are found in the orbits.
$\operatorname{Li}_s$ for $s \ge 3$ (other than at roots of unity) is
reduced by the inversion formula
$\operatorname{Li}_s(z) + (-1)^s \operatorname{Li}_s(1/z) =
-\frac{(2\pi i)^s}{s!} B_s(\frac12 + \frac{\log(-z)}{2\pi i})$
(with $\log(-z) = \log z + \pi i$ on the cut $z > 1$) to $|z| < 1$,
or $|z| = 1$ with $\operatorname{Im}(z) > 0$; $\operatorname{Li}_3(1/2)
= 7\zeta(3)/8 - \pi^2\log(2)/12 + \log^3(2)/6$. The zero test finds
the distribution relations
$\operatorname{Li}_s(y^n) = n^{s-1} \sum_{k<n} \operatorname{Li}_s(\zeta_n^k y)$
between the generators ($n = 2, 3, 4$), which are linear with no
elementary part at these arguments. Landen's relations for
$\operatorname{Li}_3$ (and the other functional equations of the
higher polylogarithms) are not used.

*Example (a relation between $\operatorname{Li}_2(1/3)$ and
$\operatorname{Li}_2(1/4)$).* The values $\operatorname{Li}_2(1/3)$ and
$\operatorname{Li}_2(1/4)$ are at canonical points, and each becomes a
generator. The zero test proves

$$
\operatorname{Li}_2\bigl(\tfrac14\bigr) = \frac{\pi^2}{6} - 2\operatorname{Li}_2\bigl(\tfrac13\bigr) - \log^2 2 - \log^2\tfrac32
$$

from the relation with $n = 2$ and $y = -1/2$, $\operatorname{Li}_2(\frac14) =
2(\operatorname{Li}_2(\frac12) + \operatorname{Li}_2(-\frac12))$, where
$\operatorname{Li}_2(\frac12) = \frac{\pi^2}{12} - \frac12\log^2 2$ and $-1/2$
lies in the orbit of $1/3$. The generator $\operatorname{Li}_2(1/4)$,
created later, is the one eliminated. The relation with $y = 1/3$ then
proves $\operatorname{Li}_2(\frac13) + \operatorname{Li}_2(-\frac13) =
\frac12\operatorname{Li}_2(\frac19)$. The dilogarithm value
$\operatorname{Li}_2((\sqrt 5 - 1)/2)$ is returned in closed form, as
$(\pi^2 - 10t^2)/10$ with $t = \log((\sqrt 5 + 1)/2)$.

## Error functions {#erf}

$\operatorname{erf}$: odd, so the argument has $\operatorname{Re}(z) > 0$;
on the imaginary axis $\operatorname{erf}(iy) = i \operatorname{erfi}(y)$
with a real generator; $\operatorname{erfc}(z) = 1 - \operatorname{erf}(z)$
and $\operatorname{erfi}(z) = -i \operatorname{erf}(iz)$.

## Complete elliptic integrals {#elliptic}

$K(m) = (\pi/2)\, {}_2F_1(1/2, 1/2; 1; m)$,
$E(m) = (\pi/2)\, {}_2F_1(-1/2, 1/2; 1; m)$ are evaluated through the
hypergeometric normal form
([](algorithms/exact/castle-hypergeometric)): $K(0) = E(0) = \pi/2$, $E(1) = 1$
(Gauss's sum), $K(1/2) = \Gamma(1/4)^2 / (4 \sqrt{\pi})$ and $E(1/2)$
(Gauss's second summation theorem, with contiguity), and the
imaginary-modulus transformation $K(m) = K(m/(m-1))/\sqrt{1-m}$,
$E(m) = \sqrt{1-m}\, E(m/(m-1))$ to a canonical argument with
$|m - 1| \le 1$ ($\operatorname{Im}(m) \ge 0$ on its boundary), which is
Pfaff's transformation. The generators $K(m)$, $E(m)$ are the closed
forms of the basis values of the coset $(1/2, 1/2; 1) + \mathbb{Z}^3$;
the singular values $m = k_r^2$ for $r = 2, 3, 4$ and their complements
(Chowla–Selberg: $K$ through $\Gamma$ at rationals and $\pi$, $E$
through the elliptic alpha function) are applied when such a generator
would be created. The zero test finds Legendre's relation
$E(m) K(1-m) + E(1-m) K(m) - K(m) K(1-m) = \pi/2$ and Landen's
transformation $K(m) = (1 + k_1) K(k_1^2)$,
$E(m) = (1 + s) E(k_1^2) - s K(m)$ ($s = \sqrt{1-m}$,
$k_1 = (1-s)/(1+s)$) between the generators present, also through the
transformation (the partner of a generator at $a$ may sit at $1 - a$,
$(a-1)/a$, $1/(1-a)$ or $1/a$), and eliminates the latest generator
involved (degree one, with square roots as root generators). Landen
pairs $(a, b)$ are detected exactly by the polynomial condition $a^2(b +
1)^2 = 4(2 - a)^2b$ and verified after adjoining $\sqrt{1 - a}$.

*Example (singular values, Landen and Legendre).* By the rules of the
canonical arguments, $K(\frac12) = \Gamma(\frac14)^2/(4\sqrt\pi)$,
$E(\frac12) = (\Gamma(\frac14)^4 + 8\pi^2)/(8\sqrt\pi\,\Gamma(\frac14)^2)$,
$K(-1) = K(\frac12)/\sqrt 2$ by the imaginary-modulus transformation,
and $K(\frac{2 - \sqrt 3}{4}) = 3^{1/4}\Gamma(\frac13)^3/(2^{7/3}\pi)$ (the
singular value for $r = 3$). The values $K(3/4)$ and $K(1/9)$ are
generators. The pair $(a, b) = (3/4, 1/9)$ satisfies the Landen condition
(both sides of $a^2(b + 1)^2 = 4(2 - a)^2b$ equal $25/36$), $\sqrt{1 - a}
= 1/2$ is rational, and the zero test proves $K(3/4) = \frac43 K(1/9)$.
With the four generators $K(\frac13)$, $E(\frac13)$, $K(\frac23)$,
$E(\frac23)$, Legendre's relation $E(\frac13)K(\frac23) + E(\frac23)K(\frac13)
- K(\frac13)K(\frac23) = \pi/2$ is proved [@BB1987].

## Hypergeometric functions and Lambert W {#other}

${}_pF_q(a_1, \ldots, a_p; b_1, \ldots, b_q; z)$ (the generators, of
kind `GR_TOWER_HYPGEOM`, are not regularized; see
[](algorithms/exact/castle-hypergeometric)).

Lambert W: $W_0(0) = 0$, $W_{0,-1}(-1/e) = -1$; the other relations are
found by Richardson's algorithm
([](algorithms/exact/zero-testing#lambertw)).

## Limitations {#limitations}

Relations beyond those above are not used: Landen's relations of the
trilogarithm and the other functional equations of the higher
polylogarithms, polylogarithms at other algebraic arguments, isogenies
of elliptic integrals beyond Landen's transformation, and relations of
the Hurwitz zeta function at non-integer $s$ on rational lines. An
identity which needs one of them is not decided (`T_UNKNOWN`, or no
termination in the unlimited mode). The generators of conjectural status
are independent as far as castle arithmetic knows, so that a decision
`T_FALSE` involving them is proved numerically, never by independence.

## Use in FLINT {#use}

- The special functions of [](gr) (`gr_gamma`, `gr_polygamma`,
  `gr_zeta`, `gr_hurwitz_zeta`, `gr_dirichlet_l`, `gr_lerch_phi`,
  `gr_polylog`, `gr_erf`, `gr_elliptic_k`, `gr_elliptic_e`,
  `gr_lambertw`, ...) of the castle fields implement the canonical
  arguments and normal forms (`src/gr_ca/special.c`,
  `src/gr_ca/functions.c`, `src/gr_ca/dirichlet.c`); the relations
  between generators are found by the zero test.
- The generators are adjoined by `gr_tower_adjoin_special_flat` and
  `gr_tower_adjoin_special_multi_flat`, with the status
  `GR_TOWER_STATUS_CONJECTURAL`.
- The options `GR_TOWER_OPT_GAMMA_LATTICE_LIMIT` and the effort limits
  bound the levels and the sizes of the special values computed exactly
  ([](gr_tower#tuning-options), [](algorithms/exact/effort-limits)).
