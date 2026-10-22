---
title: Simplification and closed forms
summary: The forms in which castle arithmetic keeps numbers (canonical elements of proven towers, structured generators for radicals, roots of unity and logarithms, relations recorded in the towers, real tangent forms), the rewriting rules of the legacy ca design, and the LLL-based methods with exact verification that recognise algebraic numbers, express them in a given field and find closed forms with radicals and trigonometric functions.
complexity: The automatic rules cost a few exact or numerical tests per operation; recognising an algebraic number of degree $d$ costs LLL reductions of dimension $d + 1$ at increasing precision, each followed by an exact verification.
order: 55
---

Simplification has no canonical meaning for exact real and complex
numbers: a "simplest" form is not defined in general, and finding one
would require at least deciding equality. FLINT therefore aims at a
weaker goal: keep representations small and recognise equal numbers as
often as possible, so that zero tests succeed and results are readable.

## Castle arithmetic {#castle}

### Canonical forms {#castle-canonical}

An element of a castle field is a reduced fraction in the flat
representation of its tower, stored in the shortest prefix of the tower
which contains it. The reduction modulo the triangular set is a normal
form ([](algorithms/exact/field-towers#prop-groebner)), so that the
representation is canonical within a tower whose algebraic steps are
proven irreducible and whose transcendental generators are independent:
then equal elements of the tower have equal representations. Elements of
a number field of one proven generator, in particular, are canonical
polynomials in that generator ([](algorithms/exact/castle-representations#dense)),
and rational results always move to $\mathbb{Q}$. Across towers, and in
towers with dynamic steps, the representation is not canonical; it
becomes so after the refinements and eliminations that the zero test
performs, which are recorded in the tower.

### Structured generators {#castle-generators}

Castle fields choose generators so that equal numbers computed in
different ways get the same generators:

- roots of rational numbers are products of rational powers of
  canonical roots of primes and roots of unity of prime power orders
  ($\sqrt 6 = \sqrt 2 \sqrt 3$, $\sqrt{-2} = i \sqrt 2$), and roots of
  monomials in such roots are monomials in roots of higher orders;
- $\exp(r \pi i)$ for rational $r$ is a root of unity, and
  $\exp(\sum_j c_j \log p_j + r\pi i)$ with rational $c_j$, $r$ and
  algebraic $p_j$ is a product of principal powers and a root of unity;
- logarithms of rational and Gaussian rational numbers are sums of
  logarithms of primes and Gaussian primes, plus a multiple of $\pi i$;
- square roots are denested exactly through quadratic steps
  ($\sqrt{5 + 2\sqrt 6} = \sqrt 2 + \sqrt 3$), and roots of perfect
  powers are read off ($\sqrt{(1 - \pi)^2} = \pi - 1$);
- the values of special functions are taken at canonical arguments
  ([](algorithms/exact/castle-special-functions#principles)).

See [](algorithms/exact/field-towers#structured) and
[](algorithms/exact/castle-roots#radicals).

### Relations recorded in the towers {#castle-relations}

When the zero test finds a relation between generators, it eliminates a
generator, which becomes algebraic over the others
([](algorithms/exact/zero-testing#elimination)); the relation is then
part of the tower and simplifies every element built on it. A relation
search also runs when an exponential or a logarithm is adjoined, so
that $\exp(\log x) = x$, $e^x e^{-x} = 1$ and $\log 8 - 3 \log 2 = 0$
are found at once. Unlike the normal form of the legacy zero test, these
eliminations change the representations of results, which become
smaller.

### Real forms {#castle-real}

In the real fields, numbers whose representation involves roots of
unity are rewritten in real terms: real cyclotomic numbers are
polynomials in $\tan(\pi/M)$, and trigonometric functions of real
arguments are rational functions of half-angle tangents
([](algorithms/exact/castle-real-fields#trig)). The printed forms use the
names of the generators, with their definitions
(`GR_TOWER_PRINT_SYMBOLIC`).

## The legacy rules (Calcium) {#legacy}

### Canonical forms where available {#canonical}

Rational numbers are stored as reduced fractions, and an element of a
field $\mathbb{Q}(a)$ with a single algebraic generator as a polynomial
in $a$ of degree less than $\deg a$; both are canonical once the field
is fixed. `ca_set_qqbar` chooses the field: a quadratic number is placed
in $\mathbb{Q}(\sqrt{N})$ with $N$ the discriminant of its minimal
polynomial divided by its small square factors (found by smooth
factorisation), so that all quadratic numbers of one field share it,
and other algebraic numbers generate their own field. After every
operation, `ca_condense_field` moves results to the smallest of these
fields that the representation allows (rational results to
$\mathbb{Q}$).

### Rewriting rules {#rules}

The functions of [](ca) apply rewriting rules when they create new
extension numbers. The most important ones are:

- *Exponentials*: $e^0 = 1$, $e^{\log z} = z$, $e^{(p/q)\log z} =
  z^{p/q}$, $e^{a \log z + s} = z^a e^{s}$ for rational $a$, and $e^{\pi i
  p/q}$ becomes an algebraic root of unity (for moderate $q$).
- *Logarithms*: $\log 1 = 0$, $\log(e^z) = z + 2\pi i k$, $\log z^a = a
  \log z + 2 \pi i k$ and $\log \sqrt{z} = \frac{1}{2} \log z + 2\pi i k$,
  with the integer $k$ determined by a numerical evaluation (which is
  exact here, since $k$ is an integer and the error is small),
  $\log(-x) = \log x + \pi i$ for $x > 0$, and logarithms of roots of
  unity become rational multiples of $\pi i$.
- *Powers and roots*: integer powers are expanded when small, $(z^a)^b =
  z^{ab}$ when the branches allow it, and $\sqrt{x}$ is computed by
  writing $x = A^2 B$ with a heuristic factorisation of the
  representation of $x$ (rational content and polynomial factors in the
  generators, with smooth factorisation of integers), giving $\pm A
  \sqrt{B}$; square roots of algebraic numbers are computed exactly when
  the degree does not grow, and $\sqrt{-x} = i \sqrt{x}$ for $x > 0$.
- *Complex parts*: $|x|$, $\operatorname{sgn} x$, $\operatorname{Re} x$,
  $\operatorname{Im} x$ and $\arg x$ are computed exactly when $x$ is
  algebraic (within the degree limit), and $\arg x$ is a rational multiple
  of $\pi$ when $\operatorname{sgn} x$ is a root of unity.
- *Trigonometric functions*: with the default exponential form, values
  at rational multiples of $\pi$ become elements of cyclotomic fields;
  inverse trigonometric functions of simple algebraic numbers become
  rational multiples of $\pi$.

Each rule either is an identity valid for all arguments or is applied
only after the conditions on branches have been checked exactly.

### Relations between generators {#relations}

When several generators meet in one field, the relations between them
are found and used to reduce elements (see
[](algorithms/exact/field-towers#legacy)). This simplifies, for example,
$e^{x} e^{-x}$ to 1, $\log 8 - 3 \log 2$ to 0 and $\sqrt{2}\sqrt{3}$ to
$\sqrt{6}$ when $\sqrt{6}$ occurs in the same field. For the zero test,
the generators are additionally rewritten in a normal form based on
complex exponentials and logarithms (see
[](algorithms/exact/zero-testing#exact)); this rewriting is not applied
to the results of ordinary operations, because it tends to make them
less readable.

## Recognising algebraic numbers {#lll}

Several operations on algebraic numbers search for an integer relation
$c_0 v_0 + \cdots + c_k v_k = 0$ between numerical values $v_j$ given as
balls, with the LLL algorithm [@LLL1982] applied to the lattice spanned by
the rows of

$$
\begin{pmatrix} I_{k+1} & \tilde v \end{pmatrix},
$$

where the columns $\tilde v$ hold the real parts (and, for nonreal
input, the imaginary parts) scaled by $2^p$ and rounded to integers
(FLINT uses `fmpz_lll`, see [](fmpz_lll)). The first vector of the
reduced basis gives candidate coefficients $c_j$; the remaining 5% of the precision is used
to check that $\sum c_j v_j$ contains zero (`_qqbar_acb_lindep`). If a
relation with coefficients of $b$ bits exists, it is found once $p$ is
somewhat larger than $(k+1) b$; but a relation found by LLL is only a
candidate, which may be spurious, so every use is followed by an exact
verification:

- `qqbar_guess` looks for an annihilating polynomial of a number $z$ with
  $v_j = z^j$, factors it and selects the factor vanishing at $z$; the
  result is not verified (the input is only a ball), so it is a guess.
- `qqbar_express_in_field` looks for $x = f(\alpha)$ with $v_j =
  \alpha^j$ for $j < \deg \alpha$ and $v_{\deg \alpha} = x$, and verifies
  the result exactly with `qqbar_equal_fmpq_poly_val`. If $x \in
  \mathbb{Q}(\alpha)$, a sufficiently high precision always finds it.
- Within arithmetic, the same searches express one operand in the field
  of the other or both in a cyclotomic field, and guess minimal
  polynomials of results, which are then proved with an explicit
  separation bound (see
  [](algorithms/exact/qqbar-arithmetic#structured)).
- In castle arithmetic, `gr_tower_express` looks for a root of a
  polynomial as an element of a tower, over the monomial basis of the
  tower, and verifies it by the complete zero test
  ([](algorithms/exact/field-towers#express)); the relation searches of
  the zero test look for $\mathbb{Q}$-linear relations between
  logarithms ([](algorithms/exact/zero-testing#richardson)).

## Closed forms {#closed-forms}

`qqbar_get_fexpr_formula` converts an algebraic number to an expression
with radicals and trigonometric functions when it can (see
[](qqbar/formulas)). Quadratic numbers are written as $(a + b
\sqrt{c})/q$ with $c$ made small by smooth factorisation
(`qqbar_get_quadratic`). For a number $x$ of degree $d$ whose field may
be contained in a cyclotomic field $\mathbb{Q}(\zeta_n)$ (by the
Kronecker–Weber theorem, these are exactly the numbers with an abelian
Galois group [@Was1997]), candidate $n$ with $d \mid \varphi(n)$ are
tried, $x$ is expressed in $\mathbb{Q}(\zeta_n)$ with LLL and the
expression is verified exactly; then $x$ is written with $e^{2\pi i k/n}$
or with $\cos$ and $\sin$ of rational multiples of $\pi$. Numbers whose
minimal polynomial is $P(t^k)$ are reduced to $x^k$ (deflation), and
nonreal numbers can be split into real and imaginary parts or modulus
and sign. General solutions of cubics and quartics by radicals are not
implemented.

## Use in FLINT {#use}

- `gr_tower_sqrt`, `gr_tower_express` and the root and exponential
  functions of the castle fields implement the structured generators
  and expressions ([](algorithms/exact/castle-roots)).
- `ca_set_qqbar`, `ca_condense_field`, `ca_exp`, `ca_log`, `ca_pow`,
  `ca_sqrt_factor`, `ca_abs`, `ca_sgn`, `ca_arg` and the trigonometric
  functions implement the legacy rules ([](#legacy)).
- `qqbar_guess`, `qqbar_express_in_field` and `_qqbar_acb_lindep` find
  integer relations; `qqbar_get_fexpr_formula` and `qqbar_get_quadratic`
  produce closed forms.
- `fexpr_expanded_normal_form` provides a purely formal simplification of
  symbolic expressions as rational functions of their non-arithmetic
  subexpressions.
