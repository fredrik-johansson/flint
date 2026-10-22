---
title: Dynamic evaluation and irreducibility proofs
summary: How castle arithmetic computes in towers whose moduli are not known to be irreducible, refining a modulus when an inversion exposes a factorization, so that the zero test of algebraic numbers is complete without irreducibility proofs or minimal polynomials; the pretend-field interface of the generic rings; and the certificates that make moduli proven — structural rules (Mordell, Eisenstein, unramified primes), Capelli's theorem and modular proofs at places, Trager norms — and how they survive the reordering of generators.
complexity: An inversion in a dynamic step costs a gcd of polynomials over the field below, recursively; a refinement replaces one modulus by a factor and strictly lowers the degree of the tower. A modular proof costs a few evaluations of the moduli at places over small finite fields; Trager's method and factorization over a tower of degree $D$ factor a norm of degree $D \deg h$ over $\mathbb{Q}$ or $\mathbb{Q}(t_1, \ldots, t_r)$.
order: 35
---

The steps of a tower ([](algorithms/exact/field-towers)) are simple
extensions $F_k = F_{k-1}[x]/(m_k)$ by a root $a_k$ of a monic
polynomial $m_k$. Castle arithmetic does not require $m_k$ to be known
to be irreducible over $F_{k-1}$: a step can be adjoined with any
polynomial that vanishes at the generator (the minimal polynomial over
$\mathbb{Q}$ of an algebraic number, which may factor over the tower; a
binomial $x^n - u$ for an $n$-th root; a Cauchy modulus in a splitting
tower). The quotient ring $F_{k-1}[x]/(m_k)$ is then not necessarily a
field, but it is treated as one, and corrected when it is found not to
be one: this is *dynamic evaluation*, the D5 principle of Della Dora,
Dicrescenzo and Duval [@DDDD1985; @Duv1994]. Steel's approach to
computing in algebraically closed fields [@Ste2002; @Ste2010] is close
in spirit: fields are towers of extensions by polynomials not known to
be irreducible, and factors are discovered when needed. The theory of
triangular sets and regular chains [@Laz1991; @Kal1993; @ALM1999]
provides the general framework, and the complexity of the D5 principle
has been studied in that setting [@DMSX2006].

## Pretend fields {#pretend}

The generic rings of [](gr) express dynamic evaluation through a small
interface: a ring may *pretend* to be a field
(`gr_ctx_is_pretend_field`, `gr_ctx_set_is_pretend_field`). Generic
algorithms then use field algorithms, and when the pretense breaks — an
attempted inverse of a zero divisor — the operation returns `GR_UNABLE`
and records the zero divisor, which the caller retrieves with
`gr_ctx_recover_zero_divisor`. This is implemented for
$\mathbb{Z}/n\mathbb{Z}$ and for the generic quotient rings $R[x]/(m)$
of [](gr_domains#poly-quotient). A pretend field keeps
`gr_ctx_is_field` false, so that algorithms which are only valid over
genuine fields (root finding, factorization) stay disabled.

For example, let $R = \mathbb{Z}/91\mathbb{Z}$ pretend to be a field and
$A = \begin{pmatrix} 14 & 3 \\ 5 & 7 \end{pmatrix}$. Its determinant
$83$ is a unit of $R$, so $A$ is invertible, but the generic inversion
eliminates with the pivot 14, which is not a unit: it returns
`GR_UNABLE`, and the recovered zero divisor gives $7 = \gcd(14, 91)$, a
factor of the modulus. A caller doing dynamic evaluation would continue
in $\mathbb{Z}/7\mathbb{Z}$ or $\mathbb{Z}/13\mathbb{Z}$; the towers
use the same mechanism with $R = F_{k-1}[x]/(m_k)$.

## Dynamic evaluation {#dynamic}

A modulus that is not known to be irreducible makes a step *dynamic*.
Arithmetic in the quotient ring is correct whatever $m_k$ is: the
evaluation map $F_{k-1}[x]/(m_k) \to \mathbb{C}$, $x \mapsto a_k$, is a
ring homomorphism as long as $m_k(a_k) = 0$. Only inversion can fail:
the inverse of $g(a_k)$ is computed by the extended Euclidean algorithm
on $g$ and $m_k$ over $F_{k-1}$, and if $\gcd(g, m_k) = h$ is not
constant, $h$ is a nontrivial factor of $m_k$. In the original
formulation of the D5 principle, the computation then splits into two
branches, one for each factor; here the generator has a definite
complex value, so only one branch is relevant.

::: lemma {#lem-split}
Let $a_k$ be a simple root of $m_k$ and let $m_k = gh$ with $g, h \in
F_{k-1}[x]$ nonconstant. Then exactly one of $g(a_k)$ and $h(a_k)$ is
zero, and replacing $m_k$ by the vanishing factor yields a ring $F_k'$
with a surjection $F_k \to F_k'$ compatible with the evaluation map.
:::

::: proof
$g(a_k) h(a_k) = m_k(a_k) = 0$, and if both vanished, $a_k$ would be a
multiple root of $m_k$. The surjection is the reduction modulo the
factor; it commutes with evaluation because the factor vanishes at
$a_k$.
:::

The vanishing factor is identified numerically: the value of the other
factor at $a_k$ is nonzero, so an enclosure at a sufficient precision
excludes zero for it and not for the right one. The modulus is
*refined* to that factor in place (`gr_tower_refine_step`), and the
inversion is retried.

::: proposition Refinement preserves the elements {#prop-refine}
Let $m'$ be a factor of $m_k$ vanishing at $a_k$. Every polynomial $p$
which represents an element of $F_{k-1}[x]/(m_k)$ represents the same
number in $F_{k-1}[x]/(m')$; reducing it modulo $m'$ does not change
the number.
:::

::: proof
Both quotient rings map to $\mathbb{C}$ by $x \mapsto a_k$, and $p(a_k)$
does not depend on the modulus; $p - (p \bmod m')$ is a multiple of $m'$
and vanishes at $a_k$.
:::

Existing elements therefore stay valid across refinements: the old
modulus is a multiple of the new one, so no element needs to be
updated, and each is reduced modulo the new modulus the next time it is
used. Refinement is transparent, apart from the degree of the tower,
which decreases. A refinement of a lower step $F_j$, $j < k$, changes
the coefficient field of the steps above it, whose moduli may then have
a leading coefficient which reduces to zero modulo the new modulus
(they become non-normalized); the arithmetic is written to tolerate
this. The same holds for the eliminations of the zero test
([](algorithms/exact/zero-testing#elimination)), which replace the
modulus of a generator, or make a transcendental generator algebraic.

::: algorithm Zero test in an algebraic tower (D5) {#alg-d5}
**Input:** an element $x$ of a tower over $\mathbb{Q}$ with algebraic
generators only.

**Output:** `T_TRUE` if $x = 0$, `T_FALSE` otherwise.

1. Reduce $x$ modulo the current moduli. If the representative is 0,
   return `T_TRUE`.
2. If a numerical enclosure of $x$ excludes 0, return `T_FALSE`.
3. If all steps are proven, return `T_FALSE` (the ring is a field).
4. Attempt to invert $x$ in the nested representation. If this
   succeeds, return `T_FALSE` (units evaluate to units).
5. Otherwise refine the step at which a zero divisor was found
   ([](#lem-split)) and go to step 1.
:::

::: theorem The zero test of algebraic numbers is complete {#thm-complete}
[](#alg-d5) terminates, and its answer is correct.
:::

::: proof
Correctness of each return: a zero representative evaluates to zero
([](#prop-refine)); an enclosure excluding zero proves a nonzero value;
if all steps are proven, the quotient ring is a field and the
evaluation map is an injective field homomorphism, so a nonzero
representative has a nonzero value; if $x$ is a unit, then its value
is a unit of $\mathbb{C}$. Termination: each refinement strictly
decreases the degree $\prod_k \deg m_k$ of the tower over $\mathbb{Q}$,
and when no refinement is possible every modulus is irreducible, in
which case the ring is a field and the inversion succeeds for every
nonzero representative.
:::

The proof shows that irreducibility is never required for correctness,
only for speed: with proven steps, the test reduces to its first three
steps. Proofs of irreducibility are therefore attempted
opportunistically ([](#certificates)), and a step that resists proof
costs only the inversions. In practice, numerical enclosures decide $x
\ne 0$ first, and the exact test is needed only for numbers that are
zero or very close to zero. The test does not compute the minimal
polynomial of anything, nor any precision beyond that which separates
the roots of a factor at refinement time; this is what makes castle
arithmetic decide every predicate of algebraic numbers
([](guide/types/castle-arithmetic#semantics)). In Calcium, by contrast,
an algebraic extension number enters as a [](qqbar) number with its
minimal polynomial, and relations between several algebraic numbers
must be found when a field containing them is created
([](algorithms/exact/field-towers#legacy)). With dynamic evaluation, a
tower may contain the stacked radical expressions of $\cos(\pi/257)$ (a
tower of nominal degree $2^{20}$) or the roots of a polynomial over an
existing tower, with relations discovered only where an operation needs
them.

*Example (refining $x^4 + 1$ over $\mathbb{Q}(\sqrt 2)$).* Over $F_1 =
\mathbb{Q}(\sqrt 2)$, adjoin the root $a_2$ of $x^4 + 1$ near $e^{\pi
i/4}$ as a dynamic step, so that the tower has nominal degree 8 over
$\mathbb{Q}$. The zero test of $y = a_2^2 + a_1 a_2 + 1$ stops at the
numerical test, since $y = 2 + 2i$, and the tower is unchanged;
inversion also works in the product of fields, for instance $1/(a_2 -
a_1) = -\tfrac15 (2a_1 + 2a_2 + a_1 a_2^2 + a_2^3)$. For $w = a_2^2 -
a_1 a_2 + 1$ the representative is nonzero and every enclosure contains
zero. The inversion of $w$ is an extended gcd of $x^4 + 1$ and $x^2 -
a_1 x + 1$ over $F_1$, which finds the proper factor $x^2 - a_1 x + 1$;
it vanishes at $a_2$ ($\zeta_8^2 - \sqrt 2 \zeta_8 + 1 = i - (1 + i) + 1
= 0$), the step is refined to it, the tower has degree 4, and $w$ now
reduces to 0. In the implementation, the inversion is preceded by an
attempt to prove the dynamic steps ([](#certificates)), which here
refines the step by Trager's method to the same factor and proves it
irreducible; with the modular and Trager proofs disabled, the inversion
returns `GR_UNABLE` with the zero divisor recorded, and the refinement
is the one described.

Before the inversion, the zero test applies structural rules which look
only at the values of the generators: a dynamic step $x^2 - c$ after
roots of unity which express $\sqrt c$ as a Gauss sum becomes a step of
degree one ($\sqrt 3$ after $\zeta_3$ and $i$), and the rules of
[](#structural) prove steps which modular proofs miss ($i$ over
$\mathbb{Q}(\zeta_{111}, 3^{1/37})$, since every prime $\equiv 1 \bmod
4$ splits $x^2 + 1$).

## Irreducibility certificates {#certificates}

A proof that $m_k$ is irreducible over $F_{k-1}$ is not needed for
correctness, but it lets the zero test skip inversions, lets root
finding and factoring trust the field structure, and makes the
representations canonical
([](algorithms/exact/castle-representations#canonical)): a nonzero
element of a tower whose steps are all proven is invertible without
computing the inverse, and a non-rational element of such a tower is
definitely not rational. The implementation uses the following
certificates, cheapest first.

### Degree and structure {#structural}

Linear moduli are irreducible. Roots of unity of prime power orders of
distinct primes generate linearly disjoint fields (their discriminants
are coprime), so a cyclotomic step $\Phi_{\ell^e}$ over a field
generated by roots of unity of orders prime to $\ell$ and radicals of
integers prime to $\ell$ is irreducible; more generally, $\Phi_{\ell^e}$
is irreducible over any field unramified at $\ell$, because
$\mathbb{Q}(\zeta_{\ell^e})$ is totally ramified at $\ell$ and the
intersection of the two fields is therefore $\mathbb{Q}$
[@Was1997; @Neu1999]. A binomial $x^n - p$ with $p$ prime is
Eisenstein at every prime above $p$ of a field unramified at $p$, hence
irreducible there. For real radicals of positive integers, a theorem of
Mordell [@Mor1953], which generalizes Besicovitch's result for distinct
primes [@Bes1940], gives the degree of $\mathbb{Q}(\sqrt[n_1]{c_1},
\ldots, \sqrt[n_s]{c_s})$: when the radicands are pairwise coprime (or
equal with coprime orders) and each $x^{n_j} - c_j$ is irreducible over
$\mathbb{Q}$, the degree is $\prod_j n_j$.

*Example.* In $\mathbb{Q}(\sqrt 2)(\sqrt[3]{2})$ the radicands are equal
and the orders 2 and 3 coprime, so Mordell's theorem gives the degree 6,
and $x^3 - 2$ is irreducible over $\mathbb{Q}(\sqrt 2)$ without
computation; the Eisenstein rule does not apply, since
$\mathbb{Q}(\sqrt 2)$ is ramified at 2. For $x^3 - 5$ over
$\mathbb{Q}(\sqrt 2)$ both rules apply. The cyclotomic polynomial
$\Phi_9 = x^6 + x^3 + 1$ is irreducible over $\mathbb{Q}(i, \sqrt 2) =
\mathbb{Q}(\zeta_8)$, which is unramified at 3; the rule for roots of
unity of coprime orders does not cover this case, because of $\sqrt 2$,
but it covers $\Phi_9$ over $\mathbb{Q}(i)$.

### Capelli's theorem {#capelli}

A binomial $x^n - a$ is irreducible over a field $K$ if and only if $a
\notin K^p$ for every prime $p \mid n$, and $a \notin -4K^4$ when $4 \mid
n$ [@Lan2002, Theorem VI.9.1]. The conditions are checked at modular
places (below), which can prove that $a$ is not a $p$-th power.

*Example.* Over $K = \mathbb{Q}(\sqrt 2)$ and at $\ell = 7$, $\sqrt 2$
maps to 3 or to 4 in $\mathbb{F}_7$, which gives two places of degree
one. At the second, $1 + \sqrt 2$ maps to 5, which is not a square in
$\mathbb{F}_7$; hence $1 + \sqrt 2$ is not a square in $K$, and $x^2 -
(1 + \sqrt 2)$ is irreducible over $K$ (at the first place the image 4
is a square, and that place proves nothing). Accordingly, the step
$\sqrt{1 + \sqrt 2}$ is proven when adjoined, while $\sqrt{3 + 2\sqrt 2}$
stays dynamic, as it must since $3 + 2\sqrt 2 = (1 + \sqrt 2)^2$; the
zero test of $\sqrt{3 + 2\sqrt 2} - 1 - \sqrt 2$ refines it to the
linear modulus $x - 1 - \sqrt 2$. The step $\sqrt[3]{2}$ over $K$ is
also proven when adjoined (2 is not a cube in $\mathbb{F}_7$, whose
nonzero cubes are $\pm 1$). The condition for $4 \mid n$ matters for
$x^4 + 4 = (x^2 - 2x + 2)(x^2 + 2x + 2)$ over $\mathbb{Q}$, with $-4 =
-4 \cdot 1^4$: the step $\sqrt[4]{-4}$ stays dynamic, and the zero test
of $a^2 - 2a + 2$ for $a = \sqrt[4]{-4} = 1 + i$ refines it to $x^2 - 2x
+ 2$. For $x^4 + 1$ over $\mathbb{Q}$ the two conditions can be
certified at different places ($-1$ is not a square in $\mathbb{F}_3$,
and $-(-1)/4 = 4$ is not a fourth power in $\mathbb{F}_5$), and the step
$\sqrt[4]{-1}$ is proven when adjoined, although $x^4 + 1$ is reducible
modulo every prime. Over $\mathbb{Q}(\sqrt 2)$, on the other hand, $-1 =
-4(1/\sqrt 2)^4$, and $x^4 + 1 = (x^2 - \sqrt 2 x + 1)(x^2 + \sqrt 2 x +
1)$.

### Modular places {#places}

A chain of simple roots $\alpha_1, \ldots, \alpha_{k-1}$ of the moduli
in a finite field $\mathbb{F}_{\ell^m}$, with $\alpha_j$ a root of the
reduction of $m_j$ at $(\alpha_1, \ldots, \alpha_{j-1})$, defines a
*place* of the tower: a homomorphism from the ring generated by the
generators (with the leading coefficients and denominators inverted) to
$\mathbb{F}_{\ell^m}$.

::: proposition Proofs at a place {#prop-place}
Let $K = F_{k-1}$ be the field defined by a tower over $\mathbb{Q}$
whose steps are all proven, let $\ell$ be a prime not dividing the
leading coefficients and denominators of the moduli, and let
$\alpha_1, \ldots, \alpha_{k-1} \in \mathbb{F}_{\ell^m}$ be such that
$\alpha_j$ is a simple root of the reduction of $m_j$ at $(\alpha_1,
\ldots, \alpha_{j-1})$. Let $f \in K[x]$ be monic with $\ell$-integral
coefficients. (a) If the reduction $\bar f \in \mathbb{F}_{\ell^m}[x]$ is
irreducible of degree $\deg f$, then $f$ is irreducible over $K$.
(b) If the image of $a \in K$ is not a $p$-th power in
$\mathbb{F}_{\ell^m}$, then $a$ is not a $p$-th power in $K$.
:::

::: proof
Let $A$ be the localization of $\mathbb{Z}_{(\ell)}[x_1, \ldots,
x_{k-1}]/(m_1, \ldots, m_{k-1})$ at the kernel $\mathfrak{p}$ of the
map to $\mathbb{F}_{\ell^m}$ sending $x_j \mapsto \alpha_j$. The
Jacobian of the triangular system at the point is $\prod_j \partial
m_j/\partial x_j(\alpha) \ne 0$ because the roots are simple, so $A$ is
a regular local ring of dimension one, that is, a discrete valuation
ring, with fraction field $K$ (here we use that $K$ is a field: the
steps are proven). A discrete valuation ring is integrally closed, so by
Gauss's lemma a factorization $f = gh$ into monic factors over $K$ has
$g, h \in A[x]$, and it reduces to a factorization of $\bar f$ over the
residue field $\mathbb{F}_{\ell^m}$ with the same degrees. The statement
about $p$-th powers follows in the same way from $x^p - a$.
:::

The hypothesis that the steps below are proven is essential: if
$F_{k-1}$ is only a product of fields, the place sees one factor and
says nothing about the others. Its omission produced a wrong answer
during development ([](algorithms/exact/castle-design#findings)): two
cube roots generating the same field over $\mathbb{Q}(\zeta_3)$ made a
step reducible (dynamic), a square root above it was proved irreducible
at a place that mapped the dynamic step to one of its factors, and kept
that status after the dynamic step was refined to the other factor,
over which the radicand became a square. Modular proofs and Capelli
tests at places therefore require the steps below to be proven.

With transcendental generators, a place also specializes each $t_j$ to
a random residue in $\mathbb{F}_\ell$; the local ring at a regular point
is again integrally closed. A modulus then has a root at the place only
for some of the values (a radical $x^n - a(t)$ for about one value of
$t$ in $n$), and when a modulus has no simple root, the values of the
generators which first occur from the earliest step sharing a generator
with it are drawn again, up to 48 times per place. Eisenstein's
criterion at a transcendental generator, which is valid in every factor
of a product of fields, is used as an exception to the proven-below
rule. The primes are drawn from arithmetic progressions $1 \bmod L$
(`n_randprime_mod`), with $L$ chosen from the orders of the roots of
unity present, so that the cyclotomic moduli split and places of degree
one exist; the residue degrees $m$ worth trying are pruned by the
factorization pattern of the first modulus. A proof
(`gr_tower_prove_step_modular`) tries a number of places proportional
to `GR_TOWER_OPT_MODULAR_TRIES`.

The method fails for moduli that are reducible modulo every prime
(such as $x^4 + 1$ over $\mathbb{Q}$) and for chains of roots that
exist only at primes of low density; Trager's method is then used.
Reducible moduli are never proven.

*Example (a modular place).* Let $f = x^3 - \sqrt 2 x - 1$ be the
modulus of a dynamic step over $K = \mathbb{Q}(\sqrt 2)$. At $\ell = 7$
and $\sqrt 2 \mapsto 3$, $\bar f = x^3 - 3x - 1$ has no root in
$\mathbb{F}_7$, so it is irreducible, and by [](#prop-place) so is $f$;
at $\sqrt 2 \mapsto 4$, $\bar f = x^3 - 4x - 1$ has the root 3, and that
place proves nothing. `gr_tower_prove_step_modular` proves the step
within the default number of places. For $x^4 + 1$ as the modulus of a
dynamic step over $\mathbb{Q}$, the same function fails at every place,
and `gr_tower_prove_step_trager`, which for a first step factors the
modulus over $\mathbb{Q}$, proves it irreducible.

The same places prove that a polynomial has no root in a tower
(`gr_tower_poly_no_roots_modular`): a root in the tower would map to a
root of the reduction. This lets merges and root finding skip lattice
searches that would fail ([](algorithms/exact/castle-roots#strategies)).

### Trager norms {#trager}

For a tower over $\mathbb{Q}$ with proven steps below $k$, the *norm*

$$
N(h) = \operatorname{Res}_{x_1}\bigl(m_1, \operatorname{Res}_{x_2}(m_2,
\cdots \operatorname{Res}_{x_{k-1}}(m_{k-1}, h) \cdots)\bigr)
$$

of $h(x) = m_k(x - \theta)$, with $\theta = \sum_j s_j a_j$ for small
integers $s_j$, is a polynomial in $\mathbb{Q}[x]$ of degree $\deg(m_k)
[F_{k-1} : \mathbb{Q}]$. If it is squarefree, the irreducible factors of
$m_k$ over $F_{k-1}$ correspond to the irreducible factors of $N(h)$ over
$\mathbb{Q}$ [@Tra1976]; in particular $m_k$ is irreducible if $N(h)$ is,
and otherwise a gcd over $F_{k-1}$ yields a proper factor, with which the
step is refined (`gr_tower_prove_step_trager`). This is exact and
complete, but its cost grows quickly with $\deg(m_k) [F_{k-1} :
\mathbb{Q}]$, so it is used within a degree limit
(`GR_TOWER_OPT_TRAGER_DEGREE_LIMIT`) after the modular places have
failed. Over $\mathbb{Q}(t_1, \ldots, t_r)$ the norm is factored as a
multivariate polynomial; this is too expensive to run automatically
inside the zero test (an attempt made the trigonometric test suite run
for minutes instead of a third of a second), so the automatic fallback
is restricted to towers without transcendental generators. In the
unlimited mode these limits grow ([](algorithms/exact/effort-limits)),
which makes the proofs, and the predicates which depend on them,
decision procedures.

*Example (Trager norms over $\mathbb{Q}(\sqrt 2)$).* For $m = x^2 -
\sqrt 2$, $N(m) = \operatorname{Res}_y(y^2 - 2, x^2 - y) = x^4 - 2$ is
squarefree and irreducible over $\mathbb{Q}$, so $m$ is irreducible over
$\mathbb{Q}(\sqrt 2)$. For $m = x^4 + 1$, the shift $\theta = 0$ gives
$N(m) = (x^4 + 1)^2$, which is not squarefree. The next shift tried is
$\theta = -\sqrt 2$, with $h = (x + \sqrt 2)^4 + 1$ and

$$
N(h) = x^8 - 8x^6 + 26x^4 - 8x^2 + 25 = (x^4 + 1)(x^4 - 8x^2 + 25).
$$

The gcds of $h$ with the two factors over $\mathbb{Q}(\sqrt 2)$ are
$x^2 + \sqrt 2 x + 1$ and $x^2 + 3\sqrt 2 x + 5$, which shifted back
($x \mapsto x - \sqrt 2$) give the factors $x^2 - \sqrt 2 x + 1$ and
$x^2 + \sqrt 2 x + 1$ of $m$. Applied to the dynamic step of the example
of [](#dynamic), `gr_tower_prove_step_trager` refines the step to $x^2 -
\sqrt 2 x + 1$, the factor that vanishes at $e^{\pi i/4}$, and marks it
proven.

The same method factors polynomials over towers
([](algorithms/exact/castle-roots#factorization)): the norm is factored
by `fmpz_poly_factor` over $\mathbb{Q}$, or over $\mathbb{Q}(t_1, \ldots,
t_r)$ with its denominators cleared as a polynomial in $t_1, \ldots,
t_r, x$ by `fmpz_mpoly_factor`; the factorization requires the steps of
the field to be proven, which is attempted first.

### Proofs and reordering {#reordering}

A proof of irreducibility is relative to the field below the step. When
generators are reordered (a merge moves a generator in front of others,
[](algorithms/exact/field-towers#merging)), the field below changes,
and the proof must be dropped unless a structural rule covers the new
prefix. A generator whose proof is dropped is, when the prefix is small,
re-expressed in the new prefix by an integer relation search: after
$\zeta_8$ moves in front of $\sqrt 2$, the step of $\sqrt 2$ becomes the
linear modulus $x - (\zeta_8 - \zeta_8^3)$. Keeping the proof was a
source of a wrong answer during development: the identity $\sqrt 2 (1 +
i)/2 = e^{\pi i/4}$ was refuted because $x^2 - 2$ had been proved
irreducible over $\mathbb{Q}(i)$, $\zeta_8$ was later moved in front of
$\sqrt 2$, the status stayed proven, but $x^2 - 2$ splits over
$\mathbb{Q}(\zeta_8)$, and the shortcut "all steps proven, so the
representative is canonical" declared a true identity false.

*Example (reordering for $\sqrt 2 (1 + i)/2 = e^{\pi i/4}$).* In a fresh
context, $\sqrt 2 (1 + i)/2$ lives in the tower $\mathbb{Q}(i)(\sqrt
2)$, both steps proven, and $e^{\pi i/4}$ in the tower
$\mathbb{Q}(\zeta_8)$. Merging the two towers represents the roots of
unity of orders 4 and 8 through one generator $\zeta_8$, which is moved
to the front; $i$ becomes the linear step $x - \zeta_8^2$, and the proof
of $x^2 - 2$ over $\mathbb{Q}(i)$ is dropped, since it says nothing over
$\mathbb{Q}(\zeta_8)$. The integer relation search over the prefix
$\mathbb{Q}(\zeta_8)$ then finds $\sqrt 2 = \zeta_8 - \zeta_8^3$, so
that the merged tower is

$$
\zeta_8 : x^4 + 1, \qquad i : x - \zeta_8^2, \qquad \sqrt 2 : x -
(\zeta_8 - \zeta_8^3),
$$

of degree 4, and $\mathbb{Q}(\zeta_8)$ after the elimination of the
linear steps. The two numbers compare equal, and $\sqrt 2 + e^{\pi i/4}$
prints as `-a2^3+2*a2 {a2 = exp(2*pi*i/8)}`.

## Use in FLINT {#use}

- `gr_ctx_is_pretend_field`, `gr_ctx_set_is_pretend_field` and
  `gr_ctx_recover_zero_divisor` are the pretend-field interface of
  [](gr).
- `gr_tower_refine`, `gr_tower_refine_step` and the complete operations
  `gr_tower_is_zero`, `gr_tower_inv`, `gr_tower_div` implement dynamic
  evaluation; the castle fields use them through the flat zero test
  (`gr_tower_flat_num_is_zero`).
- `gr_tower_prove_modular`, `gr_tower_prove_step_modular`,
  `gr_tower_binomial_irreducible_modular` and
  `gr_tower_poly_no_roots_modular` are the proofs at places;
  `gr_tower_prove_step_trager`, `gr_tower_poly_factor` and
  `gr_tower_poly_roots` use Trager's method.
- The options `GR_TOWER_OPT_MODULAR_TRIES`,
  `GR_TOWER_OPT_NO_ROOTS_TRIES`, `GR_TOWER_OPT_TRAGER_DEGREE_LIMIT` and
  `GR_TOWER_OPT_FACTOR_DEGREE_LIMIT` set the effort
  ([](gr_tower#tuning-options)).
