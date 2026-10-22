---
title: Zero testing
summary: Why deciding x = 0 is hard for exact real and complex numbers and what is known about decidability; how castle arithmetic proves x ≠ 0 numerically and x = 0 exactly, by dynamic evaluation for algebraic numbers and by Richardson's algorithm (relation searches between the generators, proved and eliminated) for exponentials, logarithms and trigonometric functions; and the zero test of the legacy Calcium design.
complexity: Proving $x \ne 0$ costs evaluations at doubling precisions up to about $-\log_2 |x|$ bits; proving $x = 0$ costs one evaluation at the precision limit of the round plus the exact tests, of which a relation search costs LLL reductions of dimension at most the number of generators involved plus one, and dynamic evaluation the inversions and refinements of the tower.
order: 40
---

The zero test is the core of exact real and complex arithmetic. Every
decision of castle arithmetic ([](gr_ca)), from comparisons and
divisions to degrees of polynomials and ranks of matrices, reduces to
deciding whether some number $x$ equals zero. This page describes the
zero test of castle arithmetic ([](#castle)) and, at the end, that of
the legacy design of [](ca) ([](#legacy)).

## Decidability {#decidability}

For algebraic numbers the problem is decidable: a number given by an
expression with rational numbers, arithmetic operations and roots can
be evaluated to its minimal polynomial and an isolating ball (see
[](algorithms/exact/qqbar-arithmetic)), which is zero if and only if
the minimal polynomial is $t$, or decided by dynamic evaluation in a
tower without computing any minimal polynomial
([](algorithms/exact/dynamic-evaluation#thm-complete)).

With transcendental functions the situation is different. Richardson
[@Ric1968] showed that it is undecidable whether an expression built
from a variable $x$, rational numbers, $\pi$, $\log 2$, $+$, $\times$,
$\sin$ and $|\cdot|$ vanishes identically on $\mathbb{R}$. For
*constants*, the *constant problem* is open: it is not known whether
equality of numbers built from rational numbers with arithmetic, $\exp$
and $\log$ is decidable. Richardson [@Ric1992; @RF1994; @Ric1997;
@Ric2007] gave algorithms for such exp-log constants that are correct
whenever they terminate and that terminate if Schanuel's conjecture
holds; see also van der Hoeven [@vdH2006]. The structure theorems of
Risch [@Ris1979] and of Rothstein and Caviness [@RC1979] are the
foundation of this approach. Even for specific numbers such as $e +
\pi$, it is not known whether they are rational.

::: conjecture Schanuel {#conj-schanuel}
If $z_1, \ldots, z_n \in \mathbb{C}$ are linearly independent over
$\mathbb{Q}$, then $\operatorname{trdeg}_{\mathbb{Q}} \mathbb{Q}(z_1,
\ldots, z_n, e^{z_1}, \ldots, e^{z_n}) \ge n$.
:::

Lang [@Lan1966] states the conjecture and derives classical
consequences; Ax proved its analogue for power series [@Ax1971]. The way
it is used in exact computation goes back to Richardson. For a tower
whose transcendental generators are exponentials $e^{u_j}$ and
logarithms $\log v_j$, define the *log side* as the multiset of numbers
$u_j$ and $\log v_j$, together with $2\pi i$ when $\pi$ and $i$ are both
present in the tower.

::: proposition {#prop-schanuel-towers}
Assume Schanuel's conjecture. Let $T$ be a tower over $\mathbb{Q}$ whose
transcendental generators are $e^{u_1}, \ldots, e^{u_a}, \log v_1,
\ldots, \log v_b$, all other generators being algebraic. If the numbers
$u_1, \ldots, u_a, \log v_1, \ldots, \log v_b$ are linearly independent
over $\mathbb{Q}$, then the transcendental generators are algebraically
independent over $\overline{\mathbb{Q}}$.
:::

::: proof
Put $z = (u_1, \ldots, u_a, \log v_1, \ldots, \log v_b)$. The field
$\mathbb{Q}(z, e^z)$ is contained in the field generated over
$\mathbb{Q}$ by the generators of $T$, since each $u_j$ and $v_j$ lies
in the field of the earlier generators and $e^{\log v_j} = v_j$. That
field is algebraic over $\mathbb{Q}(e^{u_1}, \ldots, e^{u_a}, \log v_1,
\ldots, \log v_b)$. By the conjecture its transcendence degree is at
least $a + b$, so the $a + b$ transcendental generators are
algebraically independent.
:::

Hence, under Schanuel's conjecture, an element of a tower that reduces
to a nonzero rational function of the transcendental generators (the
outcome we call *field-nonzero*) is a nonzero number *unless* the log
side satisfies a nontrivial $\mathbb{Q}$-linear relation. The zero test
looks for such a relation. When one exists (and the element is zero,
there must be one, under the conjecture), it is eventually found by
integer relation search at sufficient precision and verified exactly; it
is then used to eliminate a generator. When none exists, the element is
nonzero, and numerical evaluation at increasing precision eventually
proves it. The conjecture thus guarantees that the alternation between
the two searches terminates; it is never used to conclude.

*Example (an exponential becomes algebraic).* Let $u = \exp(\tfrac12
\log 2 + 1)$ and $v = \sqrt 2 e$. In a fresh context $u$ is a generator
$t_2 = \exp((t_1 + 2)/2)$ over $t_1 = \log 2$, both of status Schanuel.
When the towers are merged to form $u - v$, the log side $\{\log 2,
\tfrac12\log 2 + 1, 1\}$ (without $2\pi i$: neither $\pi$ nor $i$
occurs) satisfies

$$
2\bigl(\tfrac12 \log 2 + 1\bigr) - \log 2 - 2 \cdot 1 = 0, \qquad \text{that is,} \qquad t_2^2 = 2e^2,
$$

and $e$ enters the tower not as a transcendental generator but as an
algebraic one, $t_3 = e$, with the dynamic modulus $x^2 - t_2^2/2$ over
$\mathbb{Q}(t_1, t_2)(\sqrt 2)$. The difference $t_2 - \sqrt 2 t_3$ is a
zero divisor of the tower, since $(t_2 - \sqrt 2 t_3)(t_2 + \sqrt 2 t_3) =
t_2^2 - 2t_3^2 \equiv 0$; the zero test refines the modulus to the factor
$x - \tfrac{\sqrt 2}{2} t_2$ that vanishes at $e$
([](algorithms/exact/dynamic-evaluation#lem-split)), and returns
`T_TRUE`.

Some independence results are theorems and are used as such: Lindemann's
theorem that $\pi$ is transcendental, the Hermite–Lindemann theorem that
$e^\alpha$ is transcendental for algebraic $\alpha \ne 0$ (and $\log
\alpha$ for algebraic $\alpha \ne 0, 1$), and the Lindemann–Weierstrass
and Gelfond–Schneider theorems [@Bak1975]. In particular:

::: lemma Hermite–Lindemann rule {#lem-hermite-lindemann}
Let $K$ be a number field given by proven algebraic steps and let $t$ be
one of $\pi$, $e^u$ or $\log u$ with $u \in K$ (and $u \ne 0$,
respectively $u \ne 0, 1$). Then every nonzero polynomial $P \in K[X]$
satisfies $P(t) \ne 0$.
:::

This rule decides $e^{10^{-10000}} - 1 \ne 0$ at once, where numerical
evaluation would need about 33000 bits. Values of $\Gamma$ at rationals
and other special values are treated by the conjectural relations of
[](algorithms/exact/castle-special-functions); known transcendence
results such as the algebraic independence of $\pi$, $e^\pi$ and
$\Gamma(1/4)$ [@Nes1996] or of $\pi$ and $\Gamma(1/3)$ [@Chu1984] are
consistent with, but not needed by, the implementation.

For the values of other special functions, the known relations come
from functional equations, special values and transformations, and
conjectures such as that of Rohrlich and Lang for the gamma function at
rational arguments assert that these are all
([](algorithms/exact/castle-special-functions#principles)). Castle
arithmetic implements Richardson's approach throughout: *every answer is
proved; conjectures affect only completeness*. The implementation never
returns `T_TRUE` or `T_FALSE` on the strength of a conjecture or a
heuristic; when a conjecture is false (or a resource limit is reached),
the outcome is `T_UNKNOWN` (`GR_UNABLE` for operations), never a wrong
answer.

## Proving inequality {#numerical}

If $x \ne 0$, then evaluating $x$ in ball arithmetic at a sufficiently
high precision gives a ball that excludes zero: the radius of the
computed ball tends to zero as the precision grows, for the operations
that occur, as long as no function is evaluated on a branch cut or at a
singularity. A ball that excludes zero is a proof. The enclosures of the
generators of a tower are computed recursively from their definitions,
refined by interval Newton iteration for algebraic generators, and
cached ([](gr_tower#numerical-evaluation)); an element is evaluated from
them (`gr_tower_flat_get_acb`).

The required precision is about $-\log_2 |x|$ plus the precision lost in
the evaluation, so numbers that are extremely close to zero without
being zero are expensive. Ramanujan's constant $e^{\pi \sqrt{163}} -
640320^3 - 744 \approx -7.5 \cdot 10^{-13}$ needs about 100 bits and is
decided numerically before any exact stage is entered. Numerical
evaluation alone can never prove $x = 0$: a ball containing zero is
compatible with both answers.

## The zero test of castle arithmetic {#castle}

An element of a castle field is a fraction $p/q$ in the flat
representation of a tower ([](algorithms/exact/field-towers#representations)),
and $x = 0$ if and only if the reduced numerator $p$ vanishes at the
generators. [](#alg-castle-zero) summarizes the zero test as it is run
on a flat element of a tower; the castle fields add a few shortcuts
before it (rational elements and dense elements of proven number fields
are canonical; comparisons of enclosures between elements of different
towers avoid merging towers merely to learn that two numbers differ).
The sections below explain the individual stages.

::: algorithm Zero test of a flat element (outline) {#alg-castle-zero}
**Input:** an element $x$ of a tower $T$.

**Output:** `T_TRUE`, `T_FALSE` or `T_UNKNOWN`.

1. Reduce $x$ modulo the triangular set; if the numerator is 0, return
   `T_TRUE`.
2. If an enclosure of $x$ at low precision excludes 0, return `T_FALSE`.
3. Restrict attention to the generators $x$ involves, closed under the
   dependencies of their definitions and moduli ([](#involved)).
4. Try cheap irreducibility certificates for the involved dynamic steps
   ([](algorithms/exact/dynamic-evaluation#certificates)).
5. If $x$ is algebraic, the tower is large and $x$ splits into parts with
   disjoint generators, return the split zero test ([](#split)).
6. Run the field zero test in the nested representation, refining by
   dynamic evaluation ([](algorithms/exact/dynamic-evaluation#alg-d5)),
   with the outcome `T_TRUE`, `T_FALSE` or *field-nonzero*.
7. If field-nonzero: if the involved transcendental generators are proven
   independent, or [](#lem-hermite-lindemann) applies, return `T_FALSE`.
   Otherwise run the relation rounds: the Hurwitz and polygamma
   relations, gamma at rationals, gamma and polygamma on lines, the
   elliptic and dilogarithm relations
   ([](algorithms/exact/castle-special-functions)), then the logarithmic
   (and, in real towers, angular) integer relation search ([](#richardson),
   [](#angles)). If a verified relation was found, eliminate the top
   generator in place ([](#elimination)) and go to step 1.
8. Evaluate $x$ at increasing precision up to
   `GR_TOWER_OPT_NUMERIC_PREC_LIMIT` ([](#precision)); return `T_FALSE` if
   an enclosure excludes 0, else `T_UNKNOWN`.
:::

::: theorem Soundness {#thm-castle-zero}
Every `T_TRUE` or `T_FALSE` answer of [](#alg-castle-zero) is correct.
:::

::: proof
`T_TRUE` is returned only when the reduced representative is zero in a
ring all of whose relations are true relations among the values: the
moduli vanish at the generator values, and refinements and eliminations
preserve this, as each new modulus is verified to vanish at the
generator — refinements by
[](algorithms/exact/dynamic-evaluation#lem-split), eliminations by an
exact verification of the relation ([](#richardson)). `T_FALSE` is
returned only when an enclosure excludes zero, when the ring is proven
to be a field in which the evaluation map is injective on the involved
generators, when $x$ is a unit
([](algorithms/exact/dynamic-evaluation#thm-complete)), or by the
theorems of transcendence cited above.
:::

::: theorem Conditional completeness {#thm-castle-complete}
Let $T$ be a tower over $\mathbb{Q}$ whose transcendental generators are
exponentials and logarithms. Assume Schanuel's conjecture, and that no
precision or size limit is reached. Then [](#alg-castle-zero) terminates
with `T_TRUE` or `T_FALSE` for every element of $T$.
:::

::: proof
(Sketch.) Each elimination reduces the number of transcendental
generators by one (the insertion of a primitive exponential adds one
generator and turns two into algebraic ones), so there are at most $r$
eliminations. Between eliminations, the algebraic part terminates by
[](algorithms/exact/dynamic-evaluation#thm-complete). If $x$ is
field-nonzero and $x = 0$, the transcendental generators involved are
algebraically dependent over $\overline{\mathbb{Q}}$, so by
[](#prop-schanuel-towers) the log side satisfies a nontrivial
$\mathbb{Q}$-linear relation, which an integer relation algorithm finds
once the precision is large compared with the size of the relation, and
which then verifies. If $x \ne 0$, numerical evaluation at sufficient
precision excludes zero.
:::

With special function generators the corresponding statements hold for
the relations that each round knows (the distribution relations of
$\Gamma$, the Hurwitz zeta values, and so on), and completeness rests on
the conjectures that these generate all relations. In the unlimited mode
([](algorithms/exact/effort-limits)) the precision and size limits grow,
so that the hypothesis "no limit is reached" holds for every zero test
that terminates at all.

*Example (a run of the zero test).* Build a fixed tower with the
generators $a = \sqrt 2$, $t_1 = e$, $t_2 = e^{1+a}$, $t_3 = e^{1-a}$, in
this order, and test $x = e^{1+\sqrt 2}e^{1-\sqrt 2} - e^2$, that is, $x =
t_2t_3 - t_1^2$, with `gr_tower_flat_num_is_zero`.

- Step 1: the reduced numerator $t_3t_2 - t_1^2$ is not 0.
- Step 2: the enclosure of $x$ at 256 bits has radius about $2 \cdot
  10^{-76}$ and contains 0.
- Step 3: all four generators are involved ($a$ through the definitions
  of $t_2$ and $t_3$).
- Steps 4 and 6: the step of $a$ is proven; in the field $\mathbb{Q}(t_1,
  t_2, t_3)(\sqrt 2)$, $x$ is a nonzero element, and the outcome is
  field-nonzero.
- Step 7: the exponentials are not proven independent. The log side
  consists of the arguments $1$, $1 + \sqrt 2$, $1 - \sqrt 2$, and LLL at
  64 bits returns the single relation $(-2, 1, 1)$. Its top generator
  $t_3$ is an exponential, and the relation is verified by the zero test
  of $-2 + (1 + a) + (1 - a)$ in $\mathbb{Q}(a)$, whose numerator reduces
  to 0. The coefficient of $t_3$ is 1, so $t_3$ becomes algebraic with the
  linear modulus $x - t_1^2/t_2$ (with the option `GR_TOWER_OPT_VERBOSE`,
  the implementation prints `t3 = t1^2/t2`).
- Step 1 again: $x$ now reduces to 0, and the answer is `T_TRUE`.

### Involved generators {#involved}

The zero test, the relation search and root finding consider only the
generators that the element involves, together with the generators that
their definitions and moduli involve. Searching over the whole tower had
two kinds of cost. In long-running contexts a tower accumulates
unrelated generators, and searches over all of them dominate. Worse, a
relation found among unrelated generators may be eliminated in a way
that blocks the relevant one: with $\exp(i)$ in the tower, the test of
$\operatorname{atan}(\tan(260515)) - 260515 + 82924\pi$ found
$\exp(260515i) = \exp(i)^{260515}$, eliminated it, and could not
proceed. Likewise, a conjectural logarithm elsewhere in the tower must
not prevent the conclusion that an element involving only proven
generators is nonzero.

*Example (three of seven generators).* A fixed tower with the generators
$\pi$, $t_1 = e$, $t_2 = \log 3$, $a_1 = \sqrt 5$, $a_2 = \sqrt 2$, $t_3 =
\log(1 + a_2)$, $t_4 = \log(3 + 2a_2)$, in this order, contains $x = 2t_3
- t_4$. The generators that $x$ involves, closed under definitions, are
$a_2$, $t_3$ and $t_4$. The relation search therefore runs on a log side
of two entries, $t_3$ and $t_4$, rather than on the four of the whole
tower (with the argument 1 of $t_1$ and the logarithm $t_2$; there is no
entry $2\pi i$, since $i$ is absent). LLL returns $(-2, 1)$, the relation
is verified as in the example below, and the zero test returns `T_TRUE`;
afterwards $t_4$ has the modulus $x - 2t_3$, and $\pi$, $t_1$, $t_2$ and
$a_1$ are unchanged.

### Integer relations on the log side {#richardson}

The log side of the involved generators consists of the arguments $u$
of exponential generators, the logarithm generators $\log v$, the number
$2\pi i$ when $\pi$ and $i$ are both present, the Lambert $W$ generators
$w = W(z)$ (whose exponential is $z/w$, [](#lambertw)), and, in complex
towers, $iA$ for each angle $A$ of a tangent or arctangent generator ($A
= 2u$ for $\tan u$, $A = 2\operatorname{atan}(u)$), whose exponential is
the point $(1 + iq)/(1 - iq)$ of the unit circle with $q = \tan(u)$ or $q
= u$ ([](#angles)). Each entry has an *exponential side* in the field of
the earlier generators: $e^u$, $v$, $1$, $z/w$, and so on.

Values are computed with adaptive precision: the exponential of a huge
argument needs the argument to absolute rather than relative precision,
so the working precision is doubled until the requested relative
accuracy is reached. An entry that cannot be evaluated is dropped from
the lattice rather than aborting the search. The integer relation search
applies LLL [@LLL1982] to the real and imaginary parts (`acb_lindep`,
which returns all relations of the reduced basis; PSLQ [@FB1992] would
be an alternative) and keeps *every* reduced row that is numerically a
relation, not just the first. Candidates whose coefficients exceed
$\max(3, p/(2n))$ bits at precision $p$ with $n$ entries are discarded.
The remaining rows are brought into Hermite normal form with the columns
ordered by decreasing definition order, so that the leading entry of
each row is the gcd of the possible coefficients of its top generator,
giving the best elimination for that generator; candidates are tried by
increasing level.

*Verification.* Let $g$ be the latest generator occurring in a candidate
$\sum_j c_j z_j = 0$.

- If $g = e^{u_g}$, all entries lie in the field of the generators
  before $g$, and the relation is verified by a recursive zero test of
  $\sum_j c_j z_j$ restricted to that prefix.
- If $g = \log v_g$, the relation is verified multiplicatively, by a
  recursive zero test of $\prod_j \exp(z_j)^{c_j} - 1$ in the prefix,
  which proves $\sum_j c_j z_j \in 2\pi i\mathbb{Z}$; the branch is then
  fixed by checking $|\sum_j c_j z_j| < 2$ numerically (with increasing
  precision, since huge terms may cancel; [](#lem-castle-log)).
- For $W$ on top, the relation $w^* e^{w^*} = z$ is verified below, and
  the branch is certified by a Krawczyk test [@Kra1969] on $we^w - z$
  around both enclosures ([](#lambertw)).

A verification whose products would exceed a cost limit
(`GR_TOWER_OPT_RELATION_COST_LIMIT`) is skipped, and an elimination whose
image cannot be rationalized is not performed.

::: lemma {#lem-castle-log}
Let $S = \sum_j c_j \log z_j + c \cdot 2\pi i$ with integers $c_j, c$.
If $\prod_j z_j^{c_j} = 1$ and $|S| < 2\pi$, then $S = 0$.
:::

::: proof
$e^S = \prod_j z_j^{c_j} = 1$, so $S \in 2\pi i \Z$, and the only
element of $2\pi i \Z$ of absolute value less than $2\pi$ is 0.
:::

*Example (verifying logarithmic relations).* In the example of
[](#involved), the candidate $2t_3 - t_4 = 0$ has the logarithm $t_4$ on
top. It is verified multiplicatively: $(1 + \sqrt 2)^2/(3 + 2\sqrt 2) - 1
= 0$ in $\mathbb{Q}(\sqrt 2)$, which proves $2\log(1 + \sqrt 2) - \log(3 +
2\sqrt 2) \in 2\pi i\mathbb{Z}$; the enclosure of the left side lies in a
disk of radius less than 2 about 0, so the multiple is 0. With $\pi$ and
$i$ (as $\zeta_4$) present and $t = \log i$, the log side is $(2\pi i,
t)$, and LLL at 64 bits returns $(1, -4)$, that is, $2\pi i - 4t = 0$.
The top generator is $t$; the product $\exp(t)^4 \exp(2\pi i)^{-1} - 1 =
i^4 - 1$ is zero, and the numerical check of $|4t - 2\pi i| < 2$ excludes
the other branches of the logarithm. The elimination gives $t$ the
linear modulus $x - \tfrac12\pi i$ (printed `t = i*pi/2` with the
verbose option), and $2t - \pi i$ is zero.

A relation search at the default precision also runs when an exponential
or a logarithm is adjoined, so that $\exp(\log x) = x$ at once, and after
a merge which brings transcendental generators together
([](algorithms/exact/field-towers#merge-devices)).

### Elimination {#elimination}

A verified relation with top generator $g$ is used to make $g$
algebraic over the earlier generators, by changing its kind and giving
it a modulus:

- a logarithm on top: $\log v_g = -\frac{1}{c_g}\sum_j c_j z_j$, a linear
  modulus with rational coefficients;
- an exponential on top with coefficient $\pm 1$: $e^{u_g} = \prod_j
  \exp(z_j)^{\mp c_j}$, a linear modulus;
- an exponential on top with $|c_g| > 1$: the radical modulus $X^{c_g} -
  \prod_j \exp(z_j)^{-c_j}$, a dynamic step since it may factor; when
  $u = r w^c$ with $r$ rational and $w$ a monomial in the generators,
  the modulus is the factor of $x^c - u$ vanishing at $g$, $w^e
  f(x/w)$ with $f$ the factor of $y^c - r$ over $\mathbb{Q}$ (so that
  $\exp(i\pi/3)$ gets the modulus $x^2 - x + 1$);
- or, when the relation is between two exponentials, $c_g u_g + c_j u_j
  = 0$ with $|c_g| > 1$ (as for $\exp(x/2)$ adjoined after $\exp(x)$),
  the *primitive exponential* of the argument $u_j/(c_g/\gcd(c_g, c_j))$
  is inserted as a new generator before both, which then become powers
  of it with linear moduli, so that the degree of the tower does not
  grow.

Eliminated generators are kept, as algebraic generators, so that
elements in the flat representation remain valid: the reduction modulo
the new modulus performs the substitution
([](algorithms/exact/dynamic-evaluation#prop-refine)). The relations are
thus recorded in the tower, shared by all the elements built on it, and
later zero tests need no search. No element is migrated, and no
versioning of elements is needed beyond the lazy conversion between
polynomial contexts. This extends the dynamic evaluation philosophy from
algebraic factorizations to transcendental relations: the tower is
refined in place, and everything computed before remains correct. When
$2\pi i$ is the top entry, $\pi$ and $i$ (which depend on nothing) are
moved to the front of the definition order and the search restarts. The
nested API (`gr_tower_is_zero`) runs the search on a copy of the tower,
because nested elements do not survive a rebuild; the flat and castle
paths run it in place, so that a castle field *learns* relations: once
$\exp(\log 2 + \log 3)$ has been found equal to 6, every element
involving it benefits.

The primitive exponentials matter for trigonometric sums. Without them,
$\sum_{n=1}^{N} \cos(na) = \sin((N + \frac12)a)/(2\sin(a/2)) -
\frac12$ produced a tower whose degree doubled with $N$; with them, the
degree stays 4 up to $N = 16$, and the identity is decided in a few
hundredths of a second.

*Example (radical moduli).* With $t_1 = \log 3$ and $t_2 = e^{t_1/2}$ in
a fixed tower, the zero test of $t_2^2 - 3$ finds the relation $2 \cdot
\frac{t_1}{2} - \log 3 = 0$, with the exponential $t_2$ on top and the
coefficient 2; the modulus of $t_2$ becomes $x^2 - \exp(t_1) = x^2 - 3$,
a dynamic step, so that $t_2$ is now $\sqrt 3$. With $\log 4$ in place of
$\log 3$, the right side 4 is a square, $x^2 - 4$ is factored over
$\mathbb{Q}$ and $t_2$ receives the linear modulus $x - 2$; with $\log
2$ and $t_2 = e^{t_1/3}$, the modulus is $x^3 - 2$. With $\pi$ and $i$
present, $t = \exp(\pi i/3)$ satisfies $2\pi i - 6 \cdot \frac{\pi i}{3}
= 0$ (LLL returns $(1, -6)$), the radical modulus $x^6 - 1$ factors as
$(x - 1)(x + 1)(x^2 + x + 1)(x^2 - x + 1)$, and $t$ receives the factor
$x^2 - x + 1$, which vanishes at $e^{\pi i/3}$. The example of
[](#decidability) shows a radical modulus over a field with
transcendental generators.

*Example (a primitive exponential).* With $a = \sqrt 2$, $t_1 = e^{a/2}$
and $t_2 = e^{a/3}$, the log side $(a/2, a/3)$ satisfies $2 \cdot
\frac a2 - 3 \cdot \frac a3 = 0$, a relation between two exponentials
with $|c_g| = 3$. The primitive exponential $t_3 = e^{a/6}$ is inserted
before both, which become $t_1 = t_3^3$ and $t_2 = t_3^2$ with linear
moduli (the definition order is now $a, t_3, t_1, t_2$), and $t_1^2 -
t_2^3$ reduces to $t_3^6 - t_3^6 = 0$. In a castle field the same
happens as soon as $e^{\sqrt 2/3}$ is created in a context that holds
$e^{\sqrt 2/2}$: the latter then prints as
`t1 {a1 = sqrt(2); t3 = exp(a1/6); t1 = exp(a1/2)}`, still in terms of
its own generator, now of degree one over $t_3$.

Calcium also finds integer relations among logarithms and exponents by
LLL and verifies them exactly ([](#legacy),
[](algorithms/exact/field-towers#ideal)), and this part of the design is
inherited. What is new is the elimination in place, which keeps the
field consistent as relations are discovered after elements were
created, and the systematic reduction of everything to a triangular set
rather than an ideal to be completed.

### Torsion {#torsion}

Arguments with rational parts in $\pi i$, $u_j = y_j + r_j \pi i$, lead
to relations $\sum_j c_j u_j + c_k 2\pi i = 0$ which are resolved
modulo the roots of unity: when $\sum_j c_j r_j + 2c_k = 0$, the
exponentials become $\exp(y_j) \exp(r_j \pi i)$ with a root of unity of
order $\operatorname{lcm}(2 \operatorname{den}(r_j))$, placed first in
the tower. Thus $\exp((16 + 30\pi i)/225)$ and $\exp((5 + 42\pi i)/60)$
are expressed through $\exp(1/900)$ and $\zeta_{60}$ (degree 16), not
through a radical of degree 375 over a root of unity of order 125. When
the rational part is hidden ($\pi$ in a denominator of the argument, as
in the transformation factors of the theta functions), a relation in
which $c_g$ divides the $c_j$ still gives a modulus of degree one when
$\exp(-2\pi i k/c_g)$ is $\pm 1$ or $\pm i$ with $i$ present:
$\exp(u + \pi i/2) = i \exp(u)$ rather than a root of the reducible
$X^2 + \exp(u)^2$.

### Trigonometric generators {#angles}

The real trigonometric generators enter through their *angles*: $A = 2u$
for $\tan(u)$ and $A = 2\arctan(u)$, whose points on the unit circle
$(\cos A, \sin A) = ((1 - q^2)/(1 + q^2), 2q/(1 + q^2))$, with $q =
\tan(u)$ or $q = u$, are rational functions of the generators. When $i$
is present, the angles are on the log side as $iA$, with $e^{iA} = (1 +
iq)/(1 - iq)$, so that relations with the complex exponentials are found
($e^{2i} = (1 + i\tan 1)/(1 - i\tan 1)$). Without $i$, a separate *real
angular search* runs on the angles and $\pi$: a relation $\sum_j c_j
A_j = 0$ is verified directly below its top generator when that is a
tangent, and otherwise through the product of the points, $\prod_j (\cos
A_j + i\sin A_j)^{c_j} = 1$ computed with pairs of real elements, plus
$|\sum_j c_j A_j| < 1$. The elimination is linear for an arctangent
($\arctan(u) = -\sum_{j \ne g} c_j A_j/(2c_g)$) and for a tangent with
a unit coefficient (by the half-angle formula); a tangent with a larger
coefficient is first swapped with another tangent of unit coefficient
when the order allows it (so that $\tan(2x)$ becomes a rational function
of $\tan(x)$), and otherwise the lattice of the angles involved gets a
basis (Hermite normal form) of new generators $\tan(B_r/2)$ inserted
before the tangents involved, which become rational functions of them by
the addition formulas. Thus $4\arctan(1/5) - \arctan(1/239) = \pi/4$,
$\arctan(\tan 2) = 2 - \pi$ and the addition formulas are decided without
complex numbers. The real angular search, its eliminations and its
termination are described in
[](algorithms/exact/castle-real-fields#angular).

### Lambert W {#lambertw}

Lambert $W$ values $w = W_k(z)$ are on the log side as well, with $e^w =
z/w$ in place of the exponential generator, so that their relations with
exponentials and logarithms are found in the same way: $W(2\log 2) =
\log 2$, $W(3e^3) = 3$, $e^{W(z)} = z/W(z)$, $\log W(z) = \log z - W(z)$.
When $w$ is the top generator of a relation, the candidate value $w^*$
below it is verified exactly as a solution, $w^* e^{w^*} = z$, and
identified with the branch of the generator by a Krawczyk test [@Kra1969]
on a ball containing the enclosures of both. The two branches
$W_0(-\log(2)/2) = -\log 2$ and $W_{-1}(-\log(2)/2) = -2\log 2$ are thus
distinguished. Under Schanuel's conjecture the zero test is complete for
towers with $W$ generators as for $\exp$ and $\log$: for algebraic $z \ne
0$, $W_k(z)$ is transcendental by the Lindemann–Weierstrass theorem (if
$w$ were algebraic, so would be $e^w = z/w$). Relations are eliminated
as linear relations only with a unit coefficient on $w$; otherwise $e^w$
would be a radical.

*Example (relations with $W(1)$).* With $w = W(1)$, the products and
functions $we^w = 1$, $e^w = 1/w$ and $\log w = -w$ are returned
directly; in the tower, $\exp(w)$ and $\log(w)$ appear as steps of degree
1 over $\mathbb{Q}(w)$. The value $t = \exp(w/2)$ is not a new
transcendental: since $t^2 = e^w = 1/w$, it is adjoined as an algebraic
step with modulus $t^2 - 1/w$ (not known to be irreducible), and $t^2 w =
1$ is decided. Similarly $W(2\log 2) = \log 2$, $W(\pi e^\pi) = \pi$, the
two real branches at $-\log 2/2$, and $\log W(2) + W(2) - \log 2 = 0$ are
proved.

### Special functions {#special-relations}

The values of other special functions are generators of *conjectural*
status: distinct generators are assumed independent unless a relation is
known, which keeps the zero test sound (a number which is not separated
from zero and has no explaining relation gives `T_UNKNOWN`). The
relations come from identities rather than lattice searches: the zero
test computes the distribution relations between gamma values at
rationals of a level, the reflection and multiplication relations of
gamma and Hurwitz zeta values on rational lines, Legendre's relation and
Landen's transformation between complete elliptic integrals, and the
distribution and five-term relations of polylogarithms, and eliminates a
generator in terms of the others. See
[](algorithms/exact/castle-special-functions).

### The split zero test {#split}

Dynamic-evaluation inversions in a tower of very large nominal degree
can be prohibitively expensive even when the field actually generated is
moderate. The motivating example was the comparison of $\cos(\pi/257)$,
which lies in $\mathbb{Q}(\zeta_{257})$ (degree 256), with its
expression by nested square roots (twelve square roots, nominal degree
4096, actual degree 1024, as two of the radicals lie in the field below
them). The difference lives in the compositum of the two towers, of
nominal degree $2^{20}$, and the nested inversions never finished.

The split zero test applies to an algebraic element of a tower of
nominal degree above a threshold (`GR_TOWER_OPT_SPLIT_DEGREE_LIMIT`,
4096 by default), not all of whose steps are proven, whose numerator
splits as $f + g$ with $f$ and $g$ involving disjoint sets of generators
(connected components of the terms, each generator taken with the
closure of its dependencies). Then $x = 0$ if and only if $f = -g$, and
each part lies in a tower of its own generators (`gr_tower_set_subset`).
Let $P \in \mathbb{Z}[x]$ be the minimal polynomial of the part in the
smaller tower, computed from the flat representation by integer linear
algebra up to a degree limit (`GR_TOWER_OPT_MINPOLY_DEGREE_LIMIT`). Then
$f = -g$ if and only if $P(-g) = 0$ (a zero test in the other tower
alone) and $-g$ is the particular root $f$ of $P$, which is decided by
isolating the roots of $P$ numerically. A sum of three or more parts
with pairwise disjoint generators would need the minimal polynomial of a
sum and is not handled ([](algorithms/exact/castle-design#open)).

*Example (Gauss's construction of $\cos(\pi/257)$).* A variant of the
motivating input is Gauss's construction of $\cos(\pi/257)$ by square
roots with each period named once (127 square roots in 7 levels, about
15 kB of text). It evaluates in a castle field to an element $g$ of a
tower of 10 square roots, of nominal degree 1024, starting with $\sqrt
{17}$, $\sqrt 2$, $\sqrt{1289}$, $\sqrt{257}$ and $\sqrt{257 + \sqrt{257}}$.
The number $\cos(\pi/257)$ lies in $\mathbb{Q}(\zeta_{257})$, of degree
256. The difference $g - \cos(\pi/257)$ lives in the stacked tower of 11
steps and nominal degree $2^{18}$; its numerator splits into the terms
in $\zeta_{257}$ and those in the square roots, and the split test
applies. The smaller tower is $\mathbb{Q}(\zeta_{257})$; the minimal
polynomial $P$ of $\cos(\pi/257)$ has degree $\varphi(514)/2 = 128$,
$P(g) = 0$ is decided in the tower of $g$, and the root isolation of $P$
at 64 bits finds the enclosures of $g$ and $\cos(\pi/257)$ in the same
isolating disk: the answer is `T_TRUE`. The fully written-out
construction (643 kB of text with 15712 square roots) is decided the
same way.

### Precision beyond the separation stage {#precision}

When no structure theorem decides an element, more precision is the only
remaining way to prove it nonzero. The default numerical separation
stage stops at `GR_TOWER_OPT_CERTIFY_PREC_LIMIT` (4096 bits); after the
exact stages have failed, the sign computation and the zero test
continue to `GR_TOWER_OPT_NUMERIC_PREC_LIMIT` (65536 bits), returning
`T_FALSE` when an enclosure excludes zero. A true zero that no relation
explains still ends as `T_UNKNOWN`, after a few more evaluations.

*Example (how far the extra precision reaches).* Since
$\operatorname{erfc}(n) = 1 - \operatorname{erf}(n) \approx
e^{-n^2}/(n\sqrt\pi)$, the number $1 - \operatorname{erf}(n)$ needs about
$n^2 \log_2 e$ bits to be separated from zero (the castle fields
represent $\operatorname{erfc}$ as $1 - \operatorname{erf}$). In fresh
fields, $1 - \operatorname{erf}(60) \approx 3.3 \cdot 10^{-1566}$ (about
5200 bits), $1 - \operatorname{erf}(100) \approx 6.4 \cdot 10^{-4346}$
(about 14500 bits) and $1 - \operatorname{erf}(170) \approx 2.6 \cdot
10^{-12554}$ (about 41700 bits) are decided nonzero, while $1 -
\operatorname{erf}(250) \approx 8.9 \cdot 10^{-27147}$ (about 90200 bits)
is undecided. With `GR_TOWER_OPT_NUMERIC_PREC_LIMIT` set to 4096, $1 -
\operatorname{erf}(60)$ is undecided as well; in the unlimited mode, the
limits grow until it is decided
([](algorithms/exact/effort-limits#measurements)). A general treatment
of numbers that are tiny or huge for structural reasons, such as Taylor
models [@MB2003], is missing
([](algorithms/exact/castle-design#open)).

## What is decided in practice {#decided}

Castle arithmetic decides:

- **inequality of distinct numbers**, given enough precision
  ([](#numerical)), within the precision limits of the context or
  without limit in the unlimited mode;
- **equality of algebraic numbers**, always (within the effort limits:
  the precisions and degrees of refinements and proofs, which the
  unlimited mode lifts);
- **equality of exp-log numbers**, when the relations that make them
  equal are $\mathbb{Q}$-linear relations between logarithms and
  multiplicative relations between exponentials, which is always the
  case under Schanuel's conjecture: Machin-like formulas with 48 terms,
  $\log(\sqrt 2 + 1) + \log(\sqrt 2 - 1) = 0$,
  $\arccos(\cos(\sqrt 2 - 1)) = \sqrt 2 - 1$;
- **identities between special function values** that follow from the
  relations they know: $\Gamma(1/3)\Gamma(2/3) = 2\pi/\sqrt 3$,
  $\psi'(1/4) = \pi^2 + 8G$, $\operatorname{Li}_2(1/2) = \pi^2/12 -
  \log^2(2)/2$, Legendre's relation, $j((1 + \sqrt{-163})/2) =
  -640320^3$, the addition formulas of the theta functions.

An identity that needs a relation which castle arithmetic does not know
(Landen's relations of the trilogarithm, isogenies of elliptic integrals
beyond Landen's transformation, a hypergeometric identity outside its
table) is not decided: the answer is `T_UNKNOWN` in the default mode, and
the test does not terminate in the unlimited mode.

Among the open issues of Calcium's tracker: $a^{50} - a^{51}a^{-1}$ for
$a = 417/(962\pi + 80808)$ is zero (Calcium proves $n = 50$ and returns
`T_UNKNOWN` from $n = 100$); $\exp((\log 2i - \pi i/2)/2)$ is the
canonical $\sqrt 2$; and $\cos(\operatorname{acos}(\sqrt 2 - \sqrt 3)/3)$
is recognized as an algebraic number of degree 12, the root $0.807193\ldots$
of $256x^{12} - 768x^{10} + 864x^8 - 592x^6 + 321x^4 - 90x^2 + 1$: castle
arithmetic computes it through the principal cube root $w$ of $\sqrt{1 -
a^2} + ia$ for $a = \sqrt 2 - \sqrt 3$, where $1 - a^2 = 2\sqrt 6 - 4$ and
$\sqrt{1 - a^2} = \sqrt 2\sqrt{\sqrt 6 - 2}$. Calcium's identity examples
that it reports as undecidable — $\operatorname{acos}(\cos 1) = 1$, tan
and cos identities at $\sqrt{2\pi}$, $\sqrt{e^{2\sqrt 2} + e^{-2\sqrt 2} -
2} = (e^{2\sqrt 2} - 1)/e^{\sqrt 2}$, $\sin 3a = 4\sin a \sin(\pi/3 - a)
\sin(\pi/3 + a)$, and $\exp(\log M) = M$ for a $2 \times 2$ matrix — are
decided in a few milliseconds. On a benchmark of 744 zero-testing and
simplification problems drawn from the literature and from earlier
systems, castle arithmetic decided every case
([](algorithms/exact/castle-design#measurements)).

## Consequences for the API {#api}

Since a zero test can fail, every function whose result depends on one
reports this: predicates return `T_UNKNOWN`, and operations which need a
predicate (division, inversion, comparisons, the generic polynomial and
matrix algorithms of [](gr_poly) and [](gr_mat)) return `GR_UNABLE`.
Comparisons $x < y$ of real numbers compare enclosures at increasing
precision and, when $x$ and $y$ cannot be separated, test $x = y$
exactly. Functions that need many zero tests, such as Gaussian
elimination, choose the algorithm and the order of the tests to keep
them few (see [](algorithms/exact/linear-algebra)).

## Use in FLINT {#use}

- `gr_tower_flat_num_is_zero` and `gr_tower_is_zero` implement the
  castle zero test; `gr_is_zero`, `gr_equal` and all predicates of the
  castle fields are built on it.
- `gr_tower_flat_get_acb` and `gr_tower_get_acb` evaluate elements in
  ball arithmetic.
- The options `GR_TOWER_OPT_PREC_LIMIT`,
  `GR_TOWER_OPT_CERTIFY_PREC_LIMIT`, `GR_TOWER_OPT_NUMERIC_PREC_LIMIT`
  and `GR_TOWER_OPT_RELATION_COST_LIMIT` bound the effort
  ([](gr_tower#tuning-options), [](algorithms/exact/effort-limits)).

## The legacy zero test (Calcium) {#legacy}

The legacy design of [](ca) represents numbers in formal fields with an
ideal of relations built when the field is created
([](algorithms/exact/field-towers#legacy)); its zero test is
`ca_check_is_zero`, with the user-level summary [](ca/check-is-zero).
It proves inequality numerically as above (up to `CA_OPT_PREC_LIMIT`)
and equality by the following methods, in order of increasing cost
[@Joh2021calcium].

### Proving equality {#exact}

::: algorithm Legacy zero test {#alg-zero}
**Input:** a `ca_t` value $x$.

**Output:** `T_TRUE`, `T_FALSE` or `T_UNKNOWN`.

1. *Representation.* If $x$ is a special value, return `T_FALSE`
   (`T_UNKNOWN` for *Unknown*). If $x \in \mathbb{Q}$ or $x$ is an element
   of a field with one algebraic generator, return whether its canonical
   representation is zero.
2. Otherwise $x = p/q$ in a field $\mathbb{Q}(a_1, \ldots, a_n)$; replace
   $x$ by the numerator $p$, which is zero if and only if $x$ is.
3. *Numerical evaluation.* For precisions $64, 128, \ldots$ up to
   `CA_OPT_PREC_LIMIT`: evaluate $p$ in ball arithmetic and return
   `T_FALSE` if the ball excludes zero. At the first precision, if every
   generator can be evaluated as an algebraic number, evaluate $p$ exactly
   with [](qqbar) (`ca_get_qqbar`) and return the answer, unless the
   degree limit `CA_OPT_QQBAR_DEG_LIMIT` is exceeded.
4. *Normal form.* Rewrite the generators recursively in a normal form
   (see below) and evaluate $p$ on the rewritten generators; if the
   result is represented as zero, return `T_TRUE`.
5. *Factoring.* Factor $p$ as a polynomial in the generators; if it has
   at least two nonconstant factors, apply steps 1 to 4 to each factor
   and return `T_TRUE` if one of them is proved zero.
6. Return `T_UNKNOWN`.
:::

::: theorem Correctness of the legacy test {#thm-zero}
If [](#alg-zero) returns `T_TRUE`, then $x = 0$; if it returns
`T_FALSE`, then $x \ne 0$.
:::

::: proof
`T_FALSE` is returned only for a special value or when a ball computed
in ball arithmetic, which is guaranteed to contain $p(a_1, \ldots,
a_n)$, excludes zero, or when the canonical representation in
$\mathbb{Q}$, in a number field, or of an exact algebraic evaluation is
nonzero. `T_TRUE` is returned when an exact representation is zero. The
representations of steps 1, 3 and 4 are elements of formal fields
reduced modulo ideals of proved relations, so by
[](algorithms/exact/field-towers#prop-sound) a zero representation
implies a zero value. The rewriting of step 4 replaces each generator by
an expression with the same value, so the rewritten $p$ has the value of
$p$. In step 5, $p$ is the product of its factors, so it vanishes if one
of them does.
:::

The normal form of step 4 (`ca_rewrite_complex_normal_form`, internal)
makes relations visible that the original generators hide:

- $\sin z$, $\cos z$, $\tan z$ become rational functions of $e^{iz}$, and
  $\operatorname{atan} z$, $\operatorname{asin} z$ and
  $\operatorname{acos} z$ become logarithms;
- $\operatorname{Re} z$, $\operatorname{Im} z$, $\overline{z}$, $|z|$,
  $\operatorname{sgn} z$ and $\arg z$ are expressed with $z$ and its deep
  conjugate $\overline{z}$, with $|z| = \sqrt{z \overline{z}}$ and $\arg z
  = -i \log(z/|z|)$;
- algebraic generators that lie in a cyclotomic field are expressed with
  a root of unity (`qqbar_try_as_cyclotomic`);
- the arguments of $\exp$, $\log$ and $\sqrt{\cdot}$ are rewritten
  recursively.

Creating the fields of the rewritten expression triggers the relation
searches of [](algorithms/exact/field-towers#ideal), which find
multiplicative relations between exponentials and linear relations
between logarithms. In this way, for example, $\cos^2 x + \sin^2 x - 1$
is proved zero even when the trigonometric functions are represented
directly (`CA_TRIG_DIRECT`), and $16 \operatorname{atan}(1/5) - 4
\operatorname{atan}(1/239) - \pi$ (Machin's formula) is proved zero.

The legacy test is sound because the ideal consists of true relations,
and incomplete because the ideal may not contain all of them, in which
case a nonzero representative of a zero number cannot be distinguished
from zero. Typical failures are identities that need a relation
discovered only later (after the field was created), relations that
require factoring a minimal polynomial over a field containing other
algebraic numbers, and transcendental identities in which an exponential
becomes algebraic over other generators.

*Example (where the legacy test fails).* For $\sqrt{5 + 2\sqrt 6} -
\sqrt 2 - \sqrt 3$, the reduced representative is $x_1 - x_3 - x_4 \ne 0$
in the field of [](algorithms/exact/field-towers#ideal), the enclosure
contains zero, and Calcium answers `T_TRUE` through the conversion to
[](qqbar). The same situation with $\pi$ in place of an algebraic number
has no such fallback: for $\sqrt{\pi^2 + 2\sqrt 3\pi + 3} - \pi - \sqrt
3$ the field is $\mathbb{Q}(x_1, x_2, x_3)/(x_1^2 - x_2^2 - 2x_2x_3 - 3,
x_3^2 - 3)$ with $x_2 = \pi$, $x_3 = \sqrt 3$, the representative $x_1 -
x_2 - x_3$ does not reduce to zero, and the answer is `T_UNKNOWN`. For
$\log(e^2e^3) - 5$ the logarithm is an extension number of its own, in
the field $\mathbb{Q}(x_1)$ with no relations, and the answer is again
`T_UNKNOWN`. Castle arithmetic decides both: the first square root is
$\pi + \sqrt 3$ as soon as it is created
([](algorithms/exact/field-towers#express)), and in the second case $e^2$
and $e^3$ are written as powers of $e$, the logarithm is that of $e^5$,
and its generator becomes algebraic with the linear modulus $x - 5$.

### What is decided in practice {#in-practice}

With the methods above, the following questions are always decided,
within the limits set in the context (see [](ca/options)).

- **Inequality of distinct numbers**, given enough precision
  ([](#numerical)). The precision is bounded by the context option `CA_OPT_PREC_LIMIT` (4096 bits by default), so a
  nonzero number smaller than about $2^{-4096}$ in absolute value is not
  recognised as nonzero by default.
- **Equality of algebraic numbers.** A number all of whose generators
  are algebraic can be evaluated exactly as a [](qqbar) number, whose
  minimal polynomial decides equality. This is limited by
  `CA_OPT_QQBAR_DEG_LIMIT` (degree 120 by default). With [](qqbar)
  itself, all predicates are decided without limit, at a cost that grows
  quickly with the degree.
- **Elements of one field.** If $x$ and $y$ belong to the same number
  field, or more generally if their difference is reduced to zero by the
  known relations of their field, equality is recognised formally.

For numbers involving transcendental functions, equality is decided
only sometimes: `ca_check_is_zero` proves $x = 0$ when it finds the
relations that make $x$ vanish, that is, algebraic relations, linear
relations between logarithms, multiplicative relations between
exponentials, the functional equations of some special functions, and
the rewriting of trigonometric functions as exponentials. This succeeds
for many identities that occur in practice, such as
$\sqrt{5 + 2\sqrt{6}} = \sqrt{2} + \sqrt{3}$,
$\log(\sqrt 2 + 1) + \log(\sqrt 2 - 1) = 0$ or Machin's formula
$16 \operatorname{atan}(1/5) - 4 \operatorname{atan}(1/239) = \pi$, and
fails for identities that need deeper knowledge, such as
$\Gamma(1/3) \Gamma(2/3) = 2\pi/\sqrt{3}$. [](ca/check-is-zero) shows these cases in a program and explains what
to do when the answer is `T_UNKNOWN`.

In castle arithmetic, the exact algebraic evaluation by [](qqbar) is
replaced by dynamic evaluation, which has no degree limit and does not
compute minimal polynomials; the normal form of step 4 by the
representation of the trigonometric functions through their angles and
the relation searches on the log side ([](#richardson), [](#angles)),
which run on the generators of one tower rather than on a rewritten
copy; and the factoring of step 5 by the refinement of dynamic steps.
The examples above, including $\Gamma(1/3) \Gamma(2/3) = 2\pi/\sqrt{3}$,
are decided by castle arithmetic ([](#decided)).

### Use in FLINT (legacy) {#legacy-use}

- `ca_check_is_zero` implements [](#alg-zero); `ca_is_zero_check_fast`
  is its first step, and `ca_check_equal`, `ca_check_is_one` and the
  other predicates are built on it.
- `ca_get_acb_raw` and `ca_ext_get_acb_raw` evaluate numbers and
  generators in ball arithmetic; `ca_get_acb` adds precision control and
  uses the zero test to return exact zeros.
- `ca_get_qqbar` performs the exact algebraic evaluation, and
  `ca_factor` the factorisation of step 5.
- [](qqbar) decides equality of algebraic numbers completely; see
  [](algorithms/exact/qqbar-representation#predicates).
