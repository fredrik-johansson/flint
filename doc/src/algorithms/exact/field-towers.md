---
title: Towers of fields
summary: How castle arithmetic represents numbers as elements of towers of field extensions of the rationals: generators and their definitions, the nested and flat representations and the triangular set as a Gröbner basis, structured generators for roots of unity and radicals, merging towers and expressing numbers in them; and the legacy design of Calcium, formal fields with an ideal of relations.
complexity: Arithmetic in the flat representation costs multivariate rational function operations followed by a reduction modulo the triangular set, one pseudo-division per step from the top; merging two towers costs a substitution per element and, for generators not identified by their definitions, a lattice search or a modular proof that they are not in the target.
order: 30
---

A number of castle arithmetic ([](gr_ca)) is an element of a *tower*
([](gr_tower)), a chain of simple extensions of a field of rational
functions over $\mathbb{Q}$:

$$
\mathbb{Q} \subseteq F_0 = \mathbb{Q}(t_1, \ldots, t_r) \subset F_1 \subset
\cdots \subset F_n, \qquad F_k = F_{k-1}[x]/(m_k), \quad a_k = x \bmod m_k,
$$

where the $t_j$ are transcendental generators ($\pi$, $\exp(u)$,
$\log(u)$, $\tan(u)$, $\arctan(u)$, values of special functions, each
with an argument $u$ in a lower field) and the $a_k$ are algebraic
generators with monic moduli $m_k \in F_{k-1}[x]$. Every generator is a
definite complex number, fixed by its definition and a numerical
enclosure, so that the tower is a subfield of $\mathbb{C}$ and
numerical values of its elements are available to any precision
([](gr_tower#numerical-evaluation)). This page describes the towers
themselves; how they are refined and how their zero test works is
described in [](algorithms/exact/dynamic-evaluation) and
[](algorithms/exact/zero-testing#castle). The design of Calcium, which
castle arithmetic succeeds, is described at the end
([](#legacy)).

## Generators and definition order {#generators}

The generators of a tower are kept in their *definition order*, the
order in which they were adjoined, which is compatible with the
dependencies (the modulus and the argument of a generator involve only
earlier ones). The first $p$ generators form a *prefix*, itself a tower;
every number records the shortest prefix in which it lives, which is
what makes merges and the restriction of computations to the relevant
generators cheap. The algebraic generators carry the status *proven*
(the modulus is known to be irreducible over the field below) or
*dynamic* (it merely vanishes at the generator,
[](algorithms/exact/dynamic-evaluation)), and the transcendental ones
the status *independent* (known to be algebraically independent of the
earlier generators: $\pi$ at the head of a tower), *Schanuel* (so under
Schanuel's conjecture: exponentials, logarithms, Lambert $W$) or
*conjectural* (beyond it: the values of special functions).
Statuses decide which shortcuts are permitted; they never decide an
answer by themselves.

Algebraically, with $t_1, \ldots, t_r$ the transcendental generators and
$a_1, \ldots, a_n$ the algebraic ones, each family in definition order,
a tower presents the ring

$$
R_T = F_0[x_1, \ldots, x_n]/(m_1, \ldots, m_n), \qquad F_0 = \mathbb{Q}(t_1, \ldots, t_r),
$$

and the nested presentation $F_k = F_{k-1}[x_k]/(m_k)$, $F_n = R_T$,
which is well defined because the modulus of $a_k$ involves only
generators that come earlier in the definition order (the
transcendental ones among them lie in $F_0$). Every generator has a
complex value — $\pi$, $\exp(\operatorname{val}(u))$, and so on for the
transcendental generators, and the root of $m_k$ isolated by its
enclosure for the algebraic ones — and the *evaluation map*
$\operatorname{ev} : R_T \to \mathbb{C}$ is the ring homomorphism these
values determine, defined on the elements whose denominators
(polynomials in the $t_j$) do not vanish at the values; every element
the implementation constructs has this property, since divisions are
preceded by zero tests. When every step is proven and the transcendental
generators are algebraically independent, $R_T$ is a field and
$\operatorname{ev}$ is injective, so that the reduced representative is
a canonical form. The point of the design is that the zero test does not
need this to hold in advance ([](algorithms/exact/dynamic-evaluation)).

*Example (a tower over $\mathbb{Q}(\pi)$).* In a fresh castle field,
$\sqrt{\pi + \sqrt 2}$ prints as `a2 {a1 = sqrt(2); a2 = sqrt(a1+pi)}`.
Its tower has the generators $\pi$ (independent), $a_1 = \sqrt 2$ with
modulus $x_1^2 - 2$ and $a_2 = \sqrt{\pi + a_1}$ with modulus $x_2^2 -
(\pi + a_1)$, in this definition order, both steps proven, so that

$$
R_T = \mathbb{Q}(\pi)[x_1, x_2]/(x_1^2 - 2, \; x_2^2 - x_1 - \pi)
$$

is a field of degree 4 over $\mathbb{Q}(\pi)$. The definition order need
not list the transcendental generators first: $\sqrt{2/\pi}$ is
represented as $\sqrt 2 \cdot \sqrt{1/\pi}$ in a tower with the
generators $\pi$, $\sqrt{1/\pi}$, $\sqrt 2$, the step of $\sqrt{1/\pi}$
with the modulus $x^2 - 1/\pi$, whose leading coefficient after clearing
denominators is $\pi$. In a tower containing $\pi$, $\sqrt 2$,
$e^{\sqrt 2}$ and $\Gamma(1/3)$, the status of $e^{\sqrt 2}$ is
Schanuel and that of $\Gamma(1/3)$ conjectural.

A generator also carries a *definition id*, which identifies the same
number across the towers of a castle field: two towers which contain
generators with the same definition id contain the same number, and
merges identify them. Generators which a relation search eliminates are
kept, as algebraic generators of degree one (or radicals) over the
others, so that elements are never invalidated
([](algorithms/exact/zero-testing#elimination)).

## Nested and flat representations {#representations}

A tower has two coexisting representations of its elements, converted
into each other through coordinate vectors (sparsely: a dense coordinate
vector of size $\prod_k \deg m_k$ is impossible for the cyclotomic
towers that occur in practice, whose degrees can exceed $10^{18}$).

In the *nested* representation, an element of $F_k$ is a polynomial of
degree less than $\deg m_k$ over $F_{k-1}$, recursively down to $F_0$,
whose elements are rational numbers (`fmpq`) or rational functions
([](fmpz_mpoly_q)); this is the generic quotient ring of
[](gr_domains#poly-quotient) iterated. Inversion is by the extended
Euclidean algorithm when the base has cheap inverses, and otherwise
(when the base is itself a quotient ring of small degree) through the
Cayley–Hamilton theorem applied to the multiplication matrix, computed
by Berkowitz's division-free algorithm [@Ber1984], so that one
inversion at level $k$ costs one inversion at level $k - 1$. The nested
representation is the *structural* one: it carries the dynamic
evaluation machinery (zero divisors are found by nested inversions),
coordinates, multiplication matrices, characteristic and minimal
polynomials, the conversion to [](qqbar), and the gcd certificates of
Trager's method.

*Example (a nested inverse).* Let $F_1 = \mathbb{Q}(a_1)$ with $a_1 =
\sqrt 2$ and $F_2 = F_1[x]/(x^2 - 1 - a_1)$, so that $a_2 = \sqrt{1 +
\sqrt 2}$. A nested element of $F_2$ is a polynomial of degree at most 1
in $a_2$ whose coefficients are polynomials of degree at most 1 in $a_1$:
$(a_1 + a_2)^3 = (6 + 5a_1) + (7 + a_1) a_2$ and $a_2^5 = (3 + 2a_1)
a_2$. Multiplication by $B = a_1 + a_2$ has the matrix

$$
M_B = \begin{pmatrix} a_1 & 1 + a_1 \\ 1 & a_1 \end{pmatrix}, \qquad
\det M_B = a_1^2 - a_1 - 1 = 1 - a_1, \qquad
\operatorname{adj}(M_B) e_0 = \begin{pmatrix} a_1 \\ -1 \end{pmatrix}
$$

over $F_1$ in the basis $1, a_2$, that is, $(a_1 + a_2)(a_1 - a_2) = 1 -
a_1$; the extended Euclidean algorithm applied to $x^2 - 1 - a_1$ and $x
+ a_1$ meets the same remainder $1 - a_1$. Either way, the inverse costs
one inversion in $F_1$, $1/(1 - a_1) = -(1 + a_1)$, and $1/(a_1 + a_2) =
-(1 + a_1)(a_1 - a_2) = (-2 - a_1) + (1 + a_1) a_2$.

Castle fields compute in the *flat* representation, the representation
of Calcium's multivariate fields, chosen for the same reason: sparse
multivariate arithmetic in FLINT is fast. An early experiment settled
the question of which representation should carry arithmetic: for the
determinant of a $9 \times 9$ matrix over a tower of degree 512, nested
arithmetic took 1.75 s, flat arithmetic 0.42 s and Calcium 0.47 s; for
dense products at degree 4096, the nested representation took 2.8 s
against 0.43 s (timings from development, indicative only). The castle
fields therefore store flat elements, and the nested form is
constructed on demand, mostly for zero tests that need inversions. In
the flat representation,
an element is a fraction $p/q$ of polynomials with integer coefficients
in variables $X_1, \ldots, X_N$ for all the generators, transcendental
and algebraic, reduced modulo the *triangular set*
$\{M_1, \ldots, M_n\}$ of the moduli written as polynomials in these
variables ($M_k$ is $m_k$ with the denominators of its coefficients
cleared, so that its leading coefficient is a polynomial in the
transcendental variables only).

The variable of the generator with definition order $d$ is $X_{c-1-d}$
for a capacity $c$ of the polynomial context, in FLINT's lexicographic
order $X_0 > X_1 > \cdots$, so that later generators are more
significant, which makes the triangular set a lexicographic Gröbner
basis whatever the kinds of the generators are, and the tower can grow
up to the capacity without renumbering. When the capacity is exceeded,
or when generators are reordered, a new context is made. Every context,
current or superseded, records for each variable the stable identifier
of the generator it represents, so that elements created in old
contexts remain convertible after reorderings and eliminations;
conversion is a renaming of variables, done lazily when an element is
next used.

*Example (layout).* For the tower $\mathbb{Q}(\sqrt 2)(\sqrt{1 + \sqrt
2})$ above, a fixed tower made by `gr_tower_init` has capacity 8: $a_1$
(definition order 0) is the variable $X_7$ and $a_2$ is $X_6$, and the
triangular set is $\{X_7^2 - 2, \; X_6^2 - X_7 - 1\}$, whose leading
monomials $X_7^2$ and $X_6^2$ are coprime, as [](#prop-groebner)
requires. The flat forms are the nested coordinates multiplied out.
Adjoining $\pi$ afterwards gives it the definition order 2 and the
variable $X_5$, below the capacity, so that no element is renamed and
the two flat moduli are unchanged; the nested fields, on the contrary,
are rebuilt over the new base field $\mathbb{Q}(\pi)$, from the flat
forms of the moduli.

::: proposition The triangular set is a Gröbner basis {#prop-groebner}
With this ordering, the leading monomial of $M_k$ is $X_{(k)}^{d_k}$,
where $X_{(k)}$ is the variable of $a_k$ and $d_k = \deg m_k$, times a
monomial in the transcendental variables. As polynomials over the field
$\mathbb{Q}(X_{t_1}, \ldots, X_{t_r})$ of the transcendental variables,
the $M_k$ have pairwise coprime leading monomials, and thus form a
Gröbner basis of the ideal they generate for the lexicographic order.
:::

::: proof
By Buchberger's criterion, the S-polynomial of two polynomials whose
leading monomials are coprime reduces to zero modulo the two
[@CLO2015, §2.9, Proposition 4]; here every S-polynomial is of this
kind.
:::

Reduction modulo the triangular set is top-down pseudo-division by
$M_n, \ldots, M_1$ in turn, with respect to the main variable of each
(`fmpz_mpoly_univar_pseudo_rem`). Multiplying by a leading coefficient,
a polynomial in the transcendental variables, never raises a degree in
an algebraic variable, and a reduction by a modulus never raises the
degree in its own variable or in a later one, so one pass from the top
down terminates with a reduced numerator; the accumulated leading
coefficients go into the denominator, which remains a polynomial in the
transcendental variables only. (An earlier version used
`fmpz_mpoly_quasidivrem_ideal`, which is not a pseudo-division with
respect to a main variable and silently failed to reduce modulo moduli
such as $tX^2 - 1$: with $P = X^3 + 1$, no term of $P$ is divisible by
the leading monomial $tX^2$, so $P$ was returned unchanged, while the
pseudo-remainder is $t^2 P \equiv tX + t^2$.) Monic univariate integer
moduli (cyclotomic polynomials, $x^n - p$) are reduced as univariate
polynomials with `fmpz_poly_rem`.

*Example (normal forms).* In the tower over $\mathbb{Q}(\pi)$ of the
example of [](#generators), the leading monomials of $m_1 = x_1^2 - 2$
and $m_2 = x_2^2 - x_1 - \pi$ are $x_1^2$ and $x_2^2$, and the normal
forms are the $\mathbb{Q}(\pi)$-linear combinations of $1, x_1, x_2,
x_1 x_2$. Reducing $x_2^4$ by $m_2$ gives $(x_1 + \pi)^2$, and then by
$m_1$, $x_2^4 \equiv \pi^2 + 2\pi x_1 + 2$, which prints as
`2*a1*pi+pi^2+2`. The inverse of $x_2$ is $(\pi - x_1) x_2/(\pi^2 -
2)$, since $x_2^2 (\pi - x_1) \equiv (\pi + x_1)(\pi - x_1) \equiv \pi^2
- 2$; its denominator lies in $\mathbb{Q}[\pi]$. For the modulus $\pi
x_3^2 - 1$ of $\sqrt{1/\pi}$, with leading coefficient $\pi$,
pseudo-division gives $\pi x_3^3 \equiv x_3$, so that $\sqrt{1/\pi}^3 =
\sqrt{1/\pi}/\pi$.

In a reduction, a numerator whose degree in each algebraic variable is
already below that of the modulus is returned before any
canonicalization (which would cost a multivariate gcd), and sums of
reduced elements are not reduced again: each element records the
version of the ideal for which it is known to be reduced. With $a_2 =
\sqrt{1 + \sqrt 2}$ as above, the reduction of $a_2^5$ starts with the
top modulus: $a_2^5 = a_2 (a_2^2)^2 \equiv (1 + a_1)^2 a_2 = (a_1^2 + 2a_1
+ 1) a_2$, which raises the degree in $a_1$ to 2, so the reduction
continues with $x_1^2 - 2$ and ends with $(2a_1 + 3) a_2$. The product
$(a_1 + a_2)(a_1 - a_2)$ has degree 2 in both variables and needs both
moduli; a product of two elements that are linear in $a_2$ and constant
in $a_1$, such as $(1 + a_2)(2 + a_2) = 2 + 3a_2 + a_2^2$, needs only
the top one, and a sum of reduced elements none.

Reduction modulo the triangular set therefore gives a numerator which is
unique up to a factor in the transcendental variables: the reduced
numerator of an element of a tower whose moduli are irreducible is zero exactly when the element is
zero. The denominator is not canonical, but castle fields *rationalize*
it, so that it involves only the transcendental generators: a monomial
in the algebraic generators through their moduli, a polynomial in the
highest algebraic generator through the adjugate of its multiplication
matrix modulo the modulus, within a size limit
(`GR_TOWER_OPT_RATIONALIZE_LIMIT`); a larger denominator is left as it
is, which is still a valid fraction. In a proven tower, the rationalized
fraction is canonical
([](algorithms/exact/castle-representations#canonical)).

*Deferred moduli.* Splitting towers
([](algorithms/exact/castle-roots#splitting)) have moduli with very many
terms: the Cauchy modulus of step $k$ of the splitting tower of a
polynomial of degree $n$ has $\binom{n+1}{k}$ terms, its coefficients
being complete homogeneous symmetric polynomials in the earlier roots.
Such a modulus is converted to the flat ideal only when a reduction
needs it (`GR_TOWER_OPT_DEFERRED_IDEAL_SIZE`), from the top down: a
reduction by step $k$ only raises degrees at lower steps, so the set of
moduli needed is discovered as the reduction proceeds. Workloads that
only use traces never expand the large moduli.

*Example (the Cauchy moduli of a quintic).* The polynomial $f = x^5 +
2x^4 + 3x^3 + x^2 + 5x + 7$ is irreducible with Galois group $S_5$. Its
roots $r_1, \ldots, r_5$ generate a splitting tower whose step $k$ has
as modulus the $k$-th Cauchy module, of degree $6 - k$ in $a_k$ and with
$\binom{6}{k} = 6, 15, 20, 15, 6$ terms for $k = 1, \ldots, 5$; the
second is

$$
a_2^4 + (a_1 + 2) a_2^3 + (a_1^2 + 2a_1 + 3) a_2^2 + (a_1^3 + 2a_1^2 +
3a_1 + 4) a_2 + (a_1^4 + 2a_1^3 + 3a_1^2 + 4a_1 + 5).
$$

In a castle field, the sum $r_1 + r_2 + r_3 + r_4$ lives in the tower of
the first four steps, of degree 120; with the default limit of 12
coefficients, the moduli of steps 2, 3 and 4 are not converted to the
flat ideal, and the first one has constant coefficients and is a
univariate polynomial. Adding $r_5 = -2 - (r_1 + \cdots + r_4)$ gives
$-2$ without converting any of them. The element $r_4^2 + r_1 + r_2 +
r_3 + r_4$, of degree 2 in $a_4$ only, converts the modulus of step 4 and
is reduced by it to a polynomial of degree at most 2 in $a_1, a_2, a_3$,
below the degrees of their moduli, so that steps 2 and 3 stay deferred.

*Rebasing.* Adjoining a transcendental generator changes $F_0$ and
therefore every nested field; the nested chain is rebuilt from the
moduli transported through the flat representation, which is
unaffected.

In the flat representation inversion is cheap (the fraction is swapped,
and the algebraic denominator is rationalized), sparse expressions stay
sparse and products use [](fmpz_mpoly) arithmetic; the castle fields add
further representations for common cases
([](algorithms/exact/castle-representations)).

## Structured generators {#structured}

Which numbers become generators, and in which form, decides the size of
the towers, the sparsity of the representation and which relations come
for free. Castle fields apply the following conventions (*generator
policies*), some of them configurable
(`gr_ca_ctx_set_gen_flags`, `GR_TOWER_OPT_CYCLOTOMIC_DEGREE_LIMIT`).

### Roots of unity of prime power orders {#roots-of-unity}

A root of unity of order $n = q_1 q_2 \cdots q_s$ (pairwise coprime prime
powers $q_i$) is written as

$$
\zeta_n = \prod_{i=1}^{s} \zeta_{q_i}^{c_i}, \qquad c_i \equiv (n/q_i)^{-1} \pmod{q_i},
$$

which holds because $\sum_i c_i n/q_i \equiv 1 \pmod n$ by the Chinese
remainder theorem. The fields $\mathbb{Q}(\zeta_{q_i})$ are linearly
disjoint, so the steps $\Phi_{q_i}$ are irreducible over the earlier
ones (proven by a structural rule,
[](algorithms/exact/dynamic-evaluation#structural)), and an element of
$\mathbb{Q}(\zeta_n)$ is a sparse polynomial in the prime-power
generators instead of a dense polynomial of degree $\varphi(n)$; a
monomial in roots of unity is a product of one reduced power of each
$\zeta_{q_i}$, which is a monomial except when a power $\zeta_{q_i}^e$
with $q_i$ odd and $e \ge \varphi(q_i)$ is reduced modulo $\Phi_{q_i}$. A
root of unity of order $\ell^e$ is a power of an existing $\zeta_{\ell^f}$
with $f \ge e$; if only $f < e$ exists, $\zeta_{\ell^e}$ is adjoined *in
front* and the old generator becomes its power, with a linear modulus.
Fieker's sparse representation of cyclotomic fields [@Fie2007] is a
related idea. The decomposition handles towers that no dense
representation could: in the exact DFT benchmark (`examples/dft.c`,
input 3), the fields have degrees up to $7 \cdot 10^{18}$ and are
handled as sparse polynomials in 27 variables.

*Example (roots of unity of orders 15 and 16).* For $n = 15 = 3 \cdot 5$,
$c_1 \equiv 5^{-1} \equiv 2 \pmod 3$ and $c_2 \equiv 3^{-1} \equiv 2
\pmod 5$, so $\zeta_{15} = \zeta_3^2 \zeta_5^2$, which prints as
`-a1*a2^2-a2^2` ($a_1 = \zeta_3$, $a_2 = \zeta_5$) since $\zeta_3^2 = -1
- \zeta_3$ is reduced modulo $\Phi_3$. Powers are products of one
reduced power of each generator: $\zeta_{15}^4 = \zeta_3^2 \zeta_5^3 =
-(1 + \zeta_3)\zeta_5^3$ has two terms, and $\zeta_{15}^7 = \zeta_3^2
\zeta_5^4 = (1 + \zeta_3)(1 + \zeta_5 + \zeta_5^2 + \zeta_5^3)$ eight;
in the dense power basis of $\mathbb{Q}(\zeta_{15})$, $\zeta_{15}^7$
would be a monomial and $\zeta_{15}^8$ a polynomial with six terms. For
the order 16: the element $u = \pi + \sqrt 3 + i$ lives in the tower
$\mathbb{Q}(i)(\sqrt 3)(\pi)$; in $u\zeta_{16} = \pi\zeta_{16} + \sqrt
3\zeta_{16} + \zeta_{16}^5$, the tower becomes

$$
\zeta_{16} : x^8 + 1, \qquad i : x - \zeta_{16}^4, \qquad \sqrt 3 : x^2 - 3 \;(\text{dynamic}), \qquad \pi,
$$

with $\zeta_{16}$ in front, $i$ a linear step, and the proof of $x^2 - 3$
dropped because the field below it changed
([](algorithms/exact/dynamic-evaluation#reordering)).

The decomposition has a price: products in $\mathbb{Q}(\zeta_{255}) =
\mathbb{Q}(\zeta_3, \zeta_5, \zeta_{17})$ use a Kronecker array of length
$3 \cdot 7 \cdot 31 = 651$ rather than 255 (but see the cyclic layout of
[](algorithms/exact/castle-representations#number-fields)). Calcium
represents $\mathbb{Q}(\zeta_N)$ as one number field. The option
`GR_TOWER_OPT_CYCLOTOMIC_DEGREE_LIMIT` makes $\zeta_N$ for composite $N$
a single generator while $\varphi(N)$ is within the limit; measurements
during development showed it 1.7–2 times faster for
$\mathbb{Q}(\zeta_{15})$, $\mathbb{Q}(\zeta_{60})$,
$\mathbb{Q}(\zeta_{120})$, 1.3–1.9 times slower for
$\mathbb{Q}(\zeta_{105})$, $\mathbb{Q}(\zeta_{210})$,
$\mathbb{Q}(\zeta_{255})$, and prohibitive for workloads with many
orders (a DFT benchmark went from 0.4 s to 13 s), so it is off by
default. With the limit 16, $e^{2\pi i/15}$ is $\zeta_{60}^4$ in
$\mathbb{Q}(\zeta_{60})$ of degree 16 (the argument $2\pi i/15$
involves $i$, and $\varphi(\operatorname{lcm}(4, 15)) = 16$ is within the
limit); by default it is $-\zeta_3\zeta_5^2 - \zeta_5^2$ in
$\mathbb{Q}(\zeta_5)(\zeta_3)$ of degree 8.

### Radicals {#radicals}

Real radicals of positive rationals are products of principal roots of
primes: $\sqrt{a/b}$ and $\sqrt[n]{a/b}$ are factored (by trial
division, with word-size cofactors factored completely; a cofactor
beyond that remains a single opaque radicand), so that $\sqrt 6 = \sqrt
2 \sqrt 3$ and $\sqrt{12} = 2\sqrt 3$ hold without any search, and
Mordell's theorem proves the degrees. Roots of prime power orders are
likewise decomposed, as for roots of unity. The square root of a
negative rational $-AB^2$ ($A$ squarefree) is $B\sqrt{-A}$ with
$\sqrt{-A}$ a single generator, a root of $x^2 + A$, except that
$\sqrt{-1} = i$ and $\sqrt{-3} = 2\zeta_3 + 1$ are roots of unity; so
$\mathbb{Q}(\sqrt{-163})$ has degree two rather than being
$\mathbb{Q}(i, \sqrt{163})$. Generator flags restore $i\sqrt A$
(`GR_TOWER_GENS_SPLIT_IMAGINARY`) or make $\sqrt A$ for composite $A$ a
single generator (`GR_TOWER_GENS_COMPOSITE_RADICALS`, Calcium's
convention). A principal root of a monomial
in radical generators is again a radical generator ($\sqrt[3]{\sqrt 2} =
2^{1/6}$), and the roots of $c$ times a monomial in roots of unity are
roots of unity times $\sqrt[n]{c}$ ($\sqrt[3]{\zeta_5^3} = \zeta_{15}^{-2}$,
with the branch determined numerically). Rational contents are split off
($\sqrt{34 + 2\sqrt{17}} = \sqrt 2 \sqrt{17 + \sqrt{17}}$); this is not
always the best choice ([](algorithms/exact/castle-design#open)).

*Example (radical generators).* In a fresh castle field, $\sqrt{3/8} =
\sqrt 2 \sqrt 3/4$ and $\sqrt[3]{12} = \sqrt[3]{2}^2 \sqrt[3]{3}$, with the
generators $\sqrt 2$, $\sqrt 3$, $\sqrt[3]{2}$, $\sqrt[3]{3}$. The radical
$2^{1/6}$ of order $6 = 2 \cdot 3$ is decomposed like $\zeta_6$: with
$1/6 = 1/2 + 2/3 - 1$,

$$
2^{1/6} = \frac{2^{1/2} (2^{1/3})^2}{2},
$$

and $\sqrt[3]{\sqrt 2}$ gives the same element, so that $\sqrt 2$,
$\sqrt[3]{2}$ and $2^{1/6}$ share the tower $\mathbb{Q}(\sqrt
2)(\sqrt[3]{2})$ of degree 6. A radical of 2-power order of a radical
generator is a new radical generator, $\sqrt{\sqrt 2} = \sqrt[4]{2}$.

*One number field, one generator.* Algebraic numbers that lie in one
number field should enter the tower through one generator; otherwise
the compositum is never discovered except by expensive dynamic
evaluation. When the sines of rational multiples of $\pi$ in the
formulas for $\Gamma(p/q)$ were separate algebraic generators, the tower
could not see that they generate a common field: four values at level 33
produced towers of nominal degree about $10^4$ that never finished,
against 30 ms once each sine was written through one root of unity per
level. The same principle shaped the tangent normal form of real
trigonometric constants
([](algorithms/exact/castle-real-fields#tangent-form)), the
prime-power decomposition of radicals in eliminations, and the Gauss sums
below. For instance, a castle field writes $2^{1/74} = 2^{1/2}
(2^{1/37})^{19}/2$ (as above, $1/74 = 1/2 + 19/37 - 1$), so that the
three radicals $2^{1/74}$, $2^{1/2}$, $2^{1/37}$ live in $\mathbb{Q}(\sqrt
2)(\sqrt[37]{2})$ of degree 74. In a fixed tower, adjoining $2^{1/74}$ as
a separate (dynamic) root of $x^{74} - 2$ over $\mathbb{Q}(\sqrt
2)(\sqrt[37]{2})$ gives the nominal degree $2 \cdot 37 \cdot 74 = 5476$:
the zero test of $a_3^{37} - \sqrt 2$ refines the step to $x^{37} -
\sqrt 2$ (degree 2738), and only the zero test of $a_3 - \sqrt
2(2^{1/37})^{19}/2$ makes it linear (degree 74). Until such refinements
happen, the tower is a product of fields; in the eliminations of the
zero test this cost 5–40 s per test during development.

### Gauss sums {#gauss-sums}

When roots of unity are present, a square root $\sqrt c$ may already
lie in the field they generate, and adjoining it as a new step would
create a reducible step. It is written as a Gauss sum,

$$
\sqrt{p^*} = \sum_{a=1}^{p-1} \left(\frac{a}{p}\right) \zeta_p^a, \qquad
p^* = (-1)^{(p-1)/2} p, \qquad \sqrt 2 = \zeta_8 + \zeta_8^{-1}, \qquad
\sqrt{-1} = \zeta_4,
$$

with the sign chosen numerically [@Was1997]. This is used by merges, by
the exact zero test (where a dynamic step $\sqrt c$ after the roots of
unity becomes of degree one) and by the eliminations of gamma values.

*Example (square roots next to $\zeta_5$ and $\zeta_7$).* With $(a/5) =
1, -1, -1, 1$ for $a = 1, 2, 3, 4$ and $\zeta_5^4 = -1 - \zeta_5 -
\zeta_5^2 - \zeta_5^3$,

$$
\sqrt 5 = \zeta_5 - \zeta_5^2 - \zeta_5^3 + \zeta_5^4 = -1 - 2\zeta_5^2 - 2\zeta_5^3,
$$

so that $\zeta_5 + \sqrt 5$ is the element $-2\zeta_5^3 - 2\zeta_5^2 +
\zeta_5 - 1$ of $\mathbb{Q}(\zeta_5)$ (in a field without real views; see
[](algorithms/exact/castle-real-fields#real-first) for the real-first
form). Likewise $\sqrt{-7} = \zeta_7 + \zeta_7^2 + \zeta_7^4 - \zeta_7^3
- \zeta_7^5 - \zeta_7^6 = 1 + 2\zeta_7 + 2\zeta_7^2 + 2\zeta_7^4$, and
$\zeta_7 + \sqrt{-7}$ prints as `2*a1^4+2*a1^2+3*a1+1` ($a_1 = \zeta_7$).
But $\sqrt 7 \notin \mathbb{Q}(\zeta_7)$, and $\zeta_7 + \sqrt 7$ has
two generators.

### Logarithms of rational and Gaussian rational numbers {#logarithms}

The logarithm of a rational or Gaussian rational number is decomposed
over the basis

$$
\{\log p : p \text{ prime}\} \cup \{\log(a + bi) : p = a^2 + b^2 \text{ prime}\} \cup \{\pi i\},
$$

with one Gaussian prime $a + bi$ per split prime $p \equiv 1 \pmod 4$,
its conjugate being $\log(a - bi) = \log p - \log(a + bi)$, and $1 + i$
expressed through $\log 2$ and $\pi i$ since $(1 + i)^2 = 2i$. The
representation $p = a^2 + b^2$ is found by Cornacchia's algorithm
[@Cor1908; @Coh1996] (a square root of $-1$ modulo $p$ followed by a
Euclidean reduction), and the factorization in $\mathbb{Z}[i]$ uses
FLINT's Gaussian integers. The branch, a multiple of $\pi i/4$, is
determined numerically. Logarithms of pairwise non-associate Gaussian
primes are linearly independent over $\mathbb{Q}$ together with $\pi i$,
by unique factorization in $\mathbb{Z}[i]$: a relation would make a
product of powers of distinct primes a root of unity. The basis above is
obtained from such logarithms by an invertible linear change, so it is
independent too, and the relation search when towers merge
([](#merging)) is skipped when every transcendental generator is $\pi$
or such a logarithm. Calcium factors integers in logarithms with a
smoothness bound (its option `CA_OPT_SMOOTH_LIMIT`) to find such
relations; the Gaussian decomposition makes the Machin-like formulas for
$\pi$ immediate: all the formulas of up to 48 terms tabulated in FLINT's
`machin` example are proved in a few seconds (Calcium proves those with
at most six terms).

Conversely, $\exp\bigl(\sum_j c_j \log u_j + r\pi i\bigr)$ with rational
$c_j$, $r$ and algebraic $u_j$ is a product of principal powers
$u_j^{c_j}$ and the root of unity $e^{r\pi i}$; for instance $\exp((\log
2i - \pi i/2)/2)$ is the canonical $\sqrt 2$.

*Example (Gaussian logarithms and Machin's formula).* The split prime $5
= 2^2 + 1^2$ is represented by the Gaussian prime $2 + i$, and

$$
\begin{aligned}
\log 12 &= 2\log 2 + \log 3, & \log(1 + 2i) &= \log 5 - \log(2 + i) + \tfrac{\pi i}{2}, \\
\log 2i &= \log 2 + \tfrac{\pi i}{2}, & \log(\tfrac34 + i) &= 2\log(2 + i) - 2\log 2, \\
\log(1 - 2i) &= \log(2 + i) - \tfrac{\pi i}{2}, & \log(7 + 24i) &= 4\log 5 - 4\log(2 + i) + \pi i,
\end{aligned}
$$

since $1 + 2i = i(2 - i)$, $3 + 4i = (2 + i)^2$ and $7 + 24i = -(2 -
i)^4$, while $\log 5$ and $\log 13$ remain generators. For Machin's
formula, $5 + i = (1 + i)(3 - 2i)$ and $1 + 239i = (1 - i)(3 + 2i)^4$
give

$$
\operatorname{atan}\tfrac15 = \tfrac{\pi}{4} - \tfrac{i}{2}\bigl(\log 13 - 2\log(3 + 2i)\bigr), \qquad
\operatorname{atan}\tfrac{1}{239} = \tfrac{3\pi}{4} - 2i\bigl(\log 13 - 2\log(3 + 2i)\bigr),
$$

so that $4\operatorname{atan}\frac15 - \operatorname{atan}\frac{1}{239}$
reduces to $\pi/4$ by linear algebra on the coordinates, without a
relation search.

### Calcium's and castle generators compared {#generators-compared}

The following table shows the representations in fresh contexts, with
the minimal polynomial (Calcium) or the definition (castle arithmetic)
of each generator.

| Number | Calcium | Castle arithmetic (defaults) |
|---|---|---|
| $\sqrt{-163}$ | $a$, $a^2 + 163 = 0$ | $a_1 = \sqrt{-163}$ |
| $\sqrt{-12}$ | $4a + 2$, $a^2 + a + 1 = 0$ | $4a_1 + 2$, $a_1 = \zeta_3$ |
| $\sqrt 6$ | $a$, $a^2 - 6 = 0$ | $a_1 a_2$, $a_1 = \sqrt 2$, $a_2 = \sqrt 3$ |
| $\sqrt{24}$ | $2a$, $a^2 - 6 = 0$ | $2a_1 a_2$ |
| $\zeta_{12}$ | $a$, $a^4 - a^2 + 1 = 0$ | $-i a_1$, $a_1 = \zeta_3$ |
| $\zeta_5$ | $a^2$, $a^4 - a^3 + a^2 - a + 1 = 0$ | $a_1 = \zeta_5$ |
| $2^{1/3}$ | $a = $ `Pow(2, 1/3)` (symbolic) | $a_1 = \sqrt[3]{2}$, modulus $x^3 - 2$ |

Calcium represents $\zeta_5$ through $a = e^{\pi i/5}$, a generator of
$\mathbb{Q}(\zeta_{10}) = \mathbb{Q}(\zeta_5)$, and $2^{1/3}$, computed as
a power rather than through `ca_set_qqbar`, as a symbolic extension
number. Castle arithmetic splits radicals of composite integers into
radicals of primes and roots of unity into roots of unity of prime power
orders; with the generator flags and the cyclotomic degree limit it
makes Calcium's choices, $\sqrt{24} = 2\sqrt 6$ and $\zeta_{12}$ as one
generator.

## Merging towers {#merging}

### Definitions and the registry {#registry}

Every generator created in a castle field has a *definition* (kind,
parameter, argument) and a *definition id* that identifies it in all
towers containing it, together with a name unique in the context
($a_1, a_2, \ldots$ for algebraic generators, $t_1, t_2, \ldots$ for
exponentials and logarithms, $\pi$, $i$). A tower is identified by the
set of definition ids of its generators, plus *aliases*: definitions
that the tower represents by an expression in its generators (after an
absorption, $\sqrt 5$ may be represented in $\mathbb{Q}(\zeta_5)$ by
$-2\zeta_5^3 - 2\zeta_5^2 - 1$). The *registry* of the context holds the
towers; a binary operation looks for the smallest existing tower
containing the definitions of both operands (a $5 \times 5$ determinant
over several algebraic numbers went from 226 towers to 10 with this
rule during development), and otherwise merges.

*Example (an alias).* In a fresh castle field, $\sqrt 5$ is the generator
$a_1$ of $\mathbb{Q}(\sqrt 5)$ and $\zeta_5 = e^{2\pi i/5}$ the generator
$a_2$ of $\mathbb{Q}(\zeta_5)$. Their sum prints as
`-2*a2^3-2*a2^2+a2-1 {a2 = exp(2*pi*i/5)}`: the merged tower is
$\mathbb{Q}(\zeta_5)$ of degree 4, with the definition of $a_1$ as an
alias, $\sqrt 5 = -2\zeta_5^3 - 2\zeta_5^2 - 1$.

Constants are hash-consed: $e = \exp(1)$, $\log 2$, $\pi$, roots of
unity and radicals of primes each exist once per context (in registries
of canonical roots of prime power orders and of constant transcendental
generators), so that independently computed values share generators.
Each root of prime power order is the generator of a tower of its own:
an element is created in the smallest field that contains it and moves
to a larger one only when it meets elements of that field, as in
Calcium ($\sqrt 2$ computed after $2^{1/128}$ is the generator of
$\mathbb{Q}(\sqrt 2)$, and becomes $a^{64}$ only in a merged tower with
$a = 2^{1/128}$). An earlier version returned the root of a smaller order
as a power of the root of the largest order of the same kind created so
far, so that a long-running context placed later values in ever larger
fields. A value which is not a canonical constant is not reduced at
creation: in a fresh context, $\exp(2)$ is the generator of a tower of
its own and $\exp(1)$ the generator of another; in the merged tower of
$\exp(1) + \exp(2)$, the relation search finds $\exp(2) = \exp(1)^2$, and
the generator $\exp(2)$ receives the modulus $x - e^2$.

### Absorption {#absorption}

An operation on two numbers of different towers needs a tower containing
both. A castle field first looks for one: when the tower of one operand
contains the definitions of all the generators of the other, it is used
(unless generators of other computations lie between them, which the
result would otherwise carry). Otherwise the towers are *merged*
(`gr_tower_absorb`): the generators of one tower $B$ are adjoined to the
other, $U$, in definition order, and a homomorphism from $B$ to $U$ is
built, through which the operand of $B$ is mapped by substitution
(`gr_tower_map_apply_flat`).

Each generator of $B$ is, in this order:

- *identified* with a generator of $U$ with the same definition id, or
  with an existing generator by its definition (roots of unity,
  radicals), or by having the same origin polynomial and an enclosure
  meeting the same isolating disk of its roots (a purely numerical but
  rigorous rule, which matters because exact searches in targets of
  degree $10^7$ to $10^9$ are hopeless), or by a structured image
  ([](#structured));
- *substituted* when its modulus has become linear over $U$ (its
  expression is the root of that linear factor);
- *factored* over the base field of $U$ when the image of its modulus
  has coefficients there, the factor vanishing at the generator
  replacing the modulus;
- *expressed* in $U$, by a lattice search for a root of its modulus
  among the elements of $U$ (`gr_tower_express_limit`, verified
  exactly), skipped when a place of $U$ shows that the modulus has no
  root in $U$ (`gr_tower_poly_no_roots_modular`) or when $U$ is too large
  (`GR_TOWER_OPT_EXPRESS_DEGREE_LIMIT`);
- otherwise *adjoined*, as a dynamic step with its modulus mapped (and
  factored over the target when the image has coefficients in
  $\mathbb{Q}(t)$), and transcendental generators with the image of their
  argument.

A map between towers (`gr_tower_map_t`) stores the flat image of each
source generator in the target; flat images survive rebasing and
eliminations in the target.

*Example (three merges).* (1) $\sqrt 2 + \pi$ is formed in the merged
tower of $\mathbb{Q}(\sqrt 2)$ (degree 2, modulus proven irreducible) and
$\mathbb{Q}(\pi)$: a tower of degree 2 over $\mathbb{Q}(\pi)$, the step
$\sqrt 2$ keeping its proven status. (2) After $\zeta_{1024}$, $\zeta_8$
is created as the generator of $\mathbb{Q}(\zeta_8)$ of degree 4, and
$\zeta_8 \zeta_{1024}$ is computed in the tower of $\zeta_{1024}$, of
degree 512, where $\zeta_8 = \zeta_{1024}^{128}$, so that the product is
$\zeta_{1024}^{129}$. Likewise $\sqrt 2$ created after $2^{1/128}$ is the
generator of $\mathbb{Q}(\sqrt 2)$, and $\sqrt 2 \cdot 2^{1/128} = a^{65}$
with $a = 2^{1/128}$. (3) Let $a \approx 3.1463$ be given as a root of
$x^4 - 10x^2 + 1$ (so $a = \sqrt 2 + \sqrt 3$, which the context does not
know). In $a\sqrt 2$, $\mathbb{Q}(\sqrt 2)$ is absorbed into
$\mathbb{Q}(a)$: the generator $\sqrt 2$ has no definition in
$\mathbb{Q}(a)$ and no structured image, the quadratic search finds no
quadratic steps to use, and the modular test does not exclude a root of
$x^2 - 2$ in $\mathbb{Q}(a)$; the integer relation search then finds
$\sqrt 2 = (a^3 - 9a)/2$, verified exactly. The product is $(a^2 - 1)/2$,
and the tower is still $\mathbb{Q}(a)$ of degree 4.

### Conjugate roots {#conjugates}

When two conjugate roots of one integer polynomial $p$ meet, the second
is adjoined as a root of the Cauchy modulus $p(x)/(x - \alpha_1)$ over
the first, without a proof of irreducibility, so that a *splitting
tower* $\alpha_1, \ldots, \alpha_{d-1}$ emerges as roots are combined,
the last root being given by Vieta's formula
(`GR_TOWER_OPT_SAME_ORIGIN_DEGREE_LIMIT`); symmetric functions of the
roots then reduce to the coefficients
([](algorithms/exact/castle-roots#splitting)).

### Keeping merges cheap {#merge-devices}

Several devices keep merging cheap; each was motivated by a workload
that was slow without it during development.

- *In-place growth.* A tower with more than six generators grows in
  place when the algebraic generators of the absorbed prefix are already
  present, rather than being copied into a new tower; sums of many
  logarithms otherwise produced a combinatorial number of towers
  (Machin-like formulas: 148 s to 3 s).
- *Prefix copies.* A transcendental generator whose argument lies at
  level $\ell$ of a tower with at least $\ell + 4$ generators is adjoined
  to a fresh copy of the $\ell$-prefix, reusing the definition id and
  name of an equal generator in the long tower. Without this, knowledge
  accumulated in a long-running context slows down new computations: one
  trigonometric identity took 0.48 s after sixteen others had filled the
  context with tangents, and 1 ms with prefix copies. Radicals are still
  adjoined to the tower of their argument, since the expressions and
  factorizations found for them there are too costly to rediscover,
  unless the tower continues with at least twelve generators beyond the
  prefix.
- *Subset towers.* Root finding and the split zero test work over the
  tower of a dependency-closed subset of generators
  (`gr_tower_set_subset`), whose statuses carry over (irreducibility
  over a larger field implies irreducibility over a subfield).
- *Merge-time relation search.* When a merge adds transcendental
  generators to a tower that has some, a relation search at default
  precision runs at once, so that relations are not rediscovered by
  expensive zero tests over unreduced elements. It is skipped when all
  transcendental generators are $\pi$ or logarithms of rational or
  Gaussian primes, which are provably independent ([](#logarithms)).
- *Comparison by enclosures.* Equality tests and comparisons between
  elements of different towers first compare enclosures computed in each
  tower; disjoint enclosures decide without merging. Sorting the roots
  of a Mignotte polynomial had tried to build composita needing 37 GB;
  it now takes about a second.

The maps between towers survive the growth of the target, since flat
elements remain valid elements of a tower as it grows. Optionally
(`GR_TOWER_OPT_PRIMITIVE_DEGREE_LIMIT`), a merged tower of several
proven algebraic steps of small degree is replaced by a single
extension $\mathbb{Q}(\theta)$ for a primitive element $\theta$, whose
elements have the dense form of [](algorithms/exact/castle-representations#dense).
This made element products 2–6 times faster in $\mathbb{Q}(\sqrt 2, \sqrt
3, \sqrt 5, \sqrt 7)$ during development, but caused coefficient growth
with cube roots (41 ms to 24 s on one benchmark) and slow nested
inversions once transcendental generators joined, so it is off by
default.

## Expressing numbers in a tower {#express}

Whether a given algebraic number lies in a tower is the question behind
square roots, polynomial roots and merges. Castle fields answer it from
cheap to expensive.

*Square roots through quadratic steps.* If $F_k = F_{k-1}(a)$ with $a^2
= r$ (after completing the square), then $x = u + v a$ is a square in
$F_k$ if and only if $N = u^2 - r v^2$ is a square $s^2$ in $F_{k-1}$
and $(u \pm s)/2$ is a square $p^2$ there, and then $\sqrt x = p + va/(2p)$;
this recurses to the base field, where squares of rational numbers and
of rational functions are tested directly (`gr_tower_sqrt`). No lattice
reduction is involved, and nested radicals such as $\sqrt{5 + 2\sqrt 6}
= \sqrt 2 + \sqrt 3$ are denested exactly.

*Lattice search.* For a root of a polynomial $q$ isolated by an
enclosure, a candidate element of the tower is found by an integer
relation search over the monomial basis of the tower [@LLL1982], then
verified exactly ($q(x) = 0$ by the complete zero test) and identified
numerically with the root (`gr_tower_express`). The search tries at most
`GR_TOWER_OPT_EXPRESS_ROUNDS` precisions: a relation is nearly always
found at the first or second precision, while a search which finds none
runs to the limit, its last rounds being the costliest.

*Proofs of absence.* At a *place* of the tower, a homomorphism to a
finite field given by a chain of simple roots of the moduli, a root of
$q$ in the tower would map to a root of the reduction of $q$; a place at
which the reduction has no root proves that $q$ has none in the tower
(relative to the irreducibility of the moduli), and saves the search
([](algorithms/exact/dynamic-evaluation#places)).

## Use in FLINT {#use}

- `gr_tower_init`, the `gr_tower_adjoin_*` functions and
  `gr_ctx_init_tower_field` build towers and compute in their fields;
  the castle fields of [](gr_ca) build them automatically.
- `gr_tower_flat_reduce`, `gr_tower_flat_rationalize` and the flat
  contexts (`gr_tower_flat`) implement the flat representation.
- `gr_tower_absorb`, `gr_tower_merge` and the maps
  (`gr_tower_map_apply_flat`) merge towers; `gr_tower_sqrt` and
  `gr_tower_express` express numbers in towers.

## The legacy design: formal fields (Calcium) {#legacy}

The [](ca) module (formerly the library Calcium) represents a number as
an element of a field $K = \mathbb{Q}(a_1, \ldots, a_n) \subset
\mathbb{C}$ generated by finitely many *extension numbers* $a_k$:
algebraic numbers, constants such as $\pi$, and values $f(z)$ of
functions such as $\exp$, $\log$ or $\sqrt{\cdot}$ at numbers $z$ that
are themselves elements of such fields [@Joh2021calcium]. Unlike a
tower, a formal field is built at once for a given list of generators,
with its ideal of relations, and cached by that list; its algebraic
generators are algebraic numbers in canonical form ([](qqbar)). Castle
arithmetic replaces the formal fields by towers grown and refined in
place, the ideal by the triangular set of the moduli together with the
relations eliminated by the zero tests, and the canonical algebraic
numbers by dynamic evaluation. The legacy design is described in the
following subsections.

### Formal fields {#formal}

Evaluation at $(a_1, \ldots, a_n)$ is a ring homomorphism $\mathbb{Q}[X_1,
\ldots, X_n] \to \mathbb{C}$. Its kernel

$$
I = \{ g \in \mathbb{Q}[X_1, \ldots, X_n] : g(a_1, \ldots, a_n) = 0 \}
$$

is a prime ideal, the *ideal of relations* of the generators, and

$$
K \cong \operatorname{Frac}\bigl(\mathbb{Q}[X_1, \ldots, X_n] / I\bigr).
$$

If $I$ were known, $K$ could be computed with exactly, using Gröbner
basis reduction modulo $I$ [@CLO2015]. In general, however, $I$ is not
known: for $K = \mathbb{Q}(\pi, e)$ it is not even known whether $I =
\{0\}$. FLINT therefore works with a computed ideal $J \subseteq I$,
generated by relations that have been *proved*, and represents an
element of $K$ by a fraction $p/q$ of polynomials with integer
coefficients, reduced modulo $J$. Two consequences are fundamental.

::: proposition Soundness {#prop-sound}
If $p \in J$, then $p(a_1, \ldots, a_n) = 0$. Hence an element whose
reduced numerator is zero is zero, and two elements with equal reduced
representations are equal. If $J = I$, the converse holds for
numerators: $p(a_1, \ldots, a_n) = 0$ if and only if $p$ reduces to
zero modulo a Gröbner basis of $J$.
:::

::: proof
$J \subseteq I$ gives the first statement. If $J = I$ and $G$ is a
Gröbner basis of $J$, then $p \in J$ if and only if the normal form of
$p$ with respect to $G$ is zero [@CLO2015, §2.6].
:::

The denominator of every element must be a polynomial that does not
vanish at the generators; this is guaranteed because a division is only
performed after the divisor has been proved nonzero (see
[](algorithms/exact/zero-testing)). Even with $J = I$ the
fraction $p/q$ is not unique: $1/\sqrt{2}$ and $\sqrt{2}/2$ are different
fractions in $\mathbb{Q}(X)/(X^2 - 2)$. Fields with a single algebraic
generator avoid this, since FLINT stores their elements as polynomials
of degree less than the degree of the generator ([](nf_elem)), which is
canonical.

*Example (Calcium fields).* $\sqrt 2$ and $\sqrt 3$ are elements of the
number fields $\mathbb{Q}(\sqrt 2)$ and $\mathbb{Q}(\sqrt 3)$, but their
sum lies in the generic field $\mathbb{Q}(x_1, x_2)/(x_1^2 - 3, x_2^2 -
2)$ with $x_1 = \sqrt 3$, $x_2 = \sqrt 2$, and prints as
`3.14626 {a+b where a = 1.73205 [a^2-3=0], b = 1.41421 [b^2-2=0]}`. Likewise $\pi +
\sqrt 2 \in \mathbb{Q}(x_1, x_2)/(x_2^2 - 2)$ with $x_1 = \pi$. In castle
arithmetic the same numbers are elements of the towers $\mathbb{Q}(\sqrt
2)(\sqrt 3)$ and $\mathbb{Q}(\pi)(\sqrt 2)$, printed
`a2+a1 {a1 = sqrt(2); a2 = sqrt(3)}` and `a1+pi {a1 = sqrt(2)}`.

### Extension numbers and fields {#extensions}

An extension number (see [](ca_ext)) is either a `qqbar_t` algebraic
number or a function code with arguments: $\pi$, $\gamma$, $e^z$, $\log
z$, $\sqrt{z}$, $z^w$, $|z|$, $\operatorname{sgn} z$, $\operatorname{Re}
z$, $\operatorname{Im} z$, $\overline{z}$, $\lfloor z \rfloor$, $\lceil z
\rceil$, $\arg z$, the direct trigonometric functions, $\Gamma(z)$,
$\operatorname{erf}(z)$, $\operatorname{erfc}(z)$ and
$\operatorname{erfi}(z)$, where the arguments are `ca_t` numbers of other
fields. The *depth* of an extension number is 0 for algebraic numbers
and constants and one more than the maximal depth of the generators of
the arguments otherwise, so fields form towers ordered by depth.

Each context keeps a cache of extension numbers without structural
duplicates and a cache of fields keyed by their (ordered) lists of
generators. An operation on two numbers in different fields first
*merges* the fields (`ca_merge_fields`): the generator lists are merged
and the elements are mapped to the common field by substituting
variables. The generators are sorted in the reverse of the structural
order of `ca_ext_cmp_repr`, in which algebraic numbers come first and
then numbers by increasing depth. Combined with a lexicographic monomial
order (`CA_OPT_MPOLY_ORD`), the variable of a generator such as $e^{\pi}$
is then larger than the variable of $\pi$, so reduction modulo $J$
eliminates more complicated generators in favour of simpler ones.

After each operation the result is *condensed* (`ca_condense_field`): a
result that does not depend on any generator is demoted to
$\mathbb{Q}$, and one that depends on a single algebraic generator is
moved to the corresponding number field. Superfluous generators are
otherwise kept.

### Building the ideal {#ideal}

When a field with several generators is created, `ca_field_build_ideal`
collects relations of the following kinds. Each one is proved before it
is added to $J$.

*Algebraic generators.* An algebraic generator $a_k$ with minimal
polynomial $P_k$ gives the relation $P_k(X_k)$. If all $d$ conjugates of
an algebraic number of degree $d \le$ `CA_OPT_VIETA_LIMIT` occur as
generators $a_{i}, \ldots, a_{i+d-1}$, Vieta's formulas

$$
c_d \, e_j(X_i, \ldots, X_{i+d-1}) - (-1)^j c_{d-j} = 0, \qquad j = 1, \ldots, d,
$$

are added, where $e_j$ is the elementary symmetric polynomial and $P =
\sum c_j t^j$. The number of terms grows like $2^d$, hence the limit.

*Roots and powers.* A generator $u = \sqrt{z}$ or $u = z^{a/b}$ with a
rational exponent $a/b$ satisfies $u^2 = z$, respectively $u^b = z^a$
(for the principal branch, $(z^{a/b})^b = e^{a \log z} = z^a$). If $z =
p/q$ can be expressed in the current field, the relation $q^a X_u^b -
p^a$ (with $p$ and $q$ swapped for negative $a$) is added.

*Logarithms.* For the generators $\log z_1, \ldots, \log z_r$, together
with $2\pi i$ if $\pi$ and $i$ are both generators, integer relations
$\sum_j c_j \log z_j + c \cdot 2 \pi i = 0$ with small coefficients are
searched with LLL [@LLL1982] from numerical values at `CA_OPT_LLL_PREC`
bits. A candidate relation is proved as follows.

::: lemma {#lem-log}
Let $S = \sum_j c_j \log z_j + c \cdot 2 \pi i$ with integers $c_j, c$.
If $\prod_j z_j^{c_j} = 1$ and $|S| < 2 \pi$, then $S = 0$.
:::

::: proof
$e^{S} = \prod_j z_j^{c_j} = 1$, so $S \in 2 \pi i \Z$, and the only
element of $2 \pi i \Z$ of absolute value less than $2\pi$ is 0.
:::

The bound $|S| < 2$ is checked in ball arithmetic and the product
$\prod_j z_j^{c_j} = 1$ is checked exactly with `ca_check_is_one`, in a
field where the logarithms do not occur. The relation $\sum_j c_j X_{\log
z_j} + 2c X_\pi X_i$ is then added. This finds, for example, $\log 8 = 3
\log 2$ and $\log(1 + \sqrt{2}) + \log(\sqrt{2} - 1) = 0$.

*Exponentials, powers and roots.* For the generators $g_j$ of the forms
$e^{z_j}$, $\sqrt{z_j}$, $z_j^{w_j}$ and algebraic numbers, logarithms
$L_j$ with $e^{L_j} = g_j$ are formed ($z_j$, $\frac{1}{2} \log z_j$, $w_j
\log z_j$ and $\log g_j$), and integer relations $\sum_j c_j L_j + c
\pi i = 0$ are searched with LLL. Such a relation implies the
multiplicative relation $\prod_j g_j^{c_j} = (-1)^c$, which is added as
the binomial $\prod_{c_j > 0} X_j^{c_j} \mp \prod_{c_j < 0} X_j^{-c_j}$.
It is proved exactly: by [](qqbar) arithmetic if all the $g_j$ involved
are algebraic, and otherwise by proving $\sum_j c_j L_j + c \pi i = 0$
with `ca_check_is_zero` (which in turn may use relations between the
logarithms). This finds, for example, $e^{z} e^{-z} = 1$ and $\sqrt{6} =
\sqrt{2} \sqrt{3}$, so that $\sqrt{2} \sqrt{3} - \sqrt{6}$ evaluates to
zero as soon as the three square roots meet in one field.

*Special functions.* For error functions, the relations implied by
$\operatorname{erf}(-z) = -\operatorname{erf}(z)$, $\operatorname{erfc}(z)
= 1 - \operatorname{erf}(z)$ and $\operatorname{erfi}(z) = -i
\operatorname{erf}(iz)$ between generators with equal, opposite or
$i$-related arguments are added (`ca_field_build_ideal_erf`). For the
gamma function, $\Gamma(z + n) = z (z+1) \cdots (z+n-1) \Gamma(z)$ is
added for generators whose arguments differ by an integer $|n| \le 10$
(`ca_field_build_ideal_gamma`).

*Gröbner basis.* Finally, if `CA_OPT_USE_GROEBNER` is set and the
relations have moderate degree, the generators of $J$ are replaced by a
reduced Gröbner basis computed with Buchberger's algorithm [@Buc1965]
(`fmpz_mpoly_buchberger_naive_with_limits`), subject to limits on the
number, length and coefficient size of the polynomials. Reduction
modulo a Gröbner basis gives normal forms, so that elements of $J$
reduce to zero; if the limits are exceeded, the unreduced relations are
kept, which is still sound.

*Example (reduction ideals).* The field of $\log 4 + \log 2$ is
$\mathbb{Q}(x_1, x_2)/(x_1 - 2x_2)$ with $x_1 = \log 4$, $x_2 = \log 2$,
and the element reduces to $3x_2 = 3\log 2$. The field of $e^2 + e^3$ is
$\mathbb{Q}(x_1, x_2)/(x_1^2 - x_2^3)$ with $x_1 = e^3$, $x_2 = e^2$, from
the relation $2 \cdot 3 - 3 \cdot 2 = 0$ between the arguments; the
element stays $x_1 + x_2$. For $\sqrt{5 + 2\sqrt 6}$, an extension number
$x_1 = \operatorname{Sqrt}(5 + 2x_2)$ with $x_2 = \sqrt 6$, next to $x_3 =
\sqrt 3$ and $x_4 = \sqrt 2$, the ideal is

$$
I = (x_1^2 - 2x_3x_4 - 5, \; x_3^2 - 3, \; x_4^2 - 2, \; x_2 - x_3x_4),
$$

which contains the relation $\sqrt 6 = \sqrt 3\sqrt 2$ found among the
algebraic numbers, but not $x_1 - x_3 - x_4$: modulo $x_3^2 - 3$ and
$x_4^2 - 2$, the first generator is a root of $(x_1 - x_3 - x_4)(x_1 +
x_3 + x_4)$, and $I$ is not prime. In castle arithmetic, $e^2$ and $e^3$
are powers of the generator $e$, so that $e^2 + e^3$ is a polynomial in
$e$, and $\log 4 + \log 2$ is $3 \log 2$; the square root $\sqrt{5 +
2\sqrt 6}$ is $\sqrt 2 + \sqrt 3$ as soon as it is created
([](#express)).

### Arithmetic and reduction {#reduction}

In a field with several generators, arithmetic operates on fractions of
multivariate polynomials ([](fmpz_mpoly_q)). After each operation, the
numerator and the denominator are reduced modulo $J$ by multivariate
pseudo-division (`_ca_mpoly_q_reduce_ideal`), the fraction is
canonicalised by removing the polynomial GCD, and monomial factors of
the denominator that can be moved to the numerator with the help of a
relation are moved (`_ca_mpoly_q_simplify_fraction_ideal`); for
example, if $e^{x}$ and $e^{-x}$ are both generators, so that $e^{x}
e^{-x} = 1$ is known, $1/e^{x}$ becomes $e^{-x}$. Elements of
number fields use [](nf_elem) arithmetic, and elements of $\mathbb{Q}$
use [](fmpq).

### Completeness {#completeness}

For a field with one algebraic generator, $J = I$ and the
representation is canonical. With several algebraic generators, the
minimal polynomials together with the relations found as above
generate $I$ in the common cases (such as square roots of integers),
but not always; this does not affect zero tests, which evaluate
algebraic numbers exactly with [](qqbar) when needed.

For transcendental generators, completeness is out of reach in general.
Schanuel's conjecture [@Lan1966] asserts that for $\mathbb{Q}$-linearly
independent complex numbers $z_1, \ldots, z_n$, the field $\mathbb{Q}(z_1,
\ldots, z_n, e^{z_1}, \ldots, e^{z_n})$ has transcendence degree at least
$n$; it would imply, for example, that $e$ and $\pi$ are algebraically
independent, and that the relations between exponentials and logarithms
of algebraic numbers are those that follow from $\mathbb{Q}$-linear
relations between the logarithms [@CP1978]. The relation searches above
are designed to find exactly those relations, so under Schanuel's
conjecture they are complete for such fields, up to the size limits of
the LLL searches. For other functions (gamma, error functions) only the
listed functional equations are known to the implementation, so that,
for instance, the reflection formula $\Gamma(z) \Gamma(1-z) =
\pi/\sin(\pi z)$ is not used.

### Use in FLINT {#legacy-use}

- `ca_field_build_ideal`, `ca_field_build_ideal_erf` and
  `ca_field_build_ideal_gamma` build the ideal; they are called by
  `ca_field_cache_insert_ext` when a new field is created.
- `ca_merge_fields`, `ca_condense_field`, `_ca_mpoly_q_reduce_ideal` and
  `_ca_mpoly_q_simplify_fraction_ideal` implement field merging,
  condensation and reduction; every arithmetic function of [](ca) uses
  them.
- The options `CA_OPT_VIETA_LIMIT`, `CA_OPT_LLL_PREC`,
  `CA_OPT_USE_GROEBNER` and the Gröbner limits control the work done
  (see [](ca/options)); `ca_ctx_print` shows the fields of a context
  with their ideals.
