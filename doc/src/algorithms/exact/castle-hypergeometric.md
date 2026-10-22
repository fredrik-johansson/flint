---
title: Hypergeometric functions in castle arithmetic
summary: The normal form of generalized hypergeometric values pFq(a; b; z) in castle arithmetic — degenerate cases, transformations choosing a canonical member of an orbit, the contiguity module expressing every value of a coset of parameters through a basis, a declarative table of known values and summation theorems, and links by quadratic and Landen transformations to generators already present.
complexity: A value in a coset of $p + q$ parameters costs products of $r \times r$ matrices ($r = \max(p, q+1)$) of rational functions along a path from the canonical cell, plus the table entries matched in a box of shifts, each evaluated only when its relation is independent of the previous ones.
order: 51
---

Castle arithmetic ([](gr_ca)) represents a value of a generalized
hypergeometric function exactly in terms of a few *basis values* per
coset of parameters, closed forms where they are known, and generators
of the kind `GR_TOWER_HYPGEOM` otherwise. This page describes the
normal form; the other special functions are described in
[](algorithms/exact/castle-special-functions).

## The normal form {#normal-form}

The values $F(P; z) = {}_pF_q(a; b; z)$, $P = (a_1, \ldots, a_p; b_1,
\ldots, b_q)$, are brought to a normal form in three layers
(`src/gr_ca/hypgeom.c`):

1. Degenerate parameters: $z = 0$, cancellation of $a_i = b_j$,
   terminating series ($a_i$ a nonpositive integer), poles ($b_j$ a
   nonpositive integer: `GR_DOMAIN`; the regularized function is
   finite there, $\tilde F = \prod_i (a_i)_{n+1} z^{n+1} \tilde F(a + n
   + 1; n + 2, b' + n + 1; z)$ for $b_1 = -n$), divergence ($p > q +
   1$), Gauss's sum at $z = 1$, ${}_0F_0 = e^z$, ${}_1F_0 = (1-z)^{-a}$,
   and the reduction of the order when $a_i - b_j = m$ is a positive
   integer: $F = (\theta + b_j)_m F_{\mathrm{red}} / (b_j)_m$,
   $\theta = z\, d/dz$, with $F_{\mathrm{red}}$ the function without
   $a_i$, $b_j$ (expanded exactly through $\theta^k = \sum_l S(k, l) z^l
   D^l$ and $D^l F = (a)_l / (b)_l\, F(P + l e)$, $e = (1, \ldots, 1)$).

2. Transformations of one term, choosing a canonical member of the orbit
   of $(P, z)$: for ${}_2F_1$, Pfaff's transformation $F(a, b; c; z) =
   (1-z)^{-a} F(a, c-b; c; z/(z-1))$ when $|z - 1| > 1$ (or $|z - 1| = 1$
   and $\operatorname{Im}(z) < 0$), then the smaller (in a fixed order of
   the cells) of $P$ and its Euler transform $(c-a, c-b; c)$ (factor
   $(1-z)^{c-a-b}$), off the cut $z \ge 1$; for ${}_1F_1$, Kummer's
   transformation $F(a; b; z) = e^z F(b-a; b; -z)$ to
   $\operatorname{Re}(z) \le 0$. Powers with an irrational exponent are
   written $x^{e_0} x^n$, $e = e_0 + n$, $0 < \operatorname{Re}(e_0) \le
   1$, so that contiguous parameters share one generator.

3. The contiguity module of the coset $P + \mathbb{Z}^{p+q}$. With $r =
   \max(p, q+1)$ and $v(P) = (F, \theta F, \ldots, \theta^{r-1} F)$, a
   unit shift of a parameter acts as $v(P') = M v(P)$ with $M = I +
   C/a_i$ (raising $a_i$) or $I + C/(b_j - 1)$ (lowering $b_j$), where
   $C$ is the companion matrix of $\theta$ modulo the hypergeometric
   operator $\theta \prod_j (\theta + b_j - 1) - z \prod_i (\theta +
   a_i)$; $\det(I + C/a_i)$ is a multiple of $\prod_j (b_j - 1 - a_i)$, so
   the shifts are invertible unless the operator is reducible. Every
   value in the coset is thus a rational combination of the $r$ basis
   values $G_k = F(P_0 + k e)$ at the canonical cell $P_0$ (parameters
   shifted to $0 < \operatorname{Re} \le 1$). A reducible pair ($b_j -
   a_i$ a positive integer) keeps $b_j - a_i = 1$ in the cell, and the
   paths stay in the chamber $b_j - a_i \ge 1$; it gives the relation
   $(\theta + a_i) F(P_0) = a_i F_{\mathrm{red}}(P_0)$ between the basis
   values. A numerator $a_i = 1$ (a positive integer, in its cell) gives
   another: the operator factors as $\theta L_1$ with $L_1 = \prod_j
   (\theta + b_j - 1) - z \prod_{k \ne i} (\theta + a_k)$, so that $L_1 F
   = \prod_j (b_j - 1)$.

## Known values {#table}

Known values give further linear relations between the basis values:
the entries of a declarative table, each a pattern (parameters affine
in symbols, the argument free or a rational constant) with a value in a
small expression language and optional conditions, for example

```text
"3 2 | a, b, c | 1+a-b, 1+a-c | 1 | gamma(1+a/2)*gamma(1+a-b)*... | re: 1+a/2-b-c"
"0 1 | | 1/2 | z | cosh(2*sqrt(z))"
"2 1 | 1/2, 1/2 | 1 | z | 2*elliptic_k_gen(z)/pi"
```

Every entry is matched against every point of the coset (in a box of
shifts of its symbols) and of the images of the coset under the
transformations of layer 2; the relations are reduced in order
(reducibility relations, then entries with a constant argument, then by
distance from the cell), each value being computed only when its
relation is independent of the previous ones, until $F(P)$ is
determined or the relations are exhausted. The basis values which
remain undetermined become generators. Thus the two entries for
${}_0F_1(; 1/2; z)$ and ${}_0F_1(; 3/2; z)$ give every ${}_0F_1(; n +
1/2; z)$ in terms of $\cosh(2\sqrt z)$ and $\sinh(2 \sqrt z)$, the
entries for $K$ and $E$ give every ${}_2F_1$ with parameters in $(1/2,
1/2, 1) + \mathbb{Z}^3$, the error function entry gives every
${}_1F_1(n + 1/2; m + 3/2; z)$ ($m \ge n$, with the reducibility
relation $(\theta + 1/2) F = e^z/2$), and the summation theorems at $z =
\pm 1, 1/2$ apply to every contiguous function of a summable one (the
"contiguous" variants of Kummer's, Gauss's and Bailey's theorems). The
current table contains the elementary cases of ${}_0F_1$, ${}_1F_1$
and ${}_2F_1$, the error function, Kummer's second formula ${}_1F_1(a;
2a; z) = e^{z/2} {}_0F_1(; a + 1/2; z^2/16)$, $K$ and $E$, the sums of
Gauss ($z = 1$), Kummer ($z = -1$), Gauss and Bailey ($z = 1/2$), and
the sums of Dixon, Watson and Whipple for ${}_3F_2$ at $z = 1$ (the
latter only at the point itself: the contiguity module at $z = 1$,
where the rank drops to $p - 1$, is not implemented). Each entry is checked
numerically by the tests.

## Links between arguments {#links}

Some entries relate the coset to the values of another
hypergeometric function at another argument; they are used only when a
generator of the context is at that argument (an `anchor:` condition,
tested numerically up to Pfaff's transformation, before the exact
evaluation), and one at a time (not inside the evaluation of another
link), so that the generators of the context are reused rather than
duplicated, and no chain of transformations is explored. These are the
quadratic transformations of ${}_2F_1$,

$$
F(a, b; 2b; z) = (1 - z/2)^{-a} F(\tfrac{a}{2}, \tfrac{a+1}{2}; b + \tfrac12; \tfrac{z^2}{(2-z)^2}) \quad (z \notin [1, \infty)),
$$

$$
F(a, b; 1 + a - b; z) = (1 + z)^{-a} F(\tfrac{a}{2}, \tfrac{a+1}{2}; 1 + a - b; \tfrac{4z}{(1+z)^2}) \quad (|z| < 1),
$$

$$
F(a, b; \tfrac{a+b+1}{2}; z) = F(\tfrac{a}{2}, \tfrac{b}{2}; \tfrac{a+b+1}{2}; 4z(1-z)) \quad (\operatorname{Re}(z) < 1/2),
$$

and their inverses (with the principal roots $\sqrt{z}$, $\sqrt{1 - z}$,
valid off the cut); for instance $F(1/4, 3/4; 1; 4z/(1+z)^2) = \sqrt{1 +
z}\, F(1/2, 1/2; 1; z) = 2\sqrt{1 + z} K(z)/\pi$ once $K(z)$ is a
generator. The complete elliptic integrals are linked by Landen's
transformation in the same way (the zero test also finds it between two
generators already present, see
[](algorithms/exact/castle-special-functions#elliptic)): with $k' = \sqrt{1
- m}$ and $k_1 = (1 - k')/(1 + k')$,

$$
K(m) = (1 + k_1) K(k_1^2), \quad E(m) = (1 + k') E(k_1^2) - k' K(m)
\quad (m \notin [1, \infty)),
$$

and a new $K(m)$ or $E(m)$ is expressed through a generator $K$, $E$ at a
point of the Landen chain of $m$, down ($m \to k_1^2$) or up ($m \to
4s/(1+s)^2$, $s = \pm \sqrt{m}$, $|s| < 1$), within three steps, rather
than becoming a generator.

## Generators {#generators}

The generators of kind `GR_TOWER_HYPGEOM` have the argument $z$ and
the additional arguments $a_1, \ldots, a_p, b_1, \ldots, b_q$ (see
`gr_tower_adjoin_special_multi_flat`), and print as
`hypgeom_2f1(a, b, c, z)` (`hypgeom_0f1(b, z)`,
`hypgeom_1f1(a, b, z)`, `hypgeom_pfq([a...], [b...], z)`).

## Use in FLINT {#use}

- `gr_hypgeom_0f1`, `gr_hypgeom_1f1`, `gr_hypgeom_2f1` and
  `gr_hypgeom_pfq` of the castle fields implement the normal form; the
  complete elliptic integrals `gr_elliptic_k` and `gr_elliptic_e` use it
  as well.
- The table of known values and its tests are in
  `src/gr_ca/hypgeom.c` and `src/gr_ca/test/t-hypgeom.c`.
