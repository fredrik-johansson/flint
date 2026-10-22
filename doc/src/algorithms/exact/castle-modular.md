---
title: Modular forms and theta functions in castle arithmetic
summary: How castle arithmetic computes modular functions and forms, theta functions and Weierstrass functions exactly — reduction to the fundamental domain with exact multipliers, all values at a point through the modular lambda function and complete elliptic integrals, CM points and the Chowla–Selberg formula, theta functions of a variable at related points (torsion, multiples, sums) with rebasing, and points related by isogenies of degree 2^k 3^l.
complexity: The reduction to the fundamental domain costs a few exact comparisons; the values at a point cost the generators $\lambda(\tau_0)$, $K$, $E$ once, then rational operations. CM values cost the factorization of a polynomial of degree $6h$ from the Hilbert class polynomial of class number $h$. Related points cost division polynomials and addition formulas over the field of the theta constants.
order: 52
---

The modular functions and forms ($j$, $\lambda$, $\eta$, $\Delta$,
the Eisenstein series), the Jacobi theta functions and the Weierstrass
functions of the castle fields ([](gr_ca#modular)) are computed through
a small set of generators per point, so that the classical identities
between them hold by construction and the remaining ones are found by
the zero test.

## Reduction to the fundamental domain {#reduction}

The modular functions and forms of a point $\tau$ of the upper half-plane
(`src/gr_ca/modular.c`) are reduced exactly to the fundamental domain

$$
\mathcal{F} = \{ -1/2 < \operatorname{Re}(\tau) \le 1/2, |\tau| > 1 \}
\cup \{ |\tau| = 1, 0 \le \operatorname{Re}(\tau) \le 1/2 \},
$$

$\tau_0 = g \tau$ with $g \in PSL_2(\mathbb{Z})$ found numerically
(`acb_modular_fundamental_domain_approx`) and then verified and
corrected by exact comparisons, and the values at $\tau$ follow from the
transformation laws with exact multipliers (the eighth roots of unity of
`acb_modular_theta_transform`, the 24th roots of unity of
`acb_modular_epsilon_arg`, powers of $c\tau + d$, and $E_2(g\tau) =
(c\tau+d)^2 E_2(\tau) + 6c(c\tau+d)/(\pi i)$). At $\tau_0$, every value is
expressed through $\Lambda = \lambda(\tau_0)$ and the complete elliptic
integrals at $\Lambda$ (canonical arguments: for $\tau_0 \in \mathcal{F}$,
$|\Lambda - 1| \le 1$, with $\operatorname{Im}(\Lambda) > 0$ on the
boundary $\operatorname{Re}(\tau_0) = 1/2$):

$$
\theta_3^2 = 2K(\Lambda)/\pi, \quad
\theta_2 = \theta_3 \Lambda^{1/4}, \quad
\theta_4 = \theta_3 (1 - \Lambda)^{1/4}, \quad
\eta = (\theta_2 \theta_3 \theta_4 / 2)^{1/3},
$$

$$
j = 256 \frac{(1 - \Lambda + \Lambda^2)^3}{\Lambda^2 (1 - \Lambda)^2}, \quad
E_4 = \theta_3^8 (1 - \Lambda + \Lambda^2), \quad
E_6 = \theta_3^{12} \frac{(1 + \Lambda)(2 - \Lambda)(1 - 2\Lambda)}{2}, \quad
E_2 = \theta_3^4 \left(\frac{3 E(\Lambda)}{K(\Lambda)} - 2 + \Lambda\right),
$$

with principal roots (the right branches on $\mathcal{F}$), and $E_{2k}$,
$k \ge 4$, by the recurrence of the Weierstrass coefficients. The
relations between these values (Jacobi's identity, $j = E_4^3/\Delta$,
$E_4^3 - E_6^2 = 1728 \Delta$, $E_8 = E_4^2$, ...) thus hold by
construction, and the relations of $K$ and $E$ known to the zero test
(Legendre's relation, which is the transformation $\tau \to -1/\tau$,
and Landen's transformation, $\tau \to 2\tau$) apply to the modular
values as well. The generator $\lambda(\tau_0)$ has the kind
`GR_TOWER_MODULAR_LAMBDA` and prints as `modular_lambda(tau0)`.

At a CM point ($\operatorname{Re}(\tau_0)$ and $|\tau_0|^2$ rational,
$\tau_0$ a root of $A\tau^2 + B\tau + C$ of discriminant $D$, class
number $h \le 16$, $|D| \le 100000$), $\Lambda$ is the algebraic number
selected numerically among the roots of the factor of $\sum_k H_k N^k
M^{h-k}$ ($H$ the Hilbert class polynomial, $N = 256(1 - x + x^2)^3$, $M
= x^2 (1 - x)^2$) which vanishes at it, so that $j(\tau)$ and
$\lambda(\tau)$ are algebraic ($j(\sqrt{-5}) = 632000 + 282880 \sqrt 5$,
$j((1 + \sqrt{-163})/2) = -640320^3$). For the fundamental
discriminants of class number one with $|D|$ at most the option
`GR_TOWER_OPT_GAMMA_LATTICE_LIMIT` (36 by default), the formula of
Chowla and Selberg gives the periods,

$$
|\eta(\tau_0)|^2 = \frac{1}{\sqrt{2\pi |D|}} \prod_{j=1}^{|D|-1}
\Gamma(j/|D|)^{w \chi_D(j)/4},
$$

with $\eta(\tau_0)^2 = e^{\pi i \operatorname{Re}(\tau_0)/6} |\eta(\tau_0)|^2$
and $\theta_3^2 = 2^{2/3} \eta^2 (\Lambda(1 - \Lambda))^{-1/6}$ (times
the sixth root of unity selected numerically), so that $\eta(i) =
\Gamma(1/4)/(2\pi^{3/4})$, $E_2(i) = 3/\pi$, and the elliptic integrals
at the singular moduli of these discriminants are expressed through
gamma values (generalizing $K(1/2)$). Elsewhere $K(\Lambda)$ and
$E(\Lambda)$ are generators.

## Weierstrass functions {#weierstrass}

For the lattice $\mathbb{Z} + \tau \mathbb{Z}$
(the conventions of `acb_elliptic_p`), with $\theta_k =
\theta_k(0, \tau)$,

$$
e_1 = \tfrac{\pi^2}{3} (\theta_3^4 + \theta_4^4), \quad
e_2 = \tfrac{\pi^2}{3} (\theta_2^4 - \theta_4^4), \quad
e_3 = -\tfrac{\pi^2}{3} (\theta_2^4 + \theta_3^4), \quad
g_2 = \tfrac{4 \pi^4}{3} E_4, \quad g_3 = \tfrac{8 \pi^6}{27} E_6,
$$

$$
\wp(z) = e_1 + \left(\frac{\pi \theta_3 \theta_4 \theta_2(z)}{\theta_1(z)}\right)^2, \quad
\wp'(z) = -2\pi^3 (\theta_2 \theta_3 \theta_4)^2 \frac{\theta_2(z) \theta_3(z) \theta_4(z)}{\theta_1(z)^3}, \quad
\sigma(z) = \frac{\theta_1(z)}{\pi \theta_2 \theta_3 \theta_4} e^{\pi^2 E_2 z^2 / 6},
$$

so that the differential equation $\wp'^2 = 4\wp^3 - g_2 \wp - g_3 =
4(\wp - e_1)(\wp - e_2)(\wp - e_3)$, the periodicity and the values at
the half periods hold by construction, and the addition and
multiplication theorems through those of the theta functions (below).
The Weierstrass zeta function (which needs $\theta_1'(z)/\theta_1(z)$)
and the inverse of $\wp$ are not implemented.

## Theta functions of a variable {#theta}

The theta functions of a variable $z$ (`jacobi_theta(z, tau)`) are
reduced in both arguments: $\tau$ to $\tau_0$ as above ($z \to
-z/(c\tau + d)$ with the factor $\sqrt{i/(c\tau+d)}\, e^{-\pi i c
z^2/(c\tau + d)}$), and then $z$ modulo the half-lattice, $z = z_1 + (p +
q\tau_0)/2$ with $z_1 = x + y\tau_0$, $x, y \in [-1/4, 1/4)$ (exact
floors; $p$, $q$ adjusted on the boundary of that box, see below), with
the characteristics: in terms of $\theta[a, b](z) = \sum_n
\exp(\pi i (n + a/2)^2 \tau + 2\pi i (n + a/2)(z + b/2))$ ($\theta_1 =
-\theta[1,1]$, $\theta_2 = \theta[1,0]$, $\theta_3 = \theta[0,0]$,
$\theta_4 = \theta[0,1]$),

$$
\theta[a,b](z_1 + (p + q\tau)/2) = e^{-\pi i q^2 \tau/4 - \pi i q (z_1
+ (b+p)/2)} \theta[a + q, b + p](z_1),
$$

$\theta[a + 2, b] = \theta[a, b]$, $\theta[a, b + 2] = (-1)^a \theta[a, b]$,
and the parity $\theta[a, b](-z) = (-1)^{ab} \theta[a, b](z)$, to the
reduced point $z_0 = \pm z_1$ in the fundamental domain of $z \to \pm z +
(\mathbb{Z} + \tau_0\mathbb{Z})/2$ given by $0 < y < 1/4$, $-1/4 \le x <
1/4$, or $y \in \{0, 1/4\}$, $0 \le x \le 1/4$ (so that $z$ and $-z$, or a
point and its reduced point, have the same reduced point). At $\tau_0 = i$ and $\tau_0 = e^{\pi i/3}$, the
stabilizer of $\tau_0$ acts on $z$ (by $z \to iz$, respectively
$e^{\pm \pi i/3} z$), and the least reduced point of the orbit is taken.
At the reduced point $z_0 \ne 0$, $\theta_1(z_0, \tau_0)$ and
$\theta_4(z_0, \tau_0)$ are generators of the kind
`GR_TOWER_JACOBI_THETA` (two arguments, printed as
`jacobi_theta_1(z0, tau0)`, `jacobi_theta_4(z0, tau0)`), and
$\theta_2(z_0, \tau_0)$, $\theta_3(z_0, \tau_0)$ are algebraic
generators over them, with the same definition kind and Jacobi's
relations

$$
\theta_2(z)^2 \sqrt{1 - \Lambda} = \theta_4(z)^2 \sqrt{\Lambda} - \theta_1(z)^2, \quad
\theta_3(z)^2 \sqrt{1 - \Lambda} = \theta_4(z)^2 - \theta_1(z)^2 \sqrt{\Lambda}
$$

as moduli (adjoined once per point: later evaluations at an equal point
find them by their definition, whatever the representation of the
point). The quasi-periodicity, the transformation laws and the quartic
relations between the theta functions of $z$ thus hold by construction.

### Related points {#related}

Every reduced point $z_0$ met at $\tau_0$ is recorded
in the context with the method of its evaluation, chosen when it is
first met from the points recorded before (relations modulo the
half-lattice, found numerically and verified exactly), so that the
values at related points are not independent generators:

* a torsion point $z_0 = (a + b\tau_0)/N$ with $3 \le N \le 6$ (the
  reduced point: $z = 1/3$ has $z_0 = 1/6$, while $(1 + 2\tau)/5$ with
  $z_0 = -1/5 + \tau/10$ is not used): the values are algebraic over the theta constants,
  through the root $X = \tilde\wp(z_0)$ of the division polynomial
  $\psi_N$ (in the units $\tilde\wp = \wp/(\pi^2\theta_3^4)$), the
  ratios $\theta_k(z_0)/\theta_1(z_0)$ (square roots of $X - e_j$ up to
  constants) and a root of $\theta_1(z_0)^{N^2}$ given by the
  multiplication formula at the lattice point $N z_0$;
* $d z_0 = \pm n u$ with $u$ recorded, $n \le 6$, $d \le 3$: the
  multiplication formulas $\theta_k(nu) = \theta_3 (\theta_1(u)/\theta_3)^{n^2}
  S_k(n)$ with $S_k(n)$ rational in the ratios at $u$ (the recursion of
  the division polynomials), then for $d > 1$ the division as for torsion
  points (with $\tilde\wp(nu) = \tilde\wp(u) - \psi_{n-1}\psi_{n+1}/\psi_n^2$);
* $z_0 = \pm(u + v)$ with $u$, $v$ and $u - v$ recorded: rational, by
  the addition formulas $\theta_k(u+v)\theta_k(u-v)\theta_4^2 =
  \theta_k(u)^2\theta_4(v)^2 - \theta_{5-k}(u)^2\theta_1(v)^2$;
* $u$, $u + z_0$, $u - z_0$ recorded: $\theta_1(z_0)^2$ and
  $\theta_4(z_0)^2$ from two of these addition formulas (linear
  equations), the other values from Jacobi's relations;
* $z_0 = \pm(nu \pm v)$ otherwise ($n \le 3$): the ratios at $z_0$
  rational in the values at $u$ and $v$ (Jacobi's addition formulas for
  $\mathrm{sn}, \mathrm{cn}, \mathrm{dn}$ from the multiple $nu$), and
  $\theta_1(z_0)$ a generator;
* $2 z_0 = \pm(u \pm v)$ otherwise: with $P = z + w$, $Q = z - w$ and
  $S_{ab} = \theta_a(P)\theta_b(Q) + \theta_b(P)\theta_a(Q)$, Jacobi's
  formulas for $\theta_a(z+w)\theta_b(z-w)$ give $\theta_1\theta_4(z) /
  \theta_2\theta_3(z) = S_{14}/S_{23}$ and its images under the
  permutations of 2, 3, 4, so that the ratios $r_k = \theta_k/\theta_1$
  at $z = (P + Q)/2$ are monomials in the square roots of the $S_{ab}$
  ($r_2 = \sqrt{S_{23} S_{24} / (S_{14} S_{13})}$, and so on, with
  $r_2 r_3 r_4$ rational; the signs of the square roots chosen
  numerically, unique by the enclosures, and checked for consistency);
  $\theta_1(z_0)$ is a generator. The point $w = (P - Q)/2$ then has
  $r_k(w)/r_k(z)$ rational in the $S_{ab}$ and the differences $D_{ab}$
  ($r_4(w)/r_4(z) = \theta_3 S_{12} / (\theta_2 D_{13})$, say), and
  $\theta_1(z)\theta_1(w)$ a monomial in the same square roots and
  $\sqrt{D_{12} D_{13}}$;
* otherwise the generators above.

### Rebasing {#rebasing}

Where the methods above would make the values at a new
point algebraic over those at points recorded with generators of their
own (roots of polynomials, square roots), the roles are exchanged: the
new point (or a related one) gets the generators, the earlier points
become rational over it, and the context records the value of each
earlier generator as a rational function of the new ones. When
elements involving such a generator and elements involving one of the
new ones meet in an operation, the generator becomes algebraic in
their tower with the linear modulus $X - \text{value}$ (the generators
of the values are moved before it: one rebuild of the tower for all the
records that apply), so that the elements built on it, before or after,
meet those built on the new generators. (Not when a tower merely
contains both: the towers of a long session hold generators of
unrelated computations.) The cases:

* $d z_0 = \pm n u$ with $d > 1$ and $u$ of generators: with $an + cd
  = 1$, the point $b = a z_0 \pm c u$ (so that $n b = z_0$, $d b = \pm
  u$) gets the generators, $u$ and $z_0$ are its multiples (the division
  by 3 of $\theta(3z)$ would be a root of degree 9 with large
  coefficients: $\theta(3z)$, $\theta(2z)$, $\theta(z)$ cost the same
  in any order);
* $u$, $P = u + x$, $Q = u - x$ with $P$, $Q$ of generators ($z = (P + Q)/2$
  with $P = z + w$, $Q = z - w$ first, and $x = w$; or $x = z_0$ in the
  case of the linear equations above): $x$ gets the generators, $P$
  the ratios of Jacobi's formulas from $u$, $x$ (its $\theta_1$ stays a
  generator) and $Q$ the addition formulas;
* $z_0 = \pm(2u \pm v)$ with $v$ of generators: as the previous case with
  $x = -(u \pm v)$, $P = v$, $Q = z_0$.

The values at the points are then those of the order in which the
points with generators came first, whatever the actual order (the
addition formulas with the values at $z$, $w$, $z \pm w$ take 0.06 to
0.25 seconds in all 24 orders, against up to seven seconds without the
rebasing). Elements built on the earlier
generators other than the values themselves ($\wp(z + w)$ as a rational
function of $\theta_1(z + w)$ and $\theta_4(z + w)$, say) become
rational functions of the new generators, which may be large.

The values are cached with the record (those of the 32 most recently
used points; the ratios $\theta_k/\theta_3$ at the 8 most recent points
$\tau_0$ and the theta constants themselves at the last one, together
with the anchor chosen for each $\tau_0$, so that values computed later
come from the same anchor), so that repeated evaluations at a point cost
a lookup; the bound keeps the cache from holding the towers of a long
session alive. At the CM points where $\lambda$ has degree more than 2
(from a class polynomial: $\lambda(1/2 + i)$ is a root of a quartic
with coefficients of five digits), the torsion points and the divisions
are not used (the roots of the division polynomials over that field are
costly: tens of seconds for $\psi_3$), nor the other relations except
the doubling $z_0 = \pm 2u$ (the inverses of the denominators of their
formulas in that field take minutes, hours for $n = 6$): the values
there are generators, and identities between the values at related
points are not decided.
Thus the duplication formula $\theta_1(2z)\theta_2\theta_3\theta_4 =
2\theta_1(z)\theta_2(z)\theta_3(z)\theta_4(z)$, the addition formulas and
the addition theorem of $\wp$ are decided whatever the order of
evaluation of $\theta(z)$, $\theta(w)$, $\theta(z \pm w)$, as rational
identities in the values at the points with generators (the forms
$S_{ab}$, $D_{ab}$ above remain for the cases the rebasing leaves out).
Relations among more than three points ($z_0 = u + v + w$, say) are not
used.

## Commensurable points {#commensurable}

A point $\tau_0 \in \mathcal{F}$ which is not a
CM point is compared with the points $\tau_1$ whose $\lambda(\tau_1)$ is
already a generator of the context (the `anchors`): when $\tau_0 =
\gamma \tau_1$ with $\gamma \in GL_2^+(\mathbb{Q})$ of determinant
$2^k 3^l$ (an integer relation between $1, \tau_1, \tau_0, \tau_0 \tau_1$
found by LLL and verified exactly; $k \le 4$ when $l = 0$, otherwise
$l \le 2$ and $k \le 1$), the values at $\tau_0$ are computed from those
at $\tau_1$ rather than from a new generator. Two points commensurable
beyond these limits have generators of their own, and the points
commensurable with both (within the limits) could come from either: of
the anchors found, those whose values come from the most recently used
generator $\lambda$ are preferred (the `root`: the values of an
evaluation in progress, $j(\tau)$, $\eta(\tau)$ and $\eta(3\tau)$ say,
then come from one root, rather than $\eta(3\tau)$ from the root of an
earlier evaluation at $(3\tau + 1)/4$, with a cheaper chain to $3\tau$,
which would leave the identity between them undecided), and of these the
one with the shortest chain of steps. A chain with a tripling and a
further step (a halving, another tripling) gives values in towers of
large degree, where everything at $\tau_0$ is costly (the theta
functions of $z$ at $4i/\pi$ through two thirdings from $1/3 + \pi i/4$:
3.5 seconds, against 0.01 at a generator of its own; $j(\tau)$,
$\eta(3\tau)$ after $\lambda((2\tau + 1)/3)$: minutes). There
$\lambda(\tau_0)$ is a new generator `linked` to the anchor: the values
at $\tau_0$ are those of a new generator (the fourth roots of
$\lambda$ and $1 - \lambda$, $K(\lambda)$, $\theta_3 = \sqrt{2K/\pi}$,
$E(\lambda)$), and the context records their values through the chain,
computed when first needed, as rebase records (see the theta functions
of $z$ above): where values from the two meet, those at $\tau_0$ are the
values of the chain, and the identities between them hold as with the
chain (Jacobi's modular equation of degree 3 between $\lambda(\tau)$
and $\lambda(3\tau)$ with $3\tau$ from a linked $(3\tau + 1)/2$);
elsewhere they cost what they cost in a new field. With $\gamma = g M$, $g
\in SL_2(\mathbb{Z})$, $M = \begin{pmatrix} a & b \\ 0 & d
\end{pmatrix}$, the steps are the duplication formulas

$$
\theta_3(2w)^2 = \frac{\theta_3^2 + \theta_4^2}{2}, \quad
\theta_4(2w)^2 = \theta_3 \theta_4, \quad
\theta_2(2w)^2 = \frac{\theta_3^2 - \theta_4^2}{2}, \quad
2 E_2(2w) - E_2(w) = \theta_3(2w)^4 + \theta_2(2w)^4
$$

(and their inverses for $w \to w/2$; the square roots are chosen
numerically), Jacobi's modular equation of degree 3 for $w \to 3w$ and
$w \to w/3$: with $u = (\theta_2(w)/\theta_3(w))^{1/2}$ and $v$ the
same at $3w$,

$$
u^4 - v^4 - 2uv(1 - u^2 v^2) = 0, \quad
\frac{\theta_3(w)^2}{\theta_3(3w)^2} = 1 + \frac{2v^3}{u}, \quad
3 E_2(3w) - E_2(w) = \sum_{k=2}^4 \theta_k(w)^2 \theta_k(3w)^2
$$

(the root of the quartic equals $\pm v$ numerically and is adjoined as an
algebraic generator), the translations and the transformation laws of
$g$. The modular equations of these levels thus hold by construction:
Landen's $\lambda(2\tau) = ((1 - k')/(1 + k'))^2$, the modular polynomial
$\Phi_2(j(\tau), j(2\tau)) = 0$, eta quotients such as $\theta_4(2\tau) =
\eta(\tau)^2/\eta(2\tau)$, the Hauptmodul $t = (\eta(\tau)/\eta(3\tau))^{12}$
of $\Gamma_0(3)$ with $j(\tau) = (t + 27)(t + 243)^3/t^3$ and $j(3\tau) =
(t + 27)(t + 3)^3/t$, and so on. When only $\lambda$ is needed, the
computation uses the ratios $\theta_k/\theta_3$ and no elliptic integral.
Matrices $M$ with $3 \mid a$ and $3 \mid d$ ($\tau_1 + 1/3$, say), other
levels (primes $\ge 5$), and levels with a factor 3 inside a conjugation
are not used: values at points related by such transformations are
independent generators, and the corresponding identities are not
decided (`T_UNKNOWN`).

The factors $\exp(\pi i \alpha)$ of the transformation laws are written
$\exp(\pi i \operatorname{Re}(\alpha)) \exp(-\pi \operatorname{Im}(\alpha))$
when $\operatorname{Re}(\alpha)$ is rational (a root of unity times a real
exponential), which keeps the multiplicative relations between them of
small degree.

Conjugation maps $F(u)$ to $F(\bar u)$ ($W_k(u)$ to $W_{-k}(\bar u)$) off
the branch cuts, and the values at real points are real generators,
which the real fields accept.

## Limitations {#limitations}

Isogenies beyond the levels $2^k 3^l$ above (primes $\ge 5$, and the
matrices excluded above), relations between theta functions at more
than three related points, and the Weierstrass zeta function and the
inverse of $\wp$ are not available: values at points related only in
these ways are independent generators, and identities between them are
not decided (`T_UNKNOWN`, or no termination in the unlimited mode).

## Use in FLINT {#use}

- `gr_modular_j`, `gr_modular_lambda`, `gr_modular_delta`,
  `gr_dedekind_eta`, `gr_eisenstein_e`, `gr_jacobi_theta` (and its
  variants), `gr_weierstrass_p`, `gr_weierstrass_p_prime`,
  `gr_weierstrass_sigma` and `gr_elliptic_invariants` of the castle
  fields implement these methods (`src/gr_ca/modular.c`).
- The generators have the kinds `GR_TOWER_MODULAR_LAMBDA` and
  `GR_TOWER_JACOBI_THETA`; the option `GR_TOWER_OPT_GAMMA_LATTICE_LIMIT`
  bounds the discriminants for which the Chowla–Selberg formula is used.
