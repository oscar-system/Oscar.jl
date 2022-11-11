```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["groebner_bases.md"]
```

# Gröbner and Standard Bases Over Fields

Let $K[x] = K[x_1, \dots, x_n]$ be a polynomial ring over a field $K$,
and let $>$ be a monomial ordering on $\text{Mon}_n(x)$.

Given a polynomial $f\in K[x]$, write $f$ as the sum of its nonzero terms as follows:

$f = a_\alpha x^\alpha + a_{\beta_1} x^{\beta_1} + \dots + a_{\beta_s} x^{\beta_s},\quad x^\alpha > x^{\beta_1} > \dots > x^{\beta_s} .$

Then, with respect to $>$, we refer to $\text{LT}_>(f) := a_\alpha x^\alpha$, $\text{LM}_>(f) := x^\alpha$, $\text{LE}_>(f) := \alpha$, $\text{LC}_>(f) := a_\alpha$,
and $\text{tail}_>(f) := f - \text{LT}_>(f)$ as the  *leading term*, the *leading monomial*, the *leading exponent*, the *leading
coefficient*, and the *tail* of $f$, respectively.

Next note that the set

$U_>:= \{u\in K[x]\setminus \{0\} \mid {\text{LM}}_>(u)=1 \}$

is a multiplicatively closed subset of $K[x]$. Consider the localization

$K[x]_>:= K[x][U^{-1}] = \left\{ \frac{f}{u} \:\bigg|\: f \in K[x], \, u\in U_>\right\}.$

Then $K[x]\subseteq K[x]_>\subseteq K[x]_{\langle x \rangle}.$ Moreover,

- ``K[x] = K[x]_>`` iff $>$ is global, and

- ``K[x]_> = K[x]_{\langle x \rangle}`` iff $>$ is local.

Extending the notation introduced for polynomials, let now $f\in K[x]_>$. Choose $u\in U_>$ such that
$uf\in K[x]$. Then, with respect to $>$, the  *leading term*  of $f$
is defined to be $\text{LT}_>(f) := \text{LT}_>(uf)$ (this definition is independent of the choice of $u$).
The *leading monomial* $\text{LM}_>(f)$, the *leading exponent* $\text{LE}_>(f)$, the
*leading coefficient* $\text{LK}_>(f)$, and the *tail* $\text{tail}_>(f)$ of $f$ are defined similarly.

!!! note
    We use the analogous notation for elements of free modules with a given monomial ordering.

!!! note
    The OSCAR functions discussed in this section depend on a monomial ordering which may be entered as a keyword argument.
    For a polynomial ring $R$, the respective default ordering is `degrevlex` except if $R$ is $\mathbb Z$-graded with
	positive weights. Then the corresponding `wdegrevlex` ordering is used. For a free $R$-module $F$, the default
	ordering is `default_ordering (R)*lex(gens(F))`.

Here are some illustrating OSCAR examples:

##### Examples

```@repl oscar
R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
default_ordering(R)
F = free_module(R, 2)
default_ordering(F)
S, _ = grade(R, [1, 2, 3])
default_ordering(S)
```

```@repl oscar
R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
f = 3*z^3+2*x*y+1
```

## Division With Remainder

Buchberger's algorithm for computing Gröbner bases makes use of reduction steps which rely on multivariate division with remainder.
For the computation of standard bases, a variant of of the division algorithm due to Mora comes into play. In that context, with notation
as above, we use the following terminology.

First suppose that $>$ is global and let polynomials $g\in K[x]$ and $f_1, \dots, f_r\in K[x]\setminus \{0\}$ be given.
In this situation, multivariate division with remainder allows us to compute expressions

$g = q_1f_1+\dots q_rf_r + h, \; h\in K[x], \;\text{ all }\; q_i \in K[x]$

such that:

- If both sides are nonzero, then $\text{LM}_>(g) \ge \text{LM}_>(q_if_i)$.
- If $h$ is nonzero, then $\text{LM}_>(h)$ is not divisible by any $\text{LM}_>(f_i)$.

Each such expression is called a *standard representation* for $g$ with *quotients* $q_i$ and *remainder* $h$
(in terms of the $f_i$, with respect to $>$). We say that the remainder $h$ in a standard representation is
*fully reduced*  if it satisfies the stronger condition below:

- If $h$ is nonzero, then no term of $h$ is divisible by any $\text{LM}_>(f_i)$.

Without restrictions on  $>$, let elements $g\in K[x]_>$ and $f_1, \dots, f_r\in K[x]\setminus \{0\}$ be given.
In this more general situation, Mora division with remainder allows us to compute expressions

$ug = q_1f_1+\dots q_rf_r + h, \; h\in K[x]_>, \;\text{ all }\; q_i \in K[x]_>$

such that:

- ``u`` is a unit of $ K[x]_>$.
- If both sides are nonzero, then $\text{LM}_>(g) \ge \text{LM}_>(q_if_i)$.
- If $h$ is nonzero, then $\text{LM}_>(h)$ is not divisible by any $\text{LM}_>(f_i)$.

Each such expression is called a *weak standard representation* for $g$ with *quotients* $q_i$ and *remainder* $h$ (in terms of the $f_i$, with respect to $>$).
If $g\in K[x])$, we speak of a *polynomial weak standard representation* if $u$ and the $q_i$ are elements of $K[x].$ While from a theoretical point of view it makes
still sense to speak of fully reduced remainders, such remainders may not be computable (in finitely many steps).

The OSCAR functions discussed below compute standard representations and polynomial weak standard representations, respectively.

```@julia
reduce(I::IdealGens, J::IdealGens; ordering::MonomialOrdering = default_ordering(base_ring(J)))
```

```@julia
reduce(f::T, F::Vector{T}; ordering::MonomialOrdering = default_ordering(parent(f))) where {T <: MPolyElem}
```

```@julia
reduce(F::Vector{T}, G::Vector{T}; ordering::MonomialOrdering = default_ordering(parent(F[1]))) where {T <: MPolyElem}
```

```@julia
reduce_with_quotients(I::IdealGens, J::IdealGens; ordering::MonomialOrdering = default_ordering(base_ring(J)))
```

```@julia
reduce_with_quotients_and_units(I::IdealGens, J::IdealGens; ordering::MonomialOrdering = default_ordering(base_ring(J)))
```


## Gröbner and Standard Bases

Still keeping the notation introduced at the beginning of this section, let $I$ be an ideal of $K[x]_>$.
Then the *leading ideal* of $I$ is the ideal of $K[x]$ defined by

$\text{L}_>(I):=\langle \text{LT}_>(f) \mid f\in I\rangle\subset K[x].$

A *standard basis* of $I$ is a finite set of elements of $I$ such that the leading terms of these elements
generate $\text{L}_>(I)$. A finite subset of $K[x]_>$ is a *standard basis* if it is a standard basis for the
ideal it generates.

We call a standard basis $\{f_1,\dots, f_r\}\subset K[x]_>$ *minimal*  if  $\text{LM}_>(f_i)\neq \text{LM}_>(f_j)$,
for $i\neq j$. We call it *reduced* if  it is minimal and no term of $f_i$ is divisible by $\text{LM}_>(f_j)$, for $i\neq j$.

!!! note
    The above definitions deviate from the usual textbook definitions as we do not ask that the leading coefficients of the standard basis elements are 1.

!!! note
    The standard bases returned by OSCAR are always minimal in the sense above.

If $>$ is global, that is, $K[x] = K[x]_>$, then we use the expression *Gröbner basis* instead of standard basis.
Note that reduced Gröbner bases can always be computed, but reduced standard bases may not be computatable (in finitely many steps).

!!! note
    We use the analogous notation for submodules of free modules with a given monomial ordering.

```@docs
groebner_basis(I::MPolyIdeal;
	ordering::MonomialOrdering = default_ordering(base_ring(I)),
	complete_reduction::Bool = false)
```
```@docs
standard_basis(I::MPolyIdeal;
	ordering::MonomialOrdering = default_ordering(base_ring(I)),
	complete_reduction::Bool = false)
```
	
### Gröbner Bases with transformation matrix

```@docs
groebner_basis_with_transformation_matrix(I::MPolyIdeal;
	ordering::MonomialOrdering = default_ordering(base_ring(I)),
	complete_reduction::Bool=false)
```

    fglm

    Gröbner walks

    Hilbert-driven

!!! warning "Expert functions for Gröbner bases"
    The following functions are low-level implementations of various Gröbner
    basis algorithms with many adjustable arguments. Only use these
    functions directly if you know what you are doing.

```@docs
f4( I::MPolyIdeal; initial_hts::Int=17, nr_thrds::Int=1, max_nr_pairs::Int=0, la_option::Int=2, reduce_gb::Int=1, info_level::Int=0)
```

## Leading Ideals

```@docs
leading_ideal(g::Vector{T}; ordering::MonomialOrdering) where { T <: MPolyElem }
leading_ideal(I::MPolyIdeal; ordering::MonomialOrdering)
```

## Normal Forms

```@docs
normal_form(f::T, J::MPolyIdeal) where { T <: MPolyElem }
normal_form(A::Vector{T}, J::MPolyIdeal) where { T <: MPolyElem }
```

## Syzygies

```@docs
syzygy_generators(a::Vector{<:MPolyElem})
```

