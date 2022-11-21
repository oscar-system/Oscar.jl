```@meta
CurrentModule = Oscar
DocTestSetup = quote
  using Oscar
end
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["groebner_bases.md"]
```

# Gröbner/Standard Bases Over Fields

Standard bases are sets of generators for ideals of multivariate polynomial rings (submodules of free modules
over such rings) which behave well under division with remainder. Each such basis is defined with respect to
a monomial ordering. If this ordering is global (that is, it is a well-ordering), we use the name Gröbner basis
instead of standard basis. In this section, we fix our notation in this context and introduce corresponding
OSCAR functionality.

Let $K[x] = K[x_1, \dots, x_n]$ be a polynomial ring over a field $K$, and let $>$ be a monomial ordering on $\text{Mon}_n(x)$.

Given a polynomial $f\in K[x]\setminus \{0\}$, write $f$ as the sum of its nonzero terms as follows:

$f = a_\alpha x^\alpha + a_{\beta_1} x^{\beta_1} + \dots + a_{\beta_s} x^{\beta_s},\quad x^\alpha > x^{\beta_1} > \dots > x^{\beta_s} .$

Then, with respect to $>$, we refer to $\text{LT}_>(f) = a_\alpha x^\alpha$, $\text{LM}_>(f) = x^\alpha$, $\text{LE}_>(f) = \alpha$, $\text{LC}_>(f) = a_\alpha$,
and $\text{tail}_>(f) = f - \text{LT}_>(f)$ as the  *leading term*, the *leading monomial*, the *leading exponent*, the *leading
coefficient*, and the *tail* of $f$, respectively.

Next note that the set

$U_>:= \{u\in K[x]\setminus \{0\} \mid {\text{LM}}_>(u)=1 \}$

is a multiplicatively closed subset of $K[x]$. Consider the localization

$K[x]_>:= K[x][U^{-1}] = \left\{ \frac{f}{u} \:\bigg|\: f \in K[x], \, u\in U_>\right\}.$

Then $K[x]\subseteq K[x]_>\subseteq K[x]_{\langle x \rangle},$ where $K[x]_{\langle x \rangle}$
is the localization of $K[x] $ at the maximal ideal $\langle x \rangle .$ Moreover,
 
- ``K[x] = K[x]_>`` iff $>$ is global, and

- ``K[x]_> = K[x]_{\langle x \rangle}`` iff $>$ is local.

Extending the notation introduced for polynomials, let now $f\in K[x]_>\setminus \{0\}$. Choose $u\in U_>$
such that $uf\in K[x]$. Then, with respect to $>$, the  *leading term*  of $f$
is defined to be $\text{LT}_>(f) = \text{LT}_>(uf)$ (this definition is independent of the choice of $u$).
The *leading monomial* $\text{LM}_>(f)$, the *leading exponent* $\text{LE}_>(f)$, the
*leading coefficient* $\text{LC}_>(f)$, and the *tail* $\text{tail}_>(f)$ of $f$ are defined similarly.

!!! note
    Given a monomial ordering $>$ on a free $K[x]$-module $F = K[x]^p$ with basis $e_1, \dots, e_p$,
    the above notation extends naturally to elements of  $K[x]^p$ and $K[x]_>^p$, respectively. There is one particularity:
	Given an element $f = K[x]^p\setminus \{0\}$ with leading term $\text{LT}(f) = x^\alpha e_i$, we write $\text{LE}_>(f) = (\alpha, i)$.

!!! note
    The OSCAR functions discussed in this section depend on a monomial ordering which may be entered as a keyword argument.
    For a polynomial ring $R$, the respective default ordering is `degrevlex` except if $R$ is $\mathbb Z$-graded with
	positive weights. Then the corresponding `wdegrevlex` ordering is used. For a free $R$-module $F$, the default
	ordering is `default_ordering(R)*lex(gens(F))`.

Here are some illustrating OSCAR examples:

##### Examples

```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> default_ordering(R)
degrevlex([x, y, z])

julia> F = free_module(R, 2)
Free module of rank 2 over Multivariate Polynomial Ring in x, y, z over Rational Field

julia> default_ordering(F)
degrevlex([x, y, z])*lex([gen(1), gen(2)])

julia> S, _ = grade(R, [1, 2, 3])
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by 
  x -> [1]
  y -> [2]
  z -> [3], MPolyElem_dec{fmpq, fmpq_mpoly}[x, y, z])

julia> default_ordering(S)
wdegrevlex([x, y, z], [1, 2, 3])

```

```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> f = 3*z^3+2*x*y+1
2*x*y + 3*z^3 + 1

julia> terms(f)
terms iterator of 3*z^3 + 2*x*y + 1

julia> collect(terms(f))
3-element Vector{fmpq_mpoly}:
 3*z^3
 2*x*y
 1

julia> monomials(f, ordering = lex(R))
monomials iterator of 2*x*y + 3*z^3 + 1

julia> exponents(f, ordering = neglex(R))
exponents iterator of 1 + 3*z^3 + 2*x*y

julia> coefficients(f)
coefficients iterator of 3*z^3 + 2*x*y + 1

julia> leading_term(f)
3*z^3

julia> leading_monomial(f, ordering = lex(R))
x*y

julia> leading_exponent(f, ordering = neglex(R))
3-element Vector{Int64}:
 0
 0
 0

julia> leading_coefficient(f)
3

julia> tail(f)
2*x*y + 1
```

## Division With Remainder

Buchberger's algorithm for computing Gröbner bases makes use of reduction steps which rely on multivariate division with remainder.
Here, Euclidean division with remainder is extended to the multivariate case, allowing more than one divisor, and using a fixed
global monomial ordering $>$ to distinguish leading terms of the polynomials (elements of free modules) involved. At each step
of the division process, the idea is to remove the leading term of the intermediate dividend using the leading term of some divisor
by which it is divisible. Termination is guaranteed since $>$ is a well-ordering. In the more general case of standard bases we also
allow local and mixed orderings. Here, a variant of the division algorithm due to Mora is used to guarantee termination.

We introduce the relevant notation.

First suppose that $>$ is global and let polynomials $g\in K[x]$ and $f_1, \dots, f_r\in K[x]\setminus \{0\}$ be given.
In this situation, multivariate division with remainder allows us to compute expressions

$g = q_1f_1+\dots q_rf_r + h, \; h\in K[x], \;\text{ all }\; q_i \in K[x]$

such that:

- ``\text{LM}_>(g) \ge \text{LM}_>(q_if_i)`` whenever both sides are nonzero.
- If $h$ is nonzero, then $\text{LM}_>(h)$ is not divisible by any $\text{LM}_>(f_i)$.

Each such expression is called a *standard representation* for $g$ with *quotients* $q_i$ and *remainder* $h$
(in terms of the $f_i$, with respect to $>$). If, at each stage of the division process, we allow to remove some
term of the current dividend instead of just focusing on the leading term, then the algorithm will return a standard
expression in which the remainder is *fully reduced*. That is, $h$ satisfies the stronger condition below:

- If $h$ is nonzero, then no term of $h$ is divisible by any $\text{LM}_>(f_i)$.

Without restrictions on  $>$, let elements $g\in K[x]_>$ and $f_1, \dots, f_r\in K[x]\setminus \{0\}$ be given.
In this more general situation, Mora division with remainder allows us to compute expressions

$ug = q_1f_1+\dots q_rf_r + h, \; h\in K[x]_>, \;\text{ all }\; q_i \in K[x]_>$

such that:

- ``u`` is a unit of $K[x]_>$, that is, $\text{LM}_>(u)=1$.
- ``\text{LM}_>(g) \ge \text{LM}_>(q_if_i)`` whenever both sides are nonzero.
- If $h$ is nonzero, then $\text{LM}_>(h)$ is not divisible by any $\text{LM}_>(f_i)$.

Each such expression is called a *weak standard representation* for $g$ with *quotients* $q_i$ and *remainder* $h$ (in terms of the $f_i$, with respect to $>$).
If $g\in K[x]$, we speak of a *polynomial weak standard representation* if $u$ and the $q_i$ are elements of $K[x].$ Using power series expansions, it makes
still sense to speak of fully reduced remainders. However, even if we start from polynomial data, such remainders may not be computable (in finitely many steps).

!!! note
    Given a monomial ordering $>$ on a free $K[x]$-module $F = K[x]^p$ with basis $e_1, \dots, e_p$,
    the above notation and the division algorithms extend naturally to $K[x]^p$ and $K[x]_>^p$, respectively.
	
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

Still keeping the notation introduced at the beginning of this section, let $G$ be a subset of $K[x]_>$.
Then the *leading ideal* of $G$ is the ideal of $K[x]$ defined by

$\text{L}_>(G)=\langle \text{LT}_>(g) \mid g\in G\setminus\{0\}\rangle\subset K[x].$

A *standard basis* of an ideal $I\subset K[x]_>$ is a finite subset $G$ of $I$ such that $\text{L}_>(G) = \text{L}_>(I)$.
A finite subset of $K[x]_>$ is a *standard basis* if it is a standard basis for the ideal it generates.

We call a standard basis $G = \{g_1,\dots, g_r\}\subset K[x]_>\setminus \{0\}$ *minimal*  if  $\text{LM}_>(g_i)\neq \text{LM}_>(g_j)$
for $i\neq j$.

!!! note
    The definition above deviates from the definition in most textbooks as we do not ask that the leading coefficients of the standard basis elements are 1.

!!! note
    Every nonzero ideal of $K[x]_>$ has a (minimal) standard basis with respect to $>$. Such bases can be computed using Buchberger's algorithm.

!!! note
    The standard bases returned by OSCAR are always minimal in the sense above.

If $>$ is global, that is, $K[x] = K[x]_>$, then we use the word *Gröbner basis* instead of standard basis. We call a
Gröbner basis $G = \{g_1,\dots, g_r\}$ *reduced* if  it is minimal and no term of $g_i$ is divisible by $\text{LM}_>(g_j)$, for $i\neq j$. Using power
series expansions, we may also define reduced standard bases. However, while reduced Gröbner bases can always be computed, reduced
standard bases may not be computable (in finitely many steps).

!!! note
    Given a monomial ordering $>$ on a free $K[x]$-module $F = K[x]^p$ with basis $e_1, \dots, e_p$,
    the above notation and Buchberger's algorithm extend naturally to submodules of $K[x]_>^p$.
   
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

