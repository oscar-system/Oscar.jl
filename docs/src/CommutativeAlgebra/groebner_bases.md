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

We fix our notation in the context of standard (Gröbner) bases and present relevant OSCAR functions.

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
    The OSCAR functions discussed in this section depend on a monomial `ordering` which is entered as a keyword argument.
    Given a polynomial ring $R$, the `default_ordering` for this is `degrevlex` except if $R$ is $\mathbb Z$-graded with
	positive weights. Then the corresponding `wdegrevlex` ordering is used. Given a free $R$-module $F$, the
	`default_ordering` is `default_ordering(R)*lex(gens(F))`.

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

julia> collect(ans)
3-element Vector{fmpq_mpoly}:
 3*z^3
 2*x*y
 1

julia> monomials(f, ordering = lex(R))
monomials iterator of 2*x*y + 3*z^3 + 1

julia> coefficients(f)
coefficients iterator of 3*z^3 + 2*x*y + 1

julia> exponents(f, ordering = neglex(R))
exponents iterator of 1 + 3*z^3 + 2*x*y

julia> coefficients_and_exponents(f)
coefficients and exponents iterator of 3*z^3 + 2*x*y + 1

julia> collect(ans)
3-element Vector{Tuple{fmpq, Vector{Int64}}}:
 (3, [0, 0, 3])
 (2, [1, 1, 0])
 (1, [0, 0, 0])

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

The computation of Gröbner (standard) bases relies on multivariate division with remainder which is interesting
in its own right. If a monomial ordering $>$ is given, the basic idea is to mimic Euclidean division with remainder,
allowing more than one divisor: At each step of the resulting process, this amounts to removing the leading term of the
intermediate dividend, using the leading term of *some* divisor by which it is divisible. In its basic form, the process
works well if $>$ is global, but may not terminate, however, for local and mixed orderings. In the latter case, Mora's division
algorithm, which relies on a restricted selection strategy for the divisors to be used, comes to our rescue.

We dicuss this in more detail:

First suppose that $>$ is global and let polynomials $g\in K[x]$ and $f_1, \dots, f_r\in K[x]\setminus \{0\}$ be given.
In this situation, multivariate division with remainder allows us to compute expressions

$g = q_1f_1+\dots q_rf_r + h, \; h\in K[x], \;\text{ all }\; q_i \in K[x]$

such that:

- ``\text{LM}_>(g) \ge \text{LM}_>(q_if_i)`` whenever both sides are nonzero.
- If $h$ is nonzero, then $\text{LM}_>(h)$ is not divisible by any $\text{LM}_>(f_i)$.

Each such expression is called a *standard representation* for $g$ with *quotients* $q_i$ and *remainder* $h$
(on division by the $f_i$, with respect to $>$). If, at each step of the division process, we allow to remove some
term of the current dividend instead of just focusing on its leading term, then the algorithm will return a standard
expression in which the remainder is *fully reduced*. That is, $h$ satisfies the stronger condition below:

- If $h$ is nonzero, then no term of $h$ is divisible by any $\text{LM}_>(f_i)$.

Without restrictions on  $>$, let elements $g\in K[x]_>$ and $f_1, \dots, f_r\in K[x]\setminus \{0\}$ be given.
In this situation, Mora division with remainder allows us to compute expressions

$ug = q_1f_1+\dots q_rf_r + h, \; h\in K[x]_>, \;\text{ all }\; q_i \in K[x]_>$

such that:

- ``u`` is a unit of $K[x]_>$, that is, $\text{LM}_>(u)=1$.
- ``\text{LM}_>(g) \ge \text{LM}_>(q_if_i)`` whenever both sides are nonzero.
- If $h$ is nonzero, then $\text{LM}_>(h)$ is not divisible by any $\text{LM}_>(f_i)$.

Each such expression is called a *weak standard representation* for $g$ with *quotients* $q_i$ and *remainder* $h$ (on division by the $f_i$, with respect to $>$).
If $g\in K[x]$, we speak of a *polynomial weak standard representation* if $u$ and the $q_i$ are elements of $K[x].$ Using power series expansions, it makes
still sense to speak of fully reduced remainders. However, even if we start from polynomial data, such remainders may not be computable (in finitely many steps).

!!! note
    Given a monomial ordering $>$ on a free $K[x]$-module $F = K[x]^p$ with basis $e_1, \dots, e_p$,
    the above notation and the division algorithms extend naturally to $K[x]^p$ and $K[x]_>^p$, respectively.
	
The OSCAR functions discussed below compute standard representations and polynomial weak standard representations, respectively.
In the global case, they always return fully reduced remainders.

```@docs
reduce(g::T, F::Vector{T}; 
	ordering::MonomialOrdering = default_ordering(parent(F[1]))) where {T <: MPolyElem}
```

```@docs
reduce_with_quotients(g::T, F::Vector{T}; 
	ordering::MonomialOrdering = default_ordering(parent(F[1]))) where {T <: MPolyElem}
```
		  
```@docs
reduce_with_quotients_and_unit(g::T, F::Vector{T}; 
	ordering::MonomialOrdering = default_ordering(parent(F[1]))) where {T <: MPolyElem}
```

## Gröbner and Standard Bases

Still keeping the notation introduced at the beginning of this section, let $G$ be a subset of $K[x]_>$.
Then the *leading ideal* of $G$ is the ideal of $K[x]$ defined by

$\text{L}_>(G)=\langle \text{LT}_>(g) \mid g\in G\setminus\{0\}\rangle\subset K[x].$

A finite subset $G$ of an ideal $I\subset K[x]_>$ is called a *standard basis* of $I$ (with respect to $>$)
if $\text{L}_>(G) = \text{L}_>(I)$.  A finite subset of $K[x]_>$ is a *standard basis*
if it is a standard basis of the ideal it generates. A standard basis with respect to a global monomial
ordering is also called a *Gröbner basis*.

!!! note
    Every standard basis of $I$ generates $I$.

!!! note
    Gröbner bases (standard bases) can be computed using Buchberger's algorithm (Buchberger's algorithm as enhanced by Mora).

We call a standard basis $G = \{g_1,\dots, g_r\}\subset K[x]_>\setminus \{0\}$ *minimal*  if $\text{LM}_>(g_i)\neq \text{LM}_>(g_j)$ for $i\neq j$.

!!! note
    The definition of minimal above deviates from the definition in most textbooks as we do not ask that the leading coefficients of the standard basis elements are 1.

!!! note
    The standard bases returned by OSCAR are always minimal in the sense above.

We call a standard basis $G = \{g_1,\dots, g_r\}$ with respect to a global monomial ordering *reduced* if it is minimal and no term of $g_i$ i
s divisible by $\text{LM}_>(g_j)$, for $i\neq j$. Using power series expansions, we may extend this notion to local and mixed orderings.
However, while reduced standard bases can be computed in the global case, they may not be computable (in finitely many steps) in the other cases.

!!! note
    Given a monomial ordering $>$ on a free $K[x]$-module $F = K[x]^p$ with basis $e_1, \dots, e_p$,
    the above notation and results extend naturally to submodules of $K[x]_>^p$.

Here are the relevant OSCAR functions for computing Gröbner and standard bases. The elements of a
computed basis can be retrieved by using the `elements` function or and its alias `gens`.

```@docs
groebner_basis(I::MPolyIdeal;
	ord::MonomialOrdering = default_ordering(base_ring(I)),
	complete_reduction::Bool = false)
```
```@docs
standard_basis(I::MPolyIdeal;
	ord::MonomialOrdering = default_ordering(base_ring(I)),
	complete_reduction::Bool = false)
```
	
### Gröbner Bases with transformation matrix

```@docs
groebner_basis_with_transformation_matrix(I::MPolyIdeal;
	ordering::MonomialOrdering = default_ordering(base_ring(I)),
	complete_reduction::Bool=false)
```
```@docs
standard_basis_with_transformation_matrix(I::MPolyIdeal;
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
leading_ideal(G::Vector{T}; ordering::MonomialOrdering) where { T <: MPolyElem }
leading_ideal(I::MPolyIdeal; ordering::MonomialOrdering)
```

## Normal Forms

Given a polynomial $g\in K[x]$, an ideal $I\subset K[x]$, and a global monomial ordering
$>$ on the monomials in $x$, the fully reduced remainder $h$ in a standard expression on division
by the elements of a Gröbner basis of $I$ with respect to $>$ is uniquely determined by
$g$, $I$, and $>$ (and does not depend on the choice of Gröbner basis). We refer to such
a remainder as the *normal form*  of $g$ mod $I$, with respect to $>$.

```@docs
    normal_form(g::T, I::MPolyIdeal; 
      ordering::MonomialOrdering = default_ordering(base_ring(I))) where { T <: MPolyElem }
```

## Syzygies

We refer to the section on modules for more on syzygies.

```@docs
syzygy_generators(G::Vector{<:MPolyElem})
```

