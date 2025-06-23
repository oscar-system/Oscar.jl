# [Gröbner/Standard Bases Over Fields](@id gb_fields)

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

## Default Orderings

!!! note
    The OSCAR functions discussed in this section depend on a monomial `ordering` which is entered as a keyword argument.
    Given a polynomial ring $R$, the `default_ordering` for this is `degrevlex` except if $R$ is $\mathbb Z$-graded with
    positive weights. Then the corresponding `wdegrevlex` ordering is used. Given a free $R$-module $F$, the
    `default_ordering` is `default_ordering(R)*lex(gens(F))`.

```@docs
default_ordering(::MPolyRing)
```

Here are some illustrating OSCAR examples:

##### Examples

```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x, y, z])

julia> default_ordering(R)
degrevlex([x, y, z])

julia> F = free_module(R, 2)
Free module of rank 2 over R

julia> default_ordering(F)
degrevlex([x, y, z])*lex([gen(1), gen(2)])

julia> S, _ = grade(R, [1, 2, 3])
(Graded multivariate polynomial ring in 3 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y, z])

julia> default_ordering(S)
wdegrevlex([x, y, z], [1, 2, 3])
```

Expert users may temporarily choose a different default ordering for a given ring.
```@docs
with_ordering
```

## [Monomials, Terms, and More](@id monomials_terms_more)

Here are examples which indicate how to recover monomials, terms, and
more from a given polynomial.

```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x, y, z])

julia> f = 3*z^3+2*x*y+1
2*x*y + 3*z^3 + 1

julia> terms(f)
terms iterator of 3*z^3 + 2*x*y + 1

julia> collect(ans)
3-element Vector{QQMPolyRingElem}:
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
3-element Vector{Tuple{QQFieldElem, Vector{Int64}}}:
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

```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y])
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

julia> F = free_module(R, 3)
Free module of rank 3 over R

julia> f = (5*x*y^2-y^10+3)*F[1]+(4*x^3+2*y) *F[2]+16*x*F[3]
(5*x*y^2 - y^10 + 3)*e[1] + (4*x^3 + 2*y)*e[2] + 16*x*e[3]

julia> default_ordering(F)
degrevlex([x, y])*lex([gen(1), gen(2), gen(3)])

julia> collect(terms(f))
6-element Vector{FreeModElem{QQMPolyRingElem}}:
 -y^10*e[1]
 4*x^3*e[2]
 5*x*y^2*e[1]
 16*x*e[3]
 2*y*e[2]
 3*e[1]

julia> collect(terms(f, ordering = invlex(F)*lex(R)))
6-element Vector{FreeModElem{QQMPolyRingElem}}:
 5*x*y^2*e[1]
 -y^10*e[1]
 3*e[1]
 4*x^3*e[2]
 2*y*e[2]
 16*x*e[3]

julia> tail(f)
(5*x*y^2 + 3)*e[1] + (4*x^3 + 2*y)*e[2] + 16*x*e[3]

julia> leading_exponent(f)
([0, 10], 1)

julia> leading_exponent(f, ordering = invlex(F)*lex(R))
([1, 2], 1)
```

## Division With Remainder

The computation of Gröbner (standard) bases relies on multivariate division with remainder which is interesting
in its own right. If a monomial ordering $>$ is given, the basic idea is to mimic Euclidean division with remainder,
allowing more than one divisor: At each step of the resulting process, this amounts to removing the leading term of the
intermediate dividend, using the leading term of *some* divisor by which it is divisible. In its basic form, the process
works well if $>$ is global, but may not terminate for local and mixed orderings. In the latter case, Mora's division
algorithm, which relies on a more restricted selection strategy for the divisors to be used, comes to our rescue.

We discuss this in more detail:

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

```@docs
reduce(g::T, F::Vector{T}; 
    ordering::MonomialOrdering = default_ordering(parent(F[1]))) where T <: MPolyRingElem
```

```@docs
reduce_with_quotients(g::T, F::Vector{T}; 
    ordering::MonomialOrdering = default_ordering(parent(F[1]))) where T <: MPolyRingElem
```

```@docs
reduce_with_quotients_and_unit(g::T, F::Vector{T}; 
    ordering::MonomialOrdering = default_ordering(parent(F[1]))) where T <: MPolyRingElem
```

## Computing Gröbner/Standard Bases

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

We call a standard basis $G = \{g_1,\dots, g_r\}$ with respect to a global monomial ordering *reduced* if it is minimal and no term of
$g_i$ is divisible by $\text{LM}_>(g_j)$, for $i\neq j$. Using power series expansions, we may extend this notion to local and mixed orderings.
However, while reduced standard bases can be computed in the global case, they may not be computable (in finitely many steps) in the other cases.

!!! note
    Given a monomial ordering $>$ on a free $K[x]$-module $F = K[x]^p$ with basis $e_1, \dots, e_p$,
    the above notation and results extend naturally to submodules of $K[x]_>^p$.

Here are the relevant OSCAR functions for computing Gröbner and standard bases. The elements of a
computed basis can be retrieved by using the `elements` function or its alias `gens`.

```@docs
groebner_basis(I::MPolyIdeal;
    ord::MonomialOrdering = default_ordering(base_ring(I)),
    complete_reduction::Bool = false, algorithm::Symbol = :buchberger)
```
```@docs
standard_basis(I::MPolyIdeal;
    ord::MonomialOrdering = default_ordering(base_ring(I)),
    complete_reduction::Bool = false, algorithm::Symbol = :buchberger)
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

!!! note
The strategy behind the `groebner_basis` function and the strategy behind the function `groebner_basis_with_transformation_matrix` differ. As a consequence, the computed generators may differ. Even if `complete_reduction` is set to `true`, the generators might still only agree up to multiplication by units.

### Gröbner Basis Conversion Algorithms

The performance of Buchberger's Gröbner basis algorithm is sensitive to the choice of monomial ordering.
A Gröbner basis computation with respect to a less favorable ordering such as `lex` may easily run out of
time or memory even in cases where a Gröbner basis computation with respect to a more efficient 
ordering such as `degrevlex` is very well feasible. *Gröbner basis conversion algorithms*
and the *Hilbert driven Buchberger algorithm*  discussed subsequently take their cue from this observation.

Gröbner basis conversion algorithms proceed along the following lines:
- Given an ideal $I$ of a multivariate polynomial ring over a field and a slow `destination_ordering`, compute a Gröbner basis for $I$ with respect to an appropriately chosen fast `start_ordering`.
- Convert the result to a Gröbner basis with respect to the given slow `destination_ordering`.
The algorithms differ in how they perform the conversion.


#### The FGLM Algorithm


```@docs
fglm(I::MPolyIdeal; start_ordering::MonomialOrdering = default_ordering(base_ring(I)),
    destination_ordering::MonomialOrdering)
```

#### Gröbner Walk Algorithms

### The Hilbert driven Buchberger Algorithm

Calling the functions `standard_basis` and `groebner_basis` with `algorithm = :hilbert` in OSCAR
triggers a version of the Hilbert driven Gröbner basis algorithm which proceeds along the following lines.

1. Given an ideal $I$ of a multivariate polynomial ring $R$ over a field $K$ and a slow `destination_ordering`, check whether $I$ is homogeneous with respect to the standard $\mathbb Z$-grading on $R$. If so, set `start_ordering` to `degrevlex` and go to step 3.

2. Check whether there exists a $\mathbb Z$-grading on $R$ with positive weights such that $I$ is homogeneous with respect to this grading. If so, let `start_ordering` be the corresponding weight ordering. If not, go to step 5.

3. Compute a Gröbner basis of $I$ with respect to `start_ordering` and use this Gröbner basis to compute the Hilbert function of $R/I$.

4. Compute a Gröbner basis with respect to `destination_ordering`,  proceeding by increasing (weighted) degree, and skipping all further Buchberger tests in a given (weighted) degree as soon as the leading terms found so far account for the Hilbert function in that (weighted) degree. Return the computed Gröbner basis.

5. Extend $R$ to a polynomial ring $S$ by appending an extra variable,
   equip $S$ with the standard $\mathbb Z$-grading, and let
   $I^{h}\subset S$ be the homogenization of $I$ with respect to the
   extra variable. Compute a Gröbner basis of $I$ with respect to
   `degrevlex` on `R`, and homogenize its elements to obtain a Gröbner
   basis of $I^{h}$ with respect to `degrevlex` on $S$. Use the latter
   basis to compute the Hilbert function of $S/I^{h}$. Extend
   `destination_ordering` to a block ordering on `S`. Following the
   recipe in step 4, compute a Gröbner basis of $S/I^{h}$ with respect
   to the extended ordering. Return the dehomogenization of this basis
   with respect to the extra variable.

If the characteristic of $K$ is zero,  by semi-continuity of the Hilbert function,
it is sufficient to perform step 3 for the reduction of $I$ modulo a conveniently
chosen prime number rather than for $I$ itself.

!!! note
    If appropriate weights and/or the Hilbert function with respect to appropriate weights
    are already known to the user, this information can be entered when calling the Hilbert
    driven Gröbner basis algorithm as follows:

```@docs
groebner_basis_hilbert_driven(I::MPolyIdeal{P};
    destination_ordering::MonomialOrdering,
    complete_reduction::Bool = false,
    weights::Vector{Int} = ones(Int, number_of_generators(base_ring(I))),
    hilbert_numerator::Union{Nothing, ZZPolyRingElem} = nothing) where {P <: MPolyRingElem}
```

### Faugère's F4 Algorithm

!!! warning "Expert function for computing Gröbner bases"
    With many adjustable keyword arguments, the following function provides low-level
    implementations of various versions of the Gröbner basis algorithm. Use these functions
    only if you know what you are doing.

```@docs
groebner_basis_f4( I::MPolyIdeal; initial_hts::Int=17, nr_thrds::Int=1, max_nr_pairs::Int=0, la_option::Int=2, reduce_gb::Int=1, info_level::Int=0)
```

## Leading Ideals


```@docs
leading_ideal(G::Vector{T}; ordering::MonomialOrdering = default_ordering(parent(G[1])))  where T <: MPolyRingElem
leading_ideal(I::MPolyIdeal; ordering::MonomialOrdering = default_ordering(base_ring(I)))
```

## Normal Forms

Given a polynomial $g\in K[x]$, an ideal $I\subset K[x]$, and a global monomial ordering
$>$ on the monomials in $x$, the fully reduced remainder $h$ in a standard expression on division
by the elements of a Gröbner basis of $I$ with respect to $>$ is uniquely determined by
$g$, $I$, and $>$ (and does not depend on the choice of Gröbner basis). We refer to such
a remainder as the *normal form*  of $g$ mod $I$, with respect to $>$.

```@docs
normal_form(g::T, I::MPolyIdeal; ordering::MonomialOrdering = default_ordering(base_ring(I))) where T <: MPolyRingElem
```

## Syzygies

We refer to the section on [modules](@ref modules_multivariate) for more on syzygies.

```@docs
syzygy_generators(G::Vector{<:MPolyRingElem})
```

