```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

# Affine Algebraic Sets
Let $\mathbb{A}^n(k)=k^n$ be the affine space of dimension $n$ over a field $k$.
For finitely many multivariate polynomials $f_1, \dots f_r \in k[x_1,\dots x_n]$
and $I = (f_1, \dots f_r) \subseteq k[x_1,\dots x_n]$ the ideal they generate,
we denote by $X = V(I)$ the variety defined by the ideal $I$.

If $k \subseteq K$ is any field extension, we denote by
\[X(K)= \{ P \in \mathbb{A}^n(K) \mid f_1(P)=\dots = f_n(P)=0\}=\{P \in \mathbb{A}^n(K) \mid \forall f\in I : f(P)=0\}\]
the set of $K$-rational points of $X$. We call $k$ the base field of $X$.

To study the $k$-rational points, one first considers the solutions $X(K)$
over an algebraically closed field extension $K$ of $k$.
Then in a second step one studies $X(k)$ as a subset of $X(K)$.

The first step involves calculations with ideals.
For instance Hilbert's Nullstellensatz implies that $X(K)$ is empty
if and only if the ideal $I=(1)$, which can be decided by computing a
Gröbner basis of $I$.
Typically, these calculations are carried out in $k[x_1,\dots x_n]$ without taking any field extensions.

The second step involves methods from number theory (if $k$ is a number field)
and from real algebraic geometry if $k = \mathbb{R}$.

Algebraic sets in Oscar are designed for the first step.
Most of their properties should be interpreted as properties
of the set $X(K)$ of their $K$-rational over an algebraic closure $K$.
A notable exception to this rule is the property of being irreducible, which
depends on the base field $k$.

For an example let $x^2+1 \in \mathbb{Q}[x]$. The variety $X = V(x^2+1)$ has
$\mathbb{Q}$ as its base field, $X(\mathbb{Q})=\emptyset$ but $X(\mathbb{Q}) = \{i, -i\}$.
Correspondingly, $X$ is irreducible over $\mathbb{R}$ (because the polynomial $x^2+1$ is)
but reducible over $\mathbb{C}$, just like $x^2+1=(x+i)(x-i)$ is reducible.
We say that $X$ is geometrically reducible.

```jldoctest
julia> R,(x,) = polynomial_ring(QQ,[:x])
(Multivariate Polynomial Ring in x over Rational Field, QQMPolyRingElem[x])

julia> X = vanishing_locus(x^2+1)
Vanishing locus
  in Affine 1-space over Rational Field
  of ideal(x^2 + 1)

julia> is_empty(X)
false

julia> irreducible(X)
true

julia> is_geometrically_irreducible(X)
false
```
For another example consider the dimension of $X$.
It is the dimension of $X(K)$ and not of $X(k)$ which in the following example
is even empty.
```jldoctest
julia> A = affine_space(QQ,2);

julia> (x1, x2) = coordinates(A)

julia> X = vanishing_locus(x^2 + 1)
Vanishing locus
  in Affine 2-space over Rational Field
  of ideal(x^2 + 1)

julia> dim(X)
1

```

## Relation to Schemes
Alternatively one may consider the variety as a scheme over $k$.
Indeed, $X = V(x^2+1)$ viewed as a scheme over $\mathbb{Q}$, has additional
scheme-theoretic points and is therefore is not empty.

Many constructions involving varieties lead naturally to schemes.
For instance the intersection of $X = V(x^2 - y)$ and $Y = V(y)$ as
sets is the point ${(0,0)}$. As a scheme it is defined by $V(x^2, y)$ which
can be interpreted as a point of multiplicity $2$ and reflects the information
that the intersection of $X$ and $Y$ is tangential in $(0,0)$.

Therefore we have two methods
- `set_theoretic_intersection` which can be thought of $X(K)\cap Y(K)$
- `intersection` which is the scheme theoretic intersection

!!! note
    If a construction returns a scheme $Z$, but you want to ignore the scheme
    structure, call the function `affine_algebraic_set(Z)` to convert the scheme
    $Z$ to an affine algebraic set.

For example `affine_algebraic_set(intersection(X, Y))`
is equivalent to `set_theoretic_intersection(X, Y)`.

### The following is for the experts familiar with schemes.
By abuse of terminology we say that a scheme is an affine algebraic set
if it is isomorphic to one. For example a hypersurface complement is an
affine algebraic set.
In particular, affine algebraic sets are not necessarily Zariski closed in their
ambient affine space.
```@docs
AbsAffineAlgebraicSet
```

## Constructors
The recommended way to create an algebraic set is as a
vanishing locus of an ideal or a polynomial.
```@docs
vanishing_locus(I::MPolyIdeal{<:MPolyElem}; check::Bool=true)
vanishing_locus(f::MPolyElem; check::Bool=true)
```
Convert an affine scheme to an affine algebraic set in order to ignore
its (non-reduced) scheme structure.
```@docs
affine_algebraic_set(X::Spec; check::Bool=true)
```

```@docs
set_theoretic_intersection(X::AbsAffineAlgebraicSet, Y::AbsAffineAlgebraicSet)
closure(X::AbsAffineAlgebraicSet{<:Field})
```

## Attributes
In addition to the attributes inherited from [Affine schemes](@ref)
the following are available.
```@docs
irreducible_components(X::AffineAlgebraicSet)
geometric_irreducible_components(X::AffineAlgebraicSet)
vanishing_ideal(X::AbsAffineAlgebraicSet)
ideal(X::AbsAffineAlgebraicSet{<:Field})
```

## Methods
Inherited from [Affine schemes](@ref)
## Properties
Inherited from [Affine schemes](@ref)
