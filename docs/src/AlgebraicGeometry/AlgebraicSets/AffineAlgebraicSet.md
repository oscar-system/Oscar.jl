```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Affine Algebraic Sets

## Introduction
Let $\mathbb{A}^n(k)=k^n$ be the affine space of dimension $n$ over a field $k$.
For finitely many multivariate polynomials $f_1, \dots f_r \in k[x_1,\dots x_n]$
and $I = (f_1, \dots f_r) \subseteq k[x_1,\dots x_n]$ the ideal they generate,
we denote by $X = V(I)$ the (affine) algebraic set defined by the ideal $I$
and call $k$ its base field.

If $k \subseteq K$ is any field extension, we denote the set of $K$-points of $X$ by

$$\begin{aligned}X(K) &= \{ P \in \mathbb{A}^n(K) \mid f_1(P)=\dots = f_n(P)=0\}\\&=\{P \in \mathbb{A}^n(K) \mid \forall f\in I : f(P)=0\}.\end{aligned}$$

Most properties of the algebraic set $X$ refer to $X(K)$ where $K$ is an
algebraically closed field. For instance `is_empty` returns whether $X(K) = \emptyset$.

Exceptions to the rule, that we refer to $X(K)$, are documented in the respective methods.
For example the property of being irreducible depends on $k$:
The algebraic set $X = V(x^2+y^2) \subseteq \mathbb{A}^2$ is irreducible over
$k = \mathbb{R}$. But it is the union of two lines over $K = \mathbb{C}$,
i.e. $X$ is irreducible but geometrically reducible.
See
[`is_irreducible(X::AbsAffineScheme{<:Field, <:MPolyAnyRing})`](@ref) for details.

## Rational points
To study the $k$-points, also called $k$-rational points, of the algebraic set $X$
one first considers the solutions $X(K)$ over an algebraically closed field
extension $K$ of $k$. Then in a second step one studies $X(k)$ as a subset of $X(K)$.

The first step involves calculations with ideals.
For instance Hilbert's Nullstellensatz implies that $X(K)$ is empty
if and only if the ideal $I=(1)$. This is decided by an ideal membership test relying on a
Gröbner basis computation of $I$ and can be carried out in $k[x_1,\dots x_n]$
without taking any field extensions.

The second step involves methods from number theory (if $k$ is a number field)
or from real algebraic geometry (if $k = \mathbb{R}$).

Algebraic sets in Oscar are designed for the first step.
Most of their properfties should be interpreted as properties
of the set $X(K)$ of their $K$-points over an algebraic closure $K$.

## Relation to Schemes
One may view an (affine) algebraic set as a geometrically reduced (affine)
scheme over a field $k$.

Many constructions involving varieties lead naturally to schemes.
For instance the intersection of $X = V(x^2 - y)$ and $Y = V(y)$ as
sets is the point ${(0,0)}=V(x,y)$. As a scheme the intersection is defined by the ideal
$(x^2, y)$ which can be interpreted as a point of multiplicity $2$ and contains
the information that the intersection of $X$ and $Y$ is tangential in $(0,0)$.

Therefore we have two methods
- [`set_theoretic_intersection(::AbsAffineAlgebraicSet)`](@ref) which can be thought of as $X(K)\cap Y(K)$
- [`intersect(::AbsAffineAlgebraicSet)`](@ref) which is the scheme theoretic intersection

!!! note
    If a construction returns a scheme $Z$, but you want to ignore the scheme
    structure, call the function `algebraic_set(Z)` to convert the scheme
    $Z$ to an affine algebraic set.

For example `algebraic_set(intersect(X, Y))`
is equivalent to `set_theoretic_intersection(X, Y)`.

Internally an `AffineAlgebraicSet` is constructed from a possibly
non-reduced affine scheme, which we call the `fat_scheme` of `X` as opposed to the `reduced_scheme` of `X` which we refer to as the `underlying_scheme`.
```@docs
fat_ideal(X::AffineAlgebraicSet{<:Field})
fat_scheme(X::AffineAlgebraicSet)
underlying_scheme(X::AffineAlgebraicSet)
```

### More general affine algebraic sets
By abuse of terminology we say that a scheme is an affine algebraic set
if it is isomorphic to one. For example a hypersurface complement is an
affine algebraic set.
In particular, we allow affine algebraic sets which are not necessarily
Zariski closed in their ambient affine space.
```@docs
AbsAffineAlgebraicSet
```

## Constructors
One can create an algebraic set from an ideal or a multivariate polynomial.
```@docs
algebraic_set(I::MPolyIdeal{<:MPolyRingElem}; check::Bool=true)
algebraic_set(f::MPolyRingElem; check::Bool=true)
```
Convert an affine scheme to an affine algebraic set in order to ignore
its (non-reduced) scheme structure.
```@docs
algebraic_set(X::AffineScheme; check::Bool=true)
```

```@docs
set_theoretic_intersection(X::AbsAffineAlgebraicSet, Y::AbsAffineAlgebraicSet)
closure(X::AbsAffineAlgebraicSet{<:Field})
```

## Attributes
In addition to the attributes inherited from [Affine schemes](@ref)
the following are available.
```@docs
irreducible_components(X::AbsAffineAlgebraicSet)
geometric_irreducible_components(X::AbsAffineAlgebraicSet)
vanishing_ideal(X::AbsAffineAlgebraicSet)
```

## Methods
Inherited from [Affine schemes](@ref)
## Properties
Inherited from [Affine schemes](@ref)
