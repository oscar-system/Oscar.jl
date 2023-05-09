```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

# Affine Algebraic Sets
An affine algebraic set over an algebraically closed
field $\overline{k}$ is the common vanishing locus $V$ of
finitely many polynomials $f_1,\dots f_r \in \overline{k}[x_1,\dots x_n]$,
or equivalently of the ideal $I \subseteq \overline{k}[x_1,\dots x_n]$ they generate.
Structural questions about varieties can be answered by considering the corresponding ideal of vanishing.
For instance Hilbert's Nullstellensatz states that $V$ is empty if and only if
$I=(1)$.

In Oscar we work with algebraic sets over non-closed fields,
by viewing them as reduced schemes.
By abuse of terminology we say that a scheme is an affine algebraic set
if it is isomorphic to one. For example a hypersurface complement is an
affine algebraic set.
In particular, affine algebraic sets are not necessarily Zariski closed in their
ambient affine space.

More formally they are defined as follows:
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
But we can also convert an affine scheme.
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
