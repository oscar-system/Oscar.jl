```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

# Projective Algebraic Sets
An projective algebraic set over an algebraically closed
field $\overline{k}$ is the common vanishing locus
$V\subseteq \mathbb{P}^n_{\overline{k}}$ of
finitely many homogeneous polynomials $f_1,\dots f_r \in \overline{k}[x_0,\dots x_n]$,
or equivalently of a homogeneous ideal $I \subseteq \overline{k}[x_0,\dots x_n]$ they generate.
Structural questions about varieties can be answered by considering the corresponding homogeneous ideal of vanishing.
For instance the projective Nullstellensatz states that $V$ is empty if and only if
$I\supseteq (x_0,\dots, x_n)$.

In Oscar we work with projective algebraic sets over non-closed fields,
by viewing them as reduced schemes. See [Projective schemes](@ref).

More formally they are defined as follows:
```@docs
AbsProjectiveAlgebraicSet
```

## Constructors
The recommended way to create an algebraic set is as a
vanishing locus of a homogeneous ideal or polynomial.
```@docs
vanishing_locus(I::MPolyIdeal{<:MPolyDecRingElem}; check::Bool=true)
vanishing_locus(p::MPolyDecRingElem; check::Bool=true)
```
Algebraic sets can also be constructed from projective schemes.
```@docs
projective_algebraic_set(X::AbsProjectiveScheme; check::Bool=true)
```

```@docs
set_theoretic_intersection(X::AbsProjectiveAlgebraicSet, Y::AbsProjectiveAlgebraicSet)
irreducible_components(X::AbsProjectiveAlgebraicSet)
geometric_irreducible_components(X::AbsProjectiveAlgebraicSet)
```

## Attributes
In addition to the attributes inherited from [Projective schemes](@ref)
the following are available.
```@docs
vanishing_ideal(X::AbsProjectiveAlgebraicSet)
ideal(X::AbsProjectiveAlgebraicSet)
```

## Methods
Inherited from [Projective schemes](@ref)
## Properties
Inherited from [Projective schemes](@ref)
