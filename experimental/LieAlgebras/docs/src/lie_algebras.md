```@meta
CurrentModule = Oscar
```

# Lie algebras

Lie algebras in OSCAR are always finite dimensional, and represented by two different types, 
namely `LinearLieAlgebra` and `AbstractLieAlgebra`, depending on whether a matrix representation is available or not.
Both types are subtypes of `LieAlgebra`. Similar to other types in OSCAR, each Lie algebra type has a corresponding element type.

```@docs
zero(::LieAlgebra{C}) where {C <: RingElement}
iszero(::LieAlgebraElem{C}) where {C <: RingElement}
dim(::LieAlgebra{C}) where {C <: RingElement}
basis(::LieAlgebra{C}) where {C <: RingElement}
basis(::LieAlgebra{C}, ::Int) where {C <: RingElement}
coefficients(::LieAlgebraElem{C}) where {C <: RingElement}
coeff(::LieAlgebraElem{C}, ::Int) where {C <: RingElement}
getindex(::LieAlgebraElem{C}, ::Int) where {C <: RingElement}
```

## Special functions for `LinearLieAlgebra`s

```@docs
matrix_repr_basis(::LinearLieAlgebra{C}) where {C <: RingElement}
matrix_repr_basis(::LinearLieAlgebra{C}, ::Int) where {C <: RingElement}
matrix_repr(::LinearLieAlgebraElem{C}) where {C <: RingElement}
```
```
(::LinearLieAlgebra{C})(::MatElem{C}) where {C <: RingElement}

```

## Element constructors

```
(::LieAlgebra{C})() where {C<:RingElement}
(::LieAlgebra{C})(::Vector{C}) where {C<:RingElement}
(::LieAlgebra{C})(::Vector{Int}) where {C<:RingElement}
(::LieAlgebra{C})(::MatElem{C}) where {C<:RingElement}
(::LieAlgebra{C})(::LieAlgebraElem{C}) where {C<:RingElement}
```

## Arithmetics
The usual arithmetics, e.g. `+`, `-`, and `*` are defined for `LieAlgebraElem`s.

!!! warning
    Please note that `*` refers to the Lie bracket and is thus not associative.

## Lie algebra constructors

```@docs
lie_algebra
```

## Classical Lie algebras

```@docs
general_linear_lie_algebra(R::Ring, n::Int)
special_linear_lie_algebra(R::Ring, n::Int)
special_orthogonal_lie_algebra(R::Ring, n::Int)
```
