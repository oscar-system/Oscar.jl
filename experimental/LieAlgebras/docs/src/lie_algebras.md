```@meta
CurrentModule = Oscar
DocTestSetup = quote
  using Oscar
end
```

# Lie algebras

Lie algebras in OSCAR are currently always finite dimensional, and represented by two different types,
namely `LinearLieAlgebra{C}` and `AbstractLieAlgebra{C}`, depending on whether a matrix
representation is available or not.
Both types are subtypes of `LieAlgebra{C}`. Similar to other types in OSCAR, each Lie algebra
type has a corresponding element type.
The type parameter `C` is the element type of the coefficient ring. 

```@docs
zero(::LieAlgebra{C}) where {C<:RingElement}
iszero(::LieAlgebraElem{C}) where {C<:RingElement}
dim(::LieAlgebra{C}) where {C<:RingElement}
basis(::LieAlgebra{C}) where {C<:RingElement}
basis(::LieAlgebra{C}, ::Int) where {C<:RingElement}
coefficients(::LieAlgebraElem{C}) where {C<:RingElement}
coeff(::LieAlgebraElem{C}, ::Int) where {C<:RingElement}
getindex(::LieAlgebraElem{C}, ::Int) where {C<:RingElement}
symbols(::LieAlgebra{C}) where {C<:RingElement}
```

## Special functions for `LinearLieAlgebra`s

```@docs
matrix_repr_basis(::LinearLieAlgebra{C}) where {C<:RingElement}
matrix_repr_basis(::LinearLieAlgebra{C}, ::Int) where {C<:RingElement}
matrix_repr(::LinearLieAlgebraElem{C}) where {C<:RingElement}
```

## Element constructors

`(L::LieAlgebra{C})()` returns the zero element of the Lie algebra `L`.

`(L::LieAlgebra{C})(x::LieAlgebraElem{C})` returns `x` if `x` is an element of `L`,
and fails otherwise.

`(L::LieAlgebra{C})(v)` constructs the element of `L` with coefficient vector `v`.
`v` can be of type `Vector{C}`, `Vector{Int}`, `SRow{C}`,
or `MatElem{C}` (of size $1 \times \dim(L)$).

If `L` is a `LinearLieAlgebra` of `dim(L) > 1`, the call
`(L::LinearLieAlgebra{C})(m::MatElem{C})` returns the Lie algebra element whose
matrix representation corresponds to `m`.
This requires `m` to be a square matrix of size `n > 1` (the dimension of `L`), and
to lie in the Lie algebra `L` (i.e. to be in the span of `basis(L)`).
The case of `m` being a $1 \times \dim(L)$ matrix still works as explained above.


## Arithmetics
The usual arithmetics, e.g. `+`, `-`, and `*`, are defined for `LieAlgebraElem`s.

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

## Relation to GAP Lie algebras

Using `Oscar.iso_oscar_gap(L)`, one can get an isomorphism from the OSCAR Lie algebra `L`
to some isomorphic GAP Lie algebra. For more details, please refer to [`iso_oscar_gap`](@ref).
