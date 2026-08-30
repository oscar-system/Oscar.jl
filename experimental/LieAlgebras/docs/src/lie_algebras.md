# Lie algebras

Lie algebras in OSCAR are currently always finite dimensional, and represented by two different types,
namely `LinearLieAlgebra{C}` and `AbstractLieAlgebra{C}`, depending on whether a matrix
representation is available or not.
Both types are subtypes of `LieAlgebra{C}`. Similar to other types in OSCAR, each Lie algebra
type has a corresponding element type.
The type parameter `C` is the element type of the coefficient ring. 

```@docs
zero(::LieAlgebra)
iszero(::LieAlgebraElem)
dim(::LieAlgebra)
basis(::LieAlgebra)
basis(::LieAlgebra, ::Int)
coefficients(::LieAlgebraElem)
coeff(::LieAlgebraElem, ::Int)
getindex(::LieAlgebraElem, ::Int)
symbols(::LieAlgebra)
characteristic(L::LieAlgebra)
```

## Special functions for `LinearLieAlgebra`s

```@docs
coerce_to_lie_algebra_elem(::LinearLieAlgebra{C}, ::MatElem{C}) where {C<:FieldElem}
matrix_repr_basis(::LinearLieAlgebra{C}) where {C<:FieldElem}
matrix_repr_basis(::LinearLieAlgebra{C}, ::Int) where {C<:FieldElem}
matrix_repr(::LinearLieAlgebraElem{C}) where {C<:FieldElem}
```

## Element constructors

`(L::LieAlgebra{C})()` returns the zero element of the Lie algebra `L`.

`(L::LieAlgebra{C})(x::LieAlgebraElem{C})` returns `x` if `x` is an element of `L`,
and fails otherwise.

`(L::LieAlgebra{C})(v)` constructs the element of `L` with coefficient vector `v`.
`v` can be of type `Vector{C}`, `Vector{Int}`, `SRow{C}`,
or `MatElem{C}` (of size $1 \times \dim(L)$).


## Arithmetics
The usual arithmetics, e.g. `+`, `-`, and `*`, are defined for `LieAlgebraElem`s.

!!! warning
    Please note that `*` refers to the Lie bracket and is thus not associative.

## Properties

```@docs
is_abelian(L::LieAlgebra)
is_nilpotent(L::LieAlgebra)
is_perfect(L::LieAlgebra)
is_semisimple(L::LieAlgebra)
is_simple(L::LieAlgebra)
is_solvable(L::LieAlgebra)
```

## More functions

```@docs
center(L::LieAlgebra)
centralizer(L::LieAlgebra, xs::AbstractVector{<:LieAlgebraElem})
centralizer(L::LieAlgebra, x::LieAlgebraElem)
derived_algebra(L::LieAlgebra)
derived_series(L::LieAlgebra)
lower_central_series(L::LieAlgebra)
```

## Lie algebra constructors

```@docs
lie_algebra
```

## Classical Lie algebras

```@docs
abelian_lie_algebra(R::Field, n::Int)
general_linear_lie_algebra(R::Field, n::Int)
special_linear_lie_algebra(R::Field, n::Int)
special_orthogonal_lie_algebra
symplectic_lie_algebra
```

## Relation to GAP Lie algebras

Using `Oscar.iso_oscar_gap(L)`, one can get an isomorphism from the OSCAR Lie algebra `L`
to some isomorphic GAP Lie algebra. For more details, please refer to [`iso_oscar_gap`](@ref).
