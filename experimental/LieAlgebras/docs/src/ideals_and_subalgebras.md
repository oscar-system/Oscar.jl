```@meta
CurrentModule = Oscar
DocTestSetup = Oscar.doctestsetup()
```

# Ideals and Lie subalgebras

Ideals and Lie subalgebras are represented by the types `LieAlgebraIdeal` and
`LieSubalgebra` respectively.
They are used similarly in most cases.

## Functions

### Ideals

```@docs
dim(::LieAlgebraIdeal)
basis(::LieAlgebraIdeal)
basis(::LieAlgebraIdeal, ::Int)
Base.in(::LieAlgebraElem{C}, ::LieAlgebraIdeal{C}) where {C<:FieldElem}
bracket(::LieAlgebraIdeal{C,LieT}, ::LieAlgebraIdeal{C,LieT}) where {C<:FieldElem,LieT<:LieAlgebraElem{C}}
normalizer(::LieAlgebra, ::LieAlgebraIdeal)
centralizer(::LieAlgebra, ::LieAlgebraIdeal)
```

### Lie subalgebras

```@docs
dim(::LieSubalgebra)
basis(::LieSubalgebra)
basis(::LieSubalgebra, ::Int)
Base.in(::LieAlgebraElem{C}, ::LieSubalgebra{C}) where {C<:FieldElem}
bracket(::LieSubalgebra{C,LieT}, ::LieSubalgebra{C,LieT}) where {C<:FieldElem,LieT<:LieAlgebraElem{C}}
normalizer(::LieAlgebra, ::LieSubalgebra)
centralizer(::LieAlgebra, ::LieSubalgebra)
is_self_normalizing(S::LieSubalgebra)
```

## Constructors

### Ideals

```@docs
ideal(::LieAlgebra, ::Vector; is_basis::Bool=false)
ideal(::LieAlgebra{C}, ::LieAlgebraElem{C}) where {C<:FieldElem}
ideal(::LieAlgebra)
```

### Lie subalgebras

```@docs
sub(::LieAlgebra, ::Vector; is_basis::Bool=false)
sub(::LieAlgebra{C}, ::LieAlgebraElem{C}) where {C<:FieldElem}
sub(::LieAlgebra)
```

## Conversions

```@docs
lie_algebra(::LieSubalgebra)
lie_algebra(::LieAlgebraIdeal)
sub(::LieAlgebra{C}, ::LieAlgebraIdeal{C}) where {C<:FieldElem}
```
