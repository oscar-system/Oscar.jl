```@meta
CurrentModule = Oscar
DocTestSetup = quote
  using Oscar
end
```

# Lie algebra modules

Lie algebra modules in OSCAR are always finite dimensional and represented by the type
`LieAlgebraModule{C}`. Similar to other types in OSCAR, there is the corresponding
element type `LieAlgebraModuleElem{C}`.
The type parameter `C` is the element type of the coefficient ring.

```@docs
base_lie_algebra(V::LieAlgebraModule{C}) where {C<:RingElement}
zero(::LieAlgebraModule{C}) where {C<:RingElement}
iszero(::LieAlgebraModuleElem{C}) where {C<:RingElement}
dim(::LieAlgebraModule{C}) where {C<:RingElement}
basis(::LieAlgebraModule{C}) where {C<:RingElement}
basis(::LieAlgebraModule{C}, ::Int) where {C<:RingElement}
coefficients(::LieAlgebraModuleElem{C}) where {C<:RingElement}
coeff(::LieAlgebraModuleElem{C}, ::Int) where {C<:RingElement}
getindex(::LieAlgebraModuleElem{C}, ::Int) where {C<:RingElement}
symbols(::LieAlgebraModule{C}) where {C<:RingElement}
```

## Element constructors

```@docs
(::LieAlgebraModule{C})() where {C<:RingElement}
(::LieAlgebraModule{C})(::Vector{C}) where {C<:RingElement}
(::LieAlgebraModule{C})(::Vector{Int}) where {C<:RingElement}
(::LieAlgebraModule{C})(::MatElem{C}) where {C<:RingElement}
(::LieAlgebraModule{C})(::SRow{C}) where {C<:RingElement}
(::LieAlgebraModule{C})(::LieAlgebraModuleElem{C}) where {C<:RingElement}
(::LieAlgebraModule{C})(::Vector{T}) where {T<:LieAlgebraModuleElem{C}} where {C<:RingElement}
```

## Arithmetics
The usual arithmetics, e.g. `+`, `-`, and `*` (scalar multiplication), are defined for `LieAlgebraModuleElem`s.

The module action is defined as `*`.
```@docs
*(x::LieAlgebraElem{C}, v::LieAlgebraModuleElem{C}) where {C<:RingElement}
```

## Module constructors

```@docs
standard_module(::LinearLieAlgebra{C}) where {C<:RingElement}
dual(::LieAlgebraModule{C}) where {C<:RingElement}
direct_sum(::LieAlgebraModule{C}, ::LieAlgebraModule{C}) where {C<:RingElement}
tensor_product(::LieAlgebraModule{C}, ::LieAlgebraModule{C}) where {C<:RingElement}
exterior_power(::LieAlgebraModule{C}, ::Int) where {C<:RingElement}
symmetric_power(::LieAlgebraModule{C}, ::Int) where {C<:RingElement}
tensor_power(::LieAlgebraModule{C}, ::Int) where {C<:RingElement}
abstract_module(::LieAlgebra{C}, ::Int, ::Vector{<:MatElem{C}}, ::Vector{<:VarName}; ::Bool) where {C<:RingElement}
abstract_module(::LieAlgebra{C}, ::Int, ::Matrix{SRow{C}}, ::Vector{<:VarName}; ::Bool) where {C<:RingElement}
```

# Type-dependent getters

```@docs
is_standard_module(::LieAlgebraModule{C}) where {C<:RingElement}
is_dual(::LieAlgebraModule{C}) where {C<:RingElement}
is_direct_sum(::LieAlgebraModule{C}) where {C<:RingElement}
is_tensor_product(::LieAlgebraModule{C}) where {C<:RingElement}
is_exterior_power(::LieAlgebraModule{C}) where {C<:RingElement}
is_symmetric_power(::LieAlgebraModule{C}) where {C<:RingElement}
is_tensor_power(::LieAlgebraModule{C}) where {C<:RingElement}
base_module(::LieAlgebraModule{C}) where {C<:RingElement}
base_modules(::LieAlgebraModule{C}) where {C<:RingElement}
```
