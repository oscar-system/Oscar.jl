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
base_lie_algebra(V::LieAlgebraModule)
zero(::LieAlgebraModule)
iszero(::LieAlgebraModuleElem)
dim(::LieAlgebraModule)
basis(::LieAlgebraModule)
basis(::LieAlgebraModule, ::Int)
coefficients(::LieAlgebraModuleElem)
coeff(::LieAlgebraModuleElem, ::Int)
getindex(::LieAlgebraModuleElem, ::Int)
symbols(::LieAlgebraModule)
```

## Element constructors

`(V::LieAlgebraModule{C})()` returns the zero element of the Lie algebra module `V`.

`(V::LieAlgebraModule{C})(v::LieAlgebraModuleElem{C})` returns `v` if `v` is an element of `L`. If `V` is the dual module of the parent of `v`, it returns the dual of `v`. In all other cases, it fails.

`(V::LieAlgebraModule{C})(v)` constructs the element of `V` with coefficient vector `v`. `v` can be of type `Vector{C}`, `Vector{Int}`, `SRow{C}`, or `MatElem{C}` (of size $1 \times \dim(L)$).

`(V::LieAlgebraModule{C})(a::Vector{T}) where {T<:LieAlgebraModuleElem{C}})`: If `V` is a direct sum, return its element, where the $i$-th component is equal to `a[i]`.
If `V` is a tensor product, return the tensor product of the `a[i]`.
If `V` is a exterior (symmetric, tensor) power, return the wedge product
(product, tensor product) of the `a[i]`.
Requires that `a` has a suitable length, and that the `a[i]` are elements of the correct modules,
where _correct_ depends on the case above.


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
is_standard_module(::LieAlgebraModule)
is_dual(::LieAlgebraModule)
is_direct_sum(::LieAlgebraModule)
is_tensor_product(::LieAlgebraModule)
is_exterior_power(::LieAlgebraModule)
is_symmetric_power(::LieAlgebraModule)
is_tensor_power(::LieAlgebraModule)
base_module(::LieAlgebraModule{C}) where {C<:RingElement}
base_modules(::LieAlgebraModule{C}) where {C<:RingElement}
```
