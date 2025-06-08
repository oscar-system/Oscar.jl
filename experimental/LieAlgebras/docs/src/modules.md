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
*(x::LieAlgebraElem{C}, v::LieAlgebraModuleElem{C}) where {C<:FieldElem}
```

## Module constructors

```@docs
trivial_module(L::LieAlgebra, d::IntegerUnion=1)
standard_module(::LinearLieAlgebra)
dual(::LieAlgebraModule{C}) where {C<:FieldElem}
direct_sum(::LieAlgebraModule{C}, ::LieAlgebraModule{C}) where {C<:FieldElem}
tensor_product(::LieAlgebraModule{C}, ::LieAlgebraModule{C}) where {C<:FieldElem}
exterior_power(::LieAlgebraModule{C}, ::Int) where {C<:FieldElem}
symmetric_power(::LieAlgebraModule{C}, ::Int) where {C<:FieldElem}
tensor_power(::LieAlgebraModule{C}, ::Int) where {C<:FieldElem}
abstract_module(::LieAlgebra{C}, ::Int, ::Vector{<:MatElem{C}}, ::Vector{<:VarName}; ::Bool) where {C<:FieldElem}
abstract_module(::LieAlgebra{C}, ::Int, ::Matrix{SRow{C}}, ::Vector{<:VarName}; ::Bool) where {C<:FieldElem}
```

## Representation theory of semisimple Lie algebras in characteristic 0

### Functions concerning simple modules

```@docs
simple_module(::LieAlgebra, ::WeightLatticeElem)
dim_of_simple_module(::LieAlgebra, ::WeightLatticeElem)
dominant_weights(::LieAlgebra, ::WeightLatticeElem)
dominant_character(::LieAlgebra, ::WeightLatticeElem)
character(::LieAlgebra, ::WeightLatticeElem)
tensor_product_decomposition(::LieAlgebra, ::WeightLatticeElem, ::WeightLatticeElem)
```

### Functions concerning Demazure modules

```@docs
demazure_operator(::RootSpaceElem, ::Dict{WeightLatticeElem,<:IntegerUnion})
demazure_character(::LieAlgebra, ::WeightLatticeElem, ::WeylGroupElem)
```
