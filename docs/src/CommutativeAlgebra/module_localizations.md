```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["module_localizations.md"]
```

# Localizations of modules over computable rings

For localizations of modules, there exists a generic implementation of 
the common methods such as membership tests, kernel computations, etc. 
based on the work of Barakat, Posur, et. al; see [Pos18](@cite).

Let $R$ be a ring of type `<:Ring`, $U \subset R$ a multiplicative set of type `<:AbsMultSet` 
and $S = R[U^{-1}]$ the localization of $R$ at $U$. Recall that $R$ 
is *computable* if one can compute *syzygies* and *lifts* over $R$. 
The results from [Pos18](@cite), Theorem 3.9, assert that then also the localization $S$ is 
computable, provided that there exists a solution to the *localization problem* 
(Definition 3.8, [Pos18](@cite) and below). 

The user who wishes to use the generic code for 
localizations therefore has to make sure the following two 
requirements are met: 

 1) The code for finitely generated modules and ideals must be functional over ``R``, including the computation of `coordinates` and `kernel`. 

 2) The user has to solve the *localization problem* by implementing `has_nonepmty_intersection(U::MultSetType, I::IdealType)` for the type `MultSetType` of multiplicative sets and the type `IdealType` of ideals in `R` that they would like to consider.
```@docs
    has_nonepmty_intersection(U::AbsMultSet, I::Ideal)
```
**Note:** In order to clear denominators of row vectors, the generic code uses the method `lcm(v::Vector{T})` where `T = elem_type(R)`. 
If no such method already exists, this has to also be provided; in the worst case by simply returning the product of the denominators. 

As soon as the above requirements are met, the methods 
```@julia
   represents_element(u::FreeModElem{T}, M::SubQuo{T}) where {T<:AbsLocalizedRingElem}
   coordinates(u::FreeModElem{T}, M::SubQuo{T}) where {T<:AbsLocalizedRingElem}
   kernel(f::FreeModuleHom{DomType, CodType, Nothing}) where {T, DomType<:FreeMod{T}, CodType<:SubQuo{T}}
   kernel(f::SubQuoHom{DomType, CodType, Nothing}) where {T, DomType<:FreeMod{T}, CodType<:SubQuo{T}}
   iszero(a::SubQuoElem{T}) where {T<:AbsLocalizedRingElem}
```
will be available for modules over $S$, i.e. for `T = elem_type(S)`. 
As can easily be seen, having the first three of these methods
is already equivalent to $S = R[U^{-1}]$ being computable; hence all higher methods can be derived 
from these basic ones. 
    
A sample implementation for various localizations of multivariate polynomial rings 
can be found in `src/Modules/mpoly-localizations.jl`. A modified version for localizations 
of affine algebras which also overwrites some of the generic methods is in 
`src/Modules/mpolyquo-localizations.jl`
