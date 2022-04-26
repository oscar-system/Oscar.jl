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
localizations therefore has to make sure the following three 
requirements are met: 

 1) The code for finitely generated modules must be functional over ``R``. TODO: What precisely does this mean? Do we have a 'module interface' that the user needs to implement?

 2) In particular, the user has to implement the method for computing syzygies; see below. TODO: Can this be subsumed under 1)? For instance, it could be implemented generically via computation of kernels.
 3) The user has to solve the *localization problem* by implementing `has_nonepmty_intersection(U::MultSetType, I::IdealType)` for the type `MultSetType` of multiplicative sets and the type `IdealType` of ideals in `R` that they would like to consider.
```@docs
    syz(A::MatrixElem{<:AbsLocalizedRingElem})
    has_nonepmty_intersection(U::AbsMultSet, I::Ideal)
```

As soon as the above requirements are met, the methods 
```@julia
   represents_element(u::FreeModElem{T}, M::SubQuo{T}) where {T<:AbsLocalizedRingElem}
   coordinates(u::FreeModElem{T}, M::SubQuo{T}) where {T<:AbsLocalizedRingElem}
   kernel(f::FreeModuleHom{DomType, CodType, Nothing}) where {T, DomType<:FreeMod{T}, CodType<:SubQuo{T}}
   kernel(f::SubQuoHom{DomType, CodType, Nothing}) where {T, DomType<:FreeMod{T}, CodType<:SubQuo{T}}
```
will be available for modules over $S$, i.e. for `T = elem_type(S)`. 
As can easily be seen, having these methods
is equivalent to $S = R[U^{-1}]$ being computable; hence all higher methods can be derived 
from these basic ones. TODO: Is this automatic from the module code? Can/should we achieve that?
    
A sample implementation for various localizations of multivariate polynomial rings 
can be found in `src/Modules/mpoly-localizations.jl`.
