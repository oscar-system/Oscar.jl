```@meta
CurrentModule = Oscar
DocTestSetup = quote
  using Oscar
end
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["products.md"]
```

# Products of groups

## Direct products

```@docs
DirectProductGroup
direct_product(L::AbstractVector{<:GAPGroup}; morphisms=false)
inner_direct_product(L::AbstractVector{T}; morphisms=false) where T<:Union{PcGroup,FPGroup}
number_of_factors(G::DirectProductGroup)
cartesian_power(G::GAPGroup, n::Int)
inner_cartesian_power(G::T, n::Int; morphisms=false) where T<: GAPGroup
factor_of_direct_product(G::DirectProductGroup, j::Int)
as_perm_group(G::DirectProductGroup)
as_polycyclic_group(G::DirectProductGroup)
embedding(G::DirectProductGroup, j::Int)
projection(G::DirectProductGroup, j::Int)
write_as_full(G::DirectProductGroup)
isfull_direct_product(G::DirectProductGroup)
```

## Semidirect products

```@docs
SemidirectProductGroup{S<:GAPGroup, T<:GAPGroup}
semidirect_product(N::S, f::GAPGroupHomomorphism{T,AutomorphismGroup{S}}, H::T) where S <: GAPGroup where T <: GAPGroup
normal_subgroup(G::SemidirectProductGroup)
acting_subgroup(G::SemidirectProductGroup)
homomorphism_of_semidirect_product(G::SemidirectProductGroup)
isfull_semidirect_product(G::SemidirectProductGroup)
embedding(G::SemidirectProductGroup, n::Int)
projection(G::SemidirectProductGroup)
```

## Wreath products

```@docs
WreathProductGroup
wreath_product(G::T, H::PermGroup) where T<: GAPGroup
normal_subgroup(W::WreathProductGroup)
acting_subgroup(W::WreathProductGroup)
homomorphism_of_wreath_product(G::WreathProductGroup)
isfull_wreath_product(G::WreathProductGroup)
projection(W::WreathProductGroup)
embedding(W::WreathProductGroup, n::Int)
```
