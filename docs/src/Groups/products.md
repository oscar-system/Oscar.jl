```@meta
CurrentModule = Oscar
DocTestSetup = quote
  using Oscar
end
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
canonical_injection(G::DirectProductGroup, j::Int)
canonical_injections(G::DirectProductGroup)
canonical_projection(G::DirectProductGroup, j::Int)
canonical_projections(G::DirectProductGroup)
write_as_full(G::DirectProductGroup)
is_full_direct_product(G::DirectProductGroup)
```

## Semidirect products

```@docs
SemidirectProductGroup{S<:GAPGroup, T<:GAPGroup}
semidirect_product(N::S, f::GAPGroupHomomorphism{T,AutomorphismGroup{S}}, H::T) where S <: GAPGroup where T <: GAPGroup
normal_subgroup(G::SemidirectProductGroup)
acting_subgroup(G::SemidirectProductGroup)
homomorphism_of_semidirect_product(G::SemidirectProductGroup)
is_full_semidirect_product(G::SemidirectProductGroup)
canonical_injection(G::SemidirectProductGroup, n::Int)
canonical_projection(G::SemidirectProductGroup)
```

## Wreath products

```@docs
WreathProductGroup
wreath_product(G::T, H::PermGroup) where T<: GAPGroup
normal_subgroup(W::WreathProductGroup)
acting_subgroup(W::WreathProductGroup)
homomorphism_of_wreath_product(G::WreathProductGroup)
is_full_wreath_product(G::WreathProductGroup)
canonical_projection(W::WreathProductGroup)
canonical_injection(W::WreathProductGroup, n::Int)
canonical_injections(W::WreathProductGroup)
```
