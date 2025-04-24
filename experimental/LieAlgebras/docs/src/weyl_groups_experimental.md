```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Weyl groups (experimental features)

This page is an addition to the documentation of [Weyl groups](@ref) with the additional experimental features.

# Exchange lemma and Braid moves

```@docs
braid_moves(::WeylGroup, ::Vector{UInt8}, ::Vector{UInt8})
apply_braid_move!(::Vector{UInt8}, ::Tuple{Int, Int, Int})
```

## Conversion to other group types

For many computations, it may be suitable to have a `WeylGroup` as a different kind of group object, to e.g. use functionality that is only available for that other type.

The conversion functions come in pairs: one only creates an isomorphic group object, the other also computes the isomorphism.

```@docs
isomorphism(::Type{FPGroup}, ::WeylGroup)
```

```@docs
isomorphism(::Type{PermGroup}, ::WeylGroup)
```

## Parabolic subgroups

```@docs
parabolic_subgroup(::WeylGroup, vec::Vector{<:Integer}, ::WeylGroupElem=one(W))
parabolic_subgroup_with_projection(::WeylGroup, ::Vector{<:Integer})
```

## Irreducible factors

```@docs
irreducible_factors(::WeylGroup)
```

## Inner direct products

```@docs
inner_direct_product(::AbstractVector{WeylGroup})
```