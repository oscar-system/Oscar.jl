# Crystals

According to Kashiwara a crystal is a set $B$ together with the following maps

This is realised in OSCAR through the `AbstractCrystal` and `AbstractCrystalElem` interfaces and the following methods

```@docs
ealpha(::AbstractCrystalElem, i::Int)
falpha(::AbstractCrystalElem, i::Int)
weight(::AbstractCrystalElem)
```

Additionally the following methods are provided
```@docs
ealpha(::AbstractCrystalElem, i::Int, n::Int)`
ealpha(::AbstractCrystalElem, i::Vector{Int}, n::Vector{Int})
falpha(::AbstractCrystalElem, i::Int, n::Int)`
falpha(::AbstractCrystalElem, i::Vector{Int}, n::Vector{Int})
```