```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# [Introduction](@id modular_methods)


## Chinese Remainder Theorem

```@docs
crt(r1::T, m1::T, r2::T, m2::T; check::Bool=true) where T <: RingElement
crt(r::Vector{T}, m::Vector{T}; check::Bool=true) where T <: RingElement
```

```@docs
crt_with_lcm(r1::T, m1::T, r2::T, m2::T; check::Bool=true) where T <: RingElement
crt_with_lcm(r::Vector{T}, m::Vector{T}; check::Bool=true) where T <: RingElement
```

```@docs
crt(r1::S, i1::T, r2::S, i2::T) where {S<:Union{Hecke.AlgAssAbsOrdElem, Hecke.RelNumFieldOrderElem, AbsSimpleNumFieldOrderElem}, T<:Union{Hecke.AlgAssAbsOrdIdl, Hecke.RelNumFieldOrderIdeal, AbsSimpleNumFieldOrderIdeal}}
crt(a::Vector{S}, I::Vector{T}) where {S<:Union{Hecke.AlgAssAbsOrdElem, Hecke.RelNumFieldOrderElem, AbsSimpleNumFieldOrderElem}, T<:Union{Hecke.AlgAssAbsOrdIdl, Hecke.RelNumFieldOrderIdeal, AbsSimpleNumFieldOrderIdeal}}
crt(a::Vector{ZZRingElem}, b::Vector{Hecke.ZZIdl})
```

```@docs
crt(b::Vector{T}, a::crt_env{T}) where T
crt_signed(b::Vector{ZZRingElem}, a::crt_env{ZZRingElem})
crt_inv(a::T, crt_env{T}) -> Vector{T}
```
