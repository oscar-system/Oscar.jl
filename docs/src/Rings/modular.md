```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# [Introduction](@id modular_methods)


## Chinese Remainder Theorem

We have firstly the `crt` methods detailed by the `euclidean_interface`:

```@docs; canonical=false
crt(r1::T, m1::T, r2::T, m2::T; check::Bool=true) where T <: RingElement
crt(r::Vector{T}, m::Vector{T}; check::Bool=true) where T <: RingElement
```

```@docs; canonical=false
crt_with_lcm(r1::T, m1::T, r2::T, m2::T; check::Bool=true) where T <: RingElement
crt_with_lcm(r::Vector{T}, m::Vector{T}; check::Bool=true) where T <: RingElement
```


We also have these variants which use idempotents.
```@docs; canonical=false
crt(r1::AbsSimpleNumFieldOrderElem, i1::AbsNumFieldOrderIdeal{AbsSimpleNumField, AbsSimpleNumFieldElem}, r2::AbsSimpleNumFieldOrderElem, i2::AbsNumFieldOrderIdeal{AbsSimpleNumField, AbsSimpleNumFieldElem})
```

This variant is also supported in the following signatures:
`crt(r1::Hecke.AlgAssAbsOrdElem, i1::Hecke.AlgAssAbsOrdIdl, r2::Hecke.AlgAssAbsOrdElem, i2::Hecke.AlgAssAbsOrdIdl)`
`crt(r1::Hecke.RelNumFieldOrderElem, i1::Hecke.RelNumFieldOrderIdeal, r2::Hecke.RelNumFieldOrderElem, i2::Hecke.RelNumFieldOrderIdeal)`
`crt(r1::AbsSimpleNumFieldOrderElem, i1::AbsSimpleNumFieldOrderIdeal, r2::AbsSimpleNumFieldOrderElem, i2::AbsSimpleNumFieldOrderIdeal)`
`crt(a::Vector{Hecke.AlgAssAbsOrdElem}, I::Vector{Hecke.AlgAssAbsOrdIdl})`
`crt(a::Vector{Hecke.RelNumFieldOrderElem}, I::Vector{Hecke.RelNumFieldOrderIdeal})`
`crt(a::Vector{AbsSimpleNumFieldOrderElem}, I::Vector{AbsSimpleNumFieldOrderIdeal})`
`crt(a::Vector{ZZRingElem}, I::Vector{Hecke.ZZIdl})`
<!-- No docs for these -->
<!-- ```@docs
crt(r1::S, i1::T, r2::S, i2::T) where {S<:Union{Hecke.AlgAssAbsOrdElem, Hecke.RelNumFieldOrderElem, AbsSimpleNumFieldOrderElem}, T<:Union{Hecke.AlgAssAbsOrdIdl, Hecke.RelNumFieldOrderIdeal, AbsSimpleNumFieldOrderIdeal}}
crt(a::Vector{S}, I::Vector{T}) where {S<:Union{Hecke.AlgAssAbsOrdElem, Hecke.RelNumFieldOrderElem, AbsSimpleNumFieldOrderElem}, T<:Union{Hecke.AlgAssAbsOrdIdl, Hecke.RelNumFieldOrderIdeal, AbsSimpleNumFieldOrderIdeal}}
crt(a::Vector{ZZRingElem}, b::Vector{Hecke.ZZIdl})
``` -->

We can also prepare an environment through the use of `crt_env` to prepare data for application of the chinese remainder theorem for a fixed set of coprime moduli in some euclidean ring.
```@docs
crt_env(p::Vector{T}) where T
crt(b::Vector{T}, a::crt_env{T}) where T
crt_inv(a::T, c::crt_env{T}) where T
```
<!-- No docs for crt_signed -->
<!-- crt_signed(b::Vector{ZZRingElem}, a::crt_env{ZZRingElem}) -->
