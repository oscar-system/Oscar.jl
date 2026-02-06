```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# [Introduction](@id modular_methods)


## Chinese Remainder Theorem

The Chinese Remainder Theorem is a method to solve systems of linear congruences.
We can apply it whenever we are able to perform division with a remainder.
Given residues `r1`, ..., `rn` and pairwise coprime moduli `m1`, ... `mn`
we return an `x` so that `x = ri mod mi` for all i.
This solution is unique modulo `N = lcm(m1, ... mn)`.

Any ring that has a `divrem` method has access to some `crt` methods provided by the `euclidean_interface`:

```@docs; canonical=false
crt(r1::T, m1::T, r2::T, m2::T; check::Bool=true) where T <: RingElement
crt(r::Vector{T}, m::Vector{T}; check::Bool=true) where T <: RingElement
```

```@docs; canonical=false
crt_with_lcm(r1::T, m1::T, r2::T, m2::T; check::Bool=true) where T <: RingElement
crt_with_lcm(r::Vector{T}, m::Vector{T}; check::Bool=true) where T <: RingElement
```

If we need to apply the chinese remainder theorem repeatedly
for a fixed set of (coprime) moduli,
we can precompute all the necessary information and save it in a prepared environment.
```@docs
crt_env(p::Vector{T}) where T
crt(b::Vector{T}, a::crt_env{T}) where T
crt_inv(a::T, c::crt_env{T}) where T
```
<!-- No docs for crt_signed -->
<!-- crt_signed(b::Vector{ZZRingElem}, a::crt_env{ZZRingElem}) -->


We also have variants for number fields that can take as arguments residues
along with respective ideals of the corresponding ring of integers.
These are calculated using a decomposition into idempotents.
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
