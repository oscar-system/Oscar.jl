# Partitions

```@docs
Partition
getindex_safe
```

## Generating and counting

```@docs
partitions(::Oscar.IntegerUnion)
num_partitions(::Oscar.IntegerUnion)
```

### Partitions with restrictions

Often, one is interested not in all partitions but in partitions satisfying some restrictions, for example on the number of parts or on the parts themselves. A typical application is to find all ways to pay 1 Euro (=100 Cents) using coins of smaller denomination, i.e. coins worth 1, 2, 5, 10, 20, and 50 cents. A refinement is to allow only at most two of each coin. See [Knu11](@cite), Exercise 11 in Section 7.2.1.4 (page 408) for this particular problem. Instead of filtering the list of partitions of 100 (which contains about 190 million entries) by the restrictions, there are dedicated algorithms performing such tasks efficiently without generating all partitions first. 
 
```@docs
partitions(::T, ::Oscar.IntegerUnion, ::Oscar.IntegerUnion, ::Oscar.IntegerUnion) where T <: Oscar.IntegerUnion
partitions(::T, ::Oscar.IntegerUnion) where T <: Oscar.IntegerUnion
num_partitions(::Oscar.IntegerUnion, ::Oscar.IntegerUnion)
partitions(::T, ::Oscar.IntegerUnion, ::Vector{T}, ::Vector{S}) where {T <: Oscar.IntegerUnion, S<:Oscar.IntegerUnion}
partitions(::T, ::Vector{T}, ::Vector{S}) where {T <: Oscar.IntegerUnion, S<:Oscar.IntegerUnion}
partitions(::T, ::Vector{T}) where T <: Oscar.IntegerUnion
```

## Operations

```@docs
conjugate
```

## Relations

```@docs
dominates
```
