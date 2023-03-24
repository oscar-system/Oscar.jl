# Partitions

```@docs
Partition
getindex_safe
```

## Generating and counting

### Unrestricted partitions

```@docs
partitions(::Oscar.IntegerUnion)
num_partitions(::Oscar.IntegerUnion)
```

### Restricted partitions

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
