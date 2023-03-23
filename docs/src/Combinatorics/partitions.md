# Partitions

```@docs
Partition
```

## Unrestricted partitions

```@docs
partitions(::Integer)
num_partitions(::Integer)
ascending_partitions
```

## Restricted partitions

```@docs
partitions(::Integer, ::Integer, ::Integer, ::Integer)
partitions(::Integer, ::Integer)
num_partitions(::Integer, ::Integer)
partitions(::Integer, ::Integer, ::Vector{Integer}, ::Vector{Integer})
partitions(::Integer, ::Vector{Integer})
```

## Operations

```@docs
conjugate
getindex_safe
```

## Relations

```@docs
dominates
```
