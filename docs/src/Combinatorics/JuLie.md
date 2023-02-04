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
    partitions(::Array{Integer,1}, ::Integer, ::Array{Integer,1}, ::Integer)
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

