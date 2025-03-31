```@meta
CurrentModule = Oscar
DocTestSetup = Oscar.doctestsetup()
```

# Multipartitions


Multipartitions are generalizations of partitions.
An $r$-component **multipartition** of an integer ``n``
is an $r$-tuple of partitions $\lambda = (\lambda_1, \lambda_2, \ldots, \lambda_r)$
where each $\lambda_i$ is a partition of some integer $n_i \geq 0$
and the $n_i$ sum to $n$.
For an overview of multipartitions, see [And08](@cite).

A multipartition can be encoded as a 1-dimensional array whose elements are partitions.
In OSCAR, we provide the parametric type `Multipartition{T}` which is a subtype of `AbstractVector{T}`,
where `T` can by any subtype of `IntegerUnion`.
By using smaller integer types when appropriate, this allows performance to be improved.

```@docs
multipartition
```
Because `Multipartition` is a subtype of `AbstractVector`, all functions that can be used for vectors (1-dimensional arrays) can be used for multipartitions as well.
```jldoctest
julia> MP = multipartition([[3, 2], [1, 1], [3, 1]])
Partition{Int64}[[3, 2], [1, 1], [3, 1]]

julia> length(MP)
3

julia> MP[1]
[3, 2]
```
However, usually, $|\lambda| := n$ is called the **size** of a multipartition $\lambda$.
In Julia, the function `size` for arrays already exists and returns the *dimension* of an array.
Instead, one can use the Julia function `sum` to get the sum of the parts.
```jldoctest
julia> MP = multipartition([[3, 2], [1, 1], [3, 1]])
Partition{Int64}[[3, 2], [1, 1], [3, 1]]

julia> sum(MP)
11
```

## Generating and counting

```@docs
multipartitions(::Oscar.IntegerUnion, ::Oscar.IntegerUnion)
number_of_multipartitions(::Oscar.IntegerUnion, ::Oscar.IntegerUnion)
```
