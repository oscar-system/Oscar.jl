```@meta
CurrentModule = Oscar
DocTestSetup = quote
  using Oscar
end
```

# Partitions

```@docs
Partition
getindex_safe
```

## Generating and counting

```@docs
partitions(::Oscar.IntegerUnion)
number_of_partitions(::Oscar.IntegerUnion)
```

### Partitions with restrictions
> How many ways are there to pay one euro, using coins worth 1, 2, 5, 10, 20, 50, and/or 100
> cents? What if you are allowed to use at most two of each coin? 

This is Exercise 11 in [Knu11](@cite), Section 7.2.1.4 (page 408). It goes back to the famous
"Ways to change one dollar" problem, see [Pol56](@cite). Generally, the problem is to
generate and/or count partitions satisfying some restrictions. Of course, one could generate
the list of all partitions of 100 (there are about 190 million) and then filter the result
by the restrictions. But for certain types of restrictions there are much more efficient
algorithms. The functions in this section implement some of these. In combination with
Julia's [filter](https://docs.julialang.org/en/v1/base/collections/#Base.filter) function
one can also handle more general types of restrictions.

For example, there are precisely six ways for the second question in the exercise quoted
above:
```jldoctest
julia> partitions(100, [1,2,5,10,20,50], [2,2,2,2,2,2])
6-element Vector{Partition{Int64}}:
 [50, 50]
 [50, 20, 20, 10]
 [50, 20, 20, 5, 5]
 [50, 20, 10, 10, 5, 5]
 [50, 20, 20, 5, 2, 2, 1]
 [50, 20, 10, 10, 5, 2, 2, 1]
```
and there are 4562 ways for the first question in the exercise:
```jldoctest
julia> length(partitions(100, [1,2,5,10,20,50]))
4562
```
The original "Ways to change one dollar" problem has 292 solutions:
```jldoctest
julia> length(partitions(100,[1,5,10,25,50]))
292
```

```@docs
partitions(::T, ::Oscar.IntegerUnion, ::Oscar.IntegerUnion, ::Oscar.IntegerUnion) where T <: Oscar.IntegerUnion
partitions(::T, ::Oscar.IntegerUnion) where T <: Oscar.IntegerUnion
number_of_partitions(::Oscar.IntegerUnion, ::Oscar.IntegerUnion)
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
