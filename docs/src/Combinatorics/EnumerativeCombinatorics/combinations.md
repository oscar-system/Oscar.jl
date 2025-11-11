```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Combinations and multicombinations

A **combination** from a set $S$ is a selection $\lambda_1, \lambda_2 \dots \lambda_r$ of elements of $S$ taken without repetition; the order of the elements is considered not to matter. A $k$-combination is a combination consisting of $k$ elements.
A general reference on combinations is [Knu11](@cite), Section 7.2.1.3.

We represent combinations as objects of type `Combination{T}`, which is
a subtype of `AbstractVector{T}`, and thus a vector-like object.
We use a bespoke type `Combination{T}` instead of for example a `Vector{T}` to
ensure the semantics of such objects are clear.

For the most basic usage of selecting $k$ values from the set $S = \{1,\ldots,n\}$,
the entries of a combination are the values $\lambda_i$ in ascending order:
```jldoctest
julia> collect(combinations(4,2))
6-element Vector{Combination{Int64}}:
 [1, 2]
 [1, 3]
 [1, 4]
 [2, 3]
 [2, 4]
 [3, 4]
```
Here `T` will be `Int` by default, but we allow using smaller
integer types as this may improve performance or reduce memory consumption.
```jldoctest
julia> collect(combinations(Int16(4), Int16(2)))
6-element Vector{Combination{Int16}}:
 Int16[1, 2]
 Int16[1, 3]
 Int16[1, 4]
 Int16[2, 3]
 Int16[2, 4]
 Int16[3, 4]
```

We also allow using `combinations` to select `k` entries from an arbitrary vector `v`,
whose entries have any type `T`. We do not require that the elements of type `T`
have any kind of ordering, or even that the entries of `v` are unique.
This is equivalent to taking `combinations(n,k)`, where `n = length(v)`,
and using the resulting integer combinations as indices for selecting elements of `v`.
```jldoctest
julia> collect(combinations([2,2,1,"a"], 2))
6-element Vector{Combination{Any}}:
 Any[2, 2]
 Any[2, 1]
 Any[2, "a"]
 Any[2, 1]
 Any[2, "a"]
 Any[1, "a"]
```

Because `Combination` is a subtype of `AbstractVector`, many functions that can be used for vectors (1-dimensional arrays) can be used for combinations as well.
For example:
```jldoctest
julia> C = Combination([6, 4, 4, 2])
[6, 4, 4, 2]

julia> length(C)
4

julia> C[1]
6
```


## Generating and counting

```@docs
combinations
number_of_combinations
```
