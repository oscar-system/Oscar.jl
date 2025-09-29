```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Combinations

A **combination** from a set $S$ is a selection $\lambda_1, \lambda_2 \dots \lambda_r$ of elements of $S$ taken without repetition; the order of the elements is considered not to matter. A $k$-combination is a combination consisting of $k$ elements.
A general reference on combinations is [Knu11](@cite), Section 7.2.1.3.

A combination can be encoded as an array with elements $\lambda_i$.
In OSCAR, the parametric type `Combination{T}` is provided which is a subtype of `AbstractVector{T}`.
Here, `T` can be any subtype of `IntegerUnion`.
The parametric type allows one to increase performance by using smaller integer types.
```julia
julia> @time collect(combinations(20,10));
  0.010399 seconds (184.76 k allocations: 26.782 MiB)

julia> @time collect(combinations(Int8(20),10));
  0.008780 seconds (184.76 k allocations: 12.686 MiB)
```


## Generating

```@docs
combinations(n::Int, k::Int)
combinations(v::AbstractVector, k::Int)
```

Because `Combination` is a subtype of `AbstractVector`, many functions that can be used for vectors (1-dimensional arrays) can be used for combinations as well.
For example:
```jldoctest
julia> C = Combination([6, 4, 4, 2])

julia> length(C)
4

julia> C[1]
6
```
