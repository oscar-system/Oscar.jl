```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# [Partitions](@id partitions_chapter)

A **partition** of a non-negative integer $n$ is a decreasing sequence $\lambda_1 \geq \lambda_2\geq \dots \geq \lambda_r$ of positive integers $\lambda_i$ such that $n = \lambda_1 + \dots + \lambda_r$.
The $\lambda_i$ are called the **parts** of the partition and $r$ is called the **length**.
General references on partitions are [Ful97](@cite) and [Knu11](@cite), Section 7.2.1.4.

A partition can be encoded as an array with elements $\lambda_i$.
In OSCAR, the parametric type `Partition{T}` is provided which is a subtype of `AbstractVector{T}`.
Here, `T` can be any subtype of `IntegerUnion`.
The parametric type allows to increase performance by using smaller integer types.

```@docs
partition
```
Because `Partition` is a subtype of `AbstractVector`, all functions that can be used for vectors (1-dimensional arrays) can be used for partitions as well.
```jldoctest
julia> P = partition(6, 4, 4, 2)
[6, 4, 4, 2]

julia> length(P)
4

julia> P[1]
6
```
However, usually, $|\lambda| := n$ is called the **size** of $\lambda$.
In Julia, the function `size` for arrays already exists and returns the *dimension* of an array.
Instead, one can use the Julia function `sum` to get the sum of the parts.
```jldoctest
julia> P = partition(6, 4, 4, 2)
[6, 4, 4, 2]

julia> sum(P)
16
```

In algorithms involving partitions it is sometimes convenient to be able to access parts
beyond the length of the partition and then one wants to get the value zero instead of an
error. For this, OSCAR provides the function `getindex_safe`:
```@docs
getindex_safe
```
If you are sure that `P[i]` exists, use `getindex` because this will be faster.

## Generating and counting

```@docs
partitions(::Oscar.IntegerUnion)
number_of_partitions(::Oscar.IntegerUnion)
```
For counting partitions, the Hardy-Ramanujan-Rademachen formula is used, see [Joh12](@cite) for details.
See also [Knu11](@cite), Section 7.2.1.4 and [OEIS](@cite), [A000041](https://oeis.org/A000041).

### Partitions with restrictions
> How many ways are there to pay one euro, using coins worth 1, 2, 5, 10, 20, 50, and/or 100
> cents? What if you are allowed to use at most two of each coin?

This is Exercise 11 in [Knu11](@cite), Section 7.2.1.4. It goes back to the famous
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
julia> collect(partitions(100, [1, 2, 5, 10, 20, 50], [2, 2, 2, 2, 2, 2]))
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
julia> length(collect(partitions(100, [1, 2, 5, 10, 20, 50])))
4562
```
The original "Ways to change one dollar" problem has 292 solutions:
```jldoctest
julia> length(collect(partitions(100, [1, 5, 10, 25, 50])))
292
```
```@docs
number_of_partitions(::Oscar.IntegerUnion, ::Oscar.IntegerUnion)
```
For counting the partitions the recurrence relation $p_k(n) = p_{k - 1}(n - 1) + p_k(n - k)$ is used, where $p_k(n)$ denotes the number of partitions of $n$ into $k$ parts; see [Knu11](@cite), Section
7.2.1.4, Equation (39), and also [OEIS](@cite), [A008284](https://oeis.org/A008284).

```@docs
partitions(::Oscar.IntegerUnion, ::Oscar.IntegerUnion, ::Oscar.IntegerUnion, ::Oscar.IntegerUnion)
partitions(::T, ::Vector{T}) where T <: Oscar.IntegerUnion
```

## Operations

The *conjugate* of a partition $\lambda$ is obtained by considering its Young diagram
(see [Tableaux](@ref)) and then flipping it along its main diagonal, see [Ful97](@cite), page 2, and [Knu11](@cite), Section 7.2.1.4.
```@docs
conjugate
```

## Relations

The **dominance order** on partitions is the partial order $\trianglerighteq$ defined by $\lambda \trianglerighteq\mu$ if and only if $\lambda_1 + \dots + \lambda_i \geq \mu_1 + \dots + \mu_i$ for all $i$.
If $\lambda\trianglerighteq\mu$ one says that $\lambda$ **dominates** $\mu$.
See [Ful97](@cite), page 26, and [Knu11](@cite), Section 7.2.1.4, Exercise 54.

Note that whereas the lexicographic ordering is a total ordering, the dominance ordering is not.
Further, [Knu11](@cite) says **majorizes** instead of **dominates** and uses the symbol $\succeq$ instead of $\trianglerighteq$.
```@docs
dominates
```
