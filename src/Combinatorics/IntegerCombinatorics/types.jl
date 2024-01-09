################################################################################
#
#  Partition
#
################################################################################

@doc raw"""
    Partition{T<:IntegerUnion} <: AbstractVector{T}

A **partition** of a non-negative integer ``n`` is a decreasing sequence ``λ₁ ≥ λ₂ ≥ … ≥
λᵣ`` of positive integers ``λᵢ`` such that ``n = λ₁ + … + λᵣ``. The ``λᵢ`` are called the
**parts** of the partition and ``r`` is called the **length**.

A partition can be encoded as an array with elements ``λᵢ``. We provide the parametric type
`Partition{T}` which is a subtype of `AbstractVector{T}` where `T` can be any subtype of
`IntegerUnion`. All functions that can be used for vectors (1-dimensional arrays) can thus
be used for partitions as well. There is no performance impact by using an own type for
partitions rather than simply using arrays. The parametric type allows to increase
performance by using smaller integer types. For efficiency, the `partition` constructor does
not check whether the given array is indeed a decreasing sequence.

A partition can be created by either calling `partition` on an array of integers or by
calling `partition` with arguments being the sequence of parts, with the possibility
to provide the element type as the first argument.

# Examples
```jldoctest
julia> P = partition([6,4,4,2]) #The partition 6+4+4+2 of 16.
[6, 4, 4, 2]

julia> P = partition(6,4,4,2) #Same as above but less to type
[6, 4, 4, 2]

julia> length(P)
4

julia> P[1]
6
```
Usually, ``|λ| ≔ n`` is called the **size** of ``λ``. In Julia, the function `size` for
arrays already exists and returns the *dimension* of an array. Instead, you can use the
Julia function `sum` to get the sum of the parts.
```jldoctest
julia> P = partition(6,4,4,2)
[6, 4, 4, 2]

julia> sum(P)
16
```
You can create partitions with smaller integer types as follows.
```jldoctest
julia> P = partition(Int8,6,4,4,2) #Or partition(Int8[6,4,4,2])
Int8[6, 4, 4, 2]
```
There is a unique partition of 0, namely the **empty partition** (of length 0). It can be
created as follows.
```jldoctest
julia> P = partition() #Or partition([])
Int64[]
julia> sum(P)
0
julia> length(P)
0
julia> P = partition(Int8) #Or partition(Int8[])
Int8[]
```

# References
1. [Ful97](@cite)
2. [Knu11](@cite), Section 7.2.1.4 (starting on page 390).
"""
struct Partition{T<:IntegerUnion} <: AbstractVector{T}
  p::Vector{T}
end

################################################################################
#
#  Young Tableaux
#
################################################################################

@doc raw"""
    YoungTableau{T} <: AbstractVector{AbstractVector{T}}

A **Young diagram** is a diagram of finitely many empty "boxes" arranged
in left-justified rows, with the row lengths in non-increasing order. The
box in row `i` and and column `j` has the **coordinates** `(i,j)`. Listing
the number of boxes in each row gives a partition ``λ`` of a non-negative
integer `n` (the total number of boxes of the diagram). The diagram is
then said to be of **shape** ``λ``. Conversely, one can associate to any
partition ``λ`` a Young diagram in the obvious way, so Young diagrams are
just another way to look at partitions.

A **Young tableau** of shape ``λ`` is a filling of the boxes of the Young
diagram of ``λ`` with elements from some set. After relabeling we can (and
will) assume that we fill from a set of integers from ``1`` up to some number,
which in applications is often equal to `n`. We encode a tableau as an
array of arrays and we have implemented an own type `YoungTableau{T}`
as subtype of `AbstractVector{AbstractVector{T}}` to work with
tableaux. As for partitions, you may increase performance by casting
into smaller integer types, e.g.

For efficiency, we do not check whether the given array is really a
tableau, i.e. whether the structure of the array defines a partition.

# Examples
```jldoctest
julia> tab=young_tableau([[1,2,3],[4,5],[6]])
[[1, 2, 3], [4, 5], [6]]

julia> tab=young_tableau(Vector{Int8}[[2,1], [], [3,2,1]]) #Using 8 bit integers
Vector{Int8}[[2, 1], [], [3, 2, 1]]
```

# References
1. Wikipedia, [Young tableau](https://en.wikipedia.org/wiki/Young_tableau).
"""
struct YoungTableau{T} <: AbstractVector{AbstractVector{T}}
  t::Vector{Vector{T}}
end
