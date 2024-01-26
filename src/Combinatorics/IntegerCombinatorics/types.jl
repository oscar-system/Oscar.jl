################################################################################
#
#  Partition
#
################################################################################

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
box in row `i` and and column `j` has the **coordinates** `(i, j)`. Listing
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

# Examples
```jldoctest
julia> tab = young_tableau([[1, 2, 3], [4, 5], [6]])
+---+---+---+
| 1 | 2 | 3 |
+---+---+---+
| 4 | 5 |
+---+---+
| 6 |
+---+
```

# References
1. Wikipedia, [Young tableau](https://en.wikipedia.org/wiki/Young_tableau).
"""
struct YoungTableau{T} <: AbstractVector{AbstractVector{T}}
  t::Vector{Vector{T}}
end
