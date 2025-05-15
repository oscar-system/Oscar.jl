```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Wreath Macdonald polynomials

The existence, integrality and positivity of wreath Macdonald polynomials
 has been conjectured by Haiman [Hai02](@cite) and proved by Bezrukavnikov
 and Finkelberg [BF14](@cite). Here we have implemented an algorithm which
 computes the wreath Macdonald polynomials as defined in the review by
 Orr and Shimozono on this topic [OS23](@cite).

Wreath Macdonald polynomials depend on two parameters. The first one is
 a multipartition since they are indexing the set of irreducible representations of
 the complex reflection group ``G(r,1,n)``. The second one is an element
 of the affine Weyl group of type ``A_r`` which is isomorphic to the semi-direct
 product of the finite Weyl group of type ``A_r`` (the symmetric group on ``r`` letters)
 and of the coroot lattice of type ``A_r``. The element of the coroot lattice is
 given in the canonical basis. It is then the sublattice of ``\mathbb{Z}^r`` of elements summing
 up to zero.

```@docs
wreath_macs(n::Int, r::Int, wperm::PermGroupElem, coroot::Vector{Int})
wreath_mac(lbb::Multipartition, wperm::PermGroupElem, coroot::Vector{Int})
```

Compare with Example 3.15 in [OS23](@cite).

```jldoctest
julia> collect(multipartitions(1,3))
3-element Vector{Multipartition{Int64}}:
 Partition{Int64}[[], [], [1]]
 Partition{Int64}[[], [1], []]
 Partition{Int64}[[1], [], []]

julia> wreath_macs(1,3,(@perm 3 (1,2,3)), [0,1,-1])[[3, 2, 1],[3, 2, 1]]
[1   q^2     q]
[1     t     q]
[1     t   t^2]

```
