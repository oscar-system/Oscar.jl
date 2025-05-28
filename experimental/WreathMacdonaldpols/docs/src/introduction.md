```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Wreath Macdonald polynomials

The existence, integrality and positivity of wreath Macdonald polynomials
 has been conjectured by Haiman [Hai02](@cite) and proved by Bezrukavnikov
 and Finkelberg [BF14](@cite). When ``r=1``, wreath Macdonald polynomials are
 equal to the Haiman-Macdonald polynomials, used to prove the Macdonald positivity conjecture.

Here we have implemented an algorithm computing the wreath Macdonald
 polynomials as defined in the survey by Orr and Shimozono on this topic [OS23](@cite).

Wreath Macdonald polynomials depend on two parameters. The first parameter is
 an ``r``-multipartition of ``n``. The second parameter is an element of the affine Weyl group
 of type ``A^{(1)}_{r-1}`` which is isomorphic to the semi-direct product of the finite Weyl group
 of type ``A_{r-1}`` (the symmetric group on ``r`` letters) and of the coroot lattice of type ``A_{r-1}``.
 The element of the coroot lattice is given in the canonical basis. It is then the sublattice
of ``\mathbb{Z}^r`` of elements summing up to zero.

```@docs
wreath_macdonald_polynomial
wreath_macdonald_polynomials
```

Compare the following computation with Example 3.15 in [OS23](@cite).

```jldoctest
julia> collect(multipartitions(1,3))
3-element Vector{Multipartition{Int64}}:
 Partition{Int64}[[], [], [1]]
 Partition{Int64}[[], [1], []]
 Partition{Int64}[[1], [], []]

julia> wreath_macdonald_polynomials(1,3,cperm(1:3),[0,1,-1])[[3, 2, 1],[3, 2, 1]]
[1   q^2     q]
[1     t     q]
[1     t   t^2]
```
