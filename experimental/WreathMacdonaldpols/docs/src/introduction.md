```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Wreath Macdonald polynomials

The existence, integrality and positivity of wreath Macdonald polynomials
 has been conjectured by Haiman [Hai02](@cite) and proved by Bezrukavnikov
 and Finkelberg [BF14](@cite).
Here we have implemented an algorithm computing the wreath Macdonald
 polynomials as defined in the survey by Orr and Shimozono on this topic [OS23](@cite).

Take ``r`` and ``n`` two integers and consider the complex reflection group ``G(r,1,n):= (\mathbb{Z}/r\mathbb{Z})^n \rtimes \mathfrak{S}_n`` with its natural ``n`` dimensional reflection representation
``\mathfrak{h}_n``. For ``\lambda`` a partition, denote respectively by ``\mathrm{core}_r(\lambda)`` and ``\mathrm{quot}_r(\lambda)`` the ``r``-core and ``r``-quotient of ``\lambda``. Moreover, let us equip
the set of partitions with ``\leq``, the dominance order. Let ``\mathbb{K}`` denote the field ``\mathbb{Q}(q,t)``.
The simplest characterization (i.e. the one using minimal notation) of these polynomials is given in the introduction of [Wen19](@cite). We recall it now.\
For each partition ``\lambda`` such that ``|\mathrm{quot}_r(\lambda)|=n``, the wreath Macdonald polynomial ``H_{\lambda}`` is characterized by

* ``H_{\lambda} \otimes \sum_{i=0}^n(-q)^i\left[\Lambda^i\mathfrak{h}_n^*\right] \in \bigoplus_{\mu \geq \lambda, \mathrm{core}_r(\mu)=\mathrm{core}_r(\lambda)}{\mathbb{K}\left[V_{\mathrm{quot}(\mu)}\right]}``,
* ``H_{\lambda} \otimes \sum_{i=0}^n(-t)^{-i}\left[\Lambda^i\mathfrak{h}_n^*\right] \in \bigoplus_{\mu \leq \lambda, \mathrm{core}_r(\mu)=\mathrm{core}_r(\lambda)}{\mathbb{K}\left[V_{\mathrm{quot}(\mu)}\right]}``,
* ``\langle H_{\lambda},[\mathrm{triv}]\rangle = 1``.

Remark that when ``r=1``, wreath Macdonald polynomials are equal to the Haiman-Macdonald polynomials, used to prove the Macdonald positivity conjecture.

These polynomials, apart from generalizing the Haiman-Macdonald polynomials and giving access to new combinatorics, have a geometric counterpart.
Denote by ``\mathbb{T}`` the maximal diagonal torus of ``\mathrm{GL}_2(\mathbb{C})``.
To be more precise, Bezrukavnikov and Finkelberg prove that ``H_{\lambda}`` can be realized as the bigraded ``G(r,1,n)`` Frobenius character of the fiber of a wreath Procesi bundle (see [Los18](@cite)) at the ``\mathbb{T}``-fixed point associated with ``\mathrm{quot}_r(\lambda)``.

Finally, the wreath Macdonald polynomials can be interpreted as the eigenbasis of explicit vertex operators see [Wen19](@cite).

In our implementation, wreath Macdonald polynomials depend on two parameters. The first parameter is
 an ``r``-multipartition of ``n`` (the ``r``-quotient of ``\lambda``). The second parameter is an element of the affine Weyl group
 of type ``A^{(1)}_{r-1}`` which is isomorphic to the semi-direct product of the finite Weyl group
 of type ``A_{r-1}`` (the symmetric group on ``r`` letters) and of the coroot lattice of type ``A_{r-1}`` (the ``r``-core of ``\lambda``).
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
