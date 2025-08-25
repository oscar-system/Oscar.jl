```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Wreath Macdonald polynomials

The existence, integrality and positivity of wreath Macdonald polynomials
 has been conjectured by Haiman [Hai02](@cite) and proved by Bezrukavnikov and Finkelberg [BF14](@cite).
Here we have implemented an algorithm computing the wreath Macdonald
 polynomials as defined in the survey by Orr and Shimozono on this topic [OS23](@cite).

Take ``r`` and ``n`` two integers and consider the complex reflection group ``G(r,1,n):= (\mathbb{Z}/r\mathbb{Z})^n \rtimes \mathfrak{S}_n`` with its natural ``n``-dimensional reflection representation
``\mathfrak{h}_n``. Let ``\Lambda^i\mathfrak{h}_n^*`` denote the exterior power of ``\mathfrak{h}_n^*``. The irreducible representations over ``\mathbb{C}`` of ``G(r,1,n)`` are indexed by the ``r``-multipartitions of size ``n``. For ``\lambda^{\bullet}`` such a multipartition, let us denote by ``V_{\lambda^{\bullet}}`` the associated irreducible representation. If ``V`` is a representation of ``G(r,1,n)``, let ``\left[V\right]`` denote the class of ``V`` in ``K_{G(r,1,n)}``, the Grothendieck ring of ``G(r,1,n)``.\

We denote by ``\mathbb{K}`` the field ``\mathbb{Q}(q,t)``, where ``q`` and ``t`` are indeterminates over ``\mathbb{Q}``. Let ``R`` denote the ring of symmetric functions over ``\mathbb{K}``, see ``\S I.2`` in [Mac15](@cite). The ring of ``r``-multisymmetric functions is the ``r``-fold tensor product of ``R`` over ``\mathbb{K}``. If one equips ``\bigoplus_{n \in \mathbb{Z}_{\geq 0}}{K_{G(r,1,n)}}`` with the induction product, then in light of Theorem ``2.3`` in [Wen19](@cite), we identify the ring of ``r``-multisymmetric functions and ``\bigoplus_{n \in \mathbb{Z}_{\geq 0}}{K_{G(r,1,n)}}``.\

For ``\lambda`` a partition, denote respectively by ``\mathrm{core}_r(\lambda)`` and ``\mathrm{quot}_r(\lambda)`` the ``r``-core and ``r``-quotient of ``\lambda``. Moreover, let us equip the set of partitions with ``\leq``, the dominance order. For an element ``\omega`` of ``\mathfrak{S}_r`` and a partition ``\lambda``, define ``\omega.\lambda`` to be the partition with the same ``r``-core as ``\lambda`` and ``r``-quotient equal to ``\omega.\mathrm{quot}_r(\lambda)`` where ``\omega`` acts by permuting the ``r`` partitions. We define the order ``\leq_{\omega}`` as follows. If ``\lambda`` and ``\mu`` are two partitions, then ``\lambda \leq_{\omega} \mu`` if ``\omega.\lambda \leq \omega.\mu``.
We now give a simple characterization of the wreath Macdonald polynomials.\
\
For each partition ``\lambda`` such that ``|\mathrm{quot}_r(\lambda)|=n`` and each ``\omega \in \mathfrak{S}_r``, the wreath Macdonald polynomial ``H^{\omega}_{\lambda}`` is the ``r``-multisymmetric function uniquely characterized by

* ``H^{\omega}_{\lambda} \otimes \sum_{i=0}^n(-q)^i\left[\Lambda^i\mathfrak{h}_n^*\right] \in \bigoplus_{\mu \geq_{\omega} \lambda, \mathrm{core}_r(\mu)=\mathrm{core}_r(\lambda)}{\mathbb{K}\left[V_{\mathrm{quot}(\mu)}\right]}``,
* ``H^{\omega}_{\lambda} \otimes \sum_{i=0}^n(-t)^{-i}\left[\Lambda^i\mathfrak{h}_n^*\right] \in \bigoplus_{\mu \leq_{\omega} \lambda, \mathrm{core}_r(\mu)=\mathrm{core}_r(\lambda)}{\mathbb{K}\left[V_{\mathrm{quot}(\mu)}\right]}``,
* ``\langle H^{\omega}_{\lambda},[\mathrm{triv}]\rangle = 1``.

Remark that when ``r=1``, wreath Macdonald polynomials are equal to the Haiman-Macdonald polynomials, used to prove the Macdonald positivity conjecture.

These polynomials, apart from generalizing the Haiman-Macdonald polynomials and giving access to new combinatorics, have a geometric counterpart.
Denote by ``\mathbb{T}`` the maximal diagonal torus of ``\mathrm{GL}_2(\mathbb{C})``.
To be more precise, Bezrukavnikov and Finkelberg prove that ``H^{\omega}_{\lambda}`` can be realized as the bigraded ``G(r,1,n)`` Frobenius character of the fiber of a wreath Procesi bundle (see [Los18](@cite)) at the ``\mathbb{T}``-fixed point associated with ``\mathrm{quot}_r(\lambda)``.

Finally, the wreath Macdonald polynomials can be interpreted as the eigenbasis of explicit vertex operators, see [Wen19](@cite).

In our implementation, wreath Macdonald polynomials depend on two parameters. The first parameter is
 an ``r``-multipartition of ``n`` (the ``r``-quotient of ``\lambda``). The second parameter is an element of the affine Weyl group
 of type ``A^{(1)}_{r-1}`` which is isomorphic to the semi-direct product of the finite Weyl group of type ``A_{r-1}`` (the symmetric group on ``r`` letters) and of the coroot lattice of type ``A_{r-1}`` (the ``r``-core of ``\lambda``). The finite Weyl group acts by permutation. The element of the coroot lattice is given in the canonical basis. It is then the sublattice
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

# Contact

Please direct questions about this part of OSCAR to:
[Raphaël Paegelow](https://paegelow.fr/en/).
