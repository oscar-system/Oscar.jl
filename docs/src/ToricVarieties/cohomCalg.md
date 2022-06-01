```@meta
CurrentModule = Oscar
```

```@contents
Pages = ["CohomCalg.md"]
```


# Line bundle cohomology with cohomCalg

We employ the cohomCalg algorithm [cohomCalg:Implementation](@cite)
to compute the dimension of line bundle cohomologies as well as vanishing sets.


## Dimensions of line bundle cohomology

```@docs
all_cohomologies(l::ToricLineBundle)
cohomology(l::ToricLineBundle, i::Int)
```

## Toric vanishing sets

Vanishing sets describe subsets of the Picard group of toric varieties.
Their computations is based on [cohomCalg:Implementation](@cite), i.e.
this functionality is only available if the toric variety in question is
either smooth and complete or alternatively, simplicial and projective.
This approach to identify vanishing sets on toric varieties was originally
introduced in [Bie18](@cite). As described there, on a technical level,
a vanishing set is the complement of a finite family of polyhedra.

For a given toric variety `tvs`, the vanishing sets are computed as follows:
```@docs
vanishing_sets(variety::AbstractNormalToricVariety)
```
The return value is a vector of vanishing sets. This vector has length
`dim(variety) + 1`. The first vanishing set in this vector describes all line
bundles for which the zero-th cohomology class vanishes. More generally,
if a line bundle is contained in the k-th vanishing set, then its `k-1`th
cohomology class vanishes. The following method can be used to check if a line
bundle `l` is contained in a vanishing sets:
toric vanishing set:
```@docs
contains(tvs::ToricVanishingSet, l::ToricLineBundle)
```

For each vanishing set, `tvs` we support the property `isfull(tvs)`.
It returns `true` if the vanishing set covers the entire Picard group.
of the toric variety in question. Otherwise, it returns `false`.

In addition, we support the following attributes:
```@docs
toric_variety(tvs::ToricVanishingSet)
polyhedra(tvs::ToricVanishingSet)
cohomology_index(tvs::ToricVanishingSet)
```