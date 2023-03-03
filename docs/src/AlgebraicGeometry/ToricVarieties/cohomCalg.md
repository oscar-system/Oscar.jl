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

For a toric variety, all vanishing sets are computed as follows:
```@docs
vanishing_sets(variety::AbstractNormalToricVariety)
```
The return value is a vector of vanishing sets. This vector has length one
larger than the dimension of the variety in question. The first vanishing
set in this vector describes all line bundles for which the zero-th sheaf
cohomology vanishes. More generally, if a line bundle is contained in the
`n`-th vanishing set, then its `n-1`-th sheaf cohomology vanishes. The
following method checks if a line bundle is contained in a vanishing set:
```@docs
contains(tvs::ToricVanishingSet, l::ToricLineBundle)
```
A vanishing set can in principle cover the entire Picard group. This
can be checked with `isfull`. This methods returns `true` if the
vanishing set is the entire Picard group and `false` otherwise.
Beyond this, we support the following attributes for vanishing sets:
```@docs
toric_variety(tvs::ToricVanishingSet)
polyhedra(tvs::ToricVanishingSet)
cohomology_index(tvs::ToricVanishingSet)
```
