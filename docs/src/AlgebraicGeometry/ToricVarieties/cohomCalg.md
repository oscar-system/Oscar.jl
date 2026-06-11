```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Sheaf Cohomology of Toric Line Bundles

OSCAR provides functionality to compute the dimensions of sheaf cohomology
groups of toric line bundles on normal toric varieties.

The computation can be carried out using three different algorithms,
selected via a keyword argument.

- `:cohomcalg` — use the cohomCalg algorithm (cf. [BJRR10](@cite), [BJRR10*1](@cite)).
  Requires the toric variety to be simplicial and projective.

- `:chamber` — chamber counting algorithm (cf. [CLS11](@cite), p.398).
  Requires the toric variety to be simplicial and complete.

- `:local` — local cohomology method (cf. [CLS11](@cite), Section 9.5).

!!! note "Applicability of cohomCalg"
    In theory, the cohomCalg algorithm applies to all smooth complete and
    simplicial projective toric varieties [RR10](@cite), [Jow11](@cite).
    However, the implementation in [BJRR10*1](@cite) is limited to the
    simplicial, projective case due to its handling of secondary cohomologies.
    Consequently, all functions relying on cohomCalg are currently restricted
    to **simplicial, projective toric varieties**.

## Dimensions of line bundle cohomology

The following methods compute the dimensions of the sheaf cohomology
groups of a toric line bundle.
- `sheaf_cohomology(l; algorithm=...)` returns all cohomology dimensions.
- `sheaf_cohomology(l, i; algorithm=...)` returns the dimension in degree `i`.

```@docs
sheaf_cohomology(l::ToricLineBundle; algorithm::Symbol=:cohomcalg)
sheaf_cohomology(l::ToricLineBundle, i::Int; algorithm::Symbol=:cohomcalg)
```

## Toric vanishing sets

Vanishing sets describe subsets of the Picard group of a fixed toric variety
for which a given sheaf cohomology group vanishes.

Their computation is based on [BJRR10*1](@cite). Consequently, as stated above, this
functionality is currently available only for toric varieties that are simplicial
and projective.

This approach to identifying vanishing sets on toric varieties was originally
introduced in [Bie18](@cite). Technically, a vanishing set is the complement
of a finite family of polyhedra.

### Computing vanishing sets

For a toric variety, all vanishing sets are computed as follows:

```@docs
vanishing_sets(variety::NormalToricVarietyType)
```

The return value is a vector of vanishing sets of length `dim(X) + 1`.

The first vanishing set describes all line bundles for which the zeroth
sheaf cohomology vanishes. More generally, if a line bundle is contained
in the `n`-th vanishing set, then its `(n−1)`-st cohomology group vanishes.

Membership in a vanishing set can be tested with:

```@docs
contains(tvs::ToricVanishingSet, l::ToricLineBundle)
```

A vanishing set may coincide with the entire Picard group. This can be
checked with `isfull`, which returns `true` if the vanishing set equals
the entire Picard group and `false` otherwise.

Additional information associated with a vanishing set is available via:

```@docs
toric_variety(tvs::ToricVanishingSet)
polyhedra(tvs::ToricVanishingSet)
cohomology_indices(tvs::ToricVanishingSet)
```

Vanishing sets can also be used to compute immaculate line bundles:

```@docs
immaculate_line_bundles(variety::NormalToricVarietyType)
```
