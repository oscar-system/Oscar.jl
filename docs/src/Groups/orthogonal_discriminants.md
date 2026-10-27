
# Orthogonal Discriminants

This section describes Oscar code that was used for computing the
orthogonal discriminants of the absolutely irreducible
orthogonal representations of even degree that are listed in the
Atlas of Finite Groups [CCNPW85](@cite) and in the
Atlas of Brauer Characters [JLPW95](@cite).

The theoretical background is described in the chapter
*An atlas of orthogonal representations* [BNP25](@cite) of
the Oscar book [DEFHJ25](@cite).

## Access to precomputed orthogonal discriminants data

```@docs
all_od_infos
orthogonal_discriminants
orthogonal_discriminant
```

## Criteria for computing orthogonal discriminants

Direct methods compute the orthogonal discriminant from a matrix group.

```@docs
Oscar.od_from_atlas_group
```

Character-theoretical criteria do not require the representation itself,
they need only its character.

```@docs
Oscar.od_from_order
Oscar.od_from_eigenvalues
Oscar.od_for_specht_module
Oscar.od_from_p_subgroup(chi::GAPGroupClassFunction, p::Int)
```

## Miscellaneous functions related to orthogonal discriminants

```@docs
is_orthogonally_stable
show_with_ODs
show_OD_info
```

The following functions deal with the irreducible representations
of symmetric groups, which can be described combinatorially via so-called
Specht modules.
The irreducible representations of the symmetric group on $n$ points
in characteristic zero are parametrized by the partitions of $n$,
and a lot of information about these representations can be computed
directly from the corresponding partitions.

```@docs
dimension_specht_module
gram_determinant_specht_module
```
