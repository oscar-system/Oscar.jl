```@meta
CurrentModule = Oscar.OrthogonalDiscriminants
DocTestSetup = quote
  using Oscar
end
```

# Criteria for computing orthogonal discriminants

## Direct methods

```@docs
od_from_atlas_group
```

## Character-theoretical criteria

```@docs
od_from_order
od_from_eigenvalues
od_for_specht_module
od_from_p_subgroup(chi::GAPGroupClassFunction, p::Int)
```
