```@meta
CurrentModule = Oscar
DocTestSetup = quote
  using Oscar
end
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["groups.md"]
```

# Groups

Oscar supports the following types of groups:

* `PermGroup` = groups of permutations
* `MatrixGroup` = groups of matrices
* `FPGroup` = finitely presented groups
* `PcGroup` = polycyclic groups
* `DirectProductOfGroups` = direct product of two groups
* `AutomorphismGroup` = group of automorphisms over a group


If `x` is an element of the group `G` of type `T`, then the type of `x` is `GAPGroupElement{T}`.
