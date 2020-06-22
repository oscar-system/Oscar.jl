```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["groups.md"]
```

# Groups


Julia supports the following types of groups:

* `PermGroup` = groups of permutations
* `MatrixGroup` = groups of matrices
* `FPGroup` = finitely presented groups
* `PcGroup` = polycyclic groups
* `DirectProductOfGroups` = direct product of two groups
* `AutomorphismGroup` = group of automorphisms over a group

## Permutations groups

In Julia, permutations are displayed as products of disjoint cycles, as in GAP.
