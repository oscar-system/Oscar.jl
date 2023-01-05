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
Pages = ["fpgroup.md"]
```

# Finitely presented groups

```@docs
FPGroup
FPGroupElem
free_group(n::Int)
is_full_fp_group(G::FPGroup)
relators(G::FPGroup)
length(g::FPGroupElem)
map_word
```
