```@meta
CurrentModule = Oscar
DocTestSetup = quote
  using Oscar
end
```

# Finitely presented groups

```@docs
FPGroup
FPGroupElem
free_group(n::Int)
is_full_fp_group(G::FPGroup)
relators(G::FPGroup)
length(g::FPGroupElem)
map_word(g::FPGroupElem, genimgs::Vector; genimgs_inv::Vector = Vector(undef, length(genimgs)), init = nothing)
```
