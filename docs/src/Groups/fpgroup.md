```@meta
CurrentModule = Oscar
DocTestSetup = Oscar.doctestsetup()
```

# Finitely presented groups

```@docs
FPGroup
FPGroupElem
SubFPGroup
SubFPGroupElem
free_group
@free_group
full_group(G::Union{SubFPGroup, SubPcGroup})
relators(G::FPGroup)
length(g::Union{FPGroupElem, SubFPGroupElem})
map_word(g::Union{FPGroupElem, SubFPGroupElem}, genimgs::Vector; genimgs_inv::Vector = Vector(undef, length(genimgs)), init = nothing)
```
