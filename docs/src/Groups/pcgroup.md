```@meta
CurrentModule = Oscar
DocTestSetup = Oscar.doctestsetup()
```

# Polycyclic groups

```@docs
PcGroup
PcGroupElem
map_word(g::Union{PcGroupElem, SubPcGroupElem}, genimgs::Vector; genimgs_inv::Vector = Vector(undef, length(genimgs)), init = nothing)
```

Julia has the following functions that allow to generate polycyclic groups:
```@docs
abelian_group(::Type{T}, v::Vector{Int}) where T <: GAPGroup
elementary_abelian_group
cyclic_group
dihedral_group
quaternion_group
```
