```@meta
CurrentModule = Oscar
DocTestSetup = quote
  using Oscar
end
```

# Polycyclic groups

```@docs
PcGroup
PcGroupElem
```

Julia has the following functions that allow to generate polycyclic groups:
```@docs
abelian_group(::Type{T}, v::Vector{Int}) where T <: GAPGroup
cyclic_group
dihedral_group
quaternion_group
```
