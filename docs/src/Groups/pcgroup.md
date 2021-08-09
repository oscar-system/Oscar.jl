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
Pages = ["pcgroup.md"]
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

The generators of a polycyclic group are displayed as `f1`, `f2`, `f3`, etc., and every element of a polycyclic group is displayed as product of such generators.

  **Example:**
```jldoctest
julia> G=abelian_group(PcGroup, [2,4]);

julia> G[1], G[2]
(f1, f2)

julia> G[2]*G[1]
f1*f2
```

Note that this does not define Julia variables named `f1`, `f2`, etc.! To get the generators of the group `G`, use `gens(G)`; for convenience they can also be accessed as `G[1]`, `G[2]`, as shown in Section [Elements of groups](@ref elements_of_groups).
