```@meta
CurrentModule = Oscar
```

# Blowups of toric varieties

## Constructors

```@docs
blow_up(v::NormalToricVarietyType, I::MPolyIdeal; coordinate_name::String = "e", set_attributes::Bool = true)
blow_up(v::NormalToricVarietyType, new_ray::AbstractVector{<:IntegerUnion}; coordinate_name::String = "e", set_attributes::Bool = true)
blow_up(v::NormalToricVarietyType, n::Int; coordinate_name::String = "e", set_attributes::Bool = true)
```

