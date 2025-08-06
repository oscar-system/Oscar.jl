######################
# 4: Properties
######################

@doc raw"""
    is_normal(v::NormalToricVarietyType)

Check if the normal toric variety `v` is normal. (This function is somewhat tautological at this point.)

# Examples
```jldoctest
julia> is_normal(projective_space(NormalToricVariety, 2))
true
```
"""
is_normal(v::NormalToricVarietyType) = true


@doc raw"""
    is_affine(v::NormalToricVarietyType)

Check if the normal toric variety `v` is affine.

# Examples
```jldoctest
julia> is_affine(projective_space(NormalToricVariety, 2))
false
```
"""
@attr Bool is_affine(v::NormalToricVarietyType) = pm_object(v).AFFINE


@doc raw"""
    is_projective(v::NormalToricVarietyType)

Check if the normal toric variety `v` is projective, i.e. if the fan of `v` is the the normal fan of a polytope.

# Examples
```jldoctest
julia> is_projective(projective_space(NormalToricVariety, 2))
true
```
"""
@attr Bool is_projective(v::NormalToricVarietyType) = pm_object(v).PROJECTIVE


@doc raw"""
    is_projective_space(v::NormalToricVarietyType)

Decides if the normal toric varieties `v` is a projective space.

# Examples
```jldoctest
julia> F5 = hirzebruch_surface(NormalToricVariety, 5)
Normal toric variety

julia> is_projective_space(F5)
false

julia> is_projective_space(projective_space(NormalToricVariety, 2))
true
```
"""
@attr Bool function is_projective_space(v::NormalToricVarietyType)
    if is_smooth(v) == false
        return false
    end
    if is_projective(v) == false
        return false
    end
    if torsion_free_rank(class_group(v)) > 1
        return false
    end
    w = [[Int(x) for x in transpose(g.coeff)] for g in gens(class_group(v))]
    for g in gens(class_group(v))
        g = [Int(x) for x in g.coeff if !iszero(x)]
        if length(g) > 1
            return false
        end
        if g[1] != 1
            return false
        end
    end
    return irrelevant_ideal(v) == ideal(gens(cox_ring(v)))
end


@doc raw"""
    is_smooth(v::NormalToricVarietyType)

Check if the normal toric variety `v` is smooth.

# Examples
```jldoctest
julia> is_smooth(projective_space(NormalToricVariety, 2))
true
```
"""
@attr Bool is_smooth(v::NormalToricVarietyType) = pm_object(v).SMOOTH


@doc raw"""
    is_complete(v::NormalToricVarietyType)

Check if the normal toric variety `v` is complete.

# Examples
```jldoctest
julia> is_complete(projective_space(NormalToricVariety, 2))
true
```
"""
@attr Bool is_complete(v::NormalToricVarietyType) = pm_object(v).COMPLETE


@doc raw"""
    has_torusfactor(v::NormalToricVarietyType)

Check if the normal toric variety `v` has a torus factor.

# Examples
```jldoctest
julia> has_torusfactor(projective_space(NormalToricVariety, 2))
false
```
"""
@attr Bool has_torusfactor(v::NormalToricVarietyType) = Polymake.common.rank(rays(v)) < ambient_dim(v)


@doc raw"""
    is_orbifold(v::NormalToricVarietyType)

Check if the normal toric variety `v` is an orbifold.

# Examples
```jldoctest
julia> is_orbifold(projective_space(NormalToricVariety, 2))
true
```
"""
@attr Bool is_orbifold(v::NormalToricVarietyType) = pm_object(v).SIMPLICIAL


@doc raw"""
    is_simplicial(v::NormalToricVarietyType)

Check if the normal toric variety `v` is simplicial. Hence, this function works just as `is_orbifold`. It is implemented for user convenience.

# Examples
```jldoctest
julia> is_simplicial(projective_space(NormalToricVariety, 2))
true
```
"""
is_simplicial(v::NormalToricVarietyType) = is_orbifold(v)


@doc raw"""
    is_gorenstein(v::NormalToricVarietyType)

Check if the normal toric variety `v` is Gorenstein.

# Examples
```jldoctest
julia> is_gorenstein(projective_space(NormalToricVariety, 2))
true
```
"""
@attr Bool is_gorenstein(v::NormalToricVarietyType) = pm_object(v).GORENSTEIN


@doc raw"""
    is_q_gorenstein(v::NormalToricVarietyType)

Check if the normal toric variety `v` is Q-Gorenstein.

# Examples
```jldoctest
julia> is_q_gorenstein(projective_space(NormalToricVariety, 2))
true
```
"""
@attr Bool is_q_gorenstein(v::NormalToricVarietyType) = pm_object(v).Q_GORENSTEIN


@doc raw"""
    is_fano(v::NormalToricVarietyType)

Check if the normal toric variety `v` is fano.

# Examples
```jldoctest
julia> is_fano(projective_space(NormalToricVariety, 2))
true
```
"""
@attr Bool function is_fano(v::NormalToricVarietyType)
  # [CLS11] Fano implies projective, thus non-projective implies non-Fano. But
  # non-complete implies non-projective.
  if !is_complete(v)
    return false
  else
    return is_ample(anticanonical_divisor(v))
  end
end
