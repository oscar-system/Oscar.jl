###################################################################
###################################################################
# 1: Attributes that work the same tor toric and non-toric settings
###################################################################
###################################################################


#####################################################
# 1.1 Hypersurface equation
#####################################################

@doc raw"""
    hypersurface_equation(h::HypersurfaceModel)

Return the hypersurface equation.

```jldoctest
julia> h = hypersurface_model_over_projective_space(2)
Hypersurface model over a concrete base

julia> hypersurface_equation(h);
```
"""
hypersurface_equation(h::HypersurfaceModel) = h.hypersurface_equation


#####################################################
# 1.2 Base, ambient space and fiber ambient space
#####################################################

@doc raw"""
    base_space(h::HypersurfaceModel)

Return the base space of the hypersurface model.

```jldoctest
julia> h = hypersurface_model_over_projective_space(2)
Hypersurface model over a concrete base

julia> base_space(h)
Scheme of a toric variety with fan spanned by RayVector{QQFieldElem}[[1, 0], [0, 1], [-1, -1]]
```
"""
function base_space(h::HypersurfaceModel)
  base_fully_specified(h) || @vprint :HypersurfaceModel 1 "Base space was not fully specified. Returning AUXILIARY base space.\n"
  return h.base_space
end


@doc raw"""
    ambient_space(h::HypersurfaceModel)

Return the ambient space of the hypersurface model.

```jldoctest
julia> h = hypersurface_model_over_projective_space(2)
Hypersurface model over a concrete base

julia> ambient_space(h)
Scheme of a toric variety with fan spanned by RayVector{QQFieldElem}[[1, 0, 0, 3], [0, 1, 0, 0], [-1, -1, 0, 0], [0, 0, -1, 1//3], [0, 0, 1, -1//2], [0, 0, 0, 1]]
```
"""
function ambient_space(h::HypersurfaceModel)
  base_fully_specified(h) || @vprint :HypersurfaceModel 1 "Base space was not fully specified. Returning AUXILIARY ambient space.\n"
  return h.ambient_space
end


@doc raw"""
    fiber_ambient_space(HypersurfaceModel)

Return the fiber ambient space of the hypersurface model.

```jldoctest
julia> h = hypersurface_model_over_projective_space(2)
Hypersurface model over a concrete base

julia> fiber_ambient_space(h)
Scheme of a toric variety with fan spanned by RayVector{QQFieldElem}[[-1, 1//3], [1, -1//2], [0, 1]]
```
"""
fiber_ambient_space(h::HypersurfaceModel) = h.fiber_ambient_space
