@doc raw"""
    base_space(m::AbstractFTheoryModel)

Return the base space of the F-theory model.

# Examples
```jldoctest
julia> w = weierstrass_model_over_projective_space(2);

julia> dim(base_space(w))
2
```
"""
function base_space(m::AbstractFTheoryModel)
  is_base_space_fully_specified(m) || @vprint :FTheoryModelPrinter 1 "Base space was not fully specified. Returning AUXILIARY base space.\n"
  return m.base_space
end


@doc raw"""
    ambient_space(m::AbstractFTheoryModel)

Return the ambient space of the F-theory model.

# Examples
```jldoctest
julia> w = weierstrass_model_over_projective_space(2);

julia> dim(ambient_space(w))
4
```
"""
function ambient_space(m::AbstractFTheoryModel)
  is_base_space_fully_specified(m) || @vprint :FTheoryModelPrinter 1 "Base space was not fully specified. Returning AUXILIARY ambient space.\n"
  return m.ambient_space
end


@doc raw"""
    fiber_ambient_space(m::AbstractFTheoryModel)

Return the fiber ambient space of an F-theory model.

# Examples
```jldoctest
julia> w = weierstrass_model_over_projective_space(2);

julia> dim(fiber_ambient_space(w))
2
```
"""
function fiber_ambient_space(m::AbstractFTheoryModel)
  @req hasfield(typeof(m), :fiber_ambient_space) "fiber_ambient_space not supported for this F-theory model"
  return m.fiber_ambient_space
end


@doc raw"""
    model_index(m::AbstractFTheoryModel)
Return database index of a literature model. This index is a unique identifier that can be used to more conveniently construct the model. 
All models have a model_index and these will not change in the future.

# Examples
```jldoctest
julia> t = literature_model(31)
Weierstrass model over a not fully specified base -- F-theory weierstrass model dual to hypersurface model with fiber ambient space F_10 based on arXiv paper 1408.4808 Eq. (3.130)

julia> model_index(t)
31
```
"""
function model_index(m::AbstractFTheoryModel)
  directory = joinpath(dirname(dirname(@__DIR__)), "LiteratureModels/")
  model_indices = JSON.parsefile(directory * "model_indices.json")
  return parse(Int, model_indices["model" * literature_identifier(m) * ".json"])
end


@doc raw"""
    calabi_yau_hypersurface(m::AbstractFTheoryModel)

Return the Calabi–Yau hypersurface that defines the hypersurface model
as a closed subvariety of its toric ambient space.

# Examples
```jldoctest
julia> w =  weierstrass_model_over_projective_space(3)
Weierstrass model over a concrete base

julia> calabi_yau_hypersurface(w)
Closed subvariety of a normal toric variety

julia> t = global_tate_model_over_projective_space(2)
Global Tate model over a concrete base

julia> calabi_yau_hypersurface(t)
Closed subvariety of a normal toric variety

julia> b = projective_space(NormalToricVariety, 2);

julia> kb = anticanonical_divisor_class(b);

julia> fiber_ambient = weighted_projective_space(NormalToricVariety, [2, 3, 1]);

julia> set_coordinate_names(fiber_ambient, ["x", "y", "z"]);

julia> p = "x^3 - y^2 + x1^12 * x * z^4 + x2^18 * z^6 + 13 * x3^3*x*y*z";

julia> h = hypersurface_model(b, fiber_ambient, [2 * kb, 3 * kb], p; completeness_check=false)
Hypersurface model over a concrete base

julia> calabi_yau_hypersurface(h)
Closed subvariety of a normal toric variety
```
"""
@attr ClosedSubvarietyOfToricVariety function calabi_yau_hypersurface(m::AbstractFTheoryModel)
  @req (m isa HypersurfaceModel || m isa WeierstrassModel || m isa GlobalTateModel) "Calabi-Yau hypersurface currently not supported for this model"
  @req base_space(m) isa NormalToricVariety "Calabi-Yau hypersurface currently only supported for toric varieties as base space"
  is_base_space_fully_specified(m) || @vprint :FTheoryModelPrinter 1 "Base space was not fully specified. Returning hypersurface in AUXILIARY ambient space.\n"
  return closed_subvariety_of_toric_variety(ambient_space(m), [hypersurface_equation(m)])
end


######################################################################################
### Attributes for flux families (not exported, rather for serialization overhaul)
######################################################################################

@attr QQMatrix function matrix_integral_quant_transverse(m::AbstractFTheoryModel; check::Bool = true)
  return matrix_integral(special_flux_family(m, check = check))
end

@attr QQMatrix function matrix_rational_quant_transverse(m::AbstractFTheoryModel; check::Bool = true)
  return matrix_rational(special_flux_family(m, check = check))
end

@attr Vector{QQFieldElem} function offset_quant_transverse(m::AbstractFTheoryModel; check::Bool = true)
  return offset(special_flux_family(m, check = check))
end

@attr QQMatrix function matrix_integral_quant_transverse_nobreak(m::AbstractFTheoryModel; check::Bool = true)
  return matrix_integral(special_flux_family(m, not_breaking = true; check = check))
end

@attr QQMatrix function matrix_rational_quant_transverse_nobreak(m::AbstractFTheoryModel; check::Bool = true)
  return matrix_rational(special_flux_family(m, not_breaking = true; check = check))
end

@attr Vector{QQFieldElem} function offset_quant_transverse_nobreak(m::AbstractFTheoryModel; check::Bool = true)
  return offset(special_flux_family(m, not_breaking = true; check = check))
end
