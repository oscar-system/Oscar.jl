#####################################################
# 1: Setters
#####################################################

@doc raw"""
    set_weierstrass_model(h::HypersurfaceModel, w::WeierstrassModel)

Allows to define the Weierstrass model corresponding to the hypersurface model.
"""
function set_weierstrass_model(h::HypersurfaceModel, w::WeierstrassModel)
  set_attribute!(h, :weierstrass_model => w)
end

@doc raw"""
    set_global_tate_model(h::HypersurfaceModel, w::GlobalTateModel)

Allows to define the global Tate model corresponding to the hypersurface model.
"""
function set_global_tate_model(h::HypersurfaceModel, t::GlobalTateModel)
  set_attribute!(h, :global_tate_model => t)
end



##############################################################
# 2: Tune a hypersurface model to custom hypersurface equation
##############################################################

@doc raw"""
    tune(h::HypersurfaceModel, p::MPolyRingElem; completeness_check::Bool = true)

Tune a hypersurface model, by specifying a custom hypersurface equation.

# Examples
```jldoctest
julia> base = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> h = hypersurface_model(base; completeness_check = false)
Hypersurface model over a concrete base

julia> Kbar = anticanonical_bundle(ambient_space(h))
Toric line bundle on a normal toric variety

julia> special_hypersurface_equation = basis_of_global_sections(Kbar)[1]
y^2

julia> h2 = tune(h, special_hypersurface_equation)
Hypersurface model over a concrete base

julia> hypersurface_equation(h2) == special_hypersurface_equation
true
```
"""
function tune(h::HypersurfaceModel, p::MPolyRingElem; completeness_check::Bool = true)
  @req !(typeof(base_space(h)) <: FamilyOfSpaces) "Currently, tuning is only possible for models over concrete toric bases"
  @req parent(p) == parent(hypersurface_equation(h)) "Parent mismatch between specified polynomial and current hypersurface equation"
  @req degree(p) == degree(hypersurface_equation(h)) "Degree mismatch between specified polynomial and current hypersurface equation"
  p == hypersurface_equation(h) && return h
  tuned_model = HypersurfaceModel(base_space(h), ambient_space(h), fiber_ambient_space(h), p)
  set_attribute!(tuned_model, :partially_resolved, false)
  return tuned_model
end
