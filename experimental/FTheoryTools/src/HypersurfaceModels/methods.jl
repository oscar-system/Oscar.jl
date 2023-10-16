@doc raw"""
    set_hypersurface_equation(h::HypersurfaceModel, p::MPolyRingElem)

Set the hypersurface equation to a custom value.

```jldoctest
julia> h = hypersurface_model_over_projective_space(2)
Hypersurface model over a concrete base

julia> R = parent(hypersurface_equation(h));

julia> new_poly = gens(R)[4]^3
x^3

julia> set_hypersurface_equation(h, new_poly);
```
"""
function set_hypersurface_equation(h::HypersurfaceModel, p::MPolyRingElem)
  @req parent(p) == parent(hypersurface_equation(h)) "Polynomial must reside in the Cox ring of the ambient space"
  @req degree(p) == degree(hypersurface_equation(h)) "Degree mismatch between specified polynomial and current hypersurface equation"
  h.hypersurface_equation = p
end


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
