@doc raw"""
    set_weierstrass_model(h::HypersurfaceModel, w::WeierstrassModel)

Assigns a Weierstrass model to the given hypersurface model.

See [`weierstrass_model(h::HypersurfaceModel)`](@ref) for an example use.
"""
function set_weierstrass_model(h::HypersurfaceModel, w::WeierstrassModel)
  set_attribute!(h, :weierstrass_model => w)
end

@doc raw"""
    set_global_tate_model(h::HypersurfaceModel, w::GlobalTateModel)

Assigns a global Tate model to the given hypersurface model.

See [`global_tate_model(h::HypersurfaceModel)`](@ref) for an example use.
"""
function set_global_tate_model(h::HypersurfaceModel, t::GlobalTateModel)
  set_attribute!(h, :global_tate_model => t)
end
