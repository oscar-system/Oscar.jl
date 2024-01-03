#####################################################
# Setters
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
