######################
# 1: The Julia type for ToricVarieties
######################

struct toric_variety
           GapToricVariety::GapObj
end
export toric_variety


######################
# 2: Generic constructors
######################

function toric_variety( rays::Vector{Vector{Int}}, cones::Vector{Vector{Int}} )
    # construct the toric variety in GAP
    gap_rays = GapObj( rays, recursive = true )
    gap_cones = GapObj( cones, recursive = true )
    fan = GAP.Globals.Fan( gap_rays, gap_cones )
    variety = GAP.Globals.ToricVariety( fan )
    
    # wrap it into a struct and return
    return toric_variety( variety )
end
export toric_variety


######################
# 3: Special constructors
######################

function projective_space( x )
    # construct the projective space in gap
    variety = GAP.Globals.ProjectiveSpace( x )
    
    # wrap it and return
    return toric_variety(  variety )
end
export projective_space


######################
# 4: Properties
######################

function is_normal_variety( v::toric_variety )
    return true
end
export is_normal_variety


function is_affine( v::toric_variety )
    return GAP.Globals.IsAffine(v.GapToricVariety)::Bool
end
export is_affine


function is_projective( v::toric_variety )
    return GAP.Globals.IsProjective( v.GapToricVariety )::Bool
end
export is_projective


function is_smooth( v::toric_variety )
    return GAP.Globals.IsSmooth( v.GapToricVariety )::Bool
end
export is_smooth


function is_complete( v::toric_variety )
    return GAP.Globals.IsComplete( v.GapToricVariety )::Bool
end
export is_complete


function has_torusfactor( v::toric_variety )
    return GAP.Globals.HasTorusfactor( v.GapToricVariety )::Bool
end
export has_torusfactor


function has_no_torusfactor( v::toric_variety )
    return GAP.Globals.HasNoTorusfactor( v.GapToricVariety )::Bool
end
export has_no_torusfactor


function is_orbifold( v::toric_variety )
    return GAP.Globals.IsOrbifold( v.GapToricVariety )::Bool
end
export is_orbifold


function is_simplicial( v::toric_variety )
    return GAP.Globals.IsSimplicial( v.GapToricVariety )::Bool
end
export is_simplicial


function is_isomorphic_to_projective_space( v::toric_variety )
    return GAP.Globals.IsIsomorphicToProjectiveSpace( v.GapToricVariety )::Bool
end
export is_isomorphic_to_projective_space


function is_direct_product_of_projective_spaces( v::toric_variety )
    return GAP.Globals.IsDirectProductOfPNs( v.GapToricVariety )::Bool
end
export is_direct_product_of_projective_spaces
