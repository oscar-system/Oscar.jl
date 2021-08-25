######################
# 1: The Julia type for ToricVarieties
######################

struct NormalToricVariety
           GapNTV::GapObj
           polymakeNTV::Polymake.BigObject
end
export NormalToricVariety


######################
# 2: Generic constructors
######################

function NormalToricVariety( rays::Matrix{Int}, cones::Vector{Vector{Int}} )
    # construct the toric variety in GAP
    gap_rays = GapObj( rays, recursive = true )
    gap_cones = GapObj( cones, recursive = true )
    fan = GAP.Globals.Fan( gap_rays, gap_cones )
    variety = GAP.Globals.ToricVariety( fan )

    Incidence = Oscar.IncidenceMatrix(cones)
    arr = @Polymake.convert_to Array{Set{Int}} Polymake.common.rows(Incidence.pm_incidencematrix)

    pmntv = Polymake.fulton.NormalToricVariety(
        INPUT_RAYS = Oscar.matrix_for_polymake(rays),
        INPUT_CONES = arr,
    )

    # wrap it into a struct and return
    return NormalToricVariety( variety, pmntv )
end
export NormalToricVariety


######################
# 3: Special constructors
######################

function projective_space( x )
    # construct the projective space in gap
    variety = GAP.Globals.ProjectiveSpace( x )
    
    # wrap it and return
    return NormalToricVariety(  variety )
end
export projective_space


######################
# 4: Properties
######################

function is_normal_variety( v::NormalToricVariety )
    return true
end
export is_normal_variety


function is_affine( v::NormalToricVariety )
    return v.polymakeNTV.AFFINE
end
export is_affine


function is_projective( v::NormalToricVariety )
    return v.polymakeNTV.PROJECTIVE
end
export is_projective


function is_smooth( v::NormalToricVariety )
    return v.polymakeNTV.SMOOTH
end
export is_smooth


function is_complete( v::NormalToricVariety )
    return GAP.Globals.IsComplete( v.GapNTV )::Bool
end
export is_complete


function has_torusfactor( v::NormalToricVariety )
    return GAP.Globals.HasTorusfactor( v.GapNTV )::Bool
end
export has_torusfactor


function is_orbifold( v::NormalToricVariety )
    return GAP.Globals.IsOrbifold( v.GapNTV )::Bool
end
export is_orbifold


function is_simplicial( v::NormalToricVariety )
    return GAP.Globals.IsSimplicial( v.GapNTV )::Bool
end
export is_simplicial


function is_isomorphic_to_projective_space( v::NormalToricVariety )
    return GAP.Globals.IsIsomorphicToProjectiveSpace( v.GapNTV )::Bool
end
export is_isomorphic_to_projective_space


function is_direct_product_of_projective_spaces( v::NormalToricVariety )
    return GAP.Globals.IsDirectProductOfPNs( v.GapNTV )::Bool
end
export is_direct_product_of_projective_spaces
