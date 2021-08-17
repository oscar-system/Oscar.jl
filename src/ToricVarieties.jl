using GAP
using CapAndHomalg


######################
# 1: The Julia type for ToricVarieties
######################

struct toric_variety
           GapToricVariety
end
export toric_variety


######################
# 2: Generic constructors
######################

function toric_variety( rays::Vector{Vector{Int64}}, cones::Vector{Vector{Int64}} )
    # load necessary gap packages
    CapAndHomalg.LoadPackage( "JuliaInterface" )
    CapAndHomalg.LoadPackage( "JConvex" )
    CapAndHomalg.LoadPackage( "ToricV" )
    
    # construct the toric variety in GAP
    gap_rays = GAP.Globals.ConvertJuliaToGAP( rays )
    gap_cones = GAP.Globals.ConvertJuliaToGAP( cones )
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
    # load necessary gap packages
    CapAndHomalg.LoadPackage( "JuliaInterface" )
    CapAndHomalg.LoadPackage( "JConvex" )
    CapAndHomalg.LoadPackage( "ToricV" )

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
    return Bool( GAP.Globals.IsAffine( v.GapToricVariety ) )
end
export is_affine


function is_projective( v::toric_variety )
    return Bool( GAP.Globals.IsProjective( v.GapToricVariety ) )
end
export is_projective


function is_smooth( v::toric_variety )
    return Bool( GAP.Globals.IsSmooth( v.GapToricVariety ) )
end
export is_smooth


function is_complete( v::toric_variety )
    return Bool( GAP.Globals.IsComplete( v.GapToricVariety ) )
end
export is_complete


function has_torusfactor( v::toric_variety )
    return Bool( GAP.Globals.HasTorusfactor( v.GapToricVariety ) )
end
export has_torusfactor


function has_no_torusfactor( v::toric_variety )
    return Bool( GAP.Globals.HasNoTorusfactor( v.GapToricVariety ) )
end
export has_no_torusfactor


function is_orbifold( v::toric_variety )
    return Bool( GAP.Globals.IsOrbifold( v.GapToricVariety ) )
end
export is_orbifold


function is_simplicial( v::toric_variety )
    return Bool( GAP.Globals.IsSimplicial( v.GapToricVariety ) )
end
export is_simplicial


function is_isomorphic_to_projective_space( v::toric_variety )
    return Bool( GAP.Globals.IsIsomorphicToProjectiveSpace( v.GapToricVariety ) )
end
export is_isomorphic_to_projective_space


function is_direct_product_of_projective_spaces( v::toric_variety )
    return Bool( GAP.Globals.IsDirectProductOfPNs( v.GapToricVariety ) )
end
export is_direct_product_of_projective_spaces
