using GAP
using CapAndHomalg


######################
# 1: The Julia type for ToricVarieties
######################

struct JToricVariety
           GapToricVariety
end
export JToricVariety


######################
# 2: Generic constructors
######################

function JToricVariety( rays::Vector{Vector{Int64}}, cones::Vector{Vector{Int64}} )
    
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
    return JToricVariety( variety )

end
export JToricVariety


######################
# 3: Special constructors
######################

function ProjectiveSpace( x )
    
    # load necessary gap packages
    CapAndHomalg.LoadPackage( "JuliaInterface" )
    CapAndHomalg.LoadPackage( "JConvex" )
    CapAndHomalg.LoadPackage( "ToricV" )

    # construct the projective space in gap
    variety = GAP.Globals.ProjectiveSpace( x )
    
    # wrap it and return
    return JToricVariety(  variety )
    
end
export ProjectiveSpace


######################
# 4: Properties
######################

function IsNormalVariety( v::JToricVariety )
    
    return Bool( GAP.Globals.IsNormalVariety( v.GapToricVariety ) )
    
end
export IsNormalVariety


function IsAffine( v::JToricVariety )
    
    return Bool( GAP.Globals.IsAffine( v.GapToricVariety ) )
    
end
export IsAffine


function IsProjective( v::JToricVariety )
    
    return Bool( GAP.Globals.IsProjective( v.GapToricVariety ) )
    
end
export IsProjective


function IsSmooth( v::JToricVariety )
    
    return Bool( GAP.Globals.IsSmooth( v.GapToricVariety ) )
    
end
export IsSmooth


function IsComplete( v::JToricVariety )
    
    return Bool( GAP.Globals.IsComplete( v.GapToricVariety ) )
    
end
export IsComplete


function HasTorusfactor( v::JToricVariety )
    
    return Bool( GAP.Globals.HasTorusfactor( v.GapToricVariety ) )
    
end
export HasTorusfactor


function HasNoTorusfactor( v::JToricVariety )
    
    return Bool( GAP.Globals.HasNoTorusfactor( v.GapToricVariety ) )
    
end
export HasNoTorusfactor


function IsOrbifold( v::JToricVariety )
    
    return Bool( GAP.Globals.IsOrbifold( v.GapToricVariety ) )
    
end
export IsOrbifold


function IsSimplicial( v::JToricVariety )
    
    return Bool( GAP.Globals.IsSimplicial( v.GapToricVariety ) )
    
end
export IsSimplicial


function IsIsomorphicToProjectiveSpace( v::JToricVariety )
    
    return Bool( GAP.Globals.IsIsomorphicToProjectiveSpace( v.GapToricVariety ) )
    
end
export IsIsomorphicToProjectiveSpace


function IsDirectProductOfPNs( v::JToricVariety )
    
    return Bool( GAP.Globals.IsDirectProductOfPNs( v.GapToricVariety ) )
    
end
export IsDirectProductOfPNs
