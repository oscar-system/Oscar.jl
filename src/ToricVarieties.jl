using GAP
using CapAndHomalg


######################
# 1: The Julia type for ToricVarieties
######################

struct JToricVariety
           bar
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
    gap_rays = CapAndHomalg.GAP.Globals.ConvertJuliaToGAP( rays )
    gap_cones = CapAndHomalg.GAP.Globals.ConvertJuliaToGAP( cones )
    fan = CapAndHomalg.GAP.Globals.Fan( gap_rays, gap_cones )
    variety = CapAndHomalg.GAP.Globals.ToricVariety( fan )
    
    # wrap it into a struct and return
    return JToricVariety( 1, variety )

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
    variety = CapAndHomalg.GAP.Globals.ProjectiveSpace( x )
    
    # wrap it and return
    return JToricVariety(  1, variety )
    
end
export ProjectiveSpace


######################
# 4: Properties
######################

function IsNormalVariety( v::JToricVariety )
    
    return GAP.Globals.GAPToJulia( CapAndHomalg.GAP.Globals.IsNormalVariety( v.GapToricVariety ) )
    
end
export IsNormalVariety


function IsAffine( v::JToricVariety )
    
    return GAP.Globals.GAPToJulia( CapAndHomalg.GAP.Globals.IsAffine( v.GapToricVariety ) )
    
end
export IsAffine


function IsProjective( v::JToricVariety )
    
    return GAP.Globals.GAPToJulia( CapAndHomalg.GAP.Globals.IsProjective( v.GapToricVariety ) )
    
end
export IsProjective


function IsSmooth( v::JToricVariety )
    
    return GAP.Globals.GAPToJulia( CapAndHomalg.GAP.Globals.IsSmooth( v.GapToricVariety ) )
    
end
export IsSmooth


function IsComplete( v::JToricVariety )
    
    return GAP.Globals.GAPToJulia( CapAndHomalg.GAP.Globals.IsComplete( v.GapToricVariety ) )
    
end
export IsComplete


function HasTorusfactor( v::JToricVariety )
    
    return GAP.Globals.GAPToJulia( CapAndHomalg.GAP.Globals.HasTorusfactor( v.GapToricVariety ) )
    
end
export HasTorusfactor


function HasNoTorusfactor( v::JToricVariety )
    
    return GAP.Globals.GAPToJulia( CapAndHomalg.GAP.Globals.HasNoTorusfactor( v.GapToricVariety ) )
    
end
export HasNoTorusfactor


function IsOrbifold( v::JToricVariety )
    
    return GAP.Globals.GAPToJulia( CapAndHomalg.GAP.Globals.IsOrbifold( v.GapToricVariety ) )
    
end
export IsOrbifold


function IsSimplicial( v::JToricVariety )
    
    return GAP.Globals.GAPToJulia( CapAndHomalg.GAP.Globals.IsSimplicial( v.GapToricVariety ) )
    
end
export IsSimplicial


function IsIsomorphicToProjectiveSpace( v::JToricVariety )
    
    return GAP.Globals.GAPToJulia( CapAndHomalg.GAP.Globals.IsIsomorphicToProjectiveSpace( v.GapToricVariety ) )
    
end
export IsIsomorphicToProjectiveSpace


function IsDirectProductOfPNs( v::JToricVariety )
    
    return GAP.Globals.GAPToJulia( CapAndHomalg.GAP.Globals.IsDirectProductOfPNs( v.GapToricVariety ) )
    
end
export IsDirectProductOfPNs
