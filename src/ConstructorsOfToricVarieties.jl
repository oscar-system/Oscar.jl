using GAP
using CapAndHomalg

# the Julia data type for ToricVarieties
struct JToricVariety
           bar
           GapToricVariety
end

# constructor for JToricVarieties
function JToricVariety( rays::Array{Array{Int64,1},1}, cones::Array{Array{Int64,1},1} )

    # construct the toric variety in GAP
    gap_rays = CapAndHomalg.GAP.Globals.ConvertJuliaToGAP( rays )
    gap_cones = CapAndHomalg.GAP.Globals.ConvertJuliaToGAP( cones )
    fan = CapAndHomalg.GAP.Globals.Fan( gap_rays, gap_cones )
    variety = CapAndHomalg.GAP.Globals.ToricVariety( fan )
    
    # wrap it into a struct and return
    return JToricVariety( 1, variety )

end


# code that wraps projective space
function ProjectiveSpace( x )
    
    variety = CapAndHomalg.GAP.Globals.ProjectiveSpace( x )
    return JToricVariety(  1, variety )
    
end
