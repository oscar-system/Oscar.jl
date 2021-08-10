# code that wraps projective space
function ProjectiveSpace( x )
    
    return CapAndHomalg.GAP.Globals.ProjectiveSpace( x )
    
end

function Fan( rays, cones )
    
    gap_rays = CapAndHomalg.GAP.Globals.ConvertJuliaToGAP( rays )
    gap_cones = CapAndHomalg.GAP.Globals.ConvertJuliaToGAP( cones )
    return CapAndHomalg.GAP.Globals.Fan( gap_rays, gap_cones )
    
end

function IsSmooth( x )
    
    return CapAndHomalg.GAP.Globals.IsSmooth( x )
    
end
