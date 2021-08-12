using GAP
using CapAndHomalg


function IsNormalVariety( v::JToricVariety )
    
    return GAP.Globals.GAPToJulia( CapAndHomalg.GAP.Globals.IsNormalVariety( v.GapToricVariety ) )
    
end


function IsAffine( v::JToricVariety )
    
    return GAP.Globals.GAPToJulia( CapAndHomalg.GAP.Globals.IsAffine( v.GapToricVariety ) )
    
end


function IsProjective( v::JToricVariety )
    
    return GAP.Globals.GAPToJulia( CapAndHomalg.GAP.Globals.IsProjective( v.GapToricVariety ) )
    
end


function IsSmooth( v::JToricVariety )
    
    return GAP.Globals.GAPToJulia( CapAndHomalg.GAP.Globals.IsSmooth( v.GapToricVariety ) )
    
end

function IsComplete( v::JToricVariety )
    
    return GAP.Globals.GAPToJulia( CapAndHomalg.GAP.Globals.IsComplete( v.GapToricVariety ) )
    
end


function HasTorusfactor( v::JToricVariety )
    
    return GAP.Globals.GAPToJulia( CapAndHomalg.GAP.Globals.HasTorusfactor( v.GapToricVariety ) )
    
end


function HasNoTorusfactor( v::JToricVariety )
    
    return GAP.Globals.GAPToJulia( CapAndHomalg.GAP.Globals.HasNoTorusfactor( v.GapToricVariety ) )
    
end


function IsOrbifold( v::JToricVariety )
    
    return GAP.Globals.GAPToJulia( CapAndHomalg.GAP.Globals.IsOrbifold( v.GapToricVariety ) )
    
end


function IsSimplicial( v::JToricVariety )
    
    return GAP.Globals.GAPToJulia( CapAndHomalg.GAP.Globals.IsSimplicial( v.GapToricVariety ) )
    
end
