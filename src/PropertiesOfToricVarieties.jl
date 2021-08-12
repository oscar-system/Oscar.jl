using GAP
using CapAndHomalg

function IsSmooth( v::JToricVariety )
    
    return GAP.Globals.GAPToJulia( CapAndHomalg.GAP.Globals.IsSmooth( v.GapToricVariety ) )
    
end

function IsComplete( v::JToricVariety )
    
    return GAP.Globals.GAPToJulia( CapAndHomalg.GAP.Globals.IsComplete( v.GapToricVariety ) )
    
end

function IsAffine( v::JToricVariety )
    
    return GAP.Globals.GAPToJulia( CapAndHomalg.GAP.Globals.IsAffine( v.GapToricVariety ) )
    
end

function IsOrbifold( v::JToricVariety )
    
    return GAP.Globals.GAPToJulia( CapAndHomalg.GAP.Globals.IsOrbifold( v.GapToricVariety ) )
    
end
