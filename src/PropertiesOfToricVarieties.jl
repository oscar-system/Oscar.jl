using GAP
using CapAndHomalg

function IsSmooth( v::JToricVariety )
    
    return GAP.Globals.GAPToJulia( CapAndHomalg.GAP.Globals.IsSmooth( v.GapToricVariety ) )
    
end

function IsComplete( v::JToricVariety )
    
    return GAP.Globals.GAPToJulia( CapAndHomalg.GAP.Globals.IsComplete( v.GapToricVariety ) )
    
end
