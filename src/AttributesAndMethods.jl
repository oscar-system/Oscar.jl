using GAP
using CapAndHomalg

include("ToricDivisors.jl")

######################
# 1: Attributes of ToricVarieties
######################

function AffineOpenCovering( v::JToricVariety )
    
    gap_cover = CapAndHomalg.GAP.Globals.AffineOpenCovering( v.GapToricVariety )
    len = GAP.Globals.GAPToJulia( GAP.Globals.Length( gap_cover ) )
    j_cover = [ JToricVariety( 1, gap_cover[ i ] ) for i in 1:len ]
    return j_cover
    
end
export AffineOpenCovering


struct JCoxRing
           bar
           GapCoxRing
end
export JCoxRing

function CoxRing( v::JToricVariety )
    
    gap_ring = CapAndHomalg.GAP.Globals.CoxRing( v.GapToricVariety )
    return JCoxRing( 1, gap_ring )
    
end
export CoxRing


function ListOfVariablesOfCoxRing( v::JToricVariety )
    
    gap_variables = CapAndHomalg.GAP.Globals.ListOfVariablesOfCoxRing( v.GapToricVariety )
    len = GAP.Globals.GAPToJulia( GAP.Globals.Length( gap_variables ) )
    julia_variables = [  GAP.Globals.GAPToJulia( gap_variables[ i ] ) for i in 1:len ]
    return julia_variables
    
end
export ListOfVariablesOfCoxRing


struct JClassGroup
           bar
           GapClassGroup
end
export JClassGroup

function ClassGroup( v::JToricVariety )
    
    gap_class_group = CapAndHomalg.GAP.Globals.ClassGroup( v.GapToricVariety )
    return JClassGroup( 1, gap_class_group )
    
end
export ClassGroup
