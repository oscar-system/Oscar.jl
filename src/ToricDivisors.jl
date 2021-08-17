using GAP
using CapAndHomalg

include("ToricVarieties.jl")


######################
# 1: The Julia type for ToricDivisors
######################

struct JToricDivisor
           bar
           GapToricDivisor
end
export JToricDivisor

######################
# 2: Generic constructors
######################

function CreateDivisor( coeffs::Vector{Int64}, v::JToricVariety )

    # create the divisor
    gap_coeffs = CapAndHomalg.GAP.Globals.ConvertJuliaToGAP( coeffs )
    gap_divisor = CapAndHomalg.GAP.Globals.CreateDivisor( gap_coeffs, v.GapToricVariety )
    
    # wrap and return
    return JToricDivisor( 1, gap_divisor )

end
export CreateDivisor

function DivisorOfCharacter( character::Vector{Int64}, v::JToricVariety )

    # create the divisor
    gap_character = CapAndHomalg.GAP.Globals.ConvertJuliaToGAP( character )
    gap_divisor = CapAndHomalg.GAP.Globals.DivisorOfCharacter( gap_character, v.GapToricVariety )
    
    # wrap and return
    return JToricDivisor( 1, gap_divisor )

end
export DivisorOfCharacter

# potentially duplicate of the former -> clean up the gap code!
function DivisorOfGivenClass( v::JToricVariety, class::Vector{Int64} )

    # create the divisor
    gap_class = CapAndHomalg.GAP.Globals.ConvertJuliaToGAP( class )
    gap_divisor = CapAndHomalg.GAP.Globals.DivisorOfGivenClass( v.GapToricVariety, gap_class )
    
    # wrap and return
    return JToricDivisor( 1, gap_divisor )

end
export DivisorOfGivenClass


######################
# 3: Properties
######################

function IsCartier( d::JToricDivisor )
    
    return GAP.Globals.GAPToJulia( CapAndHomalg.GAP.Globals.IsCartier( d.GapToricDivisor ) )
    
end
export IsCartier


function IsPrincipal( d::JToricDivisor )
    
    return GAP.Globals.GAPToJulia( CapAndHomalg.GAP.Globals.IsPrincipal( d.GapToricDivisor ) )
    
end
export IsPrincipal


function IsPrimedivisor( d::JToricDivisor )
    
    return GAP.Globals.GAPToJulia( CapAndHomalg.GAP.Globals.IsPrimedivisor( d.GapToricDivisor ) )
    
end
export IsPrimedivisor


function IsBasepointFree( d::JToricDivisor )
    
    return GAP.Globals.GAPToJulia( CapAndHomalg.GAP.Globals.IsBasepointFree( d.GapToricDivisor ) )
    
end
export IsBasepointFree


function IsAmple( d::JToricDivisor )
    
    return GAP.Globals.GAPToJulia( CapAndHomalg.GAP.Globals.IsAmple( d.GapToricDivisor ) )
    
end
export IsAmple


function IsVeryAmple( d::JToricDivisor )
    
    if ! IsAmple( d )
        @warn "Can (current) only tell for ample toric divisors if they are very ample."
        return "fail"
    end
    
    return GAP.Globals.GAPToJulia( CapAndHomalg.GAP.Globals.IsVeryAmple( d.GapToricDivisor ) )
    
end
export IsVeryAmple


function IsNumericallyEffective( d::JToricDivisor )
    
    return GAP.Globals.GAPToJulia( CapAndHomalg.GAP.Globals.IsNumericallyEffective( d.GapToricDivisor ) )
    
end
export IsNumericallyEffective
