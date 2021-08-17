using GAP
using CapAndHomalg

include("ToricVarieties.jl")


######################
# 1: The Julia type for ToricDivisors
######################

struct JToricDivisor
           GapToricDivisor
end
export JToricDivisor

######################
# 2: Generic constructors
######################

function CreateDivisor( coeffs::Vector{Int64}, v::toric_variety )
    # create the divisor
    gap_coeffs = GapObj( coeffs )
    gap_divisor = GAP.Globals.CreateDivisor( gap_coeffs, v.GapToricVariety )
    
    # wrap and return
    return JToricDivisor( gap_divisor )
end
export CreateDivisor

function DivisorOfCharacter( character::Vector{Int64}, v::toric_variety )
    # create the divisor
    gap_character = GapObj( character )
    gap_divisor = GAP.Globals.DivisorOfCharacter( gap_character, v.GapToricVariety )
    
    # wrap and return
    return JToricDivisor( gap_divisor )
end
export DivisorOfCharacter

function DivisorOfGivenClass( v::toric_variety, class::Vector{Int64} )
    # create the divisor
    gap_class = GapObj( class )
    gap_divisor = GAP.Globals.DivisorOfGivenClass( v.GapToricVariety, gap_class )
    
    # wrap and return
    return JToricDivisor( gap_divisor )
end
export DivisorOfGivenClass


######################
# 3: Properties
######################

function IsCartier( d::JToricDivisor )
    return Bool( GAP.Globals.IsCartier( d.GapToricDivisor ) )
end
export IsCartier


function IsPrincipal( d::JToricDivisor )
    return Bool( GAP.Globals.IsPrincipal( d.GapToricDivisor ) )
end
export IsPrincipal


function IsPrimedivisor( d::JToricDivisor )
    return Bool( GAP.Globals.IsPrimedivisor( d.GapToricDivisor ) )
end
export IsPrimedivisor


function IsBasepointFree( d::JToricDivisor )
    return Bool( GAP.Globals.IsBasepointFree( d.GapToricDivisor ) )
end
export IsBasepointFree


function IsAmple( d::JToricDivisor )
    return Bool( GAP.Globals.IsAmple( d.GapToricDivisor ) )
end
export IsAmple


function IsVeryAmple( d::JToricDivisor )
    if ! IsAmple( d )
        @warn "Can (current) only tell for ample toric divisors if they are very ample."
        return "fail"
    end
    
    return Bool( GAP.Globals.IsVeryAmple( d.GapToricDivisor ) )
end
export IsVeryAmple


function IsNumericallyEffective( d::JToricDivisor )
    return Bool( GAP.Globals.IsNumericallyEffective( d.GapToricDivisor ) )
end
export IsNumericallyEffective
