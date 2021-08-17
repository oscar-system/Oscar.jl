using GAP
using CapAndHomalg

include("ToricVarieties.jl")


######################
# 1: The Julia type for ToricDivisors
######################

struct toric_divisor
           GapToricDivisor::GapObj
end
export toric_divisor

######################
# 2: Generic constructors
######################

function create_divisor( coeffs::Vector{Int64}, v::toric_variety )
    # create the divisor
    gap_coeffs = GapObj( coeffs )
    gap_divisor = GAP.Globals.CreateDivisor( gap_coeffs, v.GapToricVariety )
    
    # wrap and return
    return toric_divisor( gap_divisor )
end
export create_divisor

function divisor_of_character( character::Vector{Int64}, v::toric_variety )
    # create the divisor
    gap_character = GapObj( character )
    gap_divisor = GAP.Globals.DivisorOfCharacter( gap_character, v.GapToricVariety )
    
    # wrap and return
    return toric_divisor( gap_divisor )
end
export divisor_of_character

function divisor_of_class( v::toric_variety, class::Vector{Int64} )
    # create the divisor
    gap_class = GapObj( class )
    gap_divisor = GAP.Globals.DivisorOfGivenClass( v.GapToricVariety, gap_class )
    
    # wrap and return
    return toric_divisor( gap_divisor )
end
export divisor_of_class


######################
# 3: Properties
######################

function is_cartier( d::toric_divisor )
    return Bool( GAP.Globals.IsCartier( d.GapToricDivisor ) )
end
export is_cartier


function is_principal( d::toric_divisor )
    return Bool( GAP.Globals.IsPrincipal( d.GapToricDivisor ) )
end
export is_principal


function is_primedivisor( d::toric_divisor )
    return Bool( GAP.Globals.IsPrimedivisor( d.GapToricDivisor ) )
end
export is_primedivisor


function is_basepoint_free( d::toric_divisor )
    return Bool( GAP.Globals.IsBasepointFree( d.GapToricDivisor ) )
end
export is_basepoint_free


function is_ample( d::toric_divisor )
    return Bool( GAP.Globals.IsAmple( d.GapToricDivisor ) )
end
export is_ample


function is_very_ample( d::toric_divisor )
    if ! is_ample( d )
        @warn "Can (current) only tell for ample toric divisors if they are very ample."
        return "fail"
    end
    
    return Bool( GAP.Globals.IsVeryAmple( d.GapToricDivisor ) )
end
export is_very_ample


function is_numerically_effective( d::toric_divisor )
    return Bool( GAP.Globals.IsNumericallyEffective( d.GapToricDivisor ) )
end
export is_numerically_effective
