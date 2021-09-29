######################
# 1: The Julia type for ToricDivisors
######################

struct ToricDivisor
           GapToricDivisor::GapObj
end
export ToricDivisor

######################
# 2: Generic constructors
######################

"""
    create_divisor( c::Vector{Int}, v::AbstractNormalToricVariety )

Construct the torus invariant divisor on the normal toric variety `v` as linear combination of the torus invariant prime divisors of `v`. The coefficients of thi linear combination are passed as list of integers as first argument.

# Examples
```julia-repl
julia> create_divisor( [1,1,2], projective_space( 2 ) )
toric_divisor(GAP: <A divisor of a toric variety with coordinates ( 1, 1, 2 )>)
```
"""
function create_divisor( coeffs::Vector{Int}, v::AbstractNormalToricVariety )
    # create the divisor
    gap_coeffs = GapObj( coeffs )
    gap_divisor = GAP.Globals.CreateDivisor( gap_coeffs, v.GapNTV )
    
    # wrap and return
    return ToricDivisor( gap_divisor )
end
export create_divisor

"""
    divisor_of_character( c::Vector{Int}, v::AbstractNormalToricVariety )

Construct the torus invariant divisor on the normal toric variety `v` corresponding to the character `c`.

# Examples
```julia-repl
julia> divisor_of_character( [1,2], projective_space( 2 ) )
toric_divisor(GAP: <A principal divisor of a toric variety with coordinates ( -3, 2, 1 )>)
```
"""
function divisor_of_character( character::Vector{Int}, v::AbstractNormalToricVariety )
    # create the divisor
    gap_character = GapObj( character )
    gap_divisor = GAP.Globals.DivisorOfCharacter( gap_character, v.GapNTV )
    
    # wrap and return
    return ToricDivisor( gap_divisor )
end
export divisor_of_character

"""
    divisor_of_class( v::AbstractNormalToricVariety, c::Vector{Int} )

Construct a torus invariant divisor on the normal toric variety `v` corresponding to the divisor class `c`.

# Examples
```julia-repl
julia> divisor_of_class( [1], projective_space( 2 ) )
toric_divisor(GAP: <A divisor of a toric variety with coordinates ( 1, 0, 0 )>)
```
"""
function divisor_of_class( v::AbstractNormalToricVariety, class::Vector{Int} )
    # create the divisor
    gap_class = GapObj( class )
    gap_divisor = GAP.Globals.DivisorOfGivenClass( v.GapNTV, gap_class )
    
    # wrap and return
    return ToricDivisor( gap_divisor )
end
export divisor_of_class


######################
# 3: Properties
######################

"""
    is_cartier( d::ToricDivisor )

Checks if the divisor `d` is Cartier.
"""
function is_cartier( d::ToricDivisor )
    return GAP.Globals.IsCartier( d.GapToricDivisor )::Bool
end
export is_cartier


"""
    is_principal( d::ToricDivisor )

Checks if the divisor `d` is principal.
"""
function is_principal( d::ToricDivisor )
    return GAP.Globals.IsPrincipal( d.GapToricDivisor )::Bool
end
export is_principal


"""
    is_primedivisor( d::ToricDivisor )

Checks if the divisor `d` is prime.
"""
function is_primedivisor( d::ToricDivisor )
    return GAP.Globals.IsPrimedivisor( d.GapToricDivisor )::Bool
end
export is_primedivisor


"""
    is_basepoint_free( d::ToricDivisor )

Checks if the divisor `d` is basepoint free.
"""
function is_basepoint_free( d::ToricDivisor )
    return GAP.Globals.IsBasepointFree( d.GapToricDivisor )::Bool
end
export is_basepoint_free


"""
    is_ample( d::ToricDivisor )

Checks if the divisor `d` is ample.
"""
function is_ample( d::ToricDivisor )
    return GAP.Globals.IsAmple( d.GapToricDivisor )::Bool
end
export is_ample


"""
    is_very_ample( d::ToricDivisor )

For ample divisors `d`, this method checks if `d` is very ample.
"""
function is_very_ample( d::ToricDivisor )
    if ! is_ample( d )
        @warn "Can (current) only tell for ample toric divisors if they are very ample."
        return "fail"
    end
    
    return GAP.Globals.IsVeryAmple( d.GapToricDivisor )::Bool
end
export is_very_ample


"""
    is_nef( d::ToricDivisor )

Checks if the divisor `d` is numerically effective.
"""
function is_nef( d::ToricDivisor )
    return GAP.Globals.IsNumericallyEffective( d.GapToricDivisor )::Bool
end
export is_nef
