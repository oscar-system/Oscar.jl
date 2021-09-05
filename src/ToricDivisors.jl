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

"""
    create_divisor( c, v )

Construct the torus invariant divisor on the normal toric variety `v` as linear combination of the torus invariant prime divisors of `v`. The coefficients of thi linear combination are passed as list of integers as first argument.

# Examples
```julia-repl
julia> create_divisor( [1,1,2], projective_space( 2 ) )
toric_divisor(GAP: <A divisor of a toric variety with coordinates ( 1, 1, 2 )>)
```
"""
function create_divisor( coeffs::Vector{Int}, v::NormalToricVariety )
    # create the divisor
    gap_coeffs = GapObj( coeffs )
    gap_divisor = GAP.Globals.CreateDivisor( gap_coeffs, v.GapNTV )
    
    # wrap and return
    return toric_divisor( gap_divisor )
end
export create_divisor

"""
    divisor_of_character( c, v )

Construct the torus invariant divisor on the normal toric variety `v` corresponding to the character `c`.

# Examples
```julia-repl
julia> divisor_of_character( [1,2], projective_space( 2 ) )
toric_divisor(GAP: <A principal divisor of a toric variety with coordinates ( -3, 2, 1 )>)
```
"""
function divisor_of_character( character::Vector{Int}, v::NormalToricVariety )
    # create the divisor
    gap_character = GapObj( character )
    gap_divisor = GAP.Globals.DivisorOfCharacter( gap_character, v.GapNTV )
    
    # wrap and return
    return toric_divisor( gap_divisor )
end
export divisor_of_character

"""
    divisor_of_class( c, v )

Construct a torus invariant divisor on the normal toric variety `v` corresponding to the divisor class `c`.

# Examples
```julia-repl
julia> divisor_of_class( [1], projective_space( 2 ) )
toric_divisor(GAP: <A divisor of a toric variety with coordinates ( 1, 0, 0 )>)
```
"""
function divisor_of_class( v::NormalToricVariety, class::Vector{Int} )
    # create the divisor
    gap_class = GapObj( class )
    gap_divisor = GAP.Globals.DivisorOfGivenClass( v.GapNTV, gap_class )
    
    # wrap and return
    return toric_divisor( gap_divisor )
end
export divisor_of_class


######################
# 3: Properties
######################

"""
    is_cartier( d )

Checks if the divisor `d` is Cartier.
"""
function is_cartier( d::toric_divisor )
    return GAP.Globals.IsCartier( d.GapToricDivisor )::Bool
end
export is_cartier


"""
    is_principal( d )

Checks if the divisor `d` is principal.
"""
function is_principal( d::toric_divisor )
    return GAP.Globals.IsPrincipal( d.GapToricDivisor )::Bool
end
export is_principal


"""
    is_primedivisor( d )

Checks if the divisor `d` is prime.
"""
function is_primedivisor( d::toric_divisor )
    return GAP.Globals.IsPrimedivisor( d.GapToricDivisor )::Bool
end
export is_primedivisor


"""
    is_basepoint_free( d )

Checks if the divisor `d` is basepoint free.
"""
function is_basepoint_free( d::toric_divisor )
    return GAP.Globals.IsBasepointFree( d.GapToricDivisor )::Bool
end
export is_basepoint_free


"""
    is_ample( d )

Checks if the divisor `d` is ample.
"""
function is_ample( d::toric_divisor )
    return GAP.Globals.IsAmple( d.GapToricDivisor )::Bool
end
export is_ample


"""
    is_very_ample( d )

For ample divisors `d`, this method checks if `d` is very ample.
"""
function is_very_ample( d::toric_divisor )
    if ! is_ample( d )
        @warn "Can (current) only tell for ample toric divisors if they are very ample."
        return "fail"
    end
    
    return GAP.Globals.IsVeryAmple( d.GapToricDivisor )::Bool
end
export is_very_ample


"""
    is_numerically_effective( d )

Checks if the divisor `d` is numerically effective.
"""
function is_numerically_effective( d::toric_divisor )
    return GAP.Globals.IsNumericallyEffective( d.GapToricDivisor )::Bool
end
export is_numerically_effective
