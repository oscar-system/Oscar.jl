######################
# 1: The Julia type for ToricDivisors
######################

@attributes mutable struct ToricDivisor
           polymake_divisor::Polymake.BigObject
           toric_variety::AbstractNormalToricVariety
           coeffs::Vector{fmpz}
           ToricDivisor(polymake_divisor::Polymake.BigObject,
                        toric_variety::AbstractNormalToricVariety,
                        coeffs::Vector{T}) where {T <: IntegerUnion} = new(polymake_divisor, toric_variety, [fmpz(c) for c in coeffs])
end
export ToricDivisor

function pm_tdivisor(td::ToricDivisor)
    return td.polymake_divisor
end


######################
# 2: Generic constructors
######################

@doc Markdown.doc"""
    ToricDivisor(v::AbstractNormalToricVariety, coeffs::Vector{T}) where {T <: IntegerUnion}

Construct the torus invariant divisor on the normal toric variety `v` as linear
 combination of the torus invariant prime divisors of `v`. The coefficients of this
 linear combination are passed as list of integers as first argument.

# Examples
```jldoctest
julia> ToricDivisor(projective_space(NormalToricVariety, 2), [1,1,2])
A torus-invariant, non-prime divisor on a normal toric variety
```
"""
function ToricDivisor(v::AbstractNormalToricVariety, coeffs::Vector{T}) where {T <: IntegerUnion}
    # check input
    if length(coeffs) != pm_object(v).N_RAYS
        throw(ArgumentError("Number of coefficients needs to match number of prime divisors!"))
    end
    
    # construct the divisor
    ptd = Polymake.fulton.TDivisor(COEFFICIENTS=Polymake.Vector{Polymake.Integer}(coeffs))
    Polymake.add(pm_object(v), "DIVISOR", ptd)
    td = ToricDivisor(ptd, v, coeffs)
    
    # set attributes
    if sum(coeffs) != 1
        set_attribute!(td, :is_prime, false)
    else
        set_attribute!(td, :is_prime, all(y -> (y == 1 || y == 0), coeffs))
    end
    
    # return the result
    return td
end
export ToricDivisor


######################
# 3: Special constructors
######################

@doc Markdown.doc"""
    DivisorOfCharacter(v::AbstractNormalToricVariety, character::Vector{T}) where {T <: IntegerUnion}

Construct the torus invariant divisor associated to a character of the normal toric variety `v`.

# Examples
```jldoctest
julia> DivisorOfCharacter(projective_space(NormalToricVariety, 2), [1,2])
A torus-invariant, non-prime divisor on a normal toric variety
```
"""
function DivisorOfCharacter(v::AbstractNormalToricVariety, character::Vector{T}) where {T <: IntegerUnion}
    if length(character) != rank(character_lattice(v))
        throw(ArgumentError("Character must consist of $(rank(character_lattice(v))) integers!"))
    end
    f = map_from_character_lattice_to_torusinvariant_weil_divisor_group(v)
    char = sum(character .* gens(domain(f)))
    coeffs = [fmpz(x) for x in transpose(f(char).coeff)][:,1]
    return ToricDivisor(v, coeffs)
end
export DivisorOfCharacter


########################
# 4: Addition and scalar multiplication
########################

function Base.:+(td1::ToricDivisor, td2::ToricDivisor)
    if toric_variety(td1) !== toric_variety(td2)
        throw(ArgumentError("The toric divisors must be defined on identically the same toric variety."))
    end
    new_coeffiicients = coefficients(td1) + coefficients(td2)
    return ToricDivisor(toric_variety(td1), new_coeffiicients)
end


function Base.:-(td1::ToricDivisor, td2::ToricDivisor)
    if toric_variety(td1) !== toric_variety(td2)
        throw(ArgumentError("The toric divisors must be defined on identically the same toric variety."))
    end
    new_coeffiicients = coefficients(td1) - coefficients(td2)
    return ToricDivisor(toric_variety(td1), new_coeffiicients)
end


Base.:*(c::T, td::ToricDivisor) where {T <: IntegerUnion} = ToricDivisor(toric_variety(td), [fmpz(c)*x for x in coefficients(td)])


######################
# 5: Equality
######################s

function Base.:(==)(td1::ToricDivisor, td2::ToricDivisor)
    return toric_variety(td1) === toric_variety(td2) && coefficients(td1) == coefficients(td2)
end


######################
# 6: Display
######################s

function Base.show(io::IO, td::ToricDivisor)
    # initiate properties string
    properties_string = ["A torus-invariant"]
    
    q_car_cb!(a,b) = push_attribute_if_exists!(a, b, :is_q_cartier, "q_cartier")
    push_attribute_if_exists!(properties_string, td, :iscartier, "cartier"; callback=q_car_cb!)
    push_attribute_if_exists!(properties_string, td, :isprincipal, "principal")
    push_attribute_if_exists!(properties_string, td, :is_basepoint_free, "basepoint-free")
    push_attribute_if_exists!(properties_string, td, :iseffective, "effective")
    push_attribute_if_exists!(properties_string, td, :isintegral, "integral")
    
    ample_cb!(a,b) = push_attribute_if_exists!(a, b, :isample, "ample")
    push_attribute_if_exists!(properties_string, td, :is_very_ample, "very-ample"; callback=ample_cb!)
    
    push_attribute_if_exists!(properties_string, td, :isnef, "nef")
    push_attribute_if_exists!(properties_string, td, :is_prime, "prime")
    
    # print
    push!(properties_string, "divisor on a normal toric variety")
    join(io, properties_string, ", ", " ")
end
