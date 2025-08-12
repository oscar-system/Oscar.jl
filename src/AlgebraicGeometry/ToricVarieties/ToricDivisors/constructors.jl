######################
# 1: The Julia type for ToricDivisors
######################

@attributes mutable struct ToricDivisor <: AbsWeilDivisor{NormalToricVariety, ZZRing}
           polymake_divisor::Polymake.BigObject
           toric_variety::NormalToricVarietyType
           coeffs::Vector{ZZRingElem}
           ToricDivisor(polymake_divisor::Polymake.BigObject,
                        toric_variety::NormalToricVarietyType,
                        coeffs::Vector{T}) where {T <: IntegerUnion} = new(polymake_divisor, toric_variety, [ZZRingElem(c) for c in coeffs])
end

pm_object(td::ToricDivisor) = td.polymake_divisor


######################
# 2: Generic constructors
######################

@doc raw"""
    toric_divisor(v::NormalToricVarietyType, coeffs::Vector{T}) where {T <: IntegerUnion}

Construct the torus invariant divisor on the normal toric variety `v` as linear
 combination of the torus invariant prime divisors of `v`. The coefficients of this
 linear combination are passed as list of integers as first argument.

# Examples
```jldoctest
julia> toric_divisor(projective_space(NormalToricVariety, 2), [1, 1, 2])
Torus-invariant, non-prime divisor on a normal toric variety
```
"""
function toric_divisor(v::NormalToricVarietyType, coeffs::Vector{T}) where {T <: IntegerUnion}
    # check input
    @req length(coeffs) == pm_object(v).N_RAYS "Number of coefficients needs to match number of prime divisors"
    
    # construct the divisor
    ptd = Polymake.fulton.TDivisor(COEFFICIENTS=Polymake.Vector{Polymake.Integer}(coeffs))
    Polymake.add(pm_object(v), "DIVISOR", ptd)
    td = ToricDivisor(ptd, v, coeffs)
    
    # set attributes
    if sum(coeffs) != 1
        set_attribute!(td, :is_prime, false)
    else
        set_attribute!(td, :is_prime, all(y -> is_zero(y) || is_one(y), coeffs))
    end
    
    # return the result
    return td
end


######################
# 3: Special constructors
######################

@doc raw"""
    divisor_of_character(v::NormalToricVarietyType, character::Vector{T}) where {T <: IntegerUnion}

Construct the torus invariant divisor associated to a character of the normal toric variety `v`.

# Examples
```jldoctest
julia> divisor_of_character(projective_space(NormalToricVariety, 2), [1, 2])
Torus-invariant, non-prime divisor on a normal toric variety
```
"""
function divisor_of_character(v::NormalToricVarietyType, character::Vector{T}) where {T <: IntegerUnion}
    r = torsion_free_rank(character_lattice(v))
    @req length(character) == r "Character must consist of $r integers"
    f = map_from_character_lattice_to_torusinvariant_weil_divisor_group(v)
    char = sum(character .* gens(domain(f)))
    coeffs = [ZZRingElem(x) for x in transpose(f(char).coeff)][:, 1]
    return toric_divisor(v, coeffs)
end


########################
# 4: Addition and scalar multiplication
########################

function Base.:+(td1::ToricDivisor, td2::ToricDivisor)
    @req toric_variety(td1) === toric_variety(td2) "The toric divisors must be defined on the same toric variety"
    new_coefficients = coefficients(td1) + coefficients(td2)
    return toric_divisor(toric_variety(td1), new_coefficients)
end


function Base.:-(td1::ToricDivisor, td2::ToricDivisor)
    @req toric_variety(td1) === toric_variety(td2) "The toric divisors must be defined on the same toric variety"
    new_coefficients = coefficients(td1) - coefficients(td2)
    return toric_divisor(toric_variety(td1), new_coefficients)
end

Base.:-(td :: ToricDivisor) = toric_divisor(toric_variety(td), -coefficients(td))

Base.:*(c::T, td::ToricDivisor) where {T <: IntegerUnion} = toric_divisor(toric_variety(td), [ZZRingElem(c)*x for x in coefficients(td)])
# method ambiguity requires the implementation of the following method
Base.:*(c::ZZRingElem, td::ToricDivisor) = toric_divisor(toric_variety(td), [c*x for x in coefficients(td)])


######################
# 5: Equality and hash
######################s

function Base.:(==)(td1::ToricDivisor, td2::ToricDivisor)
    @req toric_variety(td1) === toric_variety(td2) "The toric divisors must be defined on the same toric variety"
    return coefficients(td1) == coefficients(td2)
end

function Base.hash(td::ToricDivisor, h::UInt)
    b = 0x92bd6ac4f87d834e % UInt
    h = hash(toric_variety(td), h)
    h = hash(coefficients(td), h)
    return xor(h, b)
end


######################
# 6: Display
######################s

function Base.show(io::IO, td::ToricDivisor)
  # initiate properties string
  properties_string = ["Torus-invariant"]
  
  q_car_cb!(a, b) = push_attribute_if_exists!(a, b, :is_q_cartier, "q_cartier")
  push_attribute_if_exists!(properties_string, td, :is_cartier, "cartier"; callback=q_car_cb!)
  push_attribute_if_exists!(properties_string, td, :is_principal, "principal")
  push_attribute_if_exists!(properties_string, td, :is_basepoint_free, "basepoint-free")
  push_attribute_if_exists!(properties_string, td, :is_effective, "effective")
  push_attribute_if_exists!(properties_string, td, :is_integral, "integral")
  
  ample_cb!(a, b) = push_attribute_if_exists!(a, b, :is_ample, "ample")
  push_attribute_if_exists!(properties_string, td, :is_very_ample, "very-ample"; callback=ample_cb!)
  
  push_attribute_if_exists!(properties_string, td, :is_nef, "nef")
  push_attribute_if_exists!(properties_string, td, :is_prime, "prime")
  
  # print
  push!(properties_string, "divisor on a normal toric variety")
  join(io, properties_string, ", ", " ")
end

function Base.show(io::IO, ::MIME"text/plain", td::ToricDivisor)
  show(io, td)
end
