######################
# 1: The Julia type for ToricDivisors
######################

@attributes mutable struct ToricDivisor
           polymake_divisor::Polymake.BigObject
           toricvariety::AbstractNormalToricVariety
           coeffs::Vector{fmpz}
           ToricDivisor(polymake_divisor::Polymake.BigObject, toricvariety::AbstractNormalToricVariety, coeffs::Vector{fmpz}) = new(polymake_divisor, toricvariety, coeffs)
end
export ToricDivisor

function pm_tdivisor(td::ToricDivisor)
    return td.polymake_divisor
end

######################
# 2: Generic constructors
######################

@doc Markdown.doc"""
    ToricDivisor(v::AbstractNormalToricVariety, coeffs::Vector{Int})

Construct the torus invariant divisor on the normal toric variety `v` as linear combination of the torus invariant prime divisors of `v`. The coefficients of thi linear combination are passed as list of integers as first argument.

# Examples
```jldoctest
julia> ToricDivisor(toric_projective_space(2), [1,1,2])
A torus-invariant, non-prime divisor on a normal toric variety
```
"""
ToricDivisor(v::AbstractNormalToricVariety, coeffs::Vector{Int}) = ToricDivisor(v, [fmpz(k) for k in coeffs])

function ToricDivisor(v::AbstractNormalToricVariety, coeffs::Vector{fmpz})
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
        set_attribute!(td, :is_prime_divisor, false)
    else
        set_attribute!(td, :is_prime_divisor, all(y -> (y == 1 || y == 0), coeffs))
    end
    
    # return the result
    return td
end
export ToricDivisor


######################
# 3: Special constructors
######################

@doc Markdown.doc"""
    DivisorOfCharacter(v::AbstractNormalToricVariety, character::Vector{Int})

Construct the torus invariant divisor associated to a character of the normal toric variety `v`.

# Examples
```jldoctest
julia> DivisorOfCharacter(toric_projective_space(2), [1,2])
A torus-invariant, non-prime divisor on a normal toric variety
```
"""
function DivisorOfCharacter(v::AbstractNormalToricVariety, character::Vector{Int})
    if length(character) != rank(character_lattice(v))
        throw(ArgumentError("Character must consist of $(rank(character_lattice(v))) integers!"))
    end
    f = map_from_character_to_principal_divisors(v)
    char = sum(character .* gens(domain(f)))
    coeffs = [fmpz(x) for x in transpose(f(char).coeff)][:,1]
    return ToricDivisor(v, coeffs)
end
export DivisorOfCharacter


######################
### 4: Display
######################s

function Base.show(io::IO, td::ToricDivisor)
    # initiate properties string
    properties_string = ["A torus-invariant"]
    
    # cartier?
    if has_attribute(td, :iscartier)
        if get_attribute(td, :iscartier)
            push!(properties_string, "cartier")
        else
            if has_attribute(td, :is_q_cartier)
                if get_attribute(td, :is_q_cartier)
                    push!(properties_string, "q-cartier")
                else
                    push!(properties_string, "non-q-cartier")
                end
            else
                push!(properties_string, "non-cartier")
            end
        end
    end
    
    # principal?
    if has_attribute(td, :isprincipal)
        if get_attribute(td, :isprincipal)
            push!(properties_string, "principal")
        else
            push!(properties_string, "non-principal")
        end
    end
    
    # basepoint free?
    if has_attribute(td, :is_basepoint_free)
        if get_attribute(td, :is_basepoint_free)
            push!(properties_string, "basepoint-free")
        else
            push!(properties_string, "non-basepoint-free")
        end
    end
    
    # effective?
    if has_attribute(td, :iseffective)
        if get_attribute(td, :iseffective)
            push!(properties_string, "effective")
        else
            push!(properties_string, "non-effective")
        end
    end
    
    # integral?
    if has_attribute(td, :isintegral)
        if get_attribute(td, :isintegral)
            push!(properties_string, "integral")
        else
            push!(properties_string, "non-integral")
        end
    end
    
    # (very) ample?
    if has_attribute(td, :isample)
        if get_attribute(td, :isample)
            push!(properties_string, "ample")
        else
            if has_attribute(td, :is_very_ample)
                if get_attribute(td, :is_very_ample)
                    push!(properties_string, "very-ample")
                else
                    push!(properties_string, "non-very-ample")
                end
            else
                push!(properties_string, "non-ample")
            end
        end
    end
    
    # nef?
    if has_attribute(td, :isnef)
        if get_attribute(td, :isnef)
            push!(properties_string, "nef")
        else
            push!(properties_string, "non-nef")
        end
    end
    
    # prime divisor?
    if has_attribute(td, :is_prime_divisor)
        if get_attribute(td, :is_prime_divisor)
            push!(properties_string, "prime")
        else
            push!(properties_string, "non-prime")
        end
    end
    
    # print
    push!(properties_string, "divisor on a normal toric variety")
    join(io, properties_string, ", ", " ")
end
