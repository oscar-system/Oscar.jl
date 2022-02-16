# A little helper for displaying properties of an object such that logical
# implications reduce the list of the display. I.e. if a divisor is very ample,
# it is also ample, so we do not need to display ample in that case.
#
# If the property is A, then the callback can be used to ask for property B
# with A=>B. Then if A is true, only A is displayed, and if B is false, only B
# is displayed.
function push_attribute_if_exists!(result::Vector{String}, 
        v::T,
        property::Symbol, name::String, 
        non_name::String="non-"*name; 
        callback=nothing
    ) where T
    if has_attribute(v, property)
        if get_attribute(v, property)
            push!(result, name)
            return true
        else
            cbval = true
            if !isnothing(callback)
                cbval = callback(result, v)
            end
            if isnothing(cbval) || cbval
                push!(result, non_name)
            end
            return false
        end
    else
        if !isnothing(callback)
            return callback(result, v)
        else
            return nothing
        end
    end
end

include("NormalToricVarieties/constructors.jl")
include("NormalToricVarieties/auxilliary.jl")
include("NormalToricVarieties/properties.jl")
include("NormalToricVarieties/attributes.jl")
include("NormalToricVarieties/methods.jl")

include("CyclicQuotientSingularities/CyclicQuotientSingularities.jl")

include("ToricDivisors/constructors.jl")
include("ToricDivisors/properties.jl")
include("ToricDivisors/attributes.jl")

include("ToricLineBundles/constructors.jl")
include("ToricLineBundles/properties.jl")
include("ToricLineBundles/attributes.jl")


include("cohomCalg.jl")

# deprecated functions
@deprecate map_from_character_to_principal_divisors(v::AbstractNormalToricVariety) map_from_character_lattice_to_torus_invariant_weil_divisor_group(v)
