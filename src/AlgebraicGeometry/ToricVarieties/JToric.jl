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
include("NormalToricVarieties/toric_ideal.jl")
include("NormalToricVarieties/properties.jl")
include("NormalToricVarieties/attributes.jl")
include("NormalToricVarieties/betti_numbers.jl")
include("NormalToricVarieties/standard_constructions.jl")

include("CyclicQuotientSingularities/CyclicQuotientSingularities.jl")

include("ToricDivisors/constructors.jl")
include("ToricDivisors/properties.jl")
include("ToricDivisors/attributes.jl")
include("ToricDivisors/standard_constructions.jl")

include("ToricDivisorClasses/constructors.jl")
include("ToricDivisorClasses/properties.jl")
include("ToricDivisorClasses/attributes.jl")
include("ToricDivisorClasses/standard_constructions.jl")

include("ToricLineBundles/constructors.jl")
include("ToricLineBundles/properties.jl")
include("ToricLineBundles/attributes.jl")
include("ToricLineBundles/standard_constructions.jl")

include("CohomologyClasses/constructors.jl")
include("CohomologyClasses/properties.jl")
include("CohomologyClasses/attributes.jl")
include("CohomologyClasses/special_attributes.jl")
include("CohomologyClasses/methods.jl")

include("cohomCalg/VanishingSets/constructors.jl")
include("cohomCalg/VanishingSets/attributes.jl")
include("cohomCalg/VanishingSets/methods.jl")

include("cohomCalg/auxiliary.jl")
include("cohomCalg/special_attributes.jl")

include("Subvarieties/constructors.jl")
include("Subvarieties/properties.jl")
include("Subvarieties/attributes.jl")

include("AlgebraicCycles/constructors.jl")
include("AlgebraicCycles/properties.jl")
include("AlgebraicCycles/attributes.jl")
include("AlgebraicCycles/special_attributes.jl")

include("ToricMorphisms/constructors.jl")
include("ToricMorphisms/attributes.jl")
include("ToricMorphisms/standard_constructions.jl")

include("ToricSchemes/attributes.jl")
