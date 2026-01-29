include("types.jl")

include("auxiliary.jl")

include("FamilyOfSpaces/constructors.jl")
include("FamilyOfSpaces/attributes.jl")

include("AbstractFTheoryModels/Attributes/basic_attributes.jl")
include("AbstractFTheoryModels/Attributes/section_attributes.jl")
include("AbstractFTheoryModels/Attributes/precomputed_attributes.jl")
include("AbstractFTheoryModels/Attributes/not_yet_computed_attributes.jl")
include("AbstractFTheoryModels/Attributes/computed_attributes.jl")
include("AbstractFTheoryModels/Attributes/qsm_attributes.jl")

include("AbstractFTheoryModels/Properties/properties.jl")

include("AbstractFTheoryModels/Methods/adders.jl")
include("AbstractFTheoryModels/Methods/put_over_concrete_base.jl")
include("AbstractFTheoryModels/Methods/blowups.jl")
include("AbstractFTheoryModels/Methods/analyze_fibers.jl")

include("WeierstrassModels/constructors.jl")
include("WeierstrassModels/attributes.jl")

include("TateModels/constructors.jl")
include("TateModels/attributes.jl")

include("HypersurfaceModels/constructors.jl")
include("HypersurfaceModels/attributes.jl")
include("HypersurfaceModels/methods.jl")

include("TunedModels/constructors.jl")

include("LiteratureModels/constructors.jl")
include("LiteratureModels/create_index.jl")

include("G4Fluxes/constructors.jl")
include("G4Fluxes/attributes.jl")
include("G4Fluxes/properties.jl")
include("G4Fluxes/auxiliary.jl")

include("FamilyOfG4Fluxes/constructors.jl")
include("FamilyOfG4Fluxes/attributes.jl")
include("FamilyOfG4Fluxes/properties.jl")
include("FamilyOfG4Fluxes/methods.jl")
include("FamilyOfG4Fluxes/special-intersection-theory.jl")
include("FamilyOfG4Fluxes/special_constructors.jl")

include("serialization.jl")

include("exports.jl")
