include("types.jl")

include("auxiliary.jl")

include("FamilyOfSpaces/constructors.jl")
include("FamilyOfSpaces/attributes.jl")

include("AbstractFTheoryModels/attributes.jl")
include("AbstractFTheoryModels/properties.jl")
include("AbstractFTheoryModels/methods.jl")

include("WeierstrassModels/constructors.jl")
include("WeierstrassModels/attributes.jl")
include("WeierstrassModels/methods.jl")

include("TateModels/constructors.jl")
include("TateModels/attributes.jl")
include("TateModels/methods.jl")

include("HypersurfaceModels/constructors.jl")
include("HypersurfaceModels/attributes.jl")
include("HypersurfaceModels/methods.jl")

include("standard_constructions.jl")

include("LiteratureModels/constructors.jl")
include("LiteratureModels/create_index.jl")

include("G4Fluxes/constructors.jl")
include("G4Fluxes/attributes.jl")
include("G4Fluxes/properties.jl")
include("G4Fluxes/auxiliary.jl")
include("G4Fluxes/special-intersection-theory.jl")
include("G4Fluxes/special_attributes.jl")

include("Serialization/tate_models.jl")
include("Serialization/weierstrass_models.jl")
include("Serialization/hypersurface_models.jl")
include("Serialization/qsm_models.jl")

include("exports.jl")
