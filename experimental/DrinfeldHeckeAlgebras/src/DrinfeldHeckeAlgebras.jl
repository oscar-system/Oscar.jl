# Project files
include("AlternatingBilinearForm.jl")
include("SymplecticDoubling.jl")
include("HelperFunctions.jl")
include("GlobalGenericFormGeneration.jl")
include("LocalGenericFormGeneration.jl")
include("DrinfeldHeckeFormValidation.jl")
include("DrinfeldHeckeForm.jl")
include("DrinfeldHeckeAlgebra.jl")
include("DrinfeldHeckeAlgebraMultiplication.jl")

# Exports
export symplectic_doubling
export drinfeld_hecke_algebra
export generic_drinfeld_hecke_algebra
export evaluate_parameters
export base_ring
export base_algebra
export group
export parameters
export params
