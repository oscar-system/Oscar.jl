# Project files
include("AlternatingBilinearForm.jl")
include("HelperFunctions.jl")
include("GlobalGenericFormGeneration.jl")
include("LocalGenericFormGeneration.jl")
include("DrinfeldHeckeFormValidation.jl")
include("DrinfeldHeckeForm.jl")
include("DrinfeldHeckeAlgebra.jl")
include("SymplecticDoubling.jl")

# Exports
export BilinearForm
export alternating_bilinear_form

export DrinfeldHeckeForm
export drinfeld_hecke_form
export generic_drinfeld_hecke_form
export evaluate_parameters
export set_forms
export nforms
export is_generic

export DrinfeldHeckeAlgebra
export drinfeld_hecke_algebra
export generic_drinfeld_hecke_algebra
export symmetric_algebra

export symplectic_doubling

# Aliases
