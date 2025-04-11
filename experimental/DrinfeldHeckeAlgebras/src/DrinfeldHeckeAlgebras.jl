# Project files
include("AlternatingBilinearForm.jl")
include("RelationHandler.jl")
include("DrinfeldHeckeForm.jl")
include("DrinfeldHeckeAlgebra.jl")

# Exports
export BilinearForm
export alternating_bilinear_form

export DrinfeldHeckeForm
export drinfeld_hecke_form
export parametrized_drinfeld_hecke_form
export evaluate_parameters
export set_forms
export nforms
export is_parametrized

export DrinfeldHeckeAlgebra
export drinfeld_hecke_algebra
export parametrized_drinfeld_hecke_algebra

# Aliases
