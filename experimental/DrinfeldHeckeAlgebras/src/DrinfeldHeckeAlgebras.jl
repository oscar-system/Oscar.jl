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
export parameters

export DrinfeldHeckeAlgebra
export drinfeld_hecke_algebra
export parametrized_drinfeld_hecke_algebra

export generate_forms_for_conjugacy_classes

# Aliases
