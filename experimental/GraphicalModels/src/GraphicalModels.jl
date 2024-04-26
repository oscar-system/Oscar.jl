include("PhylogeneticModels.jl")
include("PhylogeneticAuxiliary.jl")
include("PhylogeneticParametrisation.jl")

#export models
export cavender_farris_neyman_model 
export jukes_cantor_model
export kimura2_model
export kimura3_model
export general_markov_model

#export probability and fourier map
export probability_map
export fourier_map

# export functions to calculate equivalent classes
export compute_equivalent_classes #are these good names for general Oscar? Do we want both exported?
export sum_equivalent_classes #are these good names for general Oscar?

# export transformation matrices
export specialized_fourier_transform
export inverse_specialized_fourier_transform

# export structs for GroupBasedPhylogeneticModel,PhylogeneticModel
export PhylogeneticModel
export GroupBasedPhylogeneticModel
