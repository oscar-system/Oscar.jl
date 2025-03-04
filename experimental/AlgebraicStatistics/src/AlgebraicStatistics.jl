include("CI.jl")
include("Markov.jl")
include("GraphicalModels.jl")

include("PhylogeneticModels.jl")
include("PhylogeneticAuxiliary.jl")
include("PhylogeneticParametrization.jl")
include("PhylogeneticInvariants.jl")
include("LoadModels.jl")

#export models
export cavender_farris_neyman_model
export jukes_cantor_model
export kimura2_model
export kimura3_model
export general_markov_model
export affine_phylogenetic_model!

#export phylogenetic models attributes
export phylogenetic_model
export graph
export number_states
export transition_matrices
export probability_ring
export root_distribution
export fourier_parameters
export fourier_ring
export group_of_model

#export probability and fourier map
export probability_map
export fourier_map

# export functions to calculate equivalent classes
export compute_equivalent_classes
export sum_equivalent_classes

# export transformation matrices
export specialized_fourier_transform
export inverse_specialized_fourier_transform

# export structs for GroupBasedPhylogeneticModel,PhylogeneticModel
export PhylogeneticModel
export GroupBasedPhylogeneticModel

# export functions to load objects (currently only graphs of phylogenetic models)
export load_phylogenetic_model
