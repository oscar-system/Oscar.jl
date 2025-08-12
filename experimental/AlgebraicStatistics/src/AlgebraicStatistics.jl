################################################################################
# important dict type, leaving here for now
const GraphGenDict = Dict{Union{Int, Edge}, T} where T <: MPolyRingElem

function Base.getindex(D::GraphGenDict, i::Int, j::Int)
  return D[Edge(i, j)]
end

const GraphTransDict = Dict{Tuple{VarName, Edge}, T} where T <: MPolyRingElem

function Base.getindex(D::GraphTransDict, s::VarName, i::Int, j::Int)
  return D[(s, Edge(i, j))]
end

################################################################################

include("CI.jl")
include("Markov.jl")
include("GraphicalModels.jl")
include("GaussianGraphicalModels.jl")
include("DiscreteGraphicalModels.jl")

include("PhylogeneticModels.jl")
include("PhylogeneticModels-functions.jl")
# include("PhylogeneticAuxiliary.jl")
# include("PhylogeneticParametrization.jl")
# include("PhylogeneticInvariants.jl")
include("LoadModels.jl")

include("./serialization.jl")

# export Abstract Graphical Model
export GraphicalModel
export model_ring
export parameter_ring
export parametrization
export graph
export vanishing_ideal
export varnames

#export models
export cavender_farris_neyman_model
export jukes_cantor_model
export kimura2_model
export kimura3_model
export general_markov_model
export affine_phylogenetic_model!

#export phylogenetic models attributes
export phylogenetic_model
export n_states
export transition_matrix
export root_distribution
export entry_transition_matrix
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
export reduced_model_ring

#export auxiliary graph functions
export n_leaves

# export functions to load objects (currently only graphs of phylogenetic models)
export load_phylogenetic_model

# Gaussian Graphical Model Exports
export GaussianGraphicalModel
export covariance_matrix
export gaussian_graphical_model
export directed_edges_matrix, error_covariance_matrix
export concentration_matrix
export gaussian_ring, GaussianRing

# Discrete graphical models
export state_space
export marginal
export markov_ring, tensor_ring, MarkovRing

# Conditional independence
export CIStmt
export ci_stmt
export @CI_str
export ci_statements
export make_elementary
export ci_ideal
export ci_polynomial
export ci_structure
export random_variables

