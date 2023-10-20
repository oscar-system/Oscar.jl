module BasisLieHighestWeight

using ..Oscar
using ..Oscar: GAPWrap, IntegerUnion

using AbstractAlgebra.PrettyPrinting

using Polymake

import Oscar: dim, monomial_ordering, monomials

import Base: length

# TODO basis_lie_highest_weight_lustzig
# TODO basis_lie_highest_weight_string
# TODO basis_lie_highest_weight_feigin_fflv
# TODO basis_lie_highest_weight_nZ

# TODO (?) Maybe export and docstring: 
# get_dim_weightspace
# orbit_weylgroup
# get_lattice_points_of_weightspace
# convert_lattice_points_to_monomials
# convert_monomials_to_lattice_points
# tensorMatricesForOperators
# weights_for_operators

# TODO Use Oscar-Lie-Algebra type instead of LieAlgebra
# TODO Data-Type for weights of Lie-Algebras? Three types, in alpha_i, w_i and eps_i, conversion is defined in RootConversion
# w_to_eps
# eps_to_w
# alpha_to_eps
# eps_to_alpha
# w_to_aplha
# alpha_to_w

# TODO GAPWrap-wrappers are missing for 
# ChevalleyBasis
# DimensionOfHighestWeightModule
# SimpleLieAlgebra
# Rationals
# HighestWeightModule
# List
# MatrixOfAction
# RootSystem
# CartanMatrix
# WeylGroup
# DominantCharacter
# DimensionOfHighestWeightModule
# CanonicalGenerators

include("LieAlgebras.jl")
include("BirationalSequence.jl")
include("MonomialBasis.jl")
include("VectorSpaceBases.jl")
include("NewMonomial.jl")
include("TensorModels.jl")
include("MonomialOrder.jl")
include("RootConversion.jl")
include("WeylPolytope.jl")
include("DemazureOperators.jl")
include("WordCalculations.jl")
include("MainAlgorithm.jl")

end

export BasisLieHighestWeight
