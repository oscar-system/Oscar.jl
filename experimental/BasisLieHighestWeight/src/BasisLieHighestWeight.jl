module BasisLieHighestWeight

using ..Oscar
using ..Oscar: IntegerUnion
using ..Oscar: _is_weighted
using ..Oscar: _root_system_type_string

using Oscar.LieAlgebras: lie_algebra_simple_module_struct_consts_gap

using AbstractAlgebra.PrettyPrinting

import Oscar: base_lie_algebra
import Oscar: character
import Oscar: dim
import Oscar: monomial_ordering
import Oscar: monomials
import Oscar: root_system
import Oscar: vector_space_dim

import Base: length

# Long-term TODO's:
# - Test if ZZx should be a graded_polynomial_ring with weights_w as weights
# - Maybe export and docstring: 
#   - get_lattice_points_of_weightspace
#   - convert_lattice_points_to_monomials
#   - convert_monomials_to_lattice_points
#   - action_matrices_of_operators
#   - weights_for_operators
# - the list of Minkowski gens contains too many elements, only include those that give us something new

include("LieAlgebras.jl")
include("ModuleData.jl")
include("BirationalSequence.jl")
include("MonomialBasis.jl")
include("NewMonomial.jl")
include("TensorModels.jl")
include("MonomialOrder.jl")
include("WeylPolytope.jl")
include("MainAlgorithm.jl")
include("UserFunctions.jl")

export MonomialBasis

export birational_sequence
export basis_coordinate_ring_kodaira
export basis_coordinate_ring_kodaira_ffl
export basis_lie_highest_weight_operators
export basis_lie_highest_weight
export basis_lie_highest_weight_ffl
export basis_lie_highest_weight_lusztig
export basis_lie_highest_weight_nz
export basis_lie_highest_weight_string
export basis_lie_demazure
export basis_lie_demazure_ffl
export basis_lie_demazure_lusztig
export basis_lie_demazure_nz
export basis_lie_demazure_string

end

using .BasisLieHighestWeight

export BasisLieHighestWeight

export MonomialBasis

export birational_sequence
export basis_coordinate_ring_kodaira
export basis_coordinate_ring_kodaira_ffl
export basis_lie_highest_weight_operators
export basis_lie_highest_weight
export basis_lie_highest_weight_ffl
export basis_lie_highest_weight_lusztig
export basis_lie_highest_weight_nz
export basis_lie_highest_weight_string
export basis_lie_demazure
