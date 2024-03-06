module BasisLieHighestWeight

using ..Oscar
using ..Oscar: GAPWrap
using ..Oscar: IntegerUnion
using ..Oscar: _is_weighted

using AbstractAlgebra.PrettyPrinting

import Oscar: dim
import Oscar: dim_of_simple_module
import Oscar: monomial_ordering
import Oscar: monomials

import Base: length

# Long-term TODO's:
# - Use Oscar-Lie-Algebra type instead of LieAlgebraStructure
# - Test if ZZx should be a graded_polynomial_ring with weights_w as weights
# - Maybe export and docstring: 
#   - get_dim_weightspace
#   - orbit_weylgroup
#   - get_lattice_points_of_weightspace
#   - convert_lattice_points_to_monomials
#   - convert_monomials_to_lattice_points
#   - action_matrices_of_operators
#   - weights_for_operators
# - Data-Type for weights of Lie-Algebras? Two types, in alpha_i and w_i, conversion is defined in RootConversion
# - the list of Minkowski gens contains too many elements, only include those that give us something new

include("LieAlgebras.jl")
include("BirationalSequence.jl")
include("MonomialBasis.jl")
include("NewMonomial.jl")
include("TensorModels.jl")
include("MonomialOrder.jl")
include("RootConversion.jl")
include("WeylPolytope.jl")
include("MainAlgorithm.jl")
include("UserFunctions.jl")

export MonomialBasis

export basis_coordinate_ring_kodaira
export basis_coordinate_ring_kodaira_ffl
export basis_lie_highest_weight_operators
export basis_lie_highest_weight
export basis_lie_highest_weight_ffl
export basis_lie_highest_weight_lusztig
export basis_lie_highest_weight_nz
export basis_lie_highest_weight_string

end

using .BasisLieHighestWeight

export BasisLieHighestWeight

export MonomialBasis

export basis_coordinate_ring_kodaira
export basis_coordinate_ring_kodaira_ffl
export basis_lie_highest_weight_operators
export basis_lie_highest_weight
export basis_lie_highest_weight_ffl
export basis_lie_highest_weight_lusztig
export basis_lie_highest_weight_nz
export basis_lie_highest_weight_string
