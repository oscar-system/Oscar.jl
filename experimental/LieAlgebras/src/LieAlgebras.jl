module LieAlgebras

using ..Oscar

import Oscar: GAPWrap

# not importet in Oscar
using AbstractAlgebra: CacheDictType, ProductIterator, get_cached!

using Base: deepcopy_internal

# functions with new methods
import ..Oscar:
  action,
  base_ring,
  basis,
  coeff,
  coefficients,
  dim,
  direct_sum,
  dual,
  elem_type,
  expressify,
  exterior_power,
  gen,
  gens,
  ngens,
  parent_type,
  symbols,
  symmetric_power,
  tensor_product,
  ⊕,
  ⊗

import Base: getindex, iszero, parent, zero

export AbstractLieAlgebra, AbstractLieAlgebraElem
export LieAlgebra, LieAlgebraElem
export LieAlgebraModule, LieAlgebraModuleElem
export LinearLieAlgebra, LinearLieAlgebraElem

export abstract_module
export base_lie_algebra
export base_module
export base_modules
export bracket
export coefficient_vector
export combinations
export exterior_power
export general_linear_lie_algebra
export highest_weight_module
export is_direct_sum
export is_dual
export is_exterior_power
export is_standard_module
export is_symmetric_power
export is_tensor_power
export is_tensor_product
export lie_algebra
export matrix_repr_basis
export multicombinations
export permutations
export permutations_with_sign
export special_linear_lie_algebra
export special_orthogonal_lie_algebra
export standard_module
export symmetric_power
export tensor_power

include("Combinatorics.jl")
include("Util.jl")
include("LieAlgebra.jl")
include("AbstractLieAlgebra.jl")
include("LinearLieAlgebra.jl")
include("LieAlgebraModule.jl")
include("GapWrapper.jl")

end

using .LieAlgebras

export AbstractLieAlgebra, AbstractLieAlgebraElem
export LieAlgebra, LieAlgebraElem
export LieAlgebraModule, LieAlgebraModuleElem
export LinearLieAlgebra, LinearLieAlgebraElem

export abstract_module
export base_lie_algebra
export base_module
export base_modules
export bracket
export exterior_power
export general_linear_lie_algebra
export highest_weight_module
export is_direct_sum
export is_dual
export is_exterior_power
export is_standard_module
export is_symmetric_power
export is_tensor_power
export is_tensor_product
export lie_algebra
export matrix_repr_basis
export special_linear_lie_algebra
export special_orthogonal_lie_algebra
export standard_module
export symmetric_power
export tensor_power
