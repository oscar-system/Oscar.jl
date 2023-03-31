import AbstractAlgebra: ProductIterator

export AbstractLieAlgebra, AbstractLieAlgebraElem
export LieAlgebra, LieAlgebraElem
export LieAlgebraModule, LieAlgebraModuleElem
export LinearLieAlgebra, LinearLieAlgebraElem

export abstract_module
export base_lie_algebra
export bracket
export coefficient_vector
export combinations
export exterior_power
export general_linear_lie_algebra
export highest_weight_module
export is_exterior_power
export is_standard_module
export is_symmetric_power
export is_tensor_power
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

include("PBWDeformations/Combinatorics.jl")
include("PBWDeformations/Util.jl")
include("PBWDeformations/LieAlgebras/LieAlgebra.jl")
include("PBWDeformations/LieAlgebras/AbstractLieAlgebra.jl")
include("PBWDeformations/LieAlgebras/LinearLieAlgebra.jl")
include("PBWDeformations/LieAlgebraModule.jl")
include("PBWDeformations/GapWrapper.jl")
