import AbstractAlgebra: ProductIterator

import Combinatorics

export AbstractLieAlgebra, AbstractLieAlgebraElem
export LieAlgebra, LieAlgebraElem
export LieAlgebraModule, LieAlgebraModuleElem
export LinearLieAlgebra, LinearLieAlgebraElem

export abstract_module
export base_liealgebra
export bracket
export coefficient_vector
export exterior_power
export general_linear_liealgebra
export highest_weight_module
export is_exterior_power
export is_standard_module
export is_symmetric_power
export is_tensor_power
export liealgebra
export matrix_repr_basis
export special_linear_liealgebra
export special_orthogonal_liealgebra
export standard_module
export symmetric_power
export tensor_power

include("PBWDeformations/Util.jl")
include("PBWDeformations/LieAlgebras/LieAlgebra.jl")
include("PBWDeformations/LieAlgebras/AbstractLieAlgebra.jl")
include("PBWDeformations/LieAlgebras/LinearLieAlgebra.jl")
include("PBWDeformations/LieAlgebraModule.jl")
include("PBWDeformations/GapWrapper.jl")
