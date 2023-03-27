import AbstractAlgebra: FPModule, FPModuleElem
import AbstractAlgebra: ProductIterator

import Combinatorics

export AbstractLieAlgebra, AbstractLieAlgebraElem
export LieAlgebra, LieAlgebraElem
export LieAlgebraAbstractModule, LieAlgebraAbstractModuleElem
export LieAlgebraExteriorPowerModule, LieAlgebraExteriorPowerModuleElem
export LieAlgebraModule, LieAlgebraModuleElem
export LieAlgebraStdModule, LieAlgebraStdModuleElem
export LieAlgebraSymmetricPowerModule, LieAlgebraSymmetricPowerModuleElem
export LieAlgebraTensorPowerModule, LieAlgebraTensorPowerModuleElem
export LinearLieAlgebra, LinearLieAlgebraElem

export abstract_module
export base_liealgebra
export bracket
export coefficient_vector
export exterior_power
export general_linear_liealgebra
export highest_weight_module
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
include("PBWDeformations/LieAlgebraModules/LieAlgebraModule.jl")
include("PBWDeformations/LieAlgebraModules/LieAlgebraAbstractModule.jl")
include("PBWDeformations/LieAlgebraModules/LieAlgebraExteriorPowerModule.jl")
include("PBWDeformations/LieAlgebraModules/LieAlgebraStdModule.jl")
include("PBWDeformations/LieAlgebraModules/LieAlgebraSymmetricPowerModule.jl")
include("PBWDeformations/LieAlgebraModules/LieAlgebraTensorPowerModule.jl")
include("PBWDeformations/GapWrapper.jl")
