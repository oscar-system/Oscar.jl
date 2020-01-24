module Oscar

#=
  We currently only import packages which:
    * are registered
    * have BinaryBuilder or build in a few minutes
    * support Linux and OSX
=#

import AbstractAlgebra
import Nemo
import Hecke
import Polymake
#import Singular

include("OscarTypes.jl")

include("Rings/integer.jl")
include("Rings/rational.jl")
include("Rings/Hecke.jl")
#include("Rings/mpoly.jl")
#include("Rings/mpoly-graded.jl")
#include("Modules/FreeModules-graded.jl")

end # module
