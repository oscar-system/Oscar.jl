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
import Gap
import Singular

include("OscarTypes.jl")

include("Rings/integer.jl")
include("Rings/rational.jl")

end # module
