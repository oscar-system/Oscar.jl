module Oscar

#=
  We currently only import packages which:
    * are registered
    * have BinaryBuilder (including for dependencies)
    * support Linux and OSX
  More to come soon....
=#

import AbstractAlgebra
import Nemo
import Hecke
import Polymake

include("OscarTypes.jl")

include("Rings/integer.jl")
include("Rings/rational.jl")

end # module
