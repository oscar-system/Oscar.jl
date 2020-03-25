module Oscar
const global VERSION_NUMBER = v"0.1.1"
#=
  We currently only import packages which:
    * are registered
    * have BinaryBuilder or build in a few minutes
    * support Linux and OSX
=#

import AbstractAlgebra
import Nemo
import Hecke
import Singular
import Polymake
using Markdown

import AbstractAlgebra: @show_name, @show_special

function __init__()
  println(" -----    -----    -----      -      -----   ")
  println("|     |  |     |  |     |    | |    |     |  ")
  println("|     |  |        |         |   |   |     |  ")
  println("|     |   -----   |        |     |  |-----   ")
  println("|     |        |  |        |-----|  |   |    ")
  println("|     |  |     |  |     |  |     |  |    |   ")
  println(" -----    -----    -----   -     -  -     -  ")
  println()
  println("...combining (and extending) GAP, Hecke, Nemo, Polymake and Singular")
  print("Version")
  printstyled(" $VERSION_NUMBER ", color = :green)
  print("... \n ... which comes with absolutely no warranty whatsoever")
  println()
  println("(c) 2019-2020 by The Oscar Development Team")
end


include("OscarTypes.jl")

include("Rings/integer.jl")
include("Rings/rational.jl")
include("Rings/Hecke.jl")
include("Rings/mpoly.jl")
include("Rings/mpoly-graded.jl")
include("Modules/FreeModules-graded.jl")
include("Polymake/Ineq.jl")
include("Polymake/NmbThy.jl")

end # module
