@doc Markdown.doc"""
Welcome to OSCAR version 0.1.1

OSCAR is developed by a large group of international collaborators, coordinated
currently, mainly at the technical university of Kaiseslautern.

Written in Julia, it combines the well established systems
 * Singular
 * Gap
 * Polymake
 * Hecke, Nemo and AbstractAlgebra
into a comprehensive tool for computational algebra.

  For more information please visit

  `https://oscar.computeralgebra.de`

Oscar is licensed under the GLP 3 (see LICENSE.md).
"""
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
# to allow access to the cornerstones! Otherwise, not even import or using from the
# user level will work as none of them will have been "added" by the user.
# possibly all should add a doc string to the module?
export Nemo, Hecke, Singular, Polymake, AbstractAlgebra

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
  println("Type: '?Oscar' for more information")
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
