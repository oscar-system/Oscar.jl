@doc Markdown.doc"""
Welcome to OSCAR version $(VERSION_NUMBER)

OSCAR is developed by a large group of international collaborators, coordinated
currently, mainly at the Technische Universität Kaiserslautern.

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
import GAP
import Pkg
using Markdown
# to allow access to the cornerstones! Otherwise, not even import or using from the
# user level will work as none of them will have been "added" by the user.
# possibly all should add a doc string to the module?
export Nemo, Hecke, Singular, Polymake, AbstractAlgebra, GAP

import AbstractAlgebra: @show_name, @show_special, force_coerce, force_op

import Base.==
import Base.rand
import Base.conj
import Base.show
import Base.conj!
import Base.parent
import Base.eltype
import Base.iterate
import Base.collect
import GAP.GapObj

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

  append!(_gap_group_types,
    [
        (GAP.Globals.IsPermGroup, PermGroup),
        (GAP.Globals.IsPcGroup, PcGroup),
        (GAP.Globals.IsMatrixGroup, MatrixGroup),
        (GAP.Globals.IsFpGroup, FPGroup),
    ])
end

if VERSION >= v"1.4"
  deps = Pkg.dependencies()
  if Base.haskey(deps, Base.UUID("f1435218-dba5-11e9-1e4d-f1a5fab5fc13"))
    ver = Pkg.dependencies()[Base.UUID("f1435218-dba5-11e9-1e4d-f1a5fab5fc13")]
    if occursin("/dev/", ver.source)
      global VERSION_NUMBER = VersionNumber("$(ver.version)-dev")
    else
      global VERSION_NUMBER = VersionNumber("$(ver.version)")
    end
  else
    global VERSION_NUMBER = "building"
  end
else
  ver = Pkg.installed()["Oscar"]
  dir = dirname(@__DIR__)
  if occursin("/dev/", dir)
    global VERSION_NUMBER = VersionNumber("$(ver)-dev")
  else
    global VERSION_NUMBER = VersionNumber("$(ver)")
  end
end

include("OscarTypes.jl")

include("Groups/types.jl")
include("Groups/group_constructors.jl")
include("Groups/sub.jl")
include("Groups/cosets.jl")
include("Groups/GAPGroups.jl")

include("Rings/integer.jl")
include("Rings/rational.jl")
include("Rings/Hecke.jl")
include("Rings/mpoly.jl")
include("Rings/mpoly-graded.jl")
include("Modules/FreeModules-graded.jl")
include("Polymake/Ineq.jl")
include("Polymake/NmbThy.jl")
const global OSCAR = Oscar
const global oscar = Oscar
export OSCAR, oscar

end # module
