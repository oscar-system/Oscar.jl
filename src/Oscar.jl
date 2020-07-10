@doc Markdown.doc"""
Welcome to OSCAR version $(VERSION_NUMBER)

OSCAR is developed by a large group of international collaborators, coordinated
currently, mainly at the Technische Universität Kaiserslautern.

Written in Julia, it combines the well established systems
 * [`Singular`](@ref Singular)
 * [`GAP`](@ref GAP)
 * [`Polymake`](@ref Polymake)
 * [`ANTIC`](@ref ANTIC) (comprising [`Hecke`](@ref Hecke), [`Nemo`](@ref Nemo) and [`AbstractAlgebra`](@ref AbstractAlgebra))
into a comprehensive tool for computational algebra.

  For more information please visit

  `https://oscar.computeralgebra.de`

Oscar is licensed under the GPL v3+ (see LICENSE.md).
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
using Markdown, Test, Requires
# to allow access to the cornerstones! Otherwise, not even import or using from the
# user level will work as none of them will have been "added" by the user.
# possibly all should add a doc string to the module?
export Nemo, Hecke, Singular, Polymake, AbstractAlgebra, GAP

import AbstractAlgebra: @show_name, @show_special, force_coerce, force_op

function __init__()
  if isinteractive()
    println(" -----    -----    -----      -      -----   ")
    println("|     |  |     |  |     |    | |    |     |  ")
    println("|     |  |        |         |   |   |     |  ")
    println("|     |   -----   |        |     |  |-----   ")
    println("|     |        |  |        |-----|  |   |    ")
    println("|     |  |     |  |     |  |     |  |    |   ")
    println(" -----    -----    -----   -     -  -     -  ")
    println()
    println("...combining (and extending) ANTIC, GAP, Polymake and Singular")
    print("Version")
    printstyled(" $VERSION_NUMBER ", color = :green)
    print("... \n ... which comes with absolutely no warranty whatsoever")
    println()
    println("Type: '?Oscar' for more information")
    println("(c) 2019-2020 by The Oscar Development Team")
  end

  append!(_gap_group_types,
    [
        (GAP.Globals.IsPermGroup, PermGroup),
        (GAP.Globals.IsPcGroup, PcGroup),
        (GAP.Globals.IsMatrixGroup, MatrixGroup),
        (GAP.Globals.IsFpGroup, FPGroup),
    ])
end

is_dev = false

if VERSION >= v"1.4"
  deps = Pkg.dependencies()
  if Base.haskey(deps, Base.UUID("f1435218-dba5-11e9-1e4d-f1a5fab5fc13"))
    ver = Pkg.dependencies()[Base.UUID("f1435218-dba5-11e9-1e4d-f1a5fab5fc13")]
    if occursin("/dev/", ver.source)
      global VERSION_NUMBER = VersionNumber("$(ver.version)-dev")
      global is_dev = true
    else
      global VERSION_NUMBER = VersionNumber("$(ver.version)")
      if ver.git_revision !== nothing && occursin("master", ver.git_revision)
        is_dev = true
      end
    end
  else
    global VERSION_NUMBER = "not installed"
  end
else
  deps = Pkg.API.__installed(Pkg.PKGMODE_MANIFEST) #to also get installed dependencies
  if haskey(deps, "Oscar")
    ver = deps["Oscar"]
    dir = dirname(@__DIR__)
    if occursin("/dev/", dir)
      global VERSION_NUMBER = VersionNumber("$(ver)-dev")
      is_dev = true
    else
      global VERSION_NUMBER = VersionNumber("$(ver)")
    end
  else
    global VERSION_NUMBER = "not installed"
  end
end

const IJuliaMime = Union{MIME"text/latex", MIME"text/html"}

const pkgdir = joinpath(dirname(pathof(Oscar)), "..")


function example(s::String)
  Base.include(Main, joinpath(dirname(pathof(Oscar)), "..", "examples", s))
end

function data(s::String)
  Base.include(Main, joinpath(dirname(pathof(Oscar)), "..", "data", s))
end

function revise(s::String)
  s = joinpath(dirname(pathof(Oscar)), "..", "examples", s)
  Main.Revise.track(Main, s)
end

function system(s::String)
  Base.include(Main, joinpath(dirname(pathof(Oscar)), "..", "system", s))
end

function build()
  system("Build.jl")
end


function test_module(x, new::Bool = true)
   julia_exe = Base.julia_cmd()
   if x == "all"
     test_file = joinpath(pkgdir, "test/runtests.jl")
   else
     test_file = joinpath(pkgdir, "test/$x.jl")
   end

   if new
     cmd = "using Test; using Oscar; Hecke.assertions(true); include(\"$test_file\");"
     @info("spawning ", `$julia_exe -e \"$cmd\"`)
     run(`$julia_exe -e $cmd`)
   else
     @info("Running tests for $x in same session")
     try
       include(test_file)
     catch e
       if isa(e, LoadError)
         println("You need to do \"using Test\"")
       else
         rethrow(e)
       end
     end
   end
end


include("OscarTypes.jl")

include("Groups/types.jl")
include("Groups/group_constructors.jl")
include("Groups/sub.jl")
include("Groups/homomorphisms.jl")
include("Groups/cosets.jl")
include("Groups/gsets.jl")
include("Groups/libraries/libraries.jl")
include("Groups/GAPGroups.jl")
include("Groups/directproducts.jl")
include("Groups/action.jl")

include("Rings/integer.jl")
include("Rings/rational.jl")
include("Rings/Hecke.jl")
include("Rings/mpoly.jl")
include("Rings/mpoly-graded.jl")
include("Modules/FreeModules-graded.jl")
include("Polymake/Ineq.jl")
include("Polymake/NmbThy.jl")

include("../StraightLinePrograms/src/StraightLinePrograms.jl")

if is_dev
  include("../examples/ModStdNF.jl")
  include("../examples/PrimDec.jl")
  include("../examples/GaloisGrp.jl")
end

const global OSCAR = Oscar
const global oscar = Oscar

@doc Markdown.doc"""
ANTIC is the project name for the number theoretic cornerstone of OSCAR, see
  ?Nemo
  ?Hecke
  ?AbstractAlgebra
  for more information
"""
module ANTIC
using Markdown
end
export ANTIC

export OSCAR, oscar

end # module
