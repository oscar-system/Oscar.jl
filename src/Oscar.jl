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

import AbstractAlgebra: @show_name, @show_special, elem_type, force_coerce, force_op,
                        parent_type, expressify, canonical_unit

# More helpful error message for users on Windows.
windows_error() = error("""

    This package unfortunately does not run natively under Windows.
    Please install Julia using Windows subsystem for Linux and try again.
    See also https://oscar.computeralgebra.de/install/.
    """)

if Sys.iswindows()
  windows_error()
end

function __init__()
  if Sys.iswindows()
    windows_error()
  end

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
    println("(c) 2019-2021 by The Oscar Development Team")
  end

  append!(_gap_group_types,
    [
        (GAP.Globals.IsPermGroup, PermGroup),
        (GAP.Globals.IsPcGroup, PcGroup),
        (GAP.Globals.IsMatrixGroup, MatrixGroup),
        (GAP.Globals.IsFpGroup, FPGroup),
    ])
    GAP.Packages.load("forms")
end

# pkgdir was added in Julia 1.4
if VERSION < v"1.4"
   pkgdir(m::Core.Module) = abspath(Base.pathof(Base.moduleroot(m)), "..", "..")
end
pkgproject(m::Core.Module) = Pkg.Operations.read_project(Pkg.Types.projectfile_path(pkgdir(m)))
pkgversion(m::Core.Module) = pkgproject(m).version
const VERSION_NUMBER = pkgversion(@__MODULE__)

const is_dev = (function(m)
        if VERSION >= v"1.4"
          uuid = pkgproject(m).uuid
          deps = Pkg.dependencies()
          if Base.haskey(deps, uuid)
            if deps[uuid].is_tracking_path
              return true
            end
          end
        end
        return occursin("-dev", lowercase(string(VERSION_NUMBER)))
    end)(@__MODULE__)

const IJuliaMime = Union{MIME"text/latex", MIME"text/html"}

const oscardir = pkgdir(Oscar)


function example(s::String)
  Base.include(Main, joinpath(oscardir, "examples", s))
end

function build_doc()
  Base.include(Main, joinpath(oscardir, "docs", "make_local.jl"))
end
export build_doc
# This can be used in
#
# module A
#   using Oscar
#   Oscar.@example("example.jl")
# end
#
# __module__ expands to the module of the call site of the macro.
macro example(s)
  :($(esc(__module__)).include(joinpath(oscardir, "examples", $(esc(s)))))
end

function data(s::String)
  Base.include(Main, joinpath(oscardir, "data", s))
end

function revise(s::String)
  s = joinpath(oscardir, "examples", s)
  Main.Revise.track(Main, s)
end

function system(s::String)
  Base.include(Main, joinpath(oscardir, "system", s))
end

function build()
  system("Build.jl")
end


function test_module(x, new::Bool = true)
   julia_exe = Base.julia_cmd()
   if x == "all"
     test_file = joinpath(oscardir, "test/runtests.jl")
   else
     test_file = joinpath(oscardir, "test/$x.jl")
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

function weights end

include("OscarTypes.jl")

include("Groups/types.jl")

include("Rings/Hecke.jl") #does all the importing from Hecke - to define names

include("GAP/gap_to_oscar.jl")
include("GAP/oscar_to_gap.jl")

include("Groups/group_constructors.jl")
include("Groups/sub.jl")
include("Groups/homomorphisms.jl")
include("Groups/cosets.jl")
include("Groups/libraries/libraries.jl")
include("Groups/GAPGroups.jl")
include("Groups/directproducts.jl")
include("Groups/action.jl")
include("Groups/gsets.jl")
include("Groups/matrices/matrices.jl")

include("Rings/integer.jl")
include("Rings/rational.jl")
include("Rings/mpoly.jl")
include("Rings/mpoly-ideals.jl")
include("Rings/groebner.jl")
include("Rings/MPolyQuo.jl")
include("Rings/mpoly-nested.jl")
include("Rings/FractionalIdeal.jl")
include("Rings/affine-algebra-homs.jl")
include("Rings/mpoly-affine-algebras.jl")
include("Rings/mpoly-graded.jl")
include("Rings/mpoly-local.jl")
include("Rings/FinField.jl")
include("Rings/NumberField.jl")

include("Modules/FreeModules-graded.jl")

include("Geometry/basics.jl")

include("Polymake/Ineq.jl")
include("Polymake/NmbThy.jl")

include("Polytopes/Polytopes.jl")

include("../StraightLinePrograms/src/StraightLinePrograms.jl")
include("Rings/lazypolys.jl")
include("Rings/slpolys.jl")

include("../experimental/Experimental.jl")

if is_dev
#  include("../examples/ModStdNF.jl")
#  include("../examples/ModStdQ.jl")
#  include("../examples/ModStdQt.jl")
  include("../examples/PrimDec.jl")
#  include("../examples/GaloisGrp.jl")

#  include("../examples/PlaneCurve.jl")
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
