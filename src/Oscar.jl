@doc """
Welcome to OSCAR version $(VERSION_NUMBER)

OSCAR is developed by a large group of international collaborators, coordinated
mainly at the University of Kaiserslautern-Landau.

Written in Julia, it combines the well established systems
 * [`Singular`](@ref Singular)
 * [`GAP`](@ref GAP)
 * [`Polymake`](@ref Polymake)
 * [`ANTIC`](@ref ANTIC) (comprising [`Hecke`](@ref Hecke), [`Nemo`](@ref Nemo) and [`AbstractAlgebra`](@ref AbstractAlgebra))
into a comprehensive tool for computational algebra.

  For more information please visit

  `https://www.oscar-system.org`

OSCAR is licensed under the GPL v3+ (see LICENSE.md).
"""
module Oscar

using Preferences
using LazyArtifacts

include("imports.jl")

include("utils/utils.jl")

# More helpful error message for users on Windows.
windows_error() = error("""

    This package unfortunately does not run natively under Windows.
    Please install Julia using Windows subsystem for Linux and try again.
    See also https://www.oscar-system.org/install/.
    """)

if Sys.iswindows()
  windows_error()
end

function _print_banner()
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
  println("(c) 2019-2024 by The OSCAR Development Team")
end



function __init__()
  if Sys.iswindows()
    windows_error()
  end

  # initialize random seed
  set_seed!(rand(UInt32))

  if isinteractive() && Base.JLOptions().banner != 0
    _print_banner()
  end

  append!(_gap_group_types,
    [
        (GAP.Globals.IsPermGroup, PermGroup),
        (GAP.Globals.IsPcGroup, PcGroup),
        (GAP.Globals.IsMatrixGroup, MatrixGroup),
        (GAP.Globals.IsSubgroupFpGroup, FPGroup),
        (GAP.Globals.IsGroupOfAutomorphisms, AutomorphismGroup),
    ])
  __GAP_info_messages_off()
  # make Oscar module accessible from GAP (it may not be available as
  # `Julia.Oscar` if Oscar is loaded indirectly as a package dependency)
  GAP.Globals.BindGlobal(GapObj("Oscar"), Oscar)
  GAP.Globals.SetPackagePath(GAP.Obj("OscarInterface"), GAP.Obj(joinpath(@__DIR__, "..", "gap", "OscarInterface")))
  GAP.Globals.LoadPackage(GAP.Obj("OscarInterface"))
  withenv("TERMINFO_DIRS" => joinpath(GAP.GAP_jll.Readline_jll.Ncurses_jll.find_artifact_dir(), "share", "terminfo")) do
    GAP.Packages.load("browse"; install=true) # needed for all_character_table_names doctest
  end
  # We want newer versions of some GAP packages than the distributed ones.
  # (But we do not complain if the installation fails.)
  for (pkg, version) in [
     ("repsn", "3.1.1"),
     ]
    GAP.Packages.install(pkg, version, interactive = false, quiet = true)
  end
  # We need some GAP packages.
  for pkg in [
     "atlasrep",
     "ctbllib",  # character tables
     "crisp",    # faster normal subgroups, socles, p-socles for finite solvable groups
     "fga",      # dealing with free groups
     "forms",    # bilinear/sesquilinear/quadratic forms
     "primgrp",  # primitive groups library
     "repsn",    # constructing representations of finite groups
     "sla",      # computing with simple Lie algebras
     "smallgrp", # small groups library
     "transgrp", # transitive groups library
     "wedderga", # provides a function to compute Schur indices
     ]
    GAP.Packages.load(pkg) || error("cannot load the GAP package $pkg")
  end
  __init_group_libraries()

  add_verbose_scope(:K3Auto)
  add_assert_scope(:K3Auto)

  add_verbose_scope(:EllipticSurface)
  add_assert_scope(:EllipticSurface)

  add_verbose_scope(:MorphismFromRationalFunctions)
  add_assert_scope(:MorphismFromRationalFunctions)

  add_verbose_scope(:Gluing)
  add_assert_scope(:Gluing)

  add_verbose_scope(:Intersections)
  add_assert_scope(:Intersections)

  add_verbose_scope(:MaximalAssociatedPoints)
  add_assert_scope(:MaximalAssociatedPoints)

  add_verbose_scope(:Divisors)
  add_assert_scope(:Divisors)

  add_verbose_scope(:Blowup)
  add_assert_scope(:Blowup)

  add_verbose_scope(:hilbert)
  add_assert_scope(:hilbert)

  add_verbose_scope(:FTheoryModelPrinter)

  add_verbosity_scope(:LinearQuotients)

  add_assertion_scope(:ZZLatWithIsom)
  add_verbosity_scope(:ZZLatWithIsom)

  # Pkg.is_manifest_current() returns false if the manifest might be out of date
  # (but might return nothing when there is no project_hash)
  if is_dev && VERSION >= v"1.8" && false === (VERSION < v"1.11.0-DEV.1135" ?
      Pkg.is_manifest_current() :
      Pkg.is_manifest_current(dirname(Base.active_project())))
    @warn "Project dependencies might have changed, please run `]up` or `]resolve`."
  end
end

const PROJECT_TOML = Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))
const VERSION_NUMBER = VersionNumber(PROJECT_TOML["version"])
const PROJECT_UUID = UUID(PROJECT_TOML["uuid"])

const is_dev = (function()
        deps = Pkg.dependencies()
        if Base.haskey(deps, PROJECT_UUID)
          if deps[PROJECT_UUID].is_tracking_path
            return true
          end
        end
        return occursin("-dev", lowercase(string(VERSION_NUMBER)))
    end)()

const IJuliaMime = Union{MIME"text/latex", MIME"text/html"}

const oscardir = Base.pkgdir(Oscar)

function example(s::String)
  Base.include(Main, joinpath(oscardir, "examples", s))
end


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

include("assertions.jl")

include("exports.jl")

include("aliases.jl")


include("printing.jl")
include("fallbacks.jl")


include("Rings/Rings.jl")
include("forward_declarations.jl")
include("Groups/Groups.jl")

include("GAP/GAP.jl")

include("../gap/OscarInterface/julia/alnuth.jl")


include("Modules/Modules.jl")
include("Rings/ReesAlgebra.jl") # Needs ModuleFP

include("NumberTheory/NmbThy.jl")

include("Combinatorics/Graphs/structs.jl")
include("PolyhedralGeometry/PolyhedralGeometry.jl")

include("Polymake/polymake_to_oscar.jl")

include("Combinatorics/Graphs/functions.jl")
include("Combinatorics/SimplicialComplexes.jl")
include("Combinatorics/OrderedMultiIndex.jl")
include("Combinatorics/Matroids/JMatroids.jl")
include("Combinatorics/Compositions.jl")
include("Combinatorics/EnumerativeCombinatorics/EnumerativeCombinatorics.jl")
include("Combinatorics/PhylogeneticTrees.jl")

include("StraightLinePrograms/StraightLinePrograms.jl")
include("Rings/lazypolys.jl") # uses StraightLinePrograms
include("Rings/slpolys.jl") # uses StraightLinePrograms
include("NumberTheory/GalThy.jl")

include("AlgebraicGeometry/AlgebraicGeometry.jl")

include("TropicalGeometry/TropicalGeometry.jl")

include("InvariantTheory/InvariantTheory.jl")

# Serialization should always come at the end of Oscar source code
# but before experimental, any experimental serialization should
# be written inside the corresponding experimental code sub directory
include("Serialization/main.jl")

include("../experimental/Experimental.jl")

include("Rings/binomial_ideals.jl") # uses QQAbModule from experimental/Rings/QQAbAndPChars.jl

if is_dev
#  include("../examples/ModStdNF.jl")
#  include("../examples/ModStdQ.jl")
#  include("../examples/ModStdQt.jl")
  include("../examples/PrimDec.jl")
end

include("deprecations.jl")

@doc raw"""
ANTIC is the project name for the number theoretic cornerstone of OSCAR, see
  ?Nemo
  ?Hecke
  ?AbstractAlgebra
  for more information
"""
module ANTIC
end

end # module
