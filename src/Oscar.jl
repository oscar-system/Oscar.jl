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

using LazyArtifacts

include("imports.jl")

AbstractAlgebra.@include_deprecated_bindings()
Nemo.@include_deprecated_bindings()
Hecke.@include_deprecated_bindings()

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

function _print_banner(;is_dev = Oscar.is_dev)
  # lets assemble a version string for the banner
  version_string = string(VERSION_NUMBER)
  if is_dev
    gitinfo = _get_oscar_git_info()
    version_string = version_string * " #$(gitinfo[:branch]) $(gitinfo[:commit][1:7]) $(gitinfo[:date][1:10])"
  else
    version_string = "Version " * version_string
  end

  if displaysize(stdout)[2] >= 80 
    println(
      raw"""  ___   ___   ___    _    ____
             / _ \ / __\ / __\  / \  |  _ \  | Combining and extending ANTIC, GAP,
            | |_| |\__ \| |__  / ^ \ |  ´ /  | Polymake and Singular
             \___/ \___/ \___//_/ \_\|_|\_\  | Type "?Oscar" for more information""")
    printstyled(raw"""o--------o-----o-----o--------o""", color = :yellow)
    println(raw"""  | Documentation: https://docs.oscar-system.org""")
    println(raw"""  S Y M B O L I C   T O O L S    | """ * version_string)
  else
    println("OSCAR $VERSION_NUMBER  https://docs.oscar-system.org  Type \"?Oscar\" for help")
  end
end

function __init__()
  if Sys.iswindows()
    windows_error()
  end

  # initialize random seed
  set_seed!(rand(UInt32))

  if AbstractAlgebra.should_show_banner()
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
  # make Oscar module accessible from GAP (it may not be available as
  # `Julia.Oscar` if Oscar is loaded indirectly as a package dependency)
  GAP.Globals.BindGlobal(GapObj("Oscar_jl"), Oscar)

  # Add the directory `gap/` as a GAP root path, so that GAP can find the
  # OscarInterface package and we can overload specific GAP library files.
  GAP.Globals.ExtendRootDirectories(GapObj([abspath(joinpath(@__DIR__, "..", "gap"))]; recursive = true))

  # We need some GAP packages (currently with unspecified versions).
  for pkg in [
     "OscarInterface", # contains all GAP code that is part of Oscar
     "atlasrep",
     "browse",   # needed for all_character_table_names doctest
     "ctbllib",  # character tables
     "crisp",    # faster normal subgroups, socles, p-socles for finite solvable groups
     "ferret",   # backtrack in permutation groups
     "fga",      # dealing with free groups
     "forms",    # bilinear/sesquilinear/quadratic forms
     "packagemanager", # has been loaded already by GAP.jl
     "polycyclic", # needed for Oscar's pc groups
     "primgrp",  # primitive groups library
     "recog",    # group recognition
     "repsn",    # constructing representations of finite groups
     "smallgrp", # small groups library
     "transgrp", # transitive groups library
     "wedderga", # provides a function to compute Schur indices
     ]
    GAP.Packages.load(pkg) || error("cannot load the GAP package $pkg")
  end

  # Switch off GAP's info messages,
  # also those that are triggered from GAP packages.
  __GAP_info_messages_off()
  __init_group_libraries()

  add_verbosity_scope(:K3Auto)
  add_assertion_scope(:K3Auto)
  
  add_verbosity_scope(:EnriquesAuto)
  add_assertion_scope(:EnriquesAuto)

  add_verbosity_scope(:EllipticSurface)
  add_assertion_scope(:EllipticSurface)

  add_verbosity_scope(:MorphismFromRationalFunctions)
  add_assertion_scope(:MorphismFromRationalFunctions)

  add_verbosity_scope(:Gluing)
  add_assertion_scope(:Gluing)

  add_verbosity_scope(:Intersections)
  add_assertion_scope(:Intersections)

  add_verbosity_scope(:MaximalAssociatedPoints)
  add_assertion_scope(:MaximalAssociatedPoints)

  add_verbosity_scope(:Divisors)
  add_assertion_scope(:Divisors)

  add_verbosity_scope(:Blowup)
  add_assertion_scope(:Blowup)

  add_verbosity_scope(:hilbert)
  add_assertion_scope(:hilbert)

  add_verbosity_scope(:BasisLieHighestWeight)

  add_verbosity_scope(:FTheoryModelPrinter)

  add_verbosity_scope(:LinearQuotients)

  add_assertion_scope(:ZZLatWithIsom)
  add_verbosity_scope(:ZZLatWithIsom)
  
  add_assertion_scope(:Vinberg)
  add_verbosity_scope(:Vinberg)

  add_assertion_scope(:IdealSheaves)
  add_verbosity_scope(:IdealSheaves)

  add_verbosity_scope(:SchurIndices)

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

include("../gap/pkg/OscarInterface/julia/alnuth.jl")


include("Modules/Modules.jl")
include("Rings/ReesAlgebra.jl") # Needs ModuleFP

include("NumberTheory/NmbThy.jl")
include("NumberTheory/vinberg.jl")

include("Combinatorics/Graphs/structs.jl")
include("PolyhedralGeometry/PolyhedralGeometry.jl")

include("Polymake/polymake_to_oscar.jl")

include("Combinatorics/Graphs/functions.jl")
include("Combinatorics/SimplicialComplexes.jl")
include("Combinatorics/OrderedMultiIndex.jl")
include("Combinatorics/Matroids/JMatroids.jl")
include("Combinatorics/EnumerativeCombinatorics/EnumerativeCombinatorics.jl")

include("PolyhedralGeometry/visualization.jl") # needs SimplicialComplex

include("Combinatorics/PhylogeneticTrees.jl")

include("StraightLinePrograms/StraightLinePrograms.jl")
include("Rings/lazypolys.jl") # uses StraightLinePrograms
include("Rings/slpolys.jl") # uses StraightLinePrograms
include("NumberTheory/GalThy.jl")

include("AlgebraicGeometry/AlgebraicGeometry.jl")

include("TropicalGeometry/TropicalGeometry.jl")

include("InvariantTheory/InvariantTheory.jl")

include("LieTheory/LieTheory.jl")

include("Misc/Misc.jl")

# Serialization should always come at the end of Oscar source code
# but before experimental, any experimental serialization should
# be written inside the corresponding experimental code sub directory
include("Serialization/main.jl")

include("../experimental/Experimental.jl")

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
