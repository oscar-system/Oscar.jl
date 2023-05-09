@doc raw"""
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

  `https://www.oscar-system.org`

OSCAR is licensed under the GPL v3+ (see LICENSE.md).
"""
module Oscar

using Preferences

include("imports.jl")

const cornerstones = String["AbstractAlgebra", "GAP", "Hecke", "Nemo", "Polymake", "Singular"];
const jll_deps = String["Antic_jll", "Arb_jll", "Calcium_jll", "FLINT_jll", "GAP_jll",
                        "libpolymake_julia_jll", "libsingular_julia_jll",
                        "polymake_jll", "Singular_jll"];

# We read experimental and filter out all packages that follow our desired
# scheme. Remember those packages to avoid doing this all over again for docs
# and test.
# We don't want to interfere with existing stuff in experimental though.
const expdir = joinpath(@__DIR__, "../experimental")
const oldexppkgs = [
  "ExteriorAlgebra",
  "GaloisGrp",
  "GITFans",
  "GModule",
  "JuLie",
  "Matrix",
  "ModStd",
  "MPolyRingSparse",
  "Rings",
  "Schemes",
  "SymmetricIntersections",
]
# DEVELOPER OPTION:
# The following lines ensure that ToricSchemes is loaded before FTheoryTools.
# DO NOT USE THIS UNLESS YOU KNOW THE CONSEQUENCES.
# For more background, see https://github.com/oscar-system/Oscar.jl/issues/2300.
const orderedpkgs = [
  "ToricSchemes",
  "FTheoryTools",
]
exppkgs = filter(x->isdir(joinpath(expdir, x)) && !(x in oldexppkgs) && !(x in orderedpkgs), readdir(expdir))
append!(exppkgs, orderedpkgs)

# Error if something is incomplete in experimental
for pkg in exppkgs
  if !isfile(joinpath(expdir, pkg, "src", "$pkg.jl"))
    error("experimental/$pkg is incomplete: $pkg/src/$pkg.jl missing.")
  end
  if !isfile(joinpath(expdir, pkg, "test", "runtests.jl"))
    error("experimental/$pkg is incomplete: $pkg/test/runtests.jl missing.")
  end
end

# force trigger recompile when folder changes
include_dependency("../experimental")

include("utils.jl")

# More helpful error message for users on Windows.
windows_error() = error("""

    This package unfortunately does not run natively under Windows.
    Please install Julia using Windows subsystem for Linux and try again.
    See also https://www.oscar-system.org/install/.
    """)

if Sys.iswindows()
  windows_error()
end

# global seed value for oscar to allow creating deterministic random number
# generators
# we initialize this here with a random value as well to allow use during
# precompilation
const rng_seed = Ref{UInt32}(rand(UInt32))

@doc raw"""
    get_seed()

Return the current random seed that is used for calls to `Oscar.get_seeded_rng`.
"""
get_seed() = return rng_seed[]

@doc raw"""
    set_seed!(s::Integer)

Set a new global seed for all subsequent calls to `Oscar.get_seeded_rng`.
"""
function set_seed!(s::Integer)
  rng_seed[] = convert(UInt32, s)
end

@doc raw"""
    get_seeded_rng()

Return a new random number generator object of type MersenneTwister which is
seeded with the global seed value `Oscar.rng_seed`. This can be used for
the testsuite to get more stable output and running times. Using a separate RNG
object for each testset (or file) makes sure previous uses don't affect later
testcases. It could also be useful for some randomized algorithms.
The seed will be initialized with a random value during OSCAR startup but can
be set to a fixed value with `Oscar.set_seed!` as it is done in `runtests.jl`.
"""
get_seeded_rng() = return MersenneTwister([get_seed()])


function __init__()
  if Sys.iswindows()
    windows_error()
  end

  # initialize random seed
  set_seed!(rand(UInt32))

  if isinteractive() && Base.JLOptions().banner != 0
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
    println("(c) 2019-2023 by The OSCAR Development Team")
  end

  append!(_gap_group_types,
    [
        (GAP.Globals.IsPermGroup, PermGroup),
        (GAP.Globals.IsPcGroup, PcGroup),
        (GAP.Globals.IsMatrixGroup, MatrixGroup),
        (GAP.Globals.IsSubgroupFpGroup, FPGroup),
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
    for pkg in ["ctbllib",
                "forms",
                "wedderga", # provides a function to compute Schur indices
                "repsn",
               ]
      GAP.Packages.load(pkg) || error("cannot load the GAP package $pkg")
    end
    __init_group_libraries()

    add_verbose_scope(:K3Auto)
    add_assert_scope(:K3Auto)

    add_verbose_scope(:GlobalTateModel)
    add_verbose_scope(:GlobalWeierstrassModel)

    add_verbosity_scope(:LinearQuotients)
end

const PROJECT_TOML = Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))
const VERSION_NUMBER = VersionNumber(PROJECT_TOML["version"])

const is_dev = (function()
        uuid = PROJECT_TOML["uuid"]
        deps = Pkg.dependencies()
        if Base.haskey(deps, uuid)
          if deps[uuid].is_tracking_path
            return true
          end
        end
        return occursin("-dev", lowercase(string(VERSION_NUMBER)))
    end)()

const IJuliaMime = Union{MIME"text/latex", MIME"text/html"}

const oscardir = Base.pkgdir(Oscar)
const aadir = Base.pkgdir(AbstractAlgebra)
const nemodir = Base.pkgdir(Nemo)
const heckedir = Base.pkgdir(Hecke)

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


@doc raw"""
    test_module(file::AbstractString; new::Bool = true)

Run the Oscar tests in the file `test/<file>.jl` where `file` may be a path.

The optional parameter `new` takes the values `false` and `true` (default). If
`true`, then the tests are run in a new session, otherwise the currently active
session is used.

For experimental modules, use [`test_experimental_module`](@ref) instead.
"""
function test_module(file::AbstractString; new::Bool=true)
  julia_exe = Base.julia_cmd()
  rel_test_file = normpath("test", "$file.jl")
  test_file = joinpath(oscardir, rel_test_file)

  if new
    cmd = "using Test; using Oscar; Hecke.assertions(true); include(\"$test_file\");"
    @info("spawning ", `$julia_exe -e \"$cmd\"`)
    run(`$julia_exe -e $cmd`)
  else
    @req isdefined(Base.Main, :Test) "You need to do \"using Test\""
    @info("Running tests for $rel_test_file in same session")
    Base.include(Base.Main, test_file)
  end
end

@doc raw"""
    test_experimental_module(project::AbstractString; file::AbstractString="runtests", new::Bool = true)

Run the Oscar tests in the file `experimental/<project>/test/<file>.jl`
where `file` may be a path.
The default is to run the entire test suite of the module `project`.

The optional parameter `new` takes the values `false` and `true` (default). If
`true`, then the tests are run in a new session, otherwise the currently active
session is used.
"""
function test_experimental_module(
  project::AbstractString; file::AbstractString="runtests", new::Bool=true
)
  test_file = "../experimental/$project/test/$file"
  test_module(test_file; new)
end

include("exports.jl")

# HACK/FIXME: remove these aliases once we have them in AA/Nemo/Hecke
@alias characteristic_polynomial charpoly  # FIXME
@alias minimal_polynomial minpoly  # FIXME

include("printing.jl")
include("fallbacks.jl")


include("Groups/Groups.jl")

include("Rings/Rings.jl")
include("GAP/GAP.jl")

include("../gap/OscarInterface/julia/alnuth.jl")

include("Groups/group_characters.jl")  # needs some Rings functionality
include("Groups/action.jl")  # needs some PolynomialRings functionality

include("Modules/Modules.jl")
include("Rings/ReesAlgebra.jl") # Needs ModuleFP

include("NumberTheory/NmbThy.jl")

include("PolyhedralGeometry/PolyhedralGeometry.jl")

include("Polymake/polymake_to_oscar.jl")

include("Combinatorics/Graphs.jl")
include("Combinatorics/SimplicialComplexes.jl")
include("Combinatorics/Matroids/JMatroids.jl")
include("Combinatorics/Matroids/matroid_strata_grassmannian.jl")

include("StraightLinePrograms/StraightLinePrograms.jl")
include("Rings/lazypolys.jl") # uses StraightLinePrograms
include("Rings/slpolys.jl") # uses StraightLinePrograms
include("NumberTheory/GalThy.jl")

include("AlgebraicGeometry/AlgebraicGeometry.jl")

include("InvariantTheory/InvariantTheory.jl")

include("../experimental/Experimental.jl")

include("Rings/binomial_ideals.jl") # uses QQAbModule from experimental/Rings/QQAbAndPChars.jl

if is_dev
#  include("../examples/ModStdNF.jl")
#  include("../examples/ModStdQ.jl")
#  include("../examples/ModStdQt.jl")
  include("../examples/PrimDec.jl")
#  include("../examples/GaloisGrp.jl")

#  include("../examples/PlaneCurve.jl")
end

include("Serialization/main.jl")

include("aliases.jl")

include("deprecations.jl")

const global OSCAR = Oscar
const global oscar = Oscar

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
