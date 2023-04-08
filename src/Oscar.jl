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
const exppkgs = filter(x->isdir(joinpath(expdir, x)) && !(x in oldexppkgs), readdir(expdir))

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


# When a specific branch is loaded via `]add Package#branch` julia will only
# create a checkout and keep a bare git repo in a separate directory.
# In a bare repo HEAD will not point to the correct commit so we use the git
# tree-hash that Pkg.jl provides and manually map this to a corresponding
# commit.
function _lookup_commit_from_cache(url::AbstractString, tree::AbstractString)
   if Sys.which("git") != nothing
      try
         path = Pkg.Types.add_repo_cache_path(url)
         if isdir(path)
            commit = readchomp(`sh -c "git -C $path log --oneline --all --pretty='tree %T;%H' | grep \"^tree $tree\" | cut -d\; -f2 | head -n1"`)
            return readchomp(`git -C $path show -s --format=", %h -- %ci" $commit`)
         end
      catch
      end
   end
   return ""
end

function _lookup_git_branch(dir::AbstractString; commit=false)
   info = ""
   if Sys.which("git") != nothing &&
         isdir(joinpath(dir,".git"))
      try
         ref = readchomp(`git -C $dir rev-parse --abbrev-ref HEAD`)
         info = " - #$(ref)"
         if commit
            c = readchomp(`git -C $dir show -s --format="%h -- %ci" HEAD`)
            info = "$info, $c"
         end
      catch
      end
   end
   return info
end

function _deps_git_info(dep::Pkg.API.PackageInfo; commit=false)
   if dep.is_tracking_repo
      info = commit ? _lookup_commit_from_cache(dep.git_source, dep.tree_hash) : ""
      return " - #$(dep.git_revision)$info"
   elseif dep.is_tracking_path
      return _lookup_git_branch(dep.source; commit=commit)
   end
   return ""
end

function _print_dependency_versions(io::IO, deps::AbstractArray{<:AbstractString}; padding="    ", suffix="", branch=false, commit=false)
   width = maximum(length.(deps))+length(suffix)+2
   deps = filter(d->d.name in deps, collect(values(Pkg.dependencies())))
   deps = sort!(deps; by=x->x.name)
   for dep in deps
      print(io, "$(padding)$(rpad(dep.name*suffix, width, ' ')) v$(dep.version)")
      println(io, branch ? _deps_git_info(dep; commit=commit) : "")
   end
end

@doc raw"""
    Oscar.versioninfo(io::IO=stdout; branch=false, jll=false, julia=false, commit=false, full=false)

Print the versions of all Oscar-related dependencies.

# Arguments
- `branch::Bool=false`: include git branch name in the output
- `commit::Bool=false`: include git commit hash and date where applicable
- `jll::Bool=false`   : include binary packages (jll) in the output
- `julia::Bool=false` : include julia `versioninfo` output
- `full::Bool=false`  : include all of the above
"""
function versioninfo(io::IO=stdout; branch=false, jll=false, julia=false, commit=false, full=false)
   if full
      branch = jll = julia = commit = true
   end
   print(io, "OSCAR version $(VERSION_NUMBER)")
   println(io, branch ? _lookup_git_branch(dirname(@__DIR__); commit=commit) : "")
   println(io, "  combining:")
   _print_dependency_versions(io, cornerstones; suffix=".jl", branch=branch, commit=commit)
   if jll
      println(io, "  building on:")
      _print_dependency_versions(io, jll_deps; branch=branch, commit=commit)
      println(io, "See `]st -m` for a full list of dependencies.")
   end
   if julia
      println(io, "")
      Main.InteractiveUtils.versioninfo(io)
      println(io, Base.TAGGED_RELEASE_BANNER)
   end
end

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
    GAP.Packages.load("ctbllib")
    GAP.Packages.load("forms")
    GAP.Packages.load("wedderga") # provides a function to compute Schur indices
    GAP.Packages.load("repsn")
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

# use tempdir by default to ensure a clean manifest (and avoid modifying the project)
function doc_init(;path=mktempdir())
  global docsproject = path
  if !isfile(joinpath(docsproject,"Project.toml"))
    cp(joinpath(oscardir, "docs", "Project.toml"), joinpath(docsproject,"Project.toml"))
  end
  Pkg.activate(docsproject) do
    # we dev all packages with the paths from where they are currently loaded
    for dir in [aadir, nemodir, heckedir, oscardir]
      Pkg.develop(path=dir)
    end
    Pkg.instantiate()
    Base.include(Main, joinpath(oscardir, "docs", "make_work.jl"))
  end
end

#function doc_update_deps()
#  Pkg.activate(Pkg.update, joinpath(oscardir, "docs"))
#end

function open_doc()
    filename = normpath(Oscar.oscardir, "docs", "build", "index.html")
    @static if Sys.isapple()
        run(`open $(filename)`; wait = false)
    elseif Sys.islinux() || Sys.isbsd()
        run(`xdg-open $(filename)`; wait = false)
    elseif Sys.iswindows()
        cmd = get(ENV, "COMSPEC", "cmd.exe")
        run(`$(cmd) /c start $(filename)`; wait = false)
    else
        @warn("Opening files the default application is not supported on this OS.",
              KERNEL = Sys.KERNEL)
    end
end


@doc raw"""
    build_doc(; doctest=false, strict=false, open_browser=true)

Build the manual of `Oscar.jl` locally and open the front page in a
browser.

The optional parameter `doctest` can take three values:
  - `false`: Do not run the doctests (default).
  - `true`: Run the doctests and report errors.
  - `:fix`: Run the doctests and replace the output in the manual with
    the output produced by Oscar. Please use this option carefully.

In github actions the Julia version used for building the manual and
running the doctests is 1.6. Using a different Julia version will produce
errors in some parts of Oscar, so please be careful, especially when setting
`doctest=:fix`.

The optional parameter `strict` is passed on to `makedocs` of `Documenter.jl`
and if set to `true` then according to the manual of `Documenter.jl` "a
doctesting error will always make makedocs throw an error in this mode".

To prevent the opening of the browser at the end, set the optional parameter
`open_browser` to `false`.

When working on the manual the `Revise` package can significantly sped
up running `build_doc`. First, install `Revise` in the following way:
```
using Pkg ; Pkg.add("Revise")
```
Second, restart Julia and load `Revise` before Oscar:
```
using Revise, Oscar;
```
The first run of `build_doc` will take the usual few minutes, subsequently runs
will be significantly faster.
"""
function build_doc(; doctest=false, strict=false, open_browser=true)
  versioncheck = (VERSION.major == 1) && (VERSION.minor == 6)
  versionwarn = 
"The Julia reference version for the doctests is 1.6, but you are using
$(VERSION). Running the doctests will produce errors that you do not expect."
  if doctest != false && !versioncheck
    @warn versionwarn
  end
  if !isdefined(Main, :BuildDoc)
    doc_init()
  end
  Pkg.activate(docsproject) do
    Base.invokelatest(Main.BuildDoc.doit, Oscar; strict=strict, local_build=true, doctest=doctest)
  end
  if open_browser
    open_doc()
  end
  if doctest != false && !versioncheck
    @warn versionwarn
  end
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
    @info("Running tests for $rel_test_file in same session")
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

include("Exports.jl")

# HACK/FIXME: remove these aliases once we have them in AA/Nemo/Hecke
@alias characteristic_polynomial charpoly  # FIXME
@alias minimal_polynomial minpoly  # FIXME

include("printing.jl")
include("fallbacks.jl")

include("GAP/wrappers.jl")

include("Groups/types.jl")
include("Groups/perm.jl")
include("Groups/group_constructors.jl")
include("Groups/sub.jl")
include("Groups/homomorphisms.jl")
include("Groups/cosets.jl")
include("Groups/GAPGroups.jl")
include("Groups/pcgroup.jl")
include("Groups/directproducts.jl")
include("Groups/matrices/matrices.jl")
include("Groups/matrices/FiniteFormOrthogonalGroup.jl")
include("Groups/libraries/libraries.jl")
include("Groups/gsets.jl")
include("Groups/MatrixDisplay.jl")
include("Groups/abelian_aut.jl")
include("Groups/spinor_norms.jl")
include("Groups/GrpAb.jl")

include("Rings/integer.jl")
include("Rings/orderings.jl")
include("Rings/mpoly.jl")
include("Rings/mpoly_types.jl")
include("Rings/mpoly-graded.jl")
include("Rings/mpoly-ideals.jl")
include("Rings/groebner.jl")
include("Rings/solving.jl")
include("Rings/MPolyQuo.jl")
include("Rings/mpoly-nested.jl")
include("Rings/FractionalIdeal.jl")

include("Rings/mpoly-affine-algebras.jl")

include("Rings/special_ideals.jl")

include("Rings/MPolyMap/MPolyAnyMap.jl")
include("Rings/MPolyMap/MPolyRing.jl")
include("Rings/MPolyMap/MPolyQuo.jl")
include("Rings/MPolyMap/AffineAlgebras.jl")

include("Rings/mpoly-local.jl")
include("Rings/localization_interface.jl")
include("Rings/mpoly-localizations.jl")
include("Rings/mpolyquo-localizations.jl")
include("Rings/MPolyMap/flattenings.jl")
include("Rings/FinField.jl")
include("Rings/NumberField.jl")
include("Rings/FunctionField.jl")
include("Rings/AbelianClosure.jl")

include("Rings/PBWAlgebra.jl")
include("Rings/PBWAlgebraQuo.jl")
include("Rings/FreeAssAlgIdeal.jl")


include("GAP/customize.jl")
include("GAP/gap_to_oscar.jl")
include("GAP/oscar_to_gap.jl")
include("GAP/iso_gap_oscar.jl")
include("GAP/iso_oscar_gap.jl")

include("Groups/group_characters.jl")  # needs some Rings functionality
include("Groups/action.jl")  # needs some PolynomialRings functionality

include("Modules/ModuleTypes.jl")
include("Modules/UngradedModules.jl")
include("Modules/homological-algebra.jl")
include("Modules/FreeModElem-orderings.jl")
#include("Modules/FreeModules-graded.jl")
include("Modules/ModulesGraded.jl")
include("Modules/module-localizations.jl")
include("Modules/local_rings.jl")
include("Modules/mpolyquo.jl")
include("Rings/ReesAlgebra.jl")

include("NumberTheory/NmbThy.jl")

include("PolyhedralGeometry/main.jl")

include("Polymake/polymake_to_oscar.jl")

include("Combinatorics/Graphs.jl")
include("Combinatorics/SimplicialComplexes.jl")
include("Combinatorics/Matroids/JMatroids.jl")
include("Combinatorics/Matroids/matroid_strata_grassmannian.jl")

include("StraightLinePrograms/StraightLinePrograms.jl")
include("Rings/lazypolys.jl")
include("Rings/slpolys.jl")
include("NumberTheory/GalThy.jl")

include("AlgebraicGeometry/Schemes/main.jl")
include("AlgebraicGeometry/ToricVarieties/JToric.jl")
include("AlgebraicGeometry/TropicalGeometry/main.jl")
include("AlgebraicGeometry/Surfaces/K3Auto.jl")
include("AlgebraicGeometry/Surfaces/SurfacesP4.jl")
include("AlgebraicGeometry/Miscellaneous/basics.jl")

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

include("Aliases.jl")

include("Deprecations.jl")

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
