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

OSCAR is licensed under the GPL v3+ (see LICENSE.md).
"""
module Oscar

using Preferences

include("imports.jl")

# to allow access to the cornerstones! Otherwise, not even import or using from the
# user level will work as none of them will have been "added" by the user.
# possibly all should add a doc string to the module?
export Nemo, Hecke, Singular, Polymake, AbstractAlgebra, GAP

const cornerstones = String["AbstractAlgebra", "GAP", "Hecke", "Nemo", "Polymake", "Singular"];
const jll_deps = String["Antic_jll", "Arb_jll", "Calcium_jll", "FLINT_jll", "GAP_jll",
                        "libpolymake_julia_jll", "libsingular_julia_jll",
                        "polymake_jll", "Singular_jll"];

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

@doc Markdown.doc"""
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
    See also https://oscar.computeralgebra.de/install/.
    """)

if Sys.iswindows()
  windows_error()
end

function __init__()
  if Sys.iswindows()
    windows_error()
  end

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
    println("(c) 2019-2022 by The OSCAR Development Team")
  end

  append!(_gap_group_types,
    [
        (GAP.Globals.IsPermGroup, PermGroup),
        (GAP.Globals.IsPcGroup, PcGroup),
        (GAP.Globals.IsMatrixGroup, MatrixGroup),
        (GAP.Globals.IsFpGroup, FPGroup),
    ])
    __GAP_info_messages_off()
    GAP.Packages.load("browse"; install=true) # needed for all_character_table_names doctest
    GAP.Packages.load("ctbllib")
    GAP.Packages.load("forms")
    GAP.Packages.load("wedderga") # provides a function to compute Schur indices
    __init_IsoGapOscar()
    __init_group_libraries()
    __init_JuliaData()
    add_verbose_scope(:K3Auto)
    add_assert_scope(:K3Auto)
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


@doc Markdown.doc"""
    build_doc(; doctest=false, strict=false)

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
function build_doc(; doctest=false, strict=false)
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
  open_doc()
  if doctest != false && !versioncheck
    @warn versionwarn
  end
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


@doc Markdown.doc"""
    test_module(x, new::Bool = true)

Run the Oscar tests in the file `x.jl` where `x` is a string.

If `x == "all"` run the entire testsuite.

The optional parameter `new` takes the values `false` and `true` (default). If
`true`, then the tests are run in a new session, otherwise the currently active
session is used.
"""
function test_module(x::AbstractString, new::Bool = true)
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
include("Rings/rational.jl")
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
include("Rings/FinField.jl")
include("Rings/NumberField.jl")
include("Rings/FunctionField.jl")
include("Rings/AbelianClosure.jl")

include("Rings/PBWAlgebra.jl")
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
#include("Modules/FreeModules-graded.jl")
include("Modules/ModulesGraded.jl")
include("Modules/module-localizations.jl")

include("Geometry/basics.jl")
include("Geometry/K3Auto.jl")

include("NumberTheory/NmbThy.jl")

include("PolyhedralGeometry/main.jl")

include("Polymake/polymake_to_oscar.jl")
include("Combinatorics/Graphs.jl")
include("Combinatorics/SimplicialComplexes.jl")

include("Combinatorics/Matroids/JMatroids.jl")

include("StraightLinePrograms/StraightLinePrograms.jl")
include("Rings/lazypolys.jl")
include("Rings/slpolys.jl")

include("ToricVarieties/JToric.jl")

include("Schemes/main.jl")

include("TropicalGeometry/main.jl")

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
