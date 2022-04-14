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

using Preferences

include("imports.jl")

# to allow access to the cornerstones! Otherwise, not even import or using from the
# user level will work as none of them will have been "added" by the user.
# possibly all should add a doc string to the module?
export Nemo, Hecke, Singular, Polymake, AbstractAlgebra, GAP

const cornerstones = String["AbstractAlgebra", "GAP", "Hecke", "Nemo", "Polymake", "Singular"];
const jll_deps = String["Antic_jll", "Arb_jll", "Calcium_jll", "FLINT_jll", "GAP_jll",
                        "libpolymake_julia_jll", "libsingular_julia_jll", "msolve_jll",
                        "polymake_jll", "Singular_jll"];

function _lookup_git_branch(dir::AbstractString)
   if length(Sys.which("git")) != nothing &&
         isdir(joinpath(dir,".git"))
      try
         ref = cd(dir) do
            readchomp(`git rev-parse --abbrev-ref HEAD`)
         end
         return " - branch #$(ref)"
      catch
      end
   end
   return ""
end

function _deps_git_info(dep::Pkg.API.PackageInfo)
   if dep.is_tracking_repo
      return " - branch #$(dep.git_revision)"
   elseif dep.is_tracking_path
      return _lookup_git_branch(dep.source)
   end
   return ""
end

function _print_dependency_versions(io::IO, deps::AbstractArray{<:AbstractString}; padding="    ", suffix="", branch=false)
   width = maximum(length.(deps))+length(suffix)+2
   deps = filter(d->d.name in deps, collect(values(Pkg.dependencies())))
   deps = sort!(deps; by=x->x.name)
   for dep in deps
      print(io, "$(padding)$(rpad(dep.name*suffix, width, ' ')) v$(dep.version)")
      println(io, branch ? _deps_git_info(dep) : "")
   end
end

@doc Markdown.doc"""
    Oscar.versioninfo(io::IO=stdout; branch=false, jll=false, julia=false)

Print the versions of all Oscar-related dependencies.

# Arguments
- `branch::Bool=false`: include git branch name in the output
- `jll::Bool=false`   : include binary packages (jll) in the output
- `julia::Bool=false` : include julia `versioninfo` output
"""
function versioninfo(io::IO=stdout; branch=false, jll=false, julia=false)
   print(io, "OSCAR version $(VERSION_NUMBER)")
   println(io, branch ? _lookup_git_branch(dirname(@__DIR__)) : "")
   println(io, "  combining:")
   _print_dependency_versions(io, cornerstones; suffix=".jl", branch=branch)
   if jll
      println(io, "  building on:")
      _print_dependency_versions(io, jll_deps; branch=branch)
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
    println("(c) 2019-2022 by The Oscar Development Team")
  end

  append!(_gap_group_types,
    [
        (GAP.Globals.IsPermGroup, PermGroup),
        (GAP.Globals.IsPcGroup, PcGroup),
        (GAP.Globals.IsMatrixGroup, MatrixGroup),
        (GAP.Globals.IsFpGroup, FPGroup),
    ])
    GAP.Packages.load("ctbllib")
    GAP.Packages.load("forms")
    __init_IsoGapOscar()
    __GAP_info_messages_off()
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

function build_doc(; doctest=false, strict=false)
  if !isdefined(Main, :BuildDoc)
    doc_init()
  end
  Pkg.activate(docsproject) do
    Base.invokelatest(Main.BuildDoc.doit, Oscar; strict=strict, local_build=true, doctest=doctest)
  end
  open_doc()
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

function iseffective end

include("printing.jl")

include("GAP/wrappers.jl")

include("Groups/types.jl")
include("Groups/perm.jl")
include("Groups/group_constructors.jl")
include("Groups/sub.jl")
include("Groups/homomorphisms.jl")
include("Groups/cosets.jl")
include("Groups/libraries/libraries.jl")
include("Groups/GAPGroups.jl")
include("Groups/directproducts.jl")
include("Groups/matrices/matrices.jl")
include("Groups/matrices/FiniteFormOrthogonalGroup.jl")
include("Groups/action.jl")
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
include("Rings/msolve/interface.jl")
include("Rings/msolve/f4.jl")
include("Rings/msolve/msolve.jl")
include("Rings/groebner.jl")
include("Rings/MPolyQuo.jl")
include("Rings/mpoly-nested.jl")
include("Rings/FractionalIdeal.jl")

include("Rings/mpoly-affine-algebras.jl")

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

include("Rings/FreeAssAlgIdeal.jl")

include("GAP/customize.jl")
include("GAP/gap_to_oscar.jl")
include("GAP/oscar_to_gap.jl")
include("GAP/iso_gap_oscar.jl")
include("GAP/iso_oscar_gap.jl")

include("Groups/group_characters.jl")  # needs some Rings functionality

include("Modules/ModuleTypes.jl")
include("Modules/UngradedModules.jl")
include("Modules/LocalizedModules.jl")
#include("Modules/FreeModules-graded.jl")
include("Modules/ModulesGraded.jl")

include("Geometry/basics.jl")

include("NumberTheory/NmbThy.jl")

include("PolyhedralGeometry/main.jl")

include("Combinatorics/Graphs.jl")
export Graphs
include("Combinatorics/SimplicialComplexes.jl")

include("StraightLinePrograms/StraightLinePrograms.jl")
include("Rings/lazypolys.jl")
include("Rings/slpolys.jl")

include("../experimental/Experimental.jl")
include("Rings/binomial_ideals.jl")

include("ToricVarieties/JToric.jl")

include("TropicalGeometry/main.jl")

if is_dev
#  include("../examples/ModStdNF.jl")
#  include("../examples/ModStdQ.jl")
#  include("../examples/ModStdQt.jl")
  include("../examples/PrimDec.jl")
#  include("../examples/GaloisGrp.jl")

#  include("../examples/PlaneCurve.jl")
end

include("Serialization/main.jl")

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
