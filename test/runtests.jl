using Oscar
using Test
using Documenter

import Random

if haskey(ENV, "JULIA_PKGEVAL") ||
    get(ENV, "CI", "") == "true" ||
    haskey(ENV, "OSCAR_RANDOM_SEED")
  seed = parse(UInt32, get(ENV, "OSCAR_RANDOM_SEED", "42"))
  @info string(@__FILE__)*" -- fixed SEED $seed"
else
  seed = rand(UInt32)
  @info string(@__FILE__)*" -- SEED $seed"
end
Oscar.set_seed!(seed)

import Oscar.Nemo.AbstractAlgebra
include(joinpath(pathof(AbstractAlgebra), "..", "..", "test", "Rings-conformance-tests.jl"))

# some helpers
import Printf
import PrettyTables

# the current code for extracting the compile times does not work on earlier
# julia version
const compiletimes = @static VERSION >= v"1.9.0-DEV" ? true : false

const stats_dict = Dict{String,NamedTuple}()

function print_stats(io::IO; fmt=PrettyTables.tf_unicode, max=25)
  sorted = sort(collect(stats_dict), by=x->x[2].time, rev=true)
  println(io, "### Stats per file")
  println(io)
  table = hcat(first.(sorted), permutedims(reduce(hcat, collect.(values.(last.(sorted))))))
  @static if compiletimes
    header=([:Filename, Symbol("Time"), :Compilation, :Recompilation, Symbol("Alloc")], [Symbol(""), :seconds, Symbol(""), Symbol("of comilation"), :MB])
    formatters = PrettyTables.ft_printf("%.2f%%", [3,4])
  else
    header=([:Filename, Symbol("Time"), Symbol("Alloc")], [Symbol(""), :seconds, :MB])
    formatters = nothing
  end
  PrettyTables.pretty_table(io, table; tf=fmt, max_num_of_rows=max, header=header, formatters=formatters)
end

# we only want to print stats for files that run tests and not those that just
# include other files
const innermost = Ref(true)
# redefine include to print and collect some extra stats
function include(str::String)
  innermost[] = true
  # we pass the identity to avoid recursing into this function again
  @static if compiletimes
    compile_elapsedtimes = Base.cumulative_compile_time_ns()
  end
  stats = @timed Base.include(identity, Main, str)
  if innermost[]
    @static if compiletimes
      compile_elapsedtimes = Base.cumulative_compile_time_ns() .- compile_elapsedtimes
      compile_elapsedtimes = compile_elapsedtimes ./ 10^9
    end
    path = Base.source_dir()
    path = joinpath(relpath(path, joinpath(Oscar.oscardir,"test")), str)
    rtime=NaN
    @static if compiletimes
      comptime = first(compile_elapsedtimes)
      rcomptime = last(compile_elapsedtimes)
      stats_dict[path] = (time=stats.time, ctime=100*comptime/stats.time, rctime=100*rcomptime/comptime, alloc=stats.bytes/2^20)
      Printf.@printf "-> Testing %s took: %.3f seconds, %.2f%% compilation, of this %.2f%% recompilation, %.2f MB\n" path stats.time 100*comptime/stats.time 100*rcomptime/comptime stats.bytes/2^20
    else
      Printf.@printf "-> Testing %s took: %.3f seconds, %.2f MB\n" path stats.time stats.bytes/2^20
      stats_dict[path] = (time=stats.time, alloc=stats.bytes/2^20)
    end
    innermost[] = false
  end
end

@static if compiletimes
  Base.cumulative_compile_timing(true);
end
# Used in both Rings/slpolys-test.jl and StraightLinePrograms/runtests.jl
const SLP = Oscar.StraightLinePrograms
include("printing.jl")

include("PolyhedralGeometry/runtests.jl")
include("Combinatorics/runtests.jl")

include("GAP/runtests.jl")
include("Groups/runtests.jl")

include("Rings/runtests.jl")

include("NumberTheory/nmbthy-test.jl")

if Oscar.is_dev
  include("Experimental/GITFans-test.jl")
end

include("Experimental/PlaneCurve-test.jl")
include("Experimental/galois-test.jl")
include("Experimental/gmodule-test.jl")
include("Experimental/ModStdQt-test.jl")
include("Experimental/ModStdNF-test.jl")
include("Experimental/MPolyRingSparse-test.jl")
include("Experimental/MatrixGroups-test.jl")

include("Geometry/K3Auto.jl")

include("Modules/UngradedModules.jl")
include("Modules/FreeModElem-orderings-test.jl")
include("Modules/ModulesGraded.jl")
include("Modules/module-localizations.jl")
include("Modules/local_rings.jl")
include("Modules/MPolyQuo.jl")
include("Modules/homological-algebra_test.jl")
include("Rings/ReesAlgebra.jl")

include("InvariantTheory/runtests.jl")

include("ToricVarieties/runtests.jl")

include("Modules/ProjectiveModules.jl")
include("Schemes/runtests.jl")

include("TropicalGeometry/runtests.jl")
include("Serialization/runtests.jl")

include("StraightLinePrograms/runtests.jl")

@static if compiletimes
  Base.cumulative_compile_timing(false);
end

# Doctests

# We want to avoid running the doctests twice so we skip them when
# "oscar_run_doctests" is set by OscarDevTools.jl
if v"1.6.0" <= VERSION < v"1.7.0" && !haskey(ENV,"oscar_run_doctests")
  @info "Running doctests (Julia version is 1.6)"
  DocMeta.setdocmeta!(Oscar, :DocTestSetup, :(using Oscar, Oscar.Graphs); recursive = true)
  doctest(Oscar)
else
  @info "Not running doctests (Julia version must be 1.6)"
end

if haskey(ENV, "GITHUB_STEP_SUMMARY")
  open(ENV["GITHUB_STEP_SUMMARY"], "a") do io
    print_stats(io, fmt=PrettyTables.tf_markdown)
  end
else
  print_stats(stdout)
end
