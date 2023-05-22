using Oscar
using Test
using Aqua
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

function print_stats(io::IO; fmt=PrettyTables.tf_unicode, max=50)
  sorted = sort(collect(stats_dict), by=x->x[2].time, rev=true)
  println(io, "### Stats per file")
  println(io)
  table = hcat(first.(sorted), permutedims(reduce(hcat, collect.(values.(last.(sorted))))))
  formatters = nothing
  @static if compiletimes
    header=[:Filename, Symbol("Runtime in s"), Symbol("+ Compilation"), Symbol("+ Recompilation"), Symbol("Allocations in MB")]
    #formatters = PrettyTables.ft_printf("%.2f%%", [3,4])
  else
    header=[:Filename, Symbol("Time in s"), Symbol("Allocations in MB")]
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
  # skip files which just include other files and ignore
  # files outside of the oscar folder
  if innermost[] && !isabspath(str)
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
      stats_dict[path] = (time=stats.time-comptime, ctime=comptime-rcomptime, rctime=rcomptime, alloc=stats.bytes/2^20)
      Printf.@printf "-> Testing %s took: runtime %.3f seconds + compilation %.3f seconds + recompilation %.3f seconds, %.2f MB\n" path stats.time-comptime comptime-rcomptime rcomptime stats.bytes/2^20
    else
      Printf.@printf "-> Testing %s took: %.3f seconds, %.2f MB\n" path stats.time stats.bytes/2^20
      stats_dict[path] = (time=stats.time, alloc=stats.bytes/2^20)
    end
    innermost[] = false
  end
end

@static if compiletimes
  Base.cumulative_compile_timing(true)
end

# Used in both Rings/slpolys.jl and StraightLinePrograms/runtests.jl
const SLP = Oscar.StraightLinePrograms

include("Aqua.jl")

include("printing.jl")

include("PolyhedralGeometry/runtests.jl")
include("Combinatorics/runtests.jl")

include("GAP/runtests.jl")
include("Groups/runtests.jl")

include("Rings/runtests.jl")

include("NumberTheory/nmbthy.jl")

if Oscar.is_dev
  include("Experimental/GITFans.jl")
end

# Will automatically include all experimental packages following our
# guidelines.
include("../experimental/runtests.jl")

include("Experimental/galois.jl")
include("Experimental/gmodule.jl")
include("Experimental/ModStdQt.jl")
include("Experimental/ModStdNF.jl")
include("Experimental/MPolyRingSparse.jl")
include("Experimental/MatrixGroups.jl")
include("Experimental/JuLie.jl")
include("Experimental/ExteriorAlgebra.jl")

include("Rings/ReesAlgebra.jl")

include("Modules/runtests.jl")

include("InvariantTheory/runtests.jl")

include("AlgebraicGeometry/Schemes/runtests.jl")
include("AlgebraicGeometry/ToricVarieties/runtests.jl")
include("AlgebraicGeometry/TropicalGeometry/runtests.jl")
include("AlgebraicGeometry/Surfaces/K3Auto.jl")

include("Serialization/runtests.jl")

include("StraightLinePrograms/runtests.jl")

@static if compiletimes
  Base.cumulative_compile_timing(false);
end

if haskey(ENV, "GITHUB_STEP_SUMMARY") && compiletimes
  open(ENV["GITHUB_STEP_SUMMARY"], "a") do io
    print_stats(io, fmt=PrettyTables.tf_markdown)
  end
else
  print_stats(stdout; max=10)
end
