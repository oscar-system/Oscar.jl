using Oscar
using Test
using Documenter
using Distributed

import Random


if !isempty(ARGS)
  jargs = [arg for arg in ARGS if startswith(arg, "-j")]
  if !isempty(jargs)
    numprocs = parse(Int,split(jargs[end], "-j")[end])
  end
elseif haskey(ENV, "NUMPROCS")
  numprocs = parse(Int,ENV["NUMPROCS"])
else
  numprocs = 1
end

if (numprocs >= 2)
  println("Adding worker processes")
  addprocs(numprocs)
end

if haskey(ENV, "JULIA_PKGEVAL") ||
  get(ENV, "CI", "") == "true" ||
  haskey(ENV, "OSCAR_RANDOM_SEED")
seed = parse(UInt32, get(ENV, "OSCAR_RANDOM_SEED", "42"))
@info string(@__FILE__)*" -- fixed SEED $seed"
else
seed = rand(UInt32)
@info string(@__FILE__)*" -- SEED $seed"
end

@everywhere using Test
@everywhere using Oscar
@everywhere Oscar.set_seed!($seed)

# hotfix, otherwise StraightLinePrograms returns something which then leads to an error
module SLPTest
end

# some helpers
@everywhere import Printf
@everywhere import PrettyTables

# the current code for extracting the compile times does not work on earlier
# julia version
@everywhere const compiletimes = @static VERSION >= v"1.9.0-DEV" ? true : false

@everywhere const stats_dict = Dict{String,NamedTuple}()

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
@everywhere const innermost = Ref(true)
# redefine include to print and collect some extra stats
@everywhere function include(str::String)
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
  Base.cumulative_compile_timing(true);
end

println("Making test list")

testlist = []
push!(testlist, "printing.jl")
push!(testlist, "PolyhedralGeometry/runtests.jl")
push!(testlist, "Combinatorics/runtests.jl")

push!(testlist, "GAP/runtests.jl")
push!(testlist, "Groups/runtests.jl")

push!(testlist, "Rings/runtests.jl")

push!(testlist, "NumberTheory/nmbthy-test.jl")

if Oscar.is_dev
  push!(testlist, "Experimental/GITFans-test.jl")
end

# Will automatically include all experimental packages following our
# guidelines.
push!(testlist, "../experimental/runtests.jl")

push!(testlist, "Experimental/galois-test.jl")
push!(testlist, "Experimental/gmodule-test.jl")
push!(testlist, "Experimental/ModStdQt-test.jl")
push!(testlist, "Experimental/ModStdNF-test.jl")
push!(testlist, "Experimental/MPolyRingSparse-test.jl")
push!(testlist, "Experimental/MatrixGroups-test.jl")
push!(testlist, "Experimental/JuLie-test.jl")
push!(testlist, "Experimental/SymmetricIntersections-test.jl")
push!(testlist, "Experimental/ExteriorAlgebra-test.jl")

push!(testlist, "Rings/ReesAlgebra.jl")

push!(testlist, "Modules/runtests.jl")

push!(testlist, "InvariantTheory/runtests.jl")

push!(testlist, "AlgebraicGeometry/Schemes/runtests.jl")
push!(testlist, "AlgebraicGeometry/ToricVarieties/runtests.jl")
push!(testlist, "AlgebraicGeometry/TropicalGeometry/runtests.jl")
push!(testlist, "AlgebraicGeometry/Surfaces/K3Auto.jl")

push!(testlist, "Serialization/runtests.jl")

push!(testlist, "StraightLinePrograms/runtests.jl")

# if many workers, distribute tasks across them
# otherwise, is essentially a serial loop
  @time pmap(x -> include(x), testlist)

@static if compiletimes
  Base.cumulative_compile_timing(false);
end

# Doctests


#currently, print_stats will fail when running tests with external workers
#TODO: potentially rewrite include as well as print_stats 
#      to comply with parallel decisions
if (numprocs == 1)
  if haskey(ENV, "GITHUB_STEP_SUMMARY") && compiletimes
    open(ENV["GITHUB_STEP_SUMMARY"], "a") do io
      print_stats(io, fmt=PrettyTables.tf_markdown)
    end
  else
    print_stats(stdout; max=10)
  end
end
