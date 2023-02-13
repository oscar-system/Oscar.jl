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

const stats_time = Dict{String,Float64}()
const stats_alloc = Dict{String,Float64}()

function print_stats(io::IO; fmt=PrettyTables.tf_unicode, max=15)
  sorted = sort(collect(stats_time), by=x->x[2], rev=true)
  println(io, "### Timings in seconds per file")
  println(io)
  PrettyTables.pretty_table(io, hcat(first.(sorted), last.(sorted)); tf=fmt, max_num_of_rows=max, header=[:Filename, Symbol("Time in s")])
  println(io)
  sorted = sort(collect(stats_alloc), by=x->x[2], rev=true)
  println(io, "### Allocations in megabytes per file")
  println(io)
  PrettyTables.pretty_table(io, hcat(first.(sorted), last.(sorted)); tf=fmt, max_num_of_rows=max, header=[:Filename, Symbol("Alloc in MB")])
end

# we only want to print stats for files that run tests and not those that just
# include other files
const innermost = Ref(true)
# redefine include to print and collect some extra stats
function include(str::String)
  innermost[] = true
  # we pass the identity to avoid recursing into this function again
  stats = @timed Base.include(identity, Main, str)
  if innermost[]
    dir = Base.source_dir()
    dir = joinpath(relpath(dir, joinpath(Oscar.oscardir,"test")), str)
    stats_time[dir] = stats.time
    stats_alloc[dir] = stats.bytes/2^20
    Printf.@printf "Testing %s took: %.2f seconds, %d MB\n" dir stats.time stats.bytes/2^20
    innermost[] = false
  end
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
