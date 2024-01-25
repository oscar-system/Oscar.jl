using Oscar
using Test
using Distributed

import Random

numprocs_str = get(ENV, "NUMPROCS", "1")

oldWorkingDirectory = pwd()
cd(joinpath(pkgdir(Oscar), "test"))

if !isempty(ARGS)
  jargs = [arg for arg in ARGS if startswith(arg, "-j")]
  if !isempty(jargs)
    numprocs_str = split(jargs[end], "-j")[end]
  end
end

const numprocs = parse(Int, numprocs_str)

if numprocs >= 2
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
@everywhere Oscar.randseed!($seed)
# setting the global julia seed does not work for distributed processes
# the RNG is task-local and each '@everywhere' runs in a separate task...
# to make sure we seed the main process we run this again
Oscar.randseed!(seed)

if VERSION >= v"1.8.0"
  # Enable GC logging to help track down certain GC related issues.
  # Note that several test files need to temporarily disable and then
  # re-enable this. If we need to disable this globally, those files
  # need to be adjusted as well.
  @everywhere  GC.enable_logging(true)
end

@everywhere import Oscar.Nemo.AbstractAlgebra
@everywhere include(joinpath(pathof(Oscar.Nemo.AbstractAlgebra), "..", "..", "test", "Rings-conformance-tests.jl"))

# hotfix, otherwise StraightLinePrograms returns something which then leads to an error
module SLPTest
end

# some helpers
@everywhere import PrettyTables


function print_stats(io::IO, stats_dict::Dict; fmt=PrettyTables.tf_unicode, max=50)
  sorted = sort(collect(stats_dict), by=x->x[2].time, rev=true)
  println(io, "### Stats per file")
  println(io)
  table = hcat(first.(sorted), permutedims(reduce(hcat, collect.(values.(last.(sorted))))))
  formatters = nothing
  if haskey(first(values(stats_dict)), :ctime)
    header=[:Filename, Symbol("Runtime in s"), Symbol("+ Compilation"), Symbol("+ Recompilation"), Symbol("Allocations in MB")]
    #formatters = PrettyTables.ft_printf("%.2f%%", [3,4])
  else
    header=[:Filename, Symbol("Time in s"), Symbol("Allocations in MB")]
  end
  PrettyTables.pretty_table(io, table; tf=fmt, max_num_of_rows=max, header=header, formatters=formatters)
end


testlist = Oscar._gather_tests("test")

for exp in [Oscar.exppkgs; Oscar.oldexppkgs]
  path = joinpath(Oscar.oscardir, "experimental", exp, "test")
  if isdir(path)
    append!(testlist, Oscar._gather_tests(path))
  end
end

# make sure we have the same list everywhere
sort!(testlist)
Random.shuffle!(Oscar.get_seeded_rng(), testlist)

# tests with the highest number of allocations / runtime / compilation time
# more or less sorted by allocations
test_large = [
              "experimental/FTheoryTools/test/weierstrass.jl",
              "test/PolyhedralGeometry/timing.jl",
              "experimental/GITFans/test/runtests.jl",
              "test/AlgebraicGeometry/ToricVarieties/toric_schemes.jl",
              "test/AlgebraicGeometry/Schemes/WeilDivisor.jl",
              "test/Rings/NumberField.jl",
              "test/Serialization/PolynomialsSeries.jl",
              "test/AlgebraicGeometry/Schemes/K3.jl",
              "test/Groups/forms.jl",
              "test/Modules/UngradedModules.jl",
              "test/GAP/oscarinterface.jl",
              "test/AlgebraicGeometry/Schemes/CoveredProjectiveSchemes.jl",
              "test/AlgebraicGeometry/Schemes/MorphismFromRationalFunctions.jl",
              "experimental/QuadFormAndIsom/test/runtests.jl",
              "experimental/GModule/test/runtests.jl",
              "test/Modules/ModulesGraded.jl",
              "test/AlgebraicGeometry/Schemes/elliptic_surface.jl",
             ]

test_subset = get(ENV, "OSCAR_TEST_SUBSET", "")
if haskey(ENV, "JULIA_PKGEVAL")
  test_subset = "short"
end

if test_subset == "short"
  filter!(x-> !in(relpath(x, Oscar.oscardir), test_large), testlist)
elseif test_subset == "long"
  testlist = joinpath.(Oscar.oscardir, test_large)
  filter!(isfile, testlist)
end


@everywhere testlist = $testlist

# if many workers, distribute tasks across them
# otherwise, is essentially a serial loop
stats = reduce(merge, pmap(testlist) do x
                        println("Starting tests for $x")
                        Oscar.test_module(x; new=false, timed=true)
                      end)

# this needs to run here to make sure it runs on the main process
# it is in the ignore list for the other tests
#
# FIXME: disabled until #3184 is solved
#if numprocs == 1 && test_subset != "short"
#  # to debug GC corruption errors
#  GC.gc()
#  GC.gc()
#
#  println("Starting tests for Serialization/IPC.jl")
#  push!(stats, Oscar._timed_include("Serialization/IPC.jl", Main))
#end

if haskey(ENV, "GITHUB_STEP_SUMMARY")
  open(ENV["GITHUB_STEP_SUMMARY"], "a") do io
    print_stats(io, stats; fmt=PrettyTables.tf_markdown)
  end
else
  print_stats(stdout, stats; max=10)
end

cd(oldWorkingDirectory)
