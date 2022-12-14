using Oscar
using Test
using Documenter


import Oscar.Nemo.AbstractAlgebra
include(joinpath(pathof(AbstractAlgebra), "..", "..", "test", "Rings-conformance-tests.jl"))

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
include("Modules/MPolyQuo.jl")

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
