using Oscar
using Test

import Oscar.Nemo.AbstractAlgebra
include(joinpath(pathof(AbstractAlgebra), "..", "..", "test", "Rings-conformance-tests.jl"))

# Used in both Rings/slpolys-test.jl and StraightLinePrograms/runtests.jl
const SLP = Oscar.StraightLinePrograms
include("printing.jl")
include("Rings/runtests.jl")
include("TropicalGeometry/runtests.jl")
