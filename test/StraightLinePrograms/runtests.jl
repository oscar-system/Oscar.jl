# avoid polluting the main namespace to to preserve proper printing in doctests
module SLPTest
using ..Test
using ..Oscar
const SLP = Oscar.StraightLinePrograms

include("setup.jl")
include("straightline.jl")
include("gap.jl")
include("atlas.jl")

end
