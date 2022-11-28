# avoid polluting the main namespace to to preserve proper printing in doctests
module SLPTest
using ..Test
using ..Oscar
using ..Oscar: SLP

include("setup.jl")
include("straightline.jl")
include("gap.jl")
include("atlas.jl")

end
