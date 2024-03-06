@doc """
OrthogonalDiscriminants.jl is a collection of tools and data
related to the orthogonal discriminants of characters
of almost simple Atlas groups.
"""
module OrthogonalDiscriminants

# the necessary Julia packages
using Oscar
using GAP
using JSON

# as long as the code is in `experimental` ...
import Oscar.IntegerUnion
import Oscar.GAPGroupClassFunction
import Oscar.Partition
import Oscar.partition

# The following code must be loaded at runtime.
function __init__()
  # This file will not be needed after the release of GAP 4.13.
  GAP.Globals.Read(GapObj(joinpath(@__DIR__, "../gap/fix.g")))
end

# The following code can be loaded at compile time.
include("utils.jl")
include("data.jl")
include("gram_det.jl")
include("direct.jl")
include("theoretical.jl")
include("exports.jl")

end # module

using .OrthogonalDiscriminants

include("exports.jl")
