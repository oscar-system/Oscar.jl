# @doc """
# OrthogonalDiscriminants.jl is a collection of tools and data
# related to the orthogonal discriminants of characters
# of almost simple Atlas groups.
# """
# module OrthogonalDiscriminants
# 
# using Oscar
# using GAP
# import JSON
# 
# import Oscar.IntegerUnion
# import Oscar.GAPGroupClassFunction
# import Oscar.Partition
# import Oscar.partition

include("utils.jl")
include("data.jl")
include("gram_det.jl")
include("direct.jl")
include("theoretical.jl")
#include("exports.jl")

# end # module
# 
# using .OrthogonalDiscriminants

