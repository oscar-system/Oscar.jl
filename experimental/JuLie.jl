module JuLie

using ..Oscar
using Markdown
import Nemo:libflint

include("JuLie/partitions.jl")
include("JuLie/schur_polynomials.jl")

end

using .JuLie
export Partition, partitions, schur_polynomial