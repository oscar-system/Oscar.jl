module JuLie

using ..Oscar
import Oscar: 
  IntegerUnion,
  number_of_partitions

include("partitions.jl")
include("schur_polynomials.jl")
include("tableaux.jl")

end

using .JuLie
