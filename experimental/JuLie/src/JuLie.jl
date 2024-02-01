module JuLie

using ..Oscar
import Oscar: 
  IntegerUnion,
  number_of_partitions, npartitions,      # aliases do not work in experimental
  weight

include("partitions.jl")
include("schur_polynomials.jl")
include("tableaux.jl")

end

using .JuLie

export Partition
export Tableau
export ascending_partitions
export bump!
export conjugate
export dominates
export getindex_safe
export hook_length
export hook_lengths
export is_semistandard
export is_standard
export number_of_partitions, npartitions    # aliases do not work in experimental
export number_of_standard_tableaux
export partition
export partitions
export reading_word
export schensted
export schur_polynomial
export semistandard_tableaux
export shape
export standard_tableaux
export weight
