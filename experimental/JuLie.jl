module JuLie

using ..Oscar
import Oscar: IntegerUnion
using Markdown

include("JuLie/partitions.jl")
include("JuLie/schur_polynomials.jl")
include("JuLie/tableaux.jl")

include("JuLie/quantum_analogs.jl")

end

using .JuLie
export Partition, partitions, ascending_partitions, dominates, conjugate, getindex_safe, num_partitions, schur_polynomial, Tableau, shape, semistandard_tableaux, is_standard, is_semistandard, standard_tableaux, schensted, hook_length, hook_lengths, num_standard_tableaux, reading_word, weight, bump!, quantum_integer, quantum_factorial, quantum_binomial
