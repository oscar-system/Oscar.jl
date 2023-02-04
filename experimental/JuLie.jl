module JuLie

using ..Oscar
using Markdown
import Nemo:libflint

include("JuLie/partitions.jl")
include("JuLie/schur_polynomials.jl")
include("JuLie/tableaux.jl")

end

using .JuLie
export Partition, partitions, schur_polynomial, Tableau, shape, semistandard_tableaux, is_standard, is_semistandard, standard_tableaux, schensted, hook_length, hook_lengths, num_standard_tableaux, reading_word, weight, bump!