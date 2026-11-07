# ---------------------------------------------------------------------------
# Native implementation of the functionality of Singular's `primdec.lib`,
# i.e. primary decomposition, minimal associated primes, radicals and
# equidimensional decompositions of ideals in polynomial rings over fields.
# The algorithms are organized by kind in the following files.
# ---------------------------------------------------------------------------

include("exports.jl")
include("helpers.jl")
include("zero-dimensional.jl")
include("radical.jl")
include("equidimensional.jl")
include("gtz.jl")
include("shimoyama-yokoyama.jl")
include("wrappers.jl")
