module GroebnerWalk
using Oscar

exponent_vectors = f->leading_exponent_vector.(monomials(f))

include("markedGB.jl")
include("generic_walk.jl")
include("walk.jl")

include("standard_walk.jl")
include("perturbed_walk.jl")
include("fractal_walk.jl")

import Oscar: weight_ordering, ZZRingElem, MonomialOrdering, ZZMatrix, IdealGens
import Oscar.Orderings: MatrixOrdering, _support_indices

function weight_ordering(w::Vector{ZZRingElem}, o::MonomialOrdering)
    i = _support_indices(o.o)
    m = ZZMatrix(1, length(w), w)
    return MonomialOrdering(base_ring(o), MatrixOrdering(i, m, false))*o
end

export groebner_walk

function __init__()
    add_verbosity_scope(:groebner_walk)

    set_verbosity_level(:groebner_walk, 0)
end

end
