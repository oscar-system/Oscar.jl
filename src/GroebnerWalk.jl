module GroebnerWalk
using Oscar

exponent_vectors = f->exponent_vector.(monomials(f), Ref(1))

include("markedGB.jl")
include("generic_walk.jl")
include("walk.jl")

import Oscar: weight_ordering, ZZRingElem, MonomialOrdering, ZZMatrix, IdealGens
import Oscar.Orderings: MatrixOrdering, _support_indices

function weight_ordering(w::Vector{ZZRingElem}, o::MonomialOrdering)
    i = _support_indices(o.o)
    m = ZZMatrix(1, length(w), w)
    return MonomialOrdering(base_ring(o), MatrixOrdering(i, m, false))*o
end

#weight_ordering(w::Vector{Int}, o::MonomialOrdering) = weight_ordering(ZZ.(w), o)


export groebner_walk
export standard_walk

export initial_form
export initial_forms



export standard_step
export next_weight
export lift2
export difference_lead_tail

function __init__()
    add_verbosity_scope(:groebner_walk)

    set_verbosity_level(:groebner_walk, 0)
end

end
