module GroebnerWalk
using Oscar

exponent_vectors = f->exponent_vector.(monomials(f), Ref(1))

include("walk.jl")

include("facet_preorder.jl")



export groebner_walk
export standard_walk

export initial_form
export initial_forms



export standard_step
export next_weight
export lift2
export difference_lead_tail

export new_next_gamma

function __init__()
    add_verbosity_scope(:groebner_walk)

    set_verbosity_level(:groebner_walk, 0)
end

end