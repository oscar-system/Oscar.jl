module GroebnerWalk
using Oscar

include("walk.jl")

export groebner_walk

function __init__()
    add_verbosity_scope(:groebner_walk)

    set_verbosity_level(:groebner_walk, 0)
end

end