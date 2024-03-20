module GroebnerWalk
using Oscar

include("walk.jl")

export groebner_walk

function __init__()
    add_verbosity_scope(:intermediate)
    add_verbosity_scope(:detailed)

    set_verbosity_level(:intermediate, 1)
    set_verbosity_level(:detailed, 2)
end

end