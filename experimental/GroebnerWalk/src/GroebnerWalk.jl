module GroebnerWalk
using Oscar

import Oscar: 
  IdealGens,
  MonomialOrdering,
  ngens
  weight_ordering,
  ZZMatrix,
  ZZRingElem 


import Oscar.Orderings: 
  MatrixOrdering,
  _support_indices

include("markedGB.jl")
include("common.jl")

include("special-ideals.jl")

include("standard_walk.jl")
include("generic_walk.jl")
include("perturbed_walk.jl")

export groebner_walk

export newell_patch
export newell_patch_with_orderings

function __init__()
    add_verbosity_scope(:groebner_walk)
end

end

import .GroebnerWalk:
  newell_patch,
  newell_patch_with_orderings,
  groebner_walk

export groebner_walk

