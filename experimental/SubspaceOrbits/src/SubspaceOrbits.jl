module SubspaceOrbits

using ..Oscar

import ..Oscar:
  @req,
  GAPWrap,
  GapObj,
  MatGroupElem,
  _as_subgroup_bare,
  _ffe_to_int,
  set_order

import Random

export subspace_orbit_representatives_and_lengths
export subspace_orbit_representatives_and_stabilizers

function __init__()
  add_verbosity_scope(:SubspaceOrbits)
end

include("labels.jl")
include("label_actions.jl")
include("orbits.jl")
include("interface.jl")

end

using .SubspaceOrbits

export subspace_orbit_representatives_and_lengths
export subspace_orbit_representatives_and_stabilizers
