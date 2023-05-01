add_assertion_scope(:LatWithIsom)
add_verbosity_scope(:LatWithIsom)

function set_lwi_level(k::Int)
  set_assertion_level(:LatWithIsom, k)
  set_verbosity_level(:LatWithIsom, k)
end

set_lwi_level(1)

include("exports.jl")
include("types.jl")
include("lattices_with_isometry.jl")
include("hermitian_miranda_morrison.jl")
include("enumeration.jl")
include("embeddings.jl")
include("printings.jl")
