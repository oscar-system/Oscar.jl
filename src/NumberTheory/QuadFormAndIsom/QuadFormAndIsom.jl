function set_lwi_level(k::Int)
  set_assertion_level(:ZZLatWithIsom, k)
  set_verbosity_level(:ZZLatWithIsom, k)
end

include("types.jl")
include("spaces_with_isometry.jl")
include("lattices_with_isometry.jl")
include("hermitian_miranda_morrison.jl")
include("enumeration.jl")
include("embeddings.jl")
include("printings.jl")
