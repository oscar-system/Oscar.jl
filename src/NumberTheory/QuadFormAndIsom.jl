function set_lwi_level(k::Int)
  set_assertion_level(:ZZLatWithIsom, k)
  set_verbosity_level(:ZZLatWithIsom, k)
end

include("QuadFormAndIsom/types.jl")
include("QuadFormAndIsom/spaces_with_isometry.jl")
include("QuadFormAndIsom/lattices_with_isometry.jl")
include("QuadFormAndIsom/hermitian_miranda_morrison.jl")
include("QuadFormAndIsom/enumeration.jl")
include("QuadFormAndIsom/embeddings.jl")
include("QuadFormAndIsom/printings.jl")
