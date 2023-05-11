include("Schemes/Schemes.jl")
include("ToricVarieties/JToric.jl")
include("Schemes/ToricSchemes/ToricSchemes.jl") # Has to come after the toric varieties due to dependencies!
include("TropicalGeometry/main.jl")
include("Surfaces/K3Auto.jl")
include("Surfaces/SurfacesP4.jl")
include("Miscellaneous/basics.jl")

