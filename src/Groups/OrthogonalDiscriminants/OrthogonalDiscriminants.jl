include("utils.jl")
include("data.jl")
include("gram_det.jl")
include("direct.jl")
include("theoretical.jl")

# backwards compatibility for the examples from the Oscar book
module OrthogonalDiscriminants

using Oscar

character_of_entry = Oscar.character_of_entry
od_for_specht_module = Oscar.od_for_specht_module
od_from_eigenvalues = Oscar.od_from_eigenvalues
od_from_order = Oscar.od_from_order
od_from_p_subgroup = Oscar.od_from_p_subgroup
show_with_ODs = Oscar.show_with_ODs

end
