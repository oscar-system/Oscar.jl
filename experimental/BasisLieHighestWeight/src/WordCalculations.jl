using ..Oscar
using ..Oscar: GAPWrap

function calculate_betas(word, simple_roots, cartan_matrix)
    # Calculate betas given simple_roots and cartan_matrix
    betas = []

    for k = 1:length(word)
        beta = copy(simple_roots[word[k]])
        for j = k-1:-1:1  # Iterate in reverse
            GAP.Globals.ApplySimpleReflection(cartan_matrix, word[j], beta)
        end
        push!(betas, beta)
    end
    
    return betas
end

function compute_betas(type::String, rank::Int, word::Vector{Int})
    lie_algebra = GAP.Globals.SimpleLieAlgebra(GAP.Obj(type), rank, GAP.Globals.Rationals)
    root_system = GAP.Globals.RootSystem(lie_algebra)

    simple_roots = GAP.Globals.SimpleSystem(root_system)
    weyl_group = GAP.Globals.WeylGroup(root_system)
    sparse_cartan_matrix = GAP.Globals.SparseCartanMatrix(weyl_group)

    betas = calculate_betas(word, simple_roots, sparse_cartan_matrix)

    return betas
end
