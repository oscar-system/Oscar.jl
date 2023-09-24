using ..Oscar
using ..Oscar: GAPWrap

function compute_betas(type::String, rank::Int, word::Vector{Int})::Vector{Int}
    """
    Calculate betas from type, rank and a longest-word from the weylgroup.
    """
    # Construct Gap-Objects
    lie_algebra = GAP.Globals.SimpleLieAlgebra(GAP.Obj(type), rank, GAP.Globals.Rationals)
    root_system = GAP.Globals.RootSystem(lie_algebra)

    simple_roots = GAP.Globals.SimpleSystem(root_system)
    weyl_group = GAP.Globals.WeylGroup(root_system)
    sparse_cartan_matrix = GAP.Globals.SparseCartanMatrix(weyl_group)

    # Calculate betas by applying simple-reflections step-by-step.
    betas = []

    for k = 1:length(word)
        beta = copy(simple_roots[word[k]])
        for j = k-1:-1:1  # Iterate in reverse
            GAP.Globals.ApplySimpleReflection(cartan_matrix, word[j], beta)
        end
        push!(betas, beta)
    end

    betas = Vector{Int}(GAP.Globals.List(b)) # Convert to Vector{Int}
    return betas
end

function root_to_root_vector(type::String, rank::Int, root::Vector{Int})
    """
    For a given positive or negative root, return the GAP root vector.
    """
    # Construct Gap-Objects
    lie_algebra = GAP.Globals.SimpleLieAlgebra(GAP.Obj(type), rank, GAP.Globals.Rationals)
    root_system = GAP.Globals.RootSystem(lie_algebra)

    # Check if root is one of positive_roots
    positive_roots = Vector{Vector{Int}}(GAP.Globals.PositiveRoots(root_system))
    positive_root_vectors = GAP.Globals.PositiveRootVectors(root_system)
    
    for (i, root_i) in enumerate(positive_roots)
        print(root_i)
        if root == root_i
            return positive_root_vectors[i]
        end
    end

    # Check if root is one of negative roots
    negative_roots = Vector{Vector{Int}}(GAP.Globals.NegativeRoots(root_system))
    negative_root_vectors = GAP.Globals.NegativeRootVectors(root_system)
    
    for (i, root_i) in enumerate(negative_roots)
        if root == root_i
            return negative_root_vectors[i]
        end
    end

    return false
end

function_compute_operators_lustzig_nz(type::String, rank::Int, reduced_expression::Vector{Int})
    """
    Computes the operators for the lustzig and nz polytopes for a longest weyl-word 
    reduced_expression.

    \beta_k := s_{i_1} … s_{i_{k-1}} (\alpha_{i_k})

    F.e. for A, 2, [1, 2, 1], we get
    \beta_1 = \alpha_1
    \beta_2 = \alpha_1 + \alpha_2
    \beta_3 = \alpha_2
    """
    betas = compute_betas(type, rank, reduced_expression)
    operators = [root_to_root_vector(type, rank, beta) for beta in betas]
    
    return operators
end