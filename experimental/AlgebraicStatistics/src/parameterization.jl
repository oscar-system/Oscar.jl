function parameterization(M)

    #parameterization
    #=
    This function will take a phylogenetic model object M as input and output the tuple (R,S,images) where R = parameter_ring(M), 
    S = model_ring(M) and images is the parameterization in Fourier coordinates of the phylogenetic network mordel M. 
    =#

    #defines the parameter_ring and model_ring for the model M with k hybrid vertices
    R = parameter_ring(M)
    S = model_ring(M)

    #this accesses the Directed Acyclic Graph from the input M (type::GroupBasedPhylogeneticModel)
    #Also determines the level of the network
    DAG = M.phylo_model.graph.graph
    N = phylogenetic_network(DAG)
    lvl = level(N)

    #defines globally hybrid parameters so we can write down parameterisation of network model

    l = Array{QQMPolyRingElem}(undef, ntuple(_ -> 2, lvl))
    for i in 1:(2^lvl)
        l[i] = S[i]
    end

    #this will return the list of directed trees obtained from the network.
    #presumably the result is an ordered list of tuples (T, theta) where theta is an element of (Z_2)^lvl i.e. a particular hybridization event
    #and T is the dispalyed tree one gets from the hybridization event theta. This will be useful as then one could index the hybrid parameters by theta.
    trees = display_trees(DAG) # or maybe input N?

    #for all trees found before it finds the tree parameterization and then sums up all these tree
    #parameterizations to obtain the network praameterization.

    images = sum(sum([network_subtree_param(T, M) for T in trees]))

    return (R,S,images)
    
end
