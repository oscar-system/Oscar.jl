struct GraphicalModel{G, T}
    graph::Graph{G}
    ring::T
    param_ring::Ring
    param_gens::Dict
end


include("GaussianGraphicalModels.jl")
include("DiscreteGraphicalModels.jl")



# add
function graph(M::GraphicalModel)

    M.graph
end


# add
function ring(M::GraphicalModel)

    M.ring
end








# computes the vanishing ideal of a directed gaussian graphical model
function vanishing_ideal(R::GraphicalModel)

end
