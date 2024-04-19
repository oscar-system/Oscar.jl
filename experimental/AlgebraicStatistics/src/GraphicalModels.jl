import Oscar.graph
import Oscar.vertices
import Oscar.gens

###################################################################################
#
#       Additional functions for manipulating variable indices
#
###################################################################################


function var_index(x)

    str = string(x)
    index = "["*split(str, "[")[2]

    eval(Meta.parse(index))
end



###################################################################################
#
#       Additional functions for graphs
#
###################################################################################


function vertices(G::Graph)

    E = [[src(e), dst(e)] for e in edges(G)]
    
    sort(unique(reduce(vcat, E)))
end



###################################################################################
#
#       Additional functions for matrices
#
###################################################################################


function cofactor(A)
    
    matrix(base_ring(A), [[(-1)^(i+j)*det(A[filter(l -> l != i, 1:n_rows(A)), filter(l -> l != j, 1:n_columns(A))]) for j in 1:n_columns(A)] for i in 1:n_rows(A)])
end


function adjoint(A::Generic.MatSpaceElem)
    
    transpose(cofactor(A))
end



###################################################################################
#
#       Generic Graphical Models
#
###################################################################################


# G, directed or undirected 
# T, GaussianRing or MarkovRing - specifies whether it is a Gaussian or Discrete graphical model
struct GraphicalModel{G, T}
    graph::G
    ring::T
    param_ring::Ring
    param_gens
end



# todo
function graph(M::GraphicalModel)

    M.graph
end



# todo
function ring(M::GraphicalModel)

    M.ring
end



function param_ring(M::GraphicalModel)

    M.param_ring
end



function param_gens(M::GraphicalModel)

    M.param_gens
end




include("GaussianGraphicalModels.jl")
#include("DiscreteGraphicalModels.jl")








# computes the vanishing ideal of a directed gaussian graphical model
function vanishing_ideal(M::GraphicalModel)

    kernel(parameterization(M))
end
