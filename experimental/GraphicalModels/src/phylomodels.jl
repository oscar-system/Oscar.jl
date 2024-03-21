# specifying all common phylogenetic models
# as separate functions using the same structure


# 
# specify the structure
#

"""
    PhylogeneticModel

Define four variables which fully specify a phylogenetic model: the number of states `n_states` as an integer, the directed graph `graph`, the polynomial ring `R` and the transition matrices `trans_matrices` as a dictionary between edges and elements of `R`.

# Examples
```julia-repl
julia> pm=(4,graph_from_edges(Directed,[[4,1],[4,2],[4,3]]),?)???
??
```
"""
struct PhylogeneticModel
    n_states::Int
    graph::Graph{Directed}
    R::MPolyRing{QQFieldElem}
    trans_matrices::Dict{Edge, MatElem{QQMPolyRingElem}}
end

# calling elements of the structure
"""
    graph(pm)

Return the graph of a PhylogeneticModel `pm`.

# Examples
```julia-repl
julia> input???
output??
```
"""
graph(pm::PhylogeneticModel) = pm.graph

"""
    transition_matrices(pm)

Return the transition matrices of a PhylogeneticModel `pm`.

# Examples
```julia-repl
julia> input???
output??
```
"""
transition_matrices(pm::PhylogeneticModel) = pm.trans_matrices

"""
    number_states(pm)

Return the number of states of a PhylogeneticModel `pm`.

# Examples
```julia-repl
julia> input???
output??
```
"""
number_states(pm::PhylogeneticModel) = pm.n_states

#
# specify the models
#

# Jukes Cantor [ref?]
"""
    jukes_cantor_model(g)

Creates an element of type `PhylogeneticModel` based on `graph g`.

# Examples
```julia-repl
julia> jukes_cantor_model(graph_from_edges(Directed,[[4,1],[4,2],[4,3]]))
PhylogeneticModel(4, Directed graph with 4 nodes and 3 edges, Multivariate polynomial ring in 6 variables over QQ, Dict{Edge, MatElem{QQMPolyRingElem}}(Edge(4, 1) => [a[1] b[1] b[1] b[1]; b[1] a[1] b[1] b[1]; b[1] b[1] a[1] b[1]; b[1] b[1] b[1] a[1]], Edge(4, 3) => [a[3] b[3] b[3] b[3]; b[3] a[3] b[3] b[3]; b[3] b[3] a[3] b[3]; b[3] b[3] b[3] a[3]], Edge(4, 2) => [a[2] b[2] b[2] b[2]; b[2] a[2] b[2] b[2]; b[2] b[2] a[2] b[2]; b[2] b[2] b[2] a[2]]))
```
"""
function jukes_cantor_model(graph::Graph{Directed})
    ns = 4
    ne = n_edges(graph)
    R, list_a, list_b = polynomial_ring(QQ, :a => 1:ne, :b => 1:ne)
    
    matrices = Dict{Edge, MatElem}(e => matrix(R, [
      a b b b
      b a b b
      b b a b
      b b b a]) for (a,b,e) in zip(list_a, list_b, edges(graph))
    )
    return PhylogeneticModel(ns, graph, R, matrices)
end

# Kimura 2
function kimura2_model(graph::Graph{Directed})
    ns = 4
    ne = n_edges(graph)
    R, list_a, list_b, list_c = polynomial_ring(QQ, :a => 1:ne, :b => 1:ne, :c => 1:ne)
    
    matrices = Dict{Edge, MatElem}(e => matrix(R, [
      a b c c
      b a b c
      c b a b
      c c b a]) for (a,b,c,e) in zip(list_a, list_b, list_c, edges(graph))
    )
    return PhylogeneticModel(ns, graph, R, matrices)
end

# Kimura 3
function kimura3_model(graph::Graph{Directed})
    ns = 4
    ne = n_edges(graph)
    R, list_a, list_b , list_c, list_d= polynomial_ring(QQ, :a => 1:ne, :b => 1:ne, :c => 1:ne, d => 1:ne)
    
    matrices = Dict{Edge, MatElem}(e => matrix(R, [
      a b c d
      b a b c
      c b a b
      d c b a]) for (a,b,c,d,e) in zip(list_a, list_b, list_c, list_d, edges(graph))
    )
    return PhylogeneticModel(ns, graph, R, matrices)
end

# general Markov model
# work in progress
# not sure whether to use this https://docs.sciml.ai/Symbolics/stable/manual/arrays/
# or array of polynomial ring elements in OSCAR
function generalmarkov_model(graph::Graph{Directed})
    ns = 4
    ne = n_edges(graph)
    R, list_m = polynomial_ring(QQ, :m => 1:(ne^2))
    
    matrices = Dict{Edge, MatElem}(e => matrix(R, [
      m11 m12 m13 m14
      m11 m12 m13 m14
      b b a b
      b b b a]) for (a,b,e) in zip(list_a, list_b, edges(graph))
    )
    return PhylogeneticModel(ns, graph, R, matrices)
end

# more models to come, Ch Mar 21