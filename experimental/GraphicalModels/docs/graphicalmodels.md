julia> pm = jukes_cantor_model(graph_from_edges(Directed,[[4,1],[4,2],[4,3]]));
```
We can recover the information of `pm`
```jldoctest
julia> graph(pm)
Directed graph with 4 nodes and the following edges:
(4, 1)(4, 2)(4, 3)

julia> number_states(pm)
4

julia> prob_ring(pm)
fourier_ring(pm)
trans_matrices(pm)
fourier_params(pm)
group(pm)
