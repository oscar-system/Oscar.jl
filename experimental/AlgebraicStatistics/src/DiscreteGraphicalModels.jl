function bron_kerbosch(G, R, P, X)

    cliques = []

    if length(P) == 0 && length(X) == 0
        
        push!(cliques, R)
    end

    Q = P
    Y = X

    for v in Q

        new_R = vcat(R, v);
        new_Q = [p for p in Q if p in neighbors(G, v)]
        new_Y = [p for p in Y if p in neighbors(G, v)]

        cliques = vcat(cliques, bron_kerbosch(G, new_R, new_Q, new_Y))

        Q = filter(i -> i != v, P)
        Y = vcat(X, v)
    end

    return cliques
end

function maximal_cliques(G::Graph{Undirected})

    max_cliques = bron_kerbosch(G, [], vertices(G), []);

    unique(sort.(max_cliques))
end




@doc raw"""
    graphical_model(G::Graph{Undirected}, S::GaussianRing, k_var_name::String="k")

A parametric statistical model associated to an undirected graph.
It contains an undirected graph `G`, a GaussianRing `S` where the vanishing ideal of the model naturally lives, 
and a parameter ring whose variables `k[i,j]` correspond to the edges `i - j` in `G`. 

## Examples

``` jldoctest undirected_ggm
julia> M = graphical_model(graph_from_edges([[1,2], [2,3]]), gaussian_ring(3))
Gaussian graphical model on an undirected graph with edges:
(1, 2), (2, 3)
```
"""
function graphical_model(G::Graph{Undirected}, S::MarkovRing, t_var_name::String="t")

    cliques = [reduce(*, string.(C)) for C in maximal_cliques(G)]

    rvs = random_variables(R)

    R, t = polynomial_ring(QQ, vcat([t_var_name*C*string([i]) for C in cliques for i in ])

    k = Dict([(Tuple(var_index(x)), x) for x in k])

    GraphicalModel(G, S, R, k)
end





S = MarkovRing("X1" => (1:2), "X2" => (1:2), "X3" => (1:2))

G = graph_from_edges(Undirected, [[1,2], [2,3]])



maximal_cliques(G)


# creates the parameterization of the model as a ring map
function parameterization(R::MarkovRing)

end




