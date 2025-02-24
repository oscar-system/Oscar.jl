@doc raw"""
    graph_curve(G::Graph)

Return the graph curve of `G`, i.e., a union of lines whose dual graph is `G`, see [BE91](@cite).

Assumes that `G` is trivalent and 3-connected.

# Examples
```jldoctest
julia> G1 = vertex_edge_graph(simplex(3))
Undirected graph with 4 nodes and the following edges:
(2, 1)(3, 1)(3, 2)(4, 1)(4, 2)(4, 3)

julia> C1 = graph_curve(G1)
Projective curve
  in projective 2-space over QQ with coordinates [x1, x2, x3]
defined by ideal (x1^2*x2*x3 - x1*x2^2*x3 + x1*x2*x3^2)

julia> G2 = vertex_edge_graph(cube(3))
Undirected graph with 8 nodes and the following edges:
(2, 1)(3, 1)(4, 2)(4, 3)(5, 1)(6, 2)(6, 5)(7, 3)(7, 5)(8, 4)(8, 6)(8, 7)

julia> C2 = graph_curve(G2)
Projective curve
  in projective 4-space over QQ with coordinates [x1, x2, x3, x4, x5]
defined by ideal with 5 generators

```
"""
function graph_curve(G::Graph; check::Bool=true)
    @req all(isequal(3),degree(G)) "G is not trivalent"
    @req connectivity(G)==3 "G is not three-connected"

    R,x = graded_polynomial_ring(QQ, div(nv(G), 2) + 1)
    rowOfVariables = matrix(R,[x])

    cycleMatrix = kernel(matrix(QQ,signed_incidence_matrix(G)); side=:right)
    cycleMatrix = matrix(R,cycleMatrix) # converting to matrix over R for vcat below

    vertexIdeals = MPolyIdeal[]
    E = collect(edges(G))
    for v in 1:n_vertices(G)
        edgesContainingV = findall(edge->(v in edge),E)
        cycleMatrix_v = cycleMatrix[edgesContainingV,:]
        push!(vertexIdeals,ideal(minors(vcat(cycleMatrix_v,rowOfVariables),3)))
    end
    graphCurveIdeal = reduce(intersect,vertexIdeals)
    graphCurve = projective_curve(graphCurveIdeal; is_radical=true)

    return graphCurve
end
