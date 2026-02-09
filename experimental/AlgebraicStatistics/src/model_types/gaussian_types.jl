const ColoredGGM{S} where {T <: Union{Mixed, Directed, Undirected}, S <: AbstractGraph{T}} = GaussianGraphicalModel{
  Graph{T}, @NamedTuple{color::GraphMap{T}}
}

@attr Tuple{
  QQMPolyRing,
  GraphDict{QQMPolyRingElem}
} function parameter_ring(GM::ColoredGGM{Undirected})
  G = graph(GM)
  colors = unique([[G.color[e] for e in edges(G)];
                   [G.color[v] for v in vertices(G)]])
  R, x = polynomial_ring(QQ, varnames(GM)[:k] => colors)
  color_dict = Dict{String, MPolyRingElem}(
    color => x[i] for (i, color) in enumerate(colors))

  gens_dict = GraphDict{QQMPolyRingElem}(
    Dict{Union{Int, Edge}, QQMPolyRingElem}(merge(
      Dict(e => color_dict[G.color[e]] for e in edges(G)),
      Dict(v => color_dict[G.color[v]] for v in vertices(G))
      )))
  return R, gens_dict
end

@attr Tuple{
  QQMPolyRing,
  GraphDict{QQMPolyRingElem}
} function parameter_ring(GM::GaussianGraphicalModel{Graph{Directed}, T}; cached=false) where T
  G = graph(GM)
  edge_colors = unique([G.color[e] for e in edges(G)])
  vertex_colors = unique([G.color[v] for v in vertices(G)])
  R, v_x, e_x = polynomial_ring(QQ, varnames(GM)[:w] => vertex_colors, varnames(GM)[:l] => edge_colors)
  color_dict = merge(Dict{String, MPolyRingElem}(vertex_color => v_x[i] for (i, vertex_color) in enumerate(vertex_colors)),
                     Dict{String, MPolyRingElem}(edge_color => e_x[i] for (i, edge_color) in enumerate(edge_colors)))

  gens_dict = merge(Dict(e => color_dict[G.color[e]] for e in edges(G)),
                    Dict(v => color_dict[G.color[v]] for v in 1:n_vertices(G)))
  return R, GraphDict(Dict{Union{Int, Edge}, QQMPolyRingElem}(gens_dict))
end
