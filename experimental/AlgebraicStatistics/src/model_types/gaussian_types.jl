const ColoredGGM{Undirected} = GaussianGraphicalModel{
  Graph{Undirected}, @NamedTuple{color::GraphMap{Undirected}}
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
