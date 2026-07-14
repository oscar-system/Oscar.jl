function reverse_digraph(D::Digraph)
  return Digraph(DigraphWrap.DigraphReverse(GapObj(D)))
end

function digraph_dual(D::Digraph)
  return Digraph(DigraphWrap.DigraphDual(GapObj(D)))
end

function reduced_digraph(D::Digraph)
  return Digraph(DigraphWrap.ReducedDigraph(GapObj(D)))
end

function induced_subdigraph(D::Digraph, vertices::Vector{Int})
  return Digraph(DigraphWrap.InducedSubdigraph(GapObj(D), GapObj(vertices)))
end

function quotient_digraph(D::Digraph, partition::Vector{Vector{Int}})
  return Digraph(DigraphWrap.QuotientDigraph(GapObj(D), GapObj(partition)))
end

function digraph_disjoint_union(D1::Digraph, D2::Digraph)
  return Digraph(DigraphWrap.DigraphDisjointUnion(GapObj(D1), GapObj(D2)))
end

function digraph_join(D1::Digraph, D2::Digraph)
  return Digraph(DigraphWrap.DigraphJoin(GapObj(D1), GapObj(D2)))
end

function digraph_union(D1::Digraph, D2::Digraph)
  return Digraph(DigraphWrap.DigraphEdgeUnion(GapObj(D1), GapObj(D2)))
end

function digraph_cartesian_product(D1::Digraph, D2::Digraph)
  return Digraph(DigraphWrap.DigraphCartesianProduct(GapObj(D1), GapObj(D2)))
end

function digraph_direct_product(D1::Digraph, D2::Digraph)
  return Digraph(DigraphWrap.DigraphDirectProduct(GapObj(D1), GapObj(D2)))
end

function digraph_lexicographic_product(D1::Digraph, D2::Digraph)
  return Digraph(DigraphWrap.DigraphLexicographicProduct(GapObj(D1), GapObj(D2)))
end

function digraph_modular_product(D1::Digraph, D2::Digraph)
  return Digraph(DigraphWrap.ModularProduct(GapObj(D1), GapObj(D2)))
end

function digraph_strong_product(D1::Digraph, D2::Digraph)
  return Digraph(DigraphWrap.StrongProduct(GapObj(D1), GapObj(D2)))
end


