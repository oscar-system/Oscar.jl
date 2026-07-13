function reverse_digraph(D::Digraph)
  return DirectedDigraph(DigraphWrap.DigraphReverse(GapObj(D)))
end

function digraph_dual(D::Digraph)
  return DirectedDigraph(DigraphWrap.DigraphDual(GapObj(D)))
end

function reduced_digraph(D::Digraph)
  return DirectedDigraph(DigraphWrap.ReducedDigraph(GapObj(D)))
end

function induced_subdigraph(D::Digraph, vertices::Vector{Int})
  return DirectedDigraph(DigraphWrap.InducedSubdigraph(GapObj(D), GapObj(vertices)))
end

function quotient_digraph(D::Digraph, partition::Vector{Vector{Int}})
  return DirectedDigraph(DigraphWrap.QuotientDigraph(GapObj(D), GapObj(partition)))
end

function digraph_disjoint_union(D1::Digraph, D2::Digraph)
  return DirectedDigraph(DigraphWrap.DigraphDisjointUnion(GapObj(D1), GapObj(D2)))
end

function digraph_join(D1::Digraph, D2::Digraph)
  return DirectedDigraph(DigraphWrap.DigraphJoin(GapObj(D1), GapObj(D2)))
end

function digraph_union(D1::Digraph, D2::Digraph)
  return DirectedDigraph(DigraphWrap.DigraphEdgeUnion(GapObj(D1), GapObj(D2)))
end

function digraph_cartesian_product(D1::Digraph, D2::Digraph)
  return DirectedDigraph(DigraphWrap.DigraphCartesianProduct(GapObj(D1), GapObj(D2)))
end

function digraph_direct_product(D1::Digraph, D2::Digraph)
  return DirectedDigraph(DigraphWrap.DigraphDirectProduct(GapObj(D1), GapObj(D2)))
end

function digraph_lexicographic_product(D1::Digraph, D2::Digraph)
  return DirectedDigraph(DigraphWrap.DigraphLexicographicProduct(GapObj(D1), GapObj(D2)))
end

function digraph_modular_product(D1::Digraph, D2::Digraph)
  return DirectedDigraph(DigraphWrap.ModularProduct(GapObj(D1), GapObj(D2)))
end

function digraph_strong_product(D1::Digraph, D2::Digraph)
  return DirectedDigraph(DigraphWrap.StrongProduct(GapObj(D1), GapObj(D2)))
end


