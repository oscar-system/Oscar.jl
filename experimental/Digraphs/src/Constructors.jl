function Digraph(out_neighbours::Vector{Vector{Int}})
  return DirectedDigraph(DigraphWrap.DigraphCons(DigraphWrap.IsMutableDigraph(), GapObj(out_neighbours)))
end

function null_digraph()
  return DirectedDigraph(DigraphWrap.NullDigraph(0))
end

function complete_digraph(n::Int)
  return DirectedDigraph(DigraphWrap.CompleteDigraph(n))
end

function complete_bipartite_digraph(m::Int, n::Int)
  return DirectedDigraph(DigraphWrap.CompleteBipartiteDigraph(m, n))
end

function cycle_digraph(n::Int)
  return DirectedDigraph(DigraphWrap.CycleDigraph(n))
end

function digraph_from_edges(n::Int, edges::Vector{Tuple{Int,Int}})
  gap_edges = GapObj([[s, t] for (s, t) in edges])
  return DirectedDigraph(DigraphWrap.DigraphByEdges(gap_edges))
end

function digraph_from_edges(n::Int, edges::Vector{Vector{Int}})
  return DirectedDigraph(DigraphWrap.DigraphByEdges(GapObj(edges)))
end

function digraph_from_adjacency_matrix(A::Matrix{Int})
  return DirectedDigraph(DigraphWrap.DigraphByAdjacencyMatrix(GapObj(A)))
end


