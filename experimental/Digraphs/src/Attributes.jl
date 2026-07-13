function nv(D::Digraph)
  return Int(D.nv)
end

function ne(D::Digraph)
  return Int(D.ne)
end

function out_neighbours(D::Digraph)
  return Vector{Vector{Int}}(DigraphWrap.OutNeighbours(GapObj(D)))
end

function in_neighbours(D::Digraph)
  return Vector{Vector{Int}}(DigraphWrap.InNeighbours(GapObj(D)))
end

function out_degrees(D::Digraph)
  return Vector{Int}(DigraphWrap.OutDegrees(GapObj(D)))
end

function in_degrees(D::Digraph)
  return Vector{Int}(DigraphWrap.InDegrees(GapObj(D)))
end

function adjacency_matrix(D::Digraph)
  return Matrix{Int}(DigraphWrap.AdjacencyMatrix(GapObj(D)))
end

function laplacian_matrix(D::Digraph)
  return Matrix{Int}(DigraphWrap.LaplacianMatrix(GapObj(D)))
end

function chromatic_number(D::Digraph)
  return Int(DigraphWrap.ChromaticNumber(GapObj(D))::GapInt)
end

function clique_number(D::Digraph)
  return Int(DigraphWrap.CliqueNumber(GapObj(D))::GapInt)
end

function has_edge(D::Digraph, s::Int, t::Int)
  adj = adjacency_matrix(D)
  return adj[s, t] != 0
end

function has_vertex(D::Digraph, v::Int)
  return 1 <= v <= nv(D)
end

function shortest_distances(D::Digraph)
  return Matrix{Int}(DigraphWrap.DigraphShortestDistances(GapObj(D)))
end

function digraph_diameter(D::Digraph)
  return Int(DigraphWrap.DigraphDiameter(GapObj(D))::GapInt)
end


