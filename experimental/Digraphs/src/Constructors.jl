@doc raw"""
    digraph(adj::Vector{Vector{Int}}; mut::Bool=false) -> Digraph
    digraph(labels::Vector{<:AbstractString}, source::Vector{<:AbstractString}, range::Vector{<:AbstractString}; mut::Bool=false) -> Digraph
    digraph(n::Int64, source::Vector{Int}, range::Vector{Int}; mut::Bool=false) -> Digraph
    digraph(list::Vector{<:Any}, func::Function; mut::Bool=false) -> Digraph
    digraph(G::T, list::Vector{<:Any}, act::Function, rel::Function; mut::Bool=false) where T<:GAPGroup -> Digraph
    digraph(obj::GapObj) -> Digraph
    digraph(name::T) where T <: AbstractString -> Digraph

# Examples
```jldoctest
julia> d = digraph([[2, 3], [1, 3], [1, 2]])
Digraph with 3 vertices, 6 edges

julia> p = ["a", "b", "c"];

julia> d = digraph(p, p[[1, 1, 2, 2, 3, 3]], p[[2, 3, 1, 3, 1, 2]])
Digraph with 3 vertices, 6 edges

julia> d = digraph(3, [1, 1, 2, 2, 3, 3], [2, 3, 1, 3, 1, 2])
Digraph with 3 vertices, 6 edges

julia> rel(x, y) = x != y;

julia> d = digraph([1, 2, 3], rel)
Digraph with 3 vertices, 6 edges

julia> act(x, g) = x^g;

julia> d = digraph(symmetric_group(3), [1, 2, 3], act, rel)
Digraph with 3 vertices, 6 edges

julia> d = digraph("Diamond")
Digraph with 4 vertices, 10 edges

julia> d = GAP.Globals.CompleteGraph(GapObj(symmetric_group(3)))
GAP: rec( adjacencies := [ [ 2, 3 ] ], group := Sym( [ 1 .. 3 ] ), isGraph := true, isSimple := true, order := 3, representatives := [ 1 ], schreierVector := [ -1, 1, 1 ] )

julia> digraph(d)
Digraph with 3 vertices, 6 edges
```
"""
function digraph(adj::Vector{Vector{Int}}; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  d = DigraphWrap.Digraph(filter, GapObj(adj, recursive = true))
  return Digraph(d)
end

function digraph(labels::Vector{<:AbstractString}, source::Vector{<:AbstractString}, range::Vector{<:AbstractString}; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  d = DigraphWrap.Digraph(filter, GapObj(labels, recursive = true), GapObj(source, recursive = true), GapObj(range, recursive = true))
  return Digraph(d)
end

function digraph(n::Int64, source::Vector{Int}, range::Vector{Int}; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  d = DigraphWrap.Digraph(filter, GapObj(n), GapObj(source), GapObj(range))
  return Digraph(d)
end

function digraph(list::Vector{<:Any}, func::Function; mut::Bool=false)
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  d = DigraphWrap.Digraph(filter, GapObj(list, recursive = true), GapObj(func))
  return Digraph(d)
end

function digraph(G::T, list::Vector{<:Any}, act::Function, rel::Function; mut::Bool=false) where T<:GAPGroup
  filter = mut ? GAP.Globals.IsMutableDigraph : GAP.Globals.IsImmutableDigraph
  d = DigraphWrap.Digraph(filter, GapObj(G), GapObj(list, recursive = true), GapObj(act), GapObj(rel))
  return Digraph(d)
end

function digraph(obj::GapObj)
  d = DigraphWrap.Digraph(obj)
  return Digraph(d)
end

function digraph(name::T) where T <: AbstractString
  d = DigraphWrap.Digraph(GapObj(name))
  return Digraph(d)
end

function null_digraph()
  return Digraph(DigraphWrap.NullDigraph(0))
end

function complete_digraph(n::Int)
  return Digraph(DigraphWrap.CompleteDigraph(n))
end

function complete_bipartite_digraph(m::Int, n::Int)
  return Digraph(DigraphWrap.CompleteBipartiteDigraph(m, n))
end

function cycle_digraph(n::Int)
  return Digraph(DigraphWrap.CycleDigraph(n))
end

function digraph_from_edges(n::Int, edges::Vector{Tuple{Int,Int}})
  gap_edges = GapObj([[s, t] for (s, t) in edges])
  return Digraph(DigraphWrap.DigraphByEdges(gap_edges))
end

function digraph_from_edges(n::Int, edges::Vector{Vector{Int}})
  return Digraph(DigraphWrap.DigraphByEdges(GapObj(edges)))
end

function digraph_from_adjacency_matrix(A::Matrix{Int})
  return Digraph(DigraphWrap.DigraphByAdjacencyMatrix(GapObj(A)))
end


