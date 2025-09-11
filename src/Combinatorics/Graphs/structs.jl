import Oscar: Polymake
import Oscar.Polymake:
  Directed, Undirected,
  EdgeMap, NodeMap

struct Mixed end

const GraphTypes = Union{Directed, Undirected, Mixed}
abstract type AbstractGraph{T <: GraphTypes} end

@attributes mutable struct Graph{T <: Union{Directed, Undirected}} <: AbstractGraph{T}
  pm_graph::Polymake.Graph{T}
  
  function Graph{T}(nverts::Int) where T  <: Union{Directed, Undirected}
    return new{T}(Polymake.Graph{T}(nverts))
  end
  function Graph{T}(pmg::Polymake.Graph{T}) where T <: Union{Directed, Undirected}
    return new{T}(pmg)
  end
end

const GraphMapValueTypes = Union{String, Bool, Int, QQFieldElem, Float64}

struct GraphMap{T, S <: Union{Nothing, EdgeMap}, U <: Union{Nothing, NodeMap}}
  graph::Graph{T}
  edge_map::S
  vertex_map::U
end

@attributes mutable struct MixedGraph <: AbstractGraph{Mixed}
  directed_component::Graph{Directed}
  undirected_component::Graph{Undirected}

  function MixedGraph(nverts::Int)
    return new(Graph{Directed}(nverts), Graph{Undirected}(nverts))
  end
end
