import Oscar: Polymake
import Oscar.Polymake:
  Directed, Undirected,
  EdgeMap, NodeMap

@attributes mutable struct Graph{T <: Union{Directed, Undirected}}
  pm_graph::Polymake.Graph{T}
  
  function Graph{T}(nverts::Int) where T  <: Union{Directed, Undirected}
    return new{T}(Polymake.Graph{T}(nverts))
  end
  function Graph{T}(pmg::Polymake.Graph{T}) where T <: Union{Directed, Undirected}
    return new{T}(pmg)
  end
end

struct GraphMap{T, S <: Union{Nothing, EdgeMap}, U <: Union{Nothing, NodeMap}}
  graph::Graph{T}
  edge_map::S
  vertex_map::U
end

# here so we can pass Mixed to the graph constructor
struct Mixed end

@attributes mutable struct MixedGraph
  directed_component::Polymake.Graph{Directed}
  undirected_component::Polymake.Graph{Undirected}

  function MixedGraph(nverts::Int)
    return new(Graph{Directed}(nverts), Graph{Undirected}(nverts))
  end
end
