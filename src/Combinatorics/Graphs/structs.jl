import Oscar: Polymake
import Oscar.Polymake:
  Directed, Undirected,
  EdgeMap, NodeMap

@attributes mutable struct Graph{T <: Union{Directed, Undirected}}
  pm_graph::Polymake.Graph{T}
  
  function Graph{T}(nverts::Int) where T  <: Union{Directed, Undirected}
    return new{T}(Polymake.Graph{T}(nverts))
  end
end

