import Oscar: Polymake
import Oscar.Polymake: Directed, Undirected

struct Graph{T <: Union{Directed, Undirected}}
  pm_graph::Polymake.Graph{T}
end

