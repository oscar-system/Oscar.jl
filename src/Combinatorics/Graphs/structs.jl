import Oscar: Polymake
import Oscar.Polymake: Directed, Undirected

struct Graph{T <: Union{Directed, Undirected}}
    pm_graph::Polymake.Graph{T}
end

mutable struct PhylogeneticTree{T <: Union{Float64, Int, QQFieldElem}}
  pm_ptree::Polymake.LibPolymake.BigObjectAllocated
end

