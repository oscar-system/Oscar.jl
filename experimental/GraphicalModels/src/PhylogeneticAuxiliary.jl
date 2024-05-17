#### GROUP OPERATIONS #### 

# This function implements the group operation on Z2 x Z2 and is subsequently 
# used to check the condition given in the computation of the Fourier
# coordinates and to calculate the sum of the leaves beneath a given edge. 

function group_sum(pm::GroupBasedPhylogeneticModel, states::Vector{Int})
  group = group_of_model(pm)
  return sum(group[states]).%2
end

function is_zero_group_sum(pm::GroupBasedPhylogeneticModel, states::Vector{Int})
  ng = length(states)
  return group_sum(pm, [states[1]]) == group_sum(pm, states[2:ng])
end
  
function which_group_element(pm::GroupBasedPhylogeneticModel, elem::Vector{Int64})
  group = group_of_model(pm)
  return findall([all(group[i].==elem) for i in 1:length(group)])[1]
end


#### OPERATIONS ON THE TREE #### 

function interior_nodes(graph::Graph)
  big_graph = Polymake.graph.Graph(ADJACENCY = pm_object(graph))
  degrees = big_graph.NODE_DEGREES
  return findall(x -> x > 1, degrees)
end

function leaves(graph::Graph)
  big_graph = Polymake.graph.Graph(ADJACENCY = pm_object(graph))
  degrees = big_graph.NODE_DEGREES
  return findall(x -> x == 1, degrees)
end

function vertex_descendants(v::Int, gr::Graph, desc::Vector{Any})
  lvs = leaves(gr)
  outn = outneighbors(gr, v)

  if v in lvs
    return [v]
  end

  innodes = setdiff(outn, lvs)
  d = unique(append!(desc, intersect(outn, lvs)))
  
  if length(innodes) > 0
    for i in innodes
      d = vertex_descendants(i, gr, d)
    end
    return d 
  end

  return d
end

function cherries(graph::Graph)
  lvs = leaves(graph)
  cherry = []
  for l in lvs
    in_node = inneighbors(graph,l)[1]
    lvs_cherr = outneighbors(graph, inneighbors(graph,l)[1])
    if issubset(lvs_cherr, lvs) == 2
      cherry = append!(cherry, [[Edge(in_node, lvs_cherr[1]), Edge(in_node, lvs_cherr[2])]])
    end
  end
  
  return unique(cherry)
end

function root(graph::Graph)
  n_parents = [length(inneighbors(graph, v)) for v in 1:n_vertices(graph)]
  return findall(x -> x == 0, n_parents)[1]
end

function order_edges(graph::Graph)
  edgs = collect(edges(graph))
  leaves_idx = findall(edge -> dst(edge) in Oscar.leaves(graph), edgs)
  return edgs[vcat(leaves_idx, setdiff(1:length(edgs), leaves_idx))]
end

## COMPUTE EQUIVALENCE CLASSES ##
# Given the dictionary of parametrization, the functions below compute 
# the equivalence classes, one time as is, the second one as in the format 
# on the website: every polynomial is multiplied by 0.25 times the size of the equivalence class. 
@doc raw"""
    compute_equivalent_classes(pm::GroupBasedPhylogeneticModel)  

Given the parametrization of a `PhylogeneticModel`, cancel all duplicate entries and return equivalence classes of states which are attached the same probabilities.

# Examples
```jldoctest
julia> pm = jukes_cantor_model(graph_from_edges(Directed,[[4,1],[4,2],[4,3]]));

julia> compute_equivalent_classes(pm)
Dict{Vector{Tuple{Int64, Int64, Int64}}, QQMPolyRingElem} with 5 entries:
 [(1, 2, 3), (3, 2, 4), (1, 3, 2), (1, 4, 2), (3, 4, 2), (3, 1, 2), (3, 4, … => 1//4*a[1]*b[2]*b[3] + 1//4*a[2]*b[1]*b[3] + 1//4*a[3]*b[1]*b[2] + 1//4*b[1]*b[2]*b[3]
 [(2, 1, 1), (4, 3, 3), (1, 2, 2), (3, 2, 2), (2, 4, 4), (2, 3, 3), (3, 1, … => 1//4*a[1]*b[2]*b[3] + 1//4*a[2]*a[3]*b[1] + 1//2*b[1]*b[2]*b[3]
 [(4, 4, 4), (1, 1, 1), (3, 3, 3), (2, 2, 2)]                                => 1//4*a[1]*a[2]*a[3] + 3//4*b[1]*b[2]*b[3]
 [(3, 1, 3), (3, 2, 3), (4, 3, 4), (1, 3, 1), (1, 4, 1), (4, 1, 4), (4, 2, … => 1//4*a[1]*a[3]*b[2] + 1//4*a[2]*b[1]*b[3] + 1//2*b[1]*b[2]*b[3]
 [(2, 2, 1), (3, 3, 2), (4, 4, 3), (1, 1, 2), (3, 3, 1), (4, 4, 2), (2, 2, … => 1//4*a[1]*a[2]*b[3] + 1//4*a[3]*b[1]*b[2] + 1//2*b[1]*b[2]*b[3]
```
"""
function compute_equivalent_classes(parametrization::Dict{Tuple{Vararg{Int64}}, QQMPolyRingElem})
  polys = unique(collect(values(parametrization)))
  
  equivalent_keys = []
  for value in polys
      eqv_class = [key for key in keys(parametrization) if parametrization[key] == value]
      #sort!(eqv_class)
      append!(equivalent_keys, [eqv_class])
  end
  equivalenceclass_dictionary = Dict{Vector{Tuple{Vararg{Int64}}}, QQMPolyRingElem}(equivalent_keys[i] => parametrization[equivalent_keys[i][1]] for i in 1:length(equivalent_keys))
  return equivalenceclass_dictionary
end

@doc raw"""
    sum_equivalent_classes(pm::PhylogeneticModel)  

Take the output of the function `compute_equivalent_classes` for `PhylogeneticModel` and multiply by a factor to obtain probabilities as specified on the original small trees database.

# Examples
```jldoctest
julia> pm = jukes_cantor_model(graph_from_edges(Directed,[[4,1],[4,2],[4,3]]));

julia> sum_equivalent_classes(pm)
Dict{Vector{Tuple{Int64, Int64, Int64}}, QQMPolyRingElem} with 5 entries:
 [(1, 2, 3), (3, 2, 4), (1, 3, 2), (1, 4, 2), (3, 4, 2), (3, 1, 2), (3, 4, 1), (4, 1, 3… => 6*a[1]*b[2]*b[3] + 6*a[2]*b[1]*b[3] + 6*a[3]*b[1]*b[2] + 6*b[1]*b[2]*b[3]
 [(2, 1, 1), (4, 3, 3), (1, 2, 2), (3, 2, 2), (2, 4, 4), (2, 3, 3), (3, 1, 1), (4, 2, 2… => 3*a[1]*b[2]*b[3] + 3*a[2]*a[3]*b[1] + 6*b[1]*b[2]*b[3]
 [(4, 4, 4), (1, 1, 1), (3, 3, 3), (2, 2, 2)]                                            => a[1]*a[2]*a[3] + 3*b[1]*b[2]*b[3]
 [(3, 1, 3), (3, 2, 3), (4, 3, 4), (1, 3, 1), (1, 4, 1), (4, 1, 4), (4, 2, 4), (1, 2, 1… => 3*a[1]*a[3]*b[2] + 3*a[2]*b[1]*b[3] + 6*b[1]*b[2]*b[3]
 [(2, 2, 1), (3, 3, 2), (4, 4, 3), (1, 1, 2), (3, 3, 1), (4, 4, 2), (2, 2, 4), (4, 4, 1… => 3*a[1]*a[2]*b[3] + 3*a[3]*b[1]*b[2] + 6*b[1]*b[2]*b[3]
```
"""
function sum_equivalent_classes(equivalent_classes::Dict{Vector{Tuple{Vararg{Int64}}}, QQMPolyRingElem})
  return Dict(key => equivalent_classes[key]*length(vcat([key]...)) for key in keys(equivalent_classes))
end



