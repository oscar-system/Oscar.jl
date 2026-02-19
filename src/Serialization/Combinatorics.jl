using Oscar: create_gs2num

###############################################################################
## Graphs
###############################################################################
@register_serialization_type Graph{Directed} "Graph{Directed}"
@register_serialization_type Graph{Undirected} "Graph{Undirected}"

function type_params(g::Graph{T}) where T <: Union{Directed, Undirected}
  isempty(labelings(g)) && return TypeParams(Graph{T}, nothing)
  labeling_types = Dict{Symbol, TypeParams}()
  for l in labelings(g)
    gm = get_attribute(g, l)
    # at least one of edge_map and vertex_map must not be nothing
    # we use the first entry to get the type since this is always defined by polymake
    # and the type of the values of the map is not accessible from the GraphMap type
    if !isnothing(gm.edge_map)
      # currently we can only label using basic types so params are nothing
      labeling_types[l] = TypeParams(typeof(gm.edge_map[first(edges(g))]), nothing)
    else
      labeling_types[l] = TypeParams(typeof(gm.vertex_map[1]), nothing)
    end
  end
  
  return TypeParams(Graph{T}, map(x -> x.first => x.second, collect(pairs(labeling_types)))...)
end

function save_object(s::SerializerState, g::Graph{T}) where T <: Union{Directed, Undirected}
  smallobject = pm_object(g)
  serialized = Polymake.call_function(Symbol("Core::Serializer"), :serialize, smallobject)
  jsonstr = Polymake.call_function(:common, :encode_json, serialized)

  if isempty(labelings(g))
    save_data_json(s, jsonstr)
  else
    save_data_dict(s) do
      save_data_json(s, jsonstr, :graph)
      labels_d = Dict{Symbol, NamedTuple}()
      for l in labelings(g)
        gm = get_attribute(g, l)
        if !isnothing(gm.edge_map)
          em = Dict((src(e), dst(e)) => gm.edge_map[e] for e in edges(g))
        else
          em = nothing
        end
        # QQ has nothing to do with serialization, just how the pm functions was implemented
        vm = !isnothing(gm.vertex_map) ? _pmdata_for_oscar(gm.vertex_map, QQ) : nothing
        
        labels_d[l] = (edge_map = em, vertex_map = vm)
      end

      save_object(s, labels_d, :labelings)
    end
  end
end

function load_object(s::DeserializerState, G::Type{Graph{T}}) where T <: Union{Directed, Undirected}
  smallobj = Polymake.call_function(:common, :deserialize_json_string, JSON.json(s.obj))
  return G(smallobj)
end

function load_object(s::DeserializerState, G::Type{Graph{T}}, params::Dict) where T <: Union{Directed, Undirected}
  g = load_node(s, :graph) do _
    G(Polymake.call_function(:common, :deserialize_json_string, JSON.json(s.obj)))
  end

  load_node(s, :labelings) do _
    for label in keys(params)
      load_node(s, label) do _
        edge_labels = load_node(s, 1) do em
          isnothing(em) ? nothing : load_object(s, Dict{Tuple{Int, Int}, params[label]}, nothing)
        end
        vertex_labels =load_node(s, 2) do vm
          isnothing(vm) ? nothing : load_object(s, Dict{Int, params[label]}, nothing)
        end
        label!(g, edge_labels, vertex_labels; name=label)
      end
    end
  end
  return g
end

###############################################################################
## Matroid
###############################################################################
@register_serialization_type Matroid "Matroid"

function save_object(s::SerializerState, m::Matroid) 
  gs = matroid_groundset(m)
  @req gs isa Vector{Int} "Groundset must be a Vector{Int}"

  save_data_dict(s) do
    save_object(s, pm_object(m), :matroid)
    save_object(s, gs, :groundset)
    save_object(s, create_gs2num(gs), :gs2num)
  end
end


function load_object(s::DeserializerState, m::Type{<:Matroid})
  # TODO: for now only objects of type `Matroid{Int}` are supported, to preserve
  # backwards compatibility.
  mt = load_object(s, Polymake.BigObject, :matroid)
  grds = load_object(s, Vector{Int}, :groundset)
  gs2num = load_object(s, Dict{Int,Int}, :gs2num)
  return Matroid{Int}(mt, grds, gs2num)
end


###############################################################################
## IncidenceMatrix
###############################################################################
@register_serialization_type Polymake.IncidenceMatrixAllocated{Polymake.NonSymmetric}

function save_object(s::SerializerState, IM::IncidenceMatrix)
  serialized = Polymake.call_function(Symbol("Core::Serializer"), :serialize, IM)
  jsonstr = Polymake.call_function(:common, :encode_json, serialized)
  save_data_json(s, jsonstr)
end

function load_object(s::DeserializerState, ::Type{<: IncidenceMatrix})
  IM = Polymake.call_function(:common, :deserialize_json_string, JSON.json(s.obj))
  return IM
end


###############################################################################
## SimplicialComplex
###############################################################################
@register_serialization_type SimplicialComplex

function save_object(s::SerializerState, K::SimplicialComplex)
  save_object(s, pm_object(K))
end

function load_object(s::DeserializerState, K::Type{SimplicialComplex})
  bigobject = Polymake.call_function(:common, :deserialize_json_string, JSON.json(s.obj))
  return K(bigobject)
end

###############################################################################
## Phylogenetic Trees
###############################################################################
@register_serialization_type PhylogeneticTree 
type_params(::PhylogeneticTree{T}) where T <: Union{QQFieldElem, Float64} = TypeParams(
  PhylogeneticTree, T == QQFieldElem ? QQ : AbstractAlgebra.Floats{Float64}())

function save_object(s::SerializerState, PT::PhylogeneticTree)
  save_data_dict(s) do
    save_object(s, pm_object(PT), :pm_tree)
    save_object(s, PT.vertex_perm, :vertex_perm)
  end
end

function load_object(s::DeserializerState, T::Type{<:PhylogeneticTree}, params::QQField)
  inner_object = load_node(s, :pm_tree) do _
    load_from_polymake(Polymake.BigObject, Dict(s.obj))
  end
  vertex_perm = load_object(s, Vector{Int}, :vertex_perm)
  PhylogeneticTree{QQFieldElem}(inner_object, vertex_perm)
end

function load_object(s::DeserializerState, T::Type{<:PhylogeneticTree}, params::AbstractAlgebra.Floats{Float64})
  inner_object = load_node(s, :pm_tree) do _
    load_from_polymake(Polymake.BigObject, Dict(s.obj))
  end
  vertex_perm = load_object(s, Vector{Int}, :vertex_perm)
  PhylogeneticTree{Float64}(inner_object, vertex_perm)
end
