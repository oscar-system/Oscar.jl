using JSON

###############################################################################
## Graphs
###############################################################################
@register_serialization_type Graph{Directed} "Graph{Directed}"
@register_serialization_type Graph{Undirected} "Graph{Undirected}"

function save_object(s::SerializerState, g::Graph{T}) where T <: Union{Directed, Undirected}
  smallobject = pm_object(g)
  serialized = Polymake.call_function(Symbol("Core::Serializer"), :serialize, smallobject)
  jsonstr = Polymake.call_function(:common, :encode_json, serialized)
  save_data_json(s, jsonstr)
end


function load_object(s::DeserializerState, g::Type{Graph{T}}) where T <: Union{Directed, Undirected}
  smallobj = Polymake.call_function(:common, :deserialize_json_string, json(s.obj))
  return g(smallobj)
end

###############################################################################
## Matroid
###############################################################################
@register_serialization_type Matroid "Matroid"

function save_object(s::SerializerState, m::Matroid) 
  @req isa(m.groundset, Vector{Int}) "Groundset must be a Vector{Int}"

  save_data_dict(s) do
    save_object(s, pm_object(m), :matroid)
    save_object(s, m.groundset, :groundset)
  end
end


function load_object(s::DeserializerState, m::Type{Matroid})
  mt = load_object(s, Polymake.BigObject, :matroid)
  grds = load_object(s, Vector{Int}, :groundset)
  return m(mt, grds)
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
  IM = Polymake.call_function(:common, :deserialize_json_string, json(s.obj))
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
  bigobject = Polymake.call_function(:common, :deserialize_json_string, json(s.obj))
  return K(bigobject)
end

###############################################################################
## Phylogenetic Trees
###############################################################################
@register_serialization_type PhylogeneticTree 

function save_object(s::SerializerState, PT::PhylogeneticTree)
  save_object(s, pm_object(PT))
end

function load_object(s::DeserializerState, T::Type{PhylogeneticTree}, dict::Dict)
  bigobject = Polymake.call_function(:common, :deserialize_json_string, json(dict))
  return T(bigobject)
end
