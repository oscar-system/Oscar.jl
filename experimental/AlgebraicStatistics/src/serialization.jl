import Oscar.Serialization: save_object, load_object,
  type_and_params, params, load_type_and_params, decode_type

################################################################################
# Special Dict Types
@register_serialization_type GraphDict 
@register_serialization_type GraphTransDict
@register_serialization_type GenDict

function type_and_params(obj::T) where T <: Union{GraphDict, GraphTransDict}
  if isempty(obj)
    return TypeAndParams(
      T,
      nothing
    )
  end
  
  value_params = type_and_params.(collect(values(obj)))
  @req allequal(value_params) "Not all params of values in $obj are the same"
  
  return TypeAndParams(
    T,
    first(value_params)
  )
end

type_and_params(D::GenDict) = TypeAndParams(GenDict, params(type_and_params(D.d)))

function save_object(s::SerializerState, e::Edge)
  save_data_array(s) do
    save_object(s, src(e))
    save_object(s, dst(e))
  end
end

load_object(s::DeserializerState, ::Type{Edge}) = Edge(load_object(s, Vector{Int})...)

function save_object(s::SerializerState, d::T) where T <: Union{GraphDict, GraphTransDict, GenDict}
  save_data_array(s) do
    for (k, v) in d
      save_object(s, (k, v))
    end
  end
end

function load_object(s::DeserializerState, tp::TypeAndParams{GraphDict, <:Ring})
  R = params(tp)
  graph_gen_dict = Dict{Union{Int, Edge}, elem_type(R)}()
  load_array_node(s) do _
    key = load_node(s, 1) do
      if s.obj isa Oscar.Serialization.JSON3.Array || s.obj isa Vector
        return load_object(s, Edge)
      else
        return load_object(s, Int)
      end
    end
    graph_gen_dict[key] = load_object(s, TypeAndParams(MPolyRingElem, R), 2)
  end
  return GraphDict{elem_type(R)}(graph_gen_dict)
end

# might need to have more type specification in the future here
# for now we know that the params are a dict with domain and codomain
function load_object(s::DeserializerState, tp::TypeAndParams{GraphDict, <:Tuple{Vararg{Pair}}})
  cdom = tp[:codomain]
  dom = tp[:domain]
  map_type = Oscar.MPolyAnyMap{typeof(dom), typeof(cdom)}
  graph_gen_dict = Dict{Union{Int, Edge}, map_type}()
  load_array_node(s) do _
    key = load_node(s, 1) do
      if s.obj isa Oscar.Serialization.JSON3.Array || s.obj isa Vector
        return load_object(s, Edge)
      else
        return load_object(s, Int)
      end
    end
    graph_gen_dict[key] = load_object(s, TypeAndParams(Oscar.MPolyAnyMap, params(tp)...), 2)
  end
  return GraphDict{map_type}(graph_gen_dict)
end

function load_object(s::DeserializerState, tp::TypeAndParams{GraphTransDict, <:Ring})
  R = params(tp)
  graph_trans_dict = Dict{Tuple{Symbol, Edge}, elem_type(R)}()
  load_array_node(s) do _
    key = load_object(s, Tuple{Symbol, Edge}, 1)
    graph_trans_dict[key] = load_object(s, TypeAndParams(MPolyRingElem, R), 2)
  end
  return GraphTransDict{elem_type(R)}(graph_trans_dict)
end

function load_type_and_params(s::DeserializerState, T::Type{GenDict})
  return load_node(s, :params) do
    key_tp = load_node(s, :key_params) do
      s.obj isa String ? TypeAndParams(decode_type(s), nothing) : load_type_and_params(s, decode_type(s))
    end
    value_tp = load_node(s, :value_params) do
      load_type_and_params(s, decode_type(s))
    end
    return TypeAndParams(GenDict{type(key_tp)}, :key_params => key_tp, :value_params => value_tp)
  end
end

function load_object(s::DeserializerState, tp::TypeAndParams{GenDict{S}, <:Tuple{Vararg{Pair}}}) where S
  return GenDict(load_object(s, TypeAndParams(Dict{S, MPolyRingElem}, params(tp)...)))
end

################################################################################
# Model Types

@register_serialization_type GaussianGraphicalModel uses_id [:parameter_ring, :model_ring]
@register_serialization_type DiscreteGraphicalModel uses_id [:parameter_ring, :model_ring]

function type_and_params(GM::S) where {T, L, S <: GraphicalModel{T, L}}
  # this should only be the graph type not the whole graph
  # need to make adjustments to TypeAndParams functionality
  TypeAndParams(S, :graph_type => TypeAndParams(typeof(graph(GM)), nothing))
end

function save_object(s::SerializerState, M::GraphicalModel)
  save_object(s, graph(M))
end

function load_object(s::DeserializerState, tp::TypeAndParams{GaussianGraphicalModel, <:Tuple{Vararg{Pair}}})
  g = load_object(s, tp[:graph_type])
  gaussian_graphical_model(g)
end

function save_object(s::SerializerState, M::DiscreteGraphicalModel)
  save_data_dict(s) do
    save_object(s, graph(M), :graph)
    save_object(s, states(M), :states)
  end
end

function load_object(s::DeserializerState, tp::TypeAndParams{DiscreteGraphicalModel, <:Tuple{Vararg{Pair}}})
  discrete_graphical_model(
    load_object(s, tp[:graph_type], :graph),
    load_object(s, Vector{Int}, :states)
  )
end

# needs to use id to have attributes
@register_serialization_type PhylogeneticModel uses_id [:parameter_ring,
                                                        :model_ring,
                                                        :full_model_ring]

type_and_params(pm::PhylogeneticModel) = TypeAndParams(
  PhylogeneticModel,
  :base_field => base_field(pm),
  # needed until serialization can handle types as parameters
  :graph_type => TypeAndParams(typeof(graph(pm)), nothing), 
  :graph_params => type_and_params(graph(pm)),
  :model_parameter_name_type => TypeAndParams(typeof(varnames(pm)), nothing),
  :transition_matrix_entry_type => TypeAndParams(eltype(transition_matrix(pm)), nothing),
  :transition_matrix_params => type_and_params(transition_matrix(pm)),
  :root_distribution_entry_type => TypeAndParams(eltype(root_distribution(pm)), nothing),
  :root_distribution_params => type_and_params(root_distribution(pm))
)

function save_object(s::SerializerState, pm::PhylogeneticModel)
  save_data_dict(s) do
    save_object(s, graph(pm), :graph)
    save_object(s, transition_matrix(pm), :transition_matrix)
    save_object(s, varnames(pm), :model_parameter_name)
    save_object(s, n_states(pm), :n_states)
    save_object(s, root_distribution(pm), :root_distribution)
  end
end

function load_object(s::DeserializerState, tp::TypeAndParams{PhylogeneticModel, <:Tuple{Vararg{Pair}}})
  return PhylogeneticModel(
    tp[:base_field],
    load_object(s, tp[:graph_params], :graph),
    load_object(s, tp[:transition_matrix_params], :transition_matrix),
    load_object(s, tp[:root_distribution_params], :root_distribution),
    load_object(s, type(tp[:model_parameter_name_type]), :model_parameter_name)
  )
end

@register_serialization_type GroupBasedPhylogeneticModel uses_id [:parameter_ring,
                                                                  :model_ring,
                                                                  :full_model_ring]

type_and_params(pm::GroupBasedPhylogeneticModel) = TypeAndParams(
  GroupBasedPhylogeneticModel,
  :phylo_model => phylogenetic_model(pm),
  # see comment in GroupBasedPhylogeneticModel constructor about group
  :group => parent(first(group(pm))),
  :model_parameter_name_type => TypeAndParams(typeof(varnames(pm)), nothing),
)

function save_object(s::SerializerState, pm::GroupBasedPhylogeneticModel)
  save_data_dict(s) do
    save_object(s, fourier_parameters(pm), :fourier_parameters)
    save_object(s, group(pm), :group_elems)
    save_object(s, varnames(pm), :varnames_group_based)
  end
end

function load_object(s::DeserializerState, tp::TypeAndParams{GroupBasedPhylogeneticModel, <:Tuple{Vararg{Pair}}})
  GroupBasedPhylogeneticModel(tp[:phylo_model],
                              load_object(s, Vector{Symbol}, :fourier_parameters),
                              load_object(s, TypeAndParams(Vector{FinGenAbGroupElem}, tp[:group]), :group_elems),
                              load_object(s, type(tp[:model_parameter_name_type]), :varnames_group_based))
end


################################################################################
# IndexedRing
@register_serialization_type IndexedRing

function type_and_params(IR::IndexedRing)
  if isempty(IR.gen_to_index)
    return TypeAndParams(IndexedRing,
               :ring => base_ring(IR),
               :index_type => TypeAndParams(Int, nothing))
  end
  return TypeAndParams(IndexedRing,
                    :ring => base_ring(IR),
                    :index_type => type_and_params(first(IR.gen_to_index).second))
end

save_object(s::SerializerState, R::IndexedRing) = save_object(s, R.gen_to_index)

function load_object(s::DeserializerState, tp::TypeAndParams{<:IndexedRing, <:Tuple{Vararg{Pair}}})
  R = tp[:ring]
  index_type = type(tp[:index_type])
  gen_to_index = load_object(s, TypeAndParams(Dict{elem_type(R), index_type},
                                           :key_params => TypeAndParams(elem_type(R), R),
                                           :value_params => TypeAndParams(index_type, nothing)))
  return IndexedRing(R, gen_to_index)
end

################################################################################
# Phylogenetic Networks
@register_serialization_type PhylogeneticNetwork

type_and_params(::PhylogeneticNetwork) = nothing

function save_object(s::SerializerState, pn::PhylogeneticNetwork)
  save_object(s, graph(pn))
end

function load_object(s::DeserializerState, ::Type{PhylogeneticNetwork})
  return phylogenetic_network(load_object(s, Graph{Directed}))
end
