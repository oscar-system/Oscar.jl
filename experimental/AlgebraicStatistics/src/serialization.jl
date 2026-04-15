import Oscar.Serialization: save_object, load_object,
  type_params, _convert_override_params, params, load_type_params, decode_type

function _convert_override_params(tp::TypeParams{<:GraphicalModel, <:Tuple{Vararg{Pair}}})
  param_dict = Dict()
  for p in Oscar.Serialization.params(tp)
    if p.first == :graph_type 
      param_dict[p.first] = Oscar.Serialization.type(p.second)
    else
      param_dict[p.first] = p.second
    end
  end
  return param_dict
end

function _convert_override_params(tp::TypeParams{T, <:Tuple{Vararg{Pair}}}) where T <: Union{GroupBasedPhylogeneticModel, PhylogeneticModel}
  param_dict = Dict()
  type_keys = [:transition_matrix_entry_type, :graph_type,
               :model_parameter_name_type, :root_distribution_entry_type]
  for p in Oscar.Serialization.params(tp)
    if p.first in type_keys
      param_dict[p.first] = Oscar.Serialization.type(p.second)
    else
      param_dict[p.first] = _convert_override_params(p.second)
    end
  end
  return param_dict
end

################################################################################
# Special Dict Types
@register_serialization_type GraphDict 
@register_serialization_type GraphTransDict
@register_serialization_type GenDict

function type_params(obj::T) where T <: Union{GraphDict, GraphTransDict}
  if isempty(obj)
    return TypeParams(
      T,
      nothing
    )
  end
  
  value_params = type_params.(collect(values(obj)))
  @req allequal(value_params) "Not all params of values in $obj are the same"
  
  return TypeParams(
    T,
    first(value_params)
  )
end

type_params(D::GenDict) = TypeParams(GenDict, params(type_params(D.d)))

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

function load_object(s::DeserializerState, ::Type{GraphDict}, R::Ring)
  graph_gen_dict = Dict{Union{Int, Edge}, elem_type(R)}()
  load_array_node(s) do (_, (k, v))
    if k isa Oscar.Serialization.JSON3.Array
      key = load_object(s, Edge, 1)
    else
      key = load_object(s, Int, 1)
    end
    graph_gen_dict[key] = load_object(s, MPolyRingElem, R, 2)
  end
  return GraphDict{elem_type(R)}(graph_gen_dict)
end

# might need to have more type specification in the future here
# for now we know that the params are a dict with domain and codomain
function load_object(s::DeserializerState, ::Type{GraphDict}, d::Dict)
  cdom = d[:codomain]
  dom = d[:domain]
  map_type = Oscar.MPolyAnyMap{typeof(dom), typeof(cdom)}
  graph_gen_dict = Dict{Union{Int, Edge}, map_type}()
  load_array_node(s) do (_, (k, v))
    if k isa Oscar.Serialization.JSON3.Array
      key = load_object(s, Edge, 1)
    else
      key = load_object(s, Int, 1)
    end
    graph_gen_dict[key] = load_object(s, Oscar.MPolyAnyMap, d, 2)
  end
  return GraphDict{map_type}(graph_gen_dict)
end

function load_object(s::DeserializerState, ::Type{GraphTransDict}, R::Ring)
  graph_trans_dict = Dict{Tuple{Symbol, Edge}, elem_type(R)}()
  load_array_node(s) do (_, (k, v))
    key = load_object(s, Tuple{Symbol, Edge}, 1)
    graph_trans_dict[key] = load_object(s, MPolyRingElem, R, 2)
  end
  return GraphTransDict{elem_type(R)}(graph_trans_dict)
end

function load_type_params(s::DeserializerState, T::Type{GenDict})
  subtype, params = load_node(s, :params) do obj
    S, key_params = load_node(s, :key_params) do params
      params isa String && return decode_type(s), nothing
      load_type_params(s, decode_type(s))
    end

    _, value_params = load_node(s, :value_params) do _
      load_type_params(s, decode_type(s))
    end

    return S, Dict(:key_params => key_params, :value_params => value_params)
  end
  return GenDict{subtype}, params
end

function load_object(s::DeserializerState, T::Type{GenDict{S}}, params::Dict) where S
  return GenDict(load_object(s, Dict{S, MPolyRingElem}, params))
end

################################################################################
# Model Types

@register_serialization_type GaussianGraphicalModel uses_id [:parameter_ring, :model_ring]
@register_serialization_type DiscreteGraphicalModel uses_id [:parameter_ring, :model_ring]

function type_params(GM::S) where {T, L, S <: GraphicalModel{T, L}}
  # this should only be the graph type not the whole graph
  # need to make adjustments to TypeParams functionality
  TypeParams(S, :graph_type => TypeParams(typeof(graph(GM)), nothing))
end

function save_object(s::SerializerState, M::GraphicalModel)
  save_object(s, graph(M))
end

function load_object(s::DeserializerState, ::Type{GaussianGraphicalModel}, params::Dict)
  g = load_object(s, params[:graph_type])
  gaussian_graphical_model(g)
end

function save_object(s::SerializerState, M::DiscreteGraphicalModel)
  save_data_dict(s) do
    save_object(s, graph(M), :graph)
    save_object(s, states(M), :states)
  end
end

function load_object(s::DeserializerState, ::Type{DiscreteGraphicalModel}, params::Dict)
  discrete_graphical_model(
    load_object(s, params[:graph_type], :graph),
    load_object(s, Vector{Int}, :states)
  )
end

# needs to use id to have attributes
@register_serialization_type PhylogeneticModel uses_id [:parameter_ring,
                                                        :model_ring,
                                                        :full_model_ring]

type_params(pm::PhylogeneticModel) = TypeParams(
  PhylogeneticModel,
  :base_field => base_field(pm),
  # needed until serialization can handle types as parameters
  :graph_type => TypeParams(typeof(graph(pm)), nothing), 
  :graph_params => type_params(graph(pm)),
  :model_parameter_name_type => TypeParams(typeof(varnames(pm)), nothing),
  :transition_matrix_entry_type => TypeParams(eltype(transition_matrix(pm)), nothing),
  :transition_matrix_params => type_params(transition_matrix(pm)),
  :root_distribution_entry_type => TypeParams(eltype(root_distribution(pm)), nothing),
  :root_distribution_params => type_params(root_distribution(pm))
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

function load_object(s::DeserializerState, ::Type{PhylogeneticModel}, params::Dict)
  T1, p1 = params[:transition_matrix_entry_type], params[:transition_matrix_params]
  T2, p2 = params[:root_distribution_entry_type], params[:root_distribution_params]
  return PhylogeneticModel(
    params[:base_field],
    load_object(s, params[:graph_type], params[:graph_params], :graph),
    load_object(s, Matrix{T1}, p1, :transition_matrix),
    load_object(s, Vector{T2}, p2, :root_distribution),
    load_object(s, params[:model_parameter_name_type], :model_parameter_name)
  )
end

@register_serialization_type GroupBasedPhylogeneticModel uses_id [:parameter_ring,
                                                                  :model_ring,
                                                                  :full_model_ring]

type_params(pm::GroupBasedPhylogeneticModel) = TypeParams(
  GroupBasedPhylogeneticModel,
  :phylo_model => phylogenetic_model(pm),
  # see comment in GroupBasedPhylogeneticModel constructor about group
  :group => parent(first(group(pm))),
  :model_parameter_name_type => TypeParams(typeof(varnames(pm)), nothing),
)

function save_object(s::SerializerState, pm::GroupBasedPhylogeneticModel)
  save_data_dict(s) do
    save_object(s, fourier_parameters(pm), :fourier_parameters)
    save_object(s, group(pm), :group_elems)
    save_object(s, varnames(pm), :varnames_group_based)
  end
end

function load_object(s::DeserializerState, ::Type{GroupBasedPhylogeneticModel}, params::Dict)
  GroupBasedPhylogeneticModel(params[:phylo_model],
                              load_object(s, Vector{Symbol}, :fourier_parameters),
                              load_object(s, Vector{FinGenAbGroupElem}, params[:group], :group_elems),
                              load_object(s, params[:model_parameter_name_type], :varnames_group_based))
end


################################################################################
# IndexedRing
@register_serialization_type IndexedRing

function type_params(IR::IndexedRing)
  if isempty(IR.gen_to_index)
    return TypeParams(IndexedRing,
               :ring => base_ring(IR),
               :index_type => TypeParams(Int, nothing))
  end
  return TypeParams(IndexedRing,
                    :ring => base_ring(IR),
                    :index_type => type_params(first(IR.gen_to_index).second))
end

save_object(s::SerializerState, R::IndexedRing) = save_object(s, R.gen_to_index)

function load_object(s::DeserializerState, ::Type{<:IndexedRing}, params::Dict)
  R = params[:ring]
  if params[:index_type] isa Tuple
    gen_to_index = load_object(s, Dict{elem_type(R), Tuple{[Int for _ in 1:fieldcount(typeof(params[:index_type]))]...}},
                                                           Dict(:key_params => R, :value_params => nothing))
  else
    gen_to_index = load_object(s, Dict{elem_type(R), params[:index_type]}, Dict(:key_params => R))
  end
  return IndexedRing(
    R, gen_to_index
  )
end

################################################################################
# Phylogenetic Networks
@register_serialization_type PhylogeneticNetwork

type_params(::PhylogeneticNetwork) = nothing

function save_object(s::SerializerState, pn::PhylogeneticNetwork)
  save_object(s, graph(pn))
end

function load_object(s::DeserializerState, ::Type{PhylogeneticNetwork})
  return phylogenetic_network(load_object(s, Graph{Directed}))
end
