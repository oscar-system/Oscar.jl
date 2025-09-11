import Oscar.Serialization: save_object, load_object,
  type_params

@register_serialization_type GraphGenDict 
@register_serialization_type GraphTransDict

function type_params(obj::T) where T <: Union{GraphGenDict, GraphTransDict}
  if isempty(obj)
    return TypeParams(
      T,
      nothing
    )
  end
  
  value_params = type_params.(collect(values(obj)))
  @req Oscar.Serialization.params_all_equal(value_params) "Not all params of values in $obj are the same"
  
  return TypeParams(
    T,
    first(value_params)
  )
end

function save_object(s::SerializerState, e::Edge)
  save_data_array(s) do
    save_object(s, src(e))
    save_object(s, dst(e))
  end
end

load_object(s::DeserializerState, ::Type{Edge}) = Edge(load_object(s, Vector{Int})...)

function save_object(s::SerializerState, d::T) where T <: Union{GraphGenDict, GraphTransDict}
  save_data_array(s) do
    for (k, v) in d
      save_object(s, (k, v))
    end
  end
end

function load_object(s::DeserializerState, ::Type{GraphGenDict}, R::Ring)
  graph_gen_dict = Dict{Union{Int, Edge}, MPolyRingElem}()
  load_array_node(s) do (_, (k, v))
    if k isa Oscar.Serialization.JSON3.Array
      key = load_object(s, Edge, 1)
    else
      key = load_object(s, Int, 1)
    end
    graph_gen_dict[key] = load_object(s, MPolyRingElem, R, 2)
  end
  return graph_gen_dict
end

function load_object(s::DeserializerState, ::Type{GraphTransDict}, R::Ring)
  graph_trans_dict = Dict{Tuple{Symbol, Edge}, MPolyRingElem}()
  load_array_node(s) do (_, (k, v))
    key = load_object(s, Tuple{Symbol, Edge}, 1)
    graph_trans_dict[key] = load_object(s, MPolyRingElem, R, 2)
  end
  return graph_trans_dict
end

@register_serialization_type GaussianGraphicalModel uses_id [:parameter_ring, :model_ring]

function type_params(GM::S) where {T, L, S <: GraphicalModel{T, L}}
  # this should only be the graph type not the whole graph
  # need to make adjustments to TypeParams functionality
  TypeParams(S, graph(GM))
end

function save_object(s::SerializerState, M::GraphicalModel)
  save_data_dict(s) do
    # needs an empty dict here
    # this requires a known issue with current general serialization
  end
end

load_object(s::DeserializerState,
            ::Type{GaussianGraphicalModel},
            g::AbstractGraph) = gaussian_graphical_model(g)


@register_serialization_type PhylogeneticModel uses_id [:parameter_ring,
                                                        :model_ring,
                                                        :reduced_parameter_ring,
                                                        :reduced_model_ring]

function type_params(pm::S) where {L, S <: PhylogeneticModel{L}}
  # we only need the graph type
  # need to change when serialization can deal with it
  return TypeParams(S, :base_field => base_field(pm), :graph => graph(pm))
end

function save_object(s::SerializerState, pm::PhylogeneticModel)
  save_data_dict(s) do
    save_object(s, transition_matrix(pm), :transition_matrix)
    save_object(s, varname(pm), :model_parameter_name)
    save_object(s, n_states(pm), :n_states)
    save_object(s, root_distribution(pm), :root_distribution)
  end
end

function load_object(s::DeserializerState, ::Type{PhylogeneticModel}, params::Dict)
  return PhylogeneticModel(
    params[:base_field],
    params[:graph],
    # TODO need to change Symbol to VarName  at some point
    load_object(s, Matrix{Symbol}, :transition_matrix),
    load_object(s, Int, :n_states),
    load_object(s, Vector{QQFieldElem}, :root_distribution),
    #TODO also needs to be come VarName
    load_object(s, Symbol, :model_parameter_name)
  )
end
