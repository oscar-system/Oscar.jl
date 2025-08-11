import Oscar.Serialization: save_object, load_object,
  type_params

@register_serialization_type GraphGenDict 

function type_params(obj::GraphGenDict)
  if isempty(obj)
    return TypeParams(
      GraphGenDict,
      nothing
    )
  end
  
  value_params = type_params.(collect(values(obj)))
  @req Oscar.params_all_equal(value_params) "Not all params of values in $obj are the same"
  
  return TypeParams(
    GraphGenDict,
    first(value_params)
  )
end

function save_object(s::SerializerState, e::Edge)
  save_data_array(s) do
    save_object(s, src(e))
    save_object(s, dst(e))
  end
end

function save_object(s::SerializerState, d::GraphGenDict)
  save_data_array(s) do
    for (k, v) in d
      save_object(s, [k, v])
    end
  end
end

@register_serialization_type GaussianGraphicalModel uses_id [:parameter_ring, :model_ring]


function save_object(s::SerializerState, M::GraphicalModel{T, L}) where {T, L}
  save_data_dict(s) do
    save_object(s, graph(M), :graph)
  end
end
