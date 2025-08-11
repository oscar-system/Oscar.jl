import Oscar.Serialization: save_object, load_object

@register_serialization_type GaussianGraphicalModel uses_id [:parameter_ring, :model_ring]


function save_object(s::SerializerState, M::GraphicalModel{T, L}) where {T, L}
  save_data_dict(s) do
    save_object(s, graph(M), :graph)
  end
end
