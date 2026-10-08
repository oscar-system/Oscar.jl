@register_serialization_type UniformHypergraph

function save_object(s::SerializerState, K::UniformHypergraph)
  save_object(s, faces(K))
end

function load_object(s::DeserializerState, K::Type{UniformHypergraph})
  return uniform_hypergraph(load_object(s, Vector{Vector{Int}}))
end
