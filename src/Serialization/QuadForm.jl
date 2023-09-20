############################################################
# QuadSpace
@register_serialization_type Hecke.QuadSpace uses_params

function save_type_params(s::SerializerState, V::Hecke.QuadSpace)
  save_data_dict(s) do
    save_object(s, encode_type(Hecke.QuadSpace), :name)
    save_type_params(s, gram_matrix(V), :params)
  end
end

function load_type_params(s::DeserializerState, ::Type{<: Hecke.QuadSpace}, str::String)
  return load_ref(s, str)
end

function load_type_params(s::DeserializerState, ::Type{<:Hecke.QuadSpace}, dict::Dict)
  return load_type_params(s, MatElem, dict[:params])
end

function save_object(s::SerializerState, V::Hecke.QuadSpace)
  save_object(s, gram_matrix(V))
end

function load_object(s::DeserializerState, ::Type{<:Hecke.QuadSpace},
                     entries::Vector, params::MatSpace)
  gram = load_object(s, MatElem, entries, params)
  F =  base_ring(params)
  return quadratic_space(F, gram)
end

############################################################
# ZZLat
@register_serialization_type ZZLat

function save_object(s::SerializerState, L::ZZLat)
  save_data_dict(s) do
    save_typed_object(s, basis_matrix(L), :basis)
    save_typed_object(s, ambient_space(L), :ambient_space)
  end
end

function load_object(s::DeserializerState, ::Type{ZZLat}, dict::Dict)
  B = load_typed_object(s, dict[:basis])
  V = load_typed_object(s, dict[:ambient_space])
  return lattice(V, B)
end
