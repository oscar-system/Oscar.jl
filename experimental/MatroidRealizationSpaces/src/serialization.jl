import Oscar.Serialization:@register_serialization_type, save_object, load_object, type_params, TypeParams
@register_serialization_type MatroidRealizationSpace uses_id
type_params(M::MatroidRealizationSpace) = TypeParams(MatroidRealizationSpace,:ambient_ring=>ambient_ring(M),:ground_ring=>M.ground_ring)

function save_object(s::SerializerState, M::MatroidRealizationSpace) 
  save_data_dict(s) do
    save_object(s, defining_ideal(M), :defining_ideal)
    save_object(s, inequations(M), :inequations)
    if isnothing(realization_matrix(M))
         save_object(s, matrix(ambient_ring(M),0,0,[]), :realization_matrix)
    else
        save_object(s, realization_matrix(M), :realization_matrix)
    end
    save_object(s, M.one_realization, :one_realization)
  end
end


function load_object(s::DeserializerState, ::Type{<:MatroidRealizationSpace}, dict::Dict)
  R = dict[:ambient_ring];
  GR = dict[:ground_ring];
  I = load_object(s, MPolyIdeal, R, :defining_ideal)
  Ineqs = load_object(s, Vector{RingElem}, R, :inequations)
  Mat = load_object(s, Union{MatElem,Nothing}, R, :realizataion_matrix)
  if isempty(Mat)
    RMat = nothing;
  else 
    RMat = Mat
  end
  char = char(coefficient_ring(R))
  if !(char == 0)
    q = order(coefficient_ring(R))
  else
    q = nothing
  end
  RS = MatroidRealizationSpace(I, Ineqs, R, RMat, char, q, GR, false)
  one_real = load_object(s, Bool, :one_realization)
  RS.one_realization = one_real
  return RS
end