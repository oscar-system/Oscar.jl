import Oscar.Serialization: save_object, load_object, type_params

@register_serialization_type MatroidRealizationSpace uses_id
type_params(M::MatroidRealizationSpace) = TypeParams(MatroidRealizationSpace,
                                                     :matrix_space=>parent(realization_matrix(M)),
                                                     :ground_ring=>M.ground_ring)

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
  MS = dict[:matrix_space];
  R = base_ring(MS)
  GR = dict[:ground_ring];
  I = load_object(s, MPolyIdeal, R, :defining_ideal)
  Ineqs = load_object(s, Vector{MPolyRingElem}, R, :inequations)
  Mat = load_object(s, MatElem, MS, :realization_matrix)
  RMat = isempty(Mat) ? nothing : Mat
  char = characteristic(coefficient_ring(R))
  q = !iszero(char) ? order(coefficient_ring(R)) : nothing

  RS = MatroidRealizationSpace(I, Ineqs, R, RMat, char, q, GR)
  RS.one_realization = load_object(s, Bool, :one_realization)
  
  return RS
end


@register_serialization_type MatroidRealizationSpace_SelfProj uses_id
type_params(M::MatroidRealizationSpace_SelfProj) = TypeParams(MatroidRealizationSpace_SelfProj,
                                                     :matrix_space=>parent(selfproj_realization_matrix(M)),
                                                     :ground_ring=>M.ground_ring)

function save_object(s::SerializerState, M::MatroidRealizationSpace_SelfProj) 
  save_data_dict(s) do
    save_object(s, defining_ideal(M), :defining_ideal)
    save_object(s, inequations(M), :inequations)
    if isnothing(selfproj_realization_matrix(M))
      save_object(s, matrix(ambient_ring(M),0,0,[]), :selfproj_realization_matrix)
    else
      save_object(s, selfproj_realization_matrix(M), :selfproj_realization_matrix)
    end
#    save_object(s, M.one_realization, :one_realization)
  end
end


function load_object(s::DeserializerState, ::Type{<:MatroidRealizationSpace_SelfProj}, dict::Dict)
  MS = dict[:matrix_space];
  R = base_ring(MS)
  GR = dict[:ground_ring];
  I = load_object(s, MPolyIdeal, R, :defining_ideal)
  Ineqs = load_object(s, Vector{MPolyRingElem}, R, :inequations)
  Mat = load_object(s, MatElem, MS, :selfproj_realization_matrix)
  RMat = isempty(Mat) ? nothing : Mat
  char = characteristic(coefficient_ring(R))
  q = !iszero(char) ? order(coefficient_ring(R)) : nothing
 # RS.one_realization = load_object(s, Bool, :one_realization)
  return MatroidRealizationSpace_SelfProj(I, Ineqs, R, RMat, char, q, GR)
end