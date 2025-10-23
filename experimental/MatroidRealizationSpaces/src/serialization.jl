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

@register_serialization_type MatroidRealizations uses_id
type_params(M::MatroidRealizations) = TypeParams(MatroidRealizations,   
                                                 :selfproj_mrs=>selfproj_realization_space(M),
                                                 :mrs=>realization_space(M))

function save_object(s::SerializerState, M::MatroidRealizations) 
  save_data_dict(s) do
    save_object(s, M.name, :name)
    save_object(s, M.matroid, :matroid)
    save_object(s, M.rk, :rk)
    save_object(s, M.length_groundset, :length_groundset)
    save_object(s, M.dim_r, :dim_r)
    save_object(s, M.dim_s, :dim_s)
    save_object(s, M.equal, :equal)
  end
end


function load_object(s::DeserializerState, ::Type{<:MatroidRealizations}, dict::Dict)
  RS = dict[:mrs]
  RSSP = dict[:selfproj_mrs]

  str = load_object(s, String, :name)
  m = load_object(s, Matroid, :matroid)
  rk = load_object(s, Int, :rk)
  n = load_object(s, Int, :length_groundset)
  dimR = load_object(s, Int, :dim_r)
  dimS = load_object(s, Int, :dim_s)
  boo = load_object(s, Bool, :equal)
  
  return MatroidRealizations(str, m, rk, n, RS, dimR, RSSP, dimS, boo)
end
