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


@register_serialization_type MatroidRealizationSpaceSelfProjecting uses_id
type_params(M::MatroidRealizationSpaceSelfProjecting) = TypeParams(MatroidRealizationSpaceSelfProjecting,
                                                     :matrix_space=>parent(selfproj_realization_matrix(M)),
                                                     :ground_ring=>M.ground_ring)

function save_object(s::SerializerState, M::MatroidRealizationSpaceSelfProjecting) 
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


function load_object(s::DeserializerState, ::Type{<:MatroidRealizationSpaceSelfProjecting}, dict::Dict)
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
  return MatroidRealizationSpaceSelfProjecting(I, Ineqs, R, RMat, char, q, GR)
end


@register_serialization_type MatroidRealizations uses_id
type_params(M::MatroidRealizations) = TypeParams(MatroidRealizations, 
                                        :matrix_space_mrs => parent(realization_matrix(realization_space(M))), 
                                        :matrix_space_sp_mrs => parent(selfproj_realization_matrix(selfproj_realization_space(M))),
                                        :ground_ring => realization_space(M).ground_ring,
                                        :ground_ring_s => selfproj_realization_space(M).ground_ring,
)

function save_object(s::SerializerState, M::MatroidRealizations) 
  save_data_dict(s) do
    save_object(s, M.name, :name)
    save_object(s, M.matroid, :matroid)
    save_object(s, M.rk, :rk)
    save_object(s, M.length_groundset, :length_groundset)
    save_object(s, M.dim_r, :dim_r)
    save_object(s, M.dim_s, :dim_s)
    save_object(s, M.equal, :equal)
    save_object(s, M.realization_space, :realization_space)
    save_object(s, M.selfproj_realization_space, :selfproj_realization_space)
  end
end


function load_object(s::DeserializerState, ::Type{<:MatroidRealizations}, dict::Dict)
  dictionary_r = Dict(:ground_ring => dict[:ground_ring], :matrix_space => dict[:matrix_space_mrs])
  dictionary_s = Dict(:ground_ring => dict[:ground_ring_s], :matrix_space => dict[:matrix_space_sp_mrs])
#  GR = dict[:ground_ring]
#  GRS = dict[:ground_ring_s]
#  R_MRS = dict[:matrix_space_mrs]
#  RSP_MRS = dict[:matrix_space_sp_mrs]
  str = load_object(s, String, :name)
  m = load_object(s, Matroid, :matroid)
  rk = load_object(s, Int, :rk)
  n = load_object(s, Int, :length_groundset)
  dimR = load_object(s, Int, :dim_r)
  dimS = load_object(s, Int, :dim_s)
  boo = load_object(s, Bool, :equal)
  RS = load_object(s, MatroidRealizationSpace, dictionary_r, :realization_space)
  RSSP = load_object(s, MatroidRealizationSpaceSelfProjecting, dictionary_s, :selfproj_realization_space)
  return MatroidRealizations(str, m, rk, n, RS, dimR, RSSP, dimS, boo)
end
