import Oscar.Serialization: save_object, load_object, type_params

@register_serialization_type MatroidRealizationSpace uses_id
type_params(M::MatroidRealizationSpace) = TypeParams(MatroidRealizationSpace,
                                                     :matrix_space => begin
                                                      mat = realization_matrix(M)
                                                      isnothing(mat) ? nothing : parent(mat)
                                                      end,
                                                      :ideal_ring => base_ring(defining_ideal(M)),
                                                     :ground_ring=>M.ground_ring)

function save_object(s::SerializerState, M::MatroidRealizationSpace) 
  save_data_dict(s) do
    save_object(s, defining_ideal(M), :defining_ideal)
    save_object(s, inequations(M), :inequations)
    if isnothing(realization_matrix(M))
      save_object(s, nothing, :realization_matrix)
    else
      save_object(s, realization_matrix(M), :realization_matrix)
    end
    save_object(s, M.one_realization, :one_realization)
  end
end


function load_object(s::DeserializerState, ::Type{<:MatroidRealizationSpace}, dict::Dict)
  MS = dict[:matrix_space]; #this might be nothing
  R = dict[:ideal_ring];
  GR = dict[:ground_ring];
  I = load_object(s, MPolyIdeal, R, :defining_ideal)
  Ineqs = load_object(s, Vector{MPolyRingElem}, R, :inequations)
  RMat = isnothing(MS) ? nothing : load_object(s, MatElem, MS, :realization_matrix)
  if R isa MPolyRing
    char = characteristic(coefficient_ring(R))
  else #in this case the realization space is zero dimensional
    char = characteristic(R)
  end
  q = !iszero(char) ? order(coefficient_ring(R)) : nothing
  RS = MatroidRealizationSpace(I, Ineqs, R, RMat, char, q, GR)
  RS.one_realization = load_object(s, Bool, :one_realization)
  
  return RS
end


@register_serialization_type MatroidRealizationSpaceSelfProjecting uses_id
type_params(M::MatroidRealizationSpaceSelfProjecting) = TypeParams(MatroidRealizationSpaceSelfProjecting,
                                                     :matrix_space => begin
                                                      mat = selfprojecting_realization_matrix(M)
                                                      isnothing(mat) ? nothing : parent(mat)
                                                      end,
                                                      :ideal_ring => base_ring(defining_ideal(M)),
                                                     :ground_ring=>M.ground_ring)

function save_object(s::SerializerState, M::MatroidRealizationSpaceSelfProjecting) 
  save_data_dict(s) do
    save_object(s, defining_ideal(M), :defining_ideal)
    save_object(s, inequations(M), :inequations)
    if isnothing(selfprojecting_realization_matrix(M))
      save_object(s, nothing, :selfprojecting_realization_matrix)
    else
      save_object(s, selfprojecting_realization_matrix(M), :selfprojecting_realization_matrix)
    end
#    save_object(s, M.one_realization, :one_realization)
  end
end


function load_object(s::DeserializerState, ::Type{<:MatroidRealizationSpaceSelfProjecting}, dict::Dict)
  MS = dict[:matrix_space]; #this could be nothing
  R = dict[:ideal_ring];
  GR = dict[:ground_ring];
  I = load_object(s, MPolyIdeal, R, :defining_ideal)
  Ineqs = load_object(s, Vector{MPolyRingElem}, R, :inequations)
  RMat = isnothing(MS) ? nothing : load_object(s, MatElem, MS, :selfprojecting_realization_matrix)
  if R isa MPolyRing
    char = characteristic(coefficient_ring(R))
  else #in this case R = QQ
    char = characteristic(R)
  end
  q = !iszero(char) ? order(coefficient_ring(R)) : nothing
 # RS.one_realization = load_object(s, Bool, :one_realization)
  return MatroidRealizationSpaceSelfProjecting(I, Ineqs, R, RMat, char, q, GR)
end


@register_serialization_type MatroidRealizations uses_id
function type_params(M::MatroidRealizations)
  rs = realization_space(M);
  matrix_rs = isnothing(rs) ? nothing : realization_matrix(rs);
  sprs = selfprojecting_realization_space(M);
  matrix_sprs = isnothing(sprs) ? nothing : selfprojecting_realization_matrix(sprs);
  TypeParams(MatroidRealizations, 
             :matrix_space_mrs => isnothing(matrix_rs) ? nothing : parent(realization_matrix(rs)), 
             :matrix_space_sp_mrs => isnothing(matrix_sprs) ? nothing : parent(selfprojecting_realization_matrix(sprs)),
             :ideal_ring_rs => isnothing(rs) ? nothing : base_ring(defining_ideal(rs)),
             :ideal_ring_sprs => isnothing(sprs) ? nothing : base_ring(defining_ideal(sprs)),
             :ground_ring => isnothing(rs) ? nothing : rs.ground_ring,
             :ground_ring_s => isnothing(sprs) ? nothing : sprs.ground_ring,
             )
end

function save_object(s::SerializerState, M::MatroidRealizations) 
  save_data_dict(s) do
    save_object(s, M.name, :name)
    save_object(s, M.matroid, :matroid)
    save_object(s, M.rank, :rank)
    save_object(s, M.length_groundset, :length_groundset)
    save_object(s, M.dim_r, :dim_r)
    save_object(s, M.dim_s, :dim_s)
    save_object(s, M.equality_of_realizationspaces, :equality_of_realizationspaces)
    save_object(s, M.realization_space, :realization_space)
    save_object(s, M.selfprojecting_realization_space, :selfprojecting_realization_space)
  end
end


function load_object(s::DeserializerState, ::Type{<:MatroidRealizations}, dict::Dict)
  dictionary_r = Dict(:ground_ring => dict[:ground_ring], :matrix_space => dict[:matrix_space_mrs], :ideal_ring => dict[:ideal_ring_rs])
  dictionary_s = Dict(:ground_ring => dict[:ground_ring_s], :matrix_space => dict[:matrix_space_sp_mrs], :ideal_ring => dict[:ideal_ring_sprs])
  str = load_object(s, String, :name)
  m = load_object(s, Matroid, :matroid)
  rk = load_object(s, Int, :rank)
  n = load_object(s, Int, :length_groundset)
  dimR = load_object(s, Int, :dim_r)
  dimS = load_object(s, Int, :dim_s)
  boo = load_object(s, Bool, :equality_of_realizationspaces)
  RS = load_object(s, MatroidRealizationSpace, dictionary_r, :realization_space)
  RSSP = load_object(s, MatroidRealizationSpaceSelfProjecting, dictionary_s, :selfprojecting_realization_space)
  return MatroidRealizations(str, m, rk, n, RS, dimR, RSSP, dimS, boo)
end
