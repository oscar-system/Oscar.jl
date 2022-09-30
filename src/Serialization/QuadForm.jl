############################################################
# QuadSpace
@registerSerializationType(Hecke.QuadSpace)

@registerSerializationType(ZLat)

function save_internal(s::SerializerState, V::Hecke.QuadSpace)
    return Dict(
        :base_ring => save_type_dispatch(s, base_ring(V)),
        :gram_matrix => save_type_dispatch(s, gram_matrix(V))
    )
end

function load_internal(s::DeserializerState, ::Type{<: Hecke.QuadSpace}, dict::Dict)
    F = load_unknown_type(s, dict[:base_ring])
    gram = load_type_dispatch(s, MatElem, dict[:gram_matrix])
    @assert base_ring(gram)===F
    return quadratic_space(F, gram)
end

# We should move this somewhere else at some point, maybe when there is a section
# on modules
function save_internal(s::SerializerState, L::ZLat)
    return Dict(
        :basis => save_type_dispatch(s, basis_matrix(L)),
        :ambient_space => save_type_dispatch(s, ambient_space(L))
    )
end

function load_internal(s::DeserializerState, ::Type{ZLat}, dict::Dict)
    B = load_type_dispatch(s, fmpq_mat, dict[:basis])
    V = load_type_dispatch(s, Hecke.QuadSpace{FlintRationalField, fmpq_mat}, dict[:ambient_space])
    return lattice(V, B)
end
