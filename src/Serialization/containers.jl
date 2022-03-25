################################################################################
# Saving and loading vectors
function save_internal(s::SerializerState, vec::Vector{T}) where T
    return Dict(
        :vector => [save_type_dispatch(s, x) for x in vec]
    )
end

function load_internal(s::DeserializerState, ::Type{Vector{T}}, dict::Dict) where T
    return Vector{T}([load_type_dispatch(s, x; check_namespace=false) for x in dict[:vector]])
end



################################################################################
# Saving and loading tuples
function save_internal(s::SerializerState, tup::NTuple{n, T}) where {n, T}
    return Dict(
        :content => [save_type_dispatch(s, x) for x in tup]
    )
end

function load_internal(s::DeserializerState, ::Type{NTuple{n, T}}, dict::Dict) where {n, T}
    content = [load_type_dispatch(s, T, x) for x in dict[:content]]
    @assert length(content) == n "Wrong length of tuple, data may be corrupted."
    return NTuple{n, T}(content)
end


