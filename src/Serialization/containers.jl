################################################################################
# Saving and loading vectors
function save_intern(s::SerializerState, vec::Vector{T}) where T
    return Dict(
        :vector => [save(s, x) for x in vec]
    )
end

function load_intern(s::DeserializerState, ::Type{Vector{T}}, dict::Dict) where T
    return Vector{T}([load(s, x; check_namespace=false) for x in dict[:vector]])
end



################################################################################
# Saving and loading tuples
function save_intern(s::SerializerState, tup::NTuple{n, T}) where {n, T}
    return Dict(
        :content => [save(s, x) for x in tup]
    )
end

function load_intern(s::DeserializerState, ::Type{NTuple{n, T}}, dict::Dict) where {n, T}
    content = [load(s, T, x) for x in dict[:content]]
    @assert length(content) == n "Wrong length of tuple, data may be corrupted."
    return NTuple{n, T}(content)
end


