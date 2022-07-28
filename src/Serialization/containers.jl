################################################################################
# Saving and loading vectors
function save_internal(s::SerializerState, vec::Vector{T}) where T
    return Dict(
        :vector => [save_type_dispatch(s, x) for x in vec]
    )
end

function load_internal(s::DeserializerState, ::Type{Vector{T}}, dict::Dict) where T
    return Vector{T}([load_unknown_type(s, x) for x in dict[:vector]])
end



################################################################################
# Saving and loading tuples
function save_internal(s::SerializerState, tup::Tuple)
    return Dict(
        :content => [save_type_dispatch(s, x) for x in tup]
    )
end

function load_internal(s::DeserializerState, T::Type{<:Tuple}, dict::Dict)
    n = fieldcount(T)
    content = dict[:content]
    @assert length(content) == n "Wrong length of tuple, data may be corrupted."
    return T(load_type_dispatch(s, fieldtype(T, i), content[i]) for i in 1:n)
end



################################################################################
# Saving and loading matrices
function save_internal(s::SerializerState, mat::Matrix{T}) where T
    m, n = size(mat)
    return Dict(
        :matrix => [save_type_dispatch(s, [mat[i, j] for j in 1:n]) for i in 1:m]
    )
end

function load_internal(s::DeserializerState, ::Type{Matrix{T}}, dict::Dict) where T
    x = dict[:matrix]
    y = reduce(vcat, [permutedims(load_type_dispatch(s, Vector{T}, x[i])) for i in 1:length(x)])
    return Matrix{T}(y)
end


