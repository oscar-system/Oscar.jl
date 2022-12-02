################################################################################
# Saving and loading vectors

@registerSerializationType(Vector)

function save_internal(s::SerializerState, vec::Vector)
    return Dict(
        :vector => [save_type_dispatch(s, x) for x in vec]
    )
end

function save_internal(s::SerializerState, vec::Vector{T}) where T <: OscarBasicType
    return Dict(
        :vector => [save_type_dispatch(s, x) for x in vec],
        :entry_type => encodeType(T)
    )
end

# deserialize with specific content type
function load_internal(s::DeserializerState, ::Type{Vector{T}}, dict::Dict) where T
    if isconcretetype(T)
      return Vector{T}([load_type_dispatch(s, T, x) for x in dict[:vector]])
    end
    return Vector{T}([load_unknown_type(s, x) for x in dict[:vector]])
end

function load_internal_with_parent(s::DeserializerState,
                                   ::Type{Vector{T}},
                                   dict::Dict,
                                   parent) where T
    if isconcretetype(T)
        return Vector{T}([load_type_dispatch(s, T, x; parent=parent) for x in dict[:vector]])
    end
    return Vector{T}([load_unknown_type(s, x; parent=parent) for x in dict[:vector]])
end

# deserialize without specific content type
function load_internal(s::DeserializerState, ::Type{Vector}, dict::Dict)
    if haskey(dict, :entry_type)
        T = decodeType(dict[:entry_type])
        
        return [load_type_dispatch(s, T, x) for x in dict[:vector]]
    end
    return [load_unknown_type(s, x) for x in dict[:vector]]
end

function load_internal_with_parent(s::DeserializerState,
                                   ::Type{Vector},
                                   dict::Dict,
                                   parent)
    if haskey(dict, :entry_type)
        T = decodeType(dict[:entry_type])
        
        return [load_type_dispatch(s, T, x) for x in dict[:vector]]
    end
    return [load_unknown_type(s, x; parent=parent) for x in dict[:vector]]
end

################################################################################
# Saving and loading Tuple
@registerSerializationType(Tuple)

function save_internal(s::SerializerState, tup::T) where T <: Tuple
    n = fieldcount(T)
    return Dict(
        :field_types => [encodeType(fieldtype(T, i)) for i in 1:n],
        :content => [save_type_dispatch(s, x) for x in tup]
    )
end

function load_internal(s::DeserializerState, T::Type{<:Tuple}, dict::Dict)
    field_types = map(decodeType, dict[:field_types])
    n = length(field_types)
    content = dict[:content]
    @assert length(content) == n  "Wrong length of tuple, data may be corrupted."
    return T(load_type_dispatch(s, field_types[i], content[i]) for i in 1:n)
end

################################################################################
# Saving and loading NamedTuple
@registerSerializationType(NamedTuple)

function save_internal(s::SerializerState, n_tup:: NamedTuple)
    return Dict(
        :keys => save_type_dispatch(s, keys(n_tup)),
        :content => save_type_dispatch(s, values(n_tup))
    )
end

function load_internal(s::DeserializerState, ::Type{<:NamedTuple}, dict::Dict)
    tup = load_unknown_type(s, dict[:content])
    keys = load_type_dispatch(s, Tuple, dict[:keys])
    return NamedTuple{keys}(tup)
end


################################################################################
# Saving and loading matrices

@registerSerializationType(Matrix)

function save_internal(s::SerializerState, mat::Matrix{T}) where T
    m, n = size(mat)
    return Dict(
        :matrix => [save_type_dispatch(s, [mat[i, j] for j in 1:n]) for i in 1:m]
    )
end

# deserialize with specific content type
function load_internal(s::DeserializerState, ::Type{Matrix{T}}, dict::Dict) where T
    x = dict[:matrix]
    y = reduce(vcat, [permutedims(load_type_dispatch(s, Vector{T}, x[i])) for i in 1:length(x)])
    return Matrix{T}(y)
end

# deserialize without specific content type
function load_internal(s::DeserializerState, ::Type{Matrix}, dict::Dict)
    x = dict[:matrix]
    y = reduce(vcat, [permutedims(load_type_dispatch(s, Vector, x[i])) for i in 1:length(x)])
    return Matrix(y)
end

# deserialize without specific content type
function load_internal_with_parent(s::DeserializerState,
                                   ::Type{Matrix}, dict::Dict, parent)
    x = dict[:matrix]
    y = reduce(vcat, [
        permutedims(load_type_dispatch(s, Vector, x[i], parent=parent)) for i in 1:length(x)
            ])
    return Matrix(y)
end
