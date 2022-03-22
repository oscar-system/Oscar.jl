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


function encodeType(::Type{Vector{T}}) where T
    return Dict(
        :container => "Vector",
        :element => encodeType(T),
    )
end

function decodeType(input::Dict{Symbol, Any})
    container = input[:container]
    if container == "Vector"
        return Vector{decodeType(input[:element])}
    else
        error("Unknown container type \"$(container)\"")
    end
end
