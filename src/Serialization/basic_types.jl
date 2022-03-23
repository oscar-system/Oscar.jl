################################################################################
# fmpz
function load_intern(s::DeserializerState, ::Type{fmpz}, str::String)
    return fmpz(str)
end

function save_intern(s::SerializerState, z::fmpz)
    return string(z)
end


################################################################################
# Number
function save_intern(s::SerializerState, z::Number)
    return string(z)
end

function load_intern(s::DeserializerState, ::Type{T}, str::String) where {T<:Number}
    return parse(T, str)
end

