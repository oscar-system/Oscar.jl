################################################################################
# fmpz
function load_internal(s::DeserializerState, ::Type{fmpz}, str::String)
    return fmpz(str)
end

function save_internal(s::SerializerState, z::fmpz)
    return string(z)
end


################################################################################
# Number
function save_internal(s::SerializerState, z::Number)
    return string(z)
end

function load_internal(s::DeserializerState, ::Type{T}, str::String) where {T<:Number}
    return parse(T, str)
end


################################################################################
# Strings
function save_internal(s::SerializerState, str::String)
    return str
end

function load_internal(s::DeserializerState, ::Type{String}, str::String)
    return str
end

################################################################################
# Symbol
function save_internal(s::SerializerState, sym::Symbol)
    return string(sym)
end

function load_internal(s::DeserializerState, ::Type{Symbol}, str::String)
    return Symbol(str)
end
