################################################################################
# fmpz
function load_intern(s::DeserializerState, z::Type{fmpz}, dict::Dict)
    return fmpz(dict[:val])
end

function save_intern(s::SerializerState, z::fmpz)
    return Dict(
        :val => string(z)
    )
end



