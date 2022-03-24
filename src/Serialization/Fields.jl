################################################################################
# non-fmpz variant
function save_intern(s::SerializerState, F::Nemo.GaloisField)
    return Dict(
        :characteristic => UInt64(characteristic(F))
    )
end

function load_intern(s::DeserializerState, ::Type{Nemo.GaloisField}, dict::Dict)
    return Nemo.GaloisField(UInt64(dict[:characteristic]))
end

# elements
function save_intern(s::SerializerState, elem::gfp_elem)
    return Dict(
        :parent => save(s, parent(elem)),
        :data => Nemo.data(elem)
    )
end

function load_intern(s::DeserializerState, z::Type{gfp_elem}, dict::Dict)
    F = load(s, Nemo.GaloisField, dict[:parent])
    return F(UInt64(dict[:data]))
end


################################################################################
# fmpz variant
function save_intern(s::SerializerState, F::Nemo.GaloisFmpzField)
    return Dict(
        :characteristic => save(s, characteristic(F))
    )
end

function load_intern(s::DeserializerState, F::Type{Nemo.GaloisFmpzField}, dict::Dict)
    return F(load(s, fmpz, dict[:characteristic]))
end

# elements
function save_intern(s::SerializerState, elem::gfp_fmpz_elem)
    return Dict(
        :parent => save(s, parent(elem)),
        :data => save(s, Nemo.data(elem))
    )
end

function load_intern(s::DeserializerState, ::Type{gfp_fmpz_elem}, dict::Dict)
    F = load(s, Nemo.GaloisFmpzField, dict[:parent])
    return F(load(s, fmpz, dict[:data]))
end
