################################################################################
# ring of integers (singleton type)
@registerSerializationType(FlintIntegerRing)


################################################################################
#  non simpleton base rings
@registerSerializationType(Nemo.NmodRing, "Nemo.NmodRing")

function save_internal(s::SerializerState, R::Nemo.NmodRing)
    return Dict(
        :modulus => save_type_dispatch(s, modulus(R))
    )
end

function load_internal(s::DeserializerState, ::Type{Nemo.NmodRing}, dict::Dict)
    modulus = load_type_dispatch(s, UInt64, dict[:modulus])
    return Nemo.NmodRing(modulus)
end

#elements
@registerSerializationType(nmod)

function save_internal(s::SerializerState, r::nmod)
    return Dict(
        :parent => save_type_dispatch(s, parent(r)),
        :class_val => save_type_dispatch(s, fmpz(r))
    )
end

function load_internal(s::DeserializerState, ::Type{nmod}, dict::Dict)
    parent_ring = load_type_dispatch(s, Nemo.NmodRing, dict[:parent])
    class_val = load_type_dispatch(s, fmpz, dict[:class_val])
    return parent_ring(class_val)
end

function load_internal_with_parent(s::DeserializerState,
                                   ::Type{nmod},
                                   dict::Dict,
                                   parent_ring::Nemo.NmodRing)
    class_val = load_type_dispatch(s, fmpz, dict[:class_val])
    return parent_ring(class_val)
end


################################################################################
#  Polynomial Rings

encodeType(::Type{<:PolyRing}) = "PolyRing"
reverseTypeMap["PolyRing"] = PolyRing

encodeType(::Type{<:MPolyRing}) = "MPolyRing"
reverseTypeMap["MPolyRing"] = MPolyRing

@registerSerializationType(FmpqMPolyRing)
@registerSerializationType(FmpqPolyRing)
@registerSerializationType(FmpzMPolyRing)
@registerSerializationType(FmpzPolyRing)
@registerSerializationType(FqNmodMPolyRing)
@registerSerializationType(FqNmodPolyRing)
@registerSerializationType(GFPPolyRing)
@registerSerializationType(NmodMPolyRing)
@registerSerializationType(NmodPolyRing)

function save_internal(s::SerializerState, R::Union{MPolyRing, PolyRing})
    return Dict(
        :symbols => save_type_dispatch(s, symbols(R)),
        :base_ring => save_type_dispatch(s, base_ring(R)),
    )
end

function load_internal(s::DeserializerState,
                       T::Type{<: Union{MPolyRing, PolyRing}},
                       dict::Dict)
    base_ring = load_unknown_type(s, dict[:base_ring])
    symbols = load_type_dispatch(s, Vector{Symbol}, dict[:symbols]) 

    if T <: PolyRing
        return PolynomialRing(base_ring, symbols..., cached=false)
    end

    return PolynomialRing(base_ring, symbols, cached=false)
end

################################################################################
# Multivariate Polynomials
@registerSerializationType(fmpq_mpoly)
@registerSerializationType(fmpz_mpoly)
@registerSerializationType(fq_nmod_mpoly)
@registerSerializationType(nmod_mpoly)

encodeType(::Type{<:MPolyElem}) = "MPolyElem"
reverseTypeMap["MPolyElem"] = MPolyElem

function save_internal(s::SerializerState, p::MPolyElem)
    parent_ring = parent(p)
    parent_ring = save_type_dispatch(s, parent_ring)
    terms = []

    for i in 1:length(p)
        term = Dict(
            :exponent => save_type_dispatch(s, exponent_vector(p, i)),
            :coeff => save_type_dispatch(s, coeff(p, i))
        )
        push!(terms, term)
    end

    return Dict(
        :terms => terms,
        :parent => parent_ring,
    )
end

function load_internal(s::DeserializerState, ::Type{<: MPolyElem}, dict::Dict)
    R, symbols = load_unknown_type(s, dict[:parent])
    coeff_ring = coefficient_ring(R)
    coeff_type = elem_type(coeff_ring)
    polynomial = MPolyBuildCtx(R)

    for term in dict[:terms]
        c = load_type_dispatch(s, coeff_type, term[:coeff])
        e = load_type_dispatch(s, Vector{Int}, term[:exponent])
        push_term!(polynomial, c, e)
    end
    return finish(polynomial)
end

function load_internal_with_parent(s::DeserializerState,
                                   ::Type{<: MPolyElem},
                                   dict::Dict,
                                   parent_ring::MPolyRing)
    # load parent in case serialized parent needs to be checked against given parent
    _, _ = load_unknown_type(s, dict[:parent])

    coeff_ring = coefficient_ring(parent_ring)
    coeff_type = elem_type(coeff_ring)
    polynomial = MPolyBuildCtx(parent_ring)
    
    for term in dict[:terms]
        c = load_type_dispatch(s, coeff_type, term[:coeff]; parent=coeff_ring)
        e = load_type_dispatch(s, Vector{Int}, term[:exponent])
        push_term!(polynomial, c, e)
    end
    return finish(polynomial)
end

################################################################################
# Univariate Polynomials

encodeType(::Type{<:PolyElem}) = "PolyElem"
reverseTypeMap["PolyElem"] = PolyElem

@registerSerializationType(fmpq_poly)
@registerSerializationType(fmpz_poly)
@registerSerializationType(fq_nmod_poly)
@registerSerializationType(gfp_poly)
@registerSerializationType(nmod_poly)

function save_internal(s::SerializerState, p::PolyElem)
    parent_ring = parent(p)
    parent_ring = save_type_dispatch(s, parent_ring)

    return Dict(
        :parent => parent_ring,
        :coeffs => save_type_dispatch(s, collect(coefficients(p)))
    )
end

function load_internal(s::DeserializerState, ::Type{<: PolyElem}, dict::Dict)
    R, y = load_unknown_type(s, dict[:parent])
    coeff_ring = coefficient_ring(R)
    coeff_type = elem_type(coeff_ring)
    coeffs = load_type_dispatch(s, Vector{coeff_type}, dict[:coeffs])

    return R(coeffs)
end

function load_internal_with_parent(s::DeserializerState,
                                   ::Type{<: PolyElem},
                                   dict::Dict,
                                   parent_ring::PolyRing)
    # cache parent inside serializer state in case parent needs
    # to be checked against the passed parent
    _, _ = load_unknown_type(s, dict[:parent])
    
    coeff_ring = coefficient_ring(parent_ring)
    coeff_type = elem_type(coeff_ring)
    coeffs = load_type_dispatch(s, Vector{coeff_type}, dict[:coeffs]; parent=coeff_ring)

    return parent_ring(coeffs)
end

################################################################################
# Polynomial Ideals

encodeType(::Type{<:MPolyIdeal}) = "MPolyIdeal"
reverseTypeMap["MPolyIdeal"] = MPolyIdeal

function save_internal(s::SerializerState, i::MPolyIdeal)
    generators = gens(i)
    parent_ring = save_type_dispatch(s, parent(generators[1]))

    return Dict(
        :parent => parent_ring,
        :gens => save_type_dispatch(s, generators),
    )
end

function load_internal(s::DeserializerState, ::Type{<: MPolyIdeal}, dict::Dict)
    parent_ring, _ = load_unknown_type(s, dict[:parent])
    gens = load_type_dispatch(s, Vector{elem_type(parent_ring)}, dict[:gens])

    return ideal(parent_ring, gens)
end


function load_internal_with_parent(s::DeserializerState,
                                   ::Type{<: MPolyIdeal},
                                   dict::Dict,
                                   parent_ring::MPolyRing)
    gens = load_type_dispatch(s, Vector{elem_type(parent_ring)}, dict[:gens])
    return ideal(parent_ring, gens)
end

################################################################################
# Matrices

encodeType(::Type{<:MatElem}) = "MatElem"
reverseTypeMap["MatElem"] = MatElem

@registerSerializationType(fmpz_mat)
@registerSerializationType(fmpq_mat)
@registerSerializationType(fq_nmod_mat)
@registerSerializationType(nmod)
@registerSerializationType(nmod_mat)

function save_internal(s::SerializerState, m::MatrixElem)
    return Dict(
        :matrix => save_type_dispatch(s, Array(m)),
    )
end

function load_internal(s::DeserializerState,
                       ::Type{<: MatElem},
                       dict::Dict) 
    mat = load_type_dispatch(s, Matrix, dict[:matrix])
    entries_ring = parent(mat[1])
    return matrix(entries_ring, mat)
end
