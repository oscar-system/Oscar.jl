################################################################################
# non-fmpz variant
function save_internal(s::SerializerState, F::Nemo.GaloisField)
    return Dict(
        :characteristic => UInt64(characteristic(F))
    )
end

function load_internal(s::DeserializerState, ::Type{Nemo.GaloisField}, dict::Dict)
    return Nemo.GaloisField(UInt64(dict[:characteristic]))
end

# elements
function save_internal(s::SerializerState, elem::gfp_elem)
    return Dict(
        :parent => save_type_dispatch(s, parent(elem)),
        :data => Nemo.data(elem)
    )
end

function load_internal(s::DeserializerState, z::Type{gfp_elem}, dict::Dict)
    F = load_type_dispatch(s, Nemo.GaloisField, dict[:parent])
    return F(UInt64(dict[:data]))
end


################################################################################
# fmpz variant
function save_internal(s::SerializerState, F::Nemo.GaloisFmpzField)
    return Dict(
        :characteristic => save_type_dispatch(s, characteristic(F))
    )
end

function load_internal(s::DeserializerState, F::Type{Nemo.GaloisFmpzField}, dict::Dict)
    return F(load_type_dispatch(s, fmpz, dict[:characteristic]))
end

# elements
function save_internal(s::SerializerState, elem::gfp_fmpz_elem)
    return Dict(
        :parent => save_type_dispatch(s, parent(elem)),
        :data => save_type_dispatch(s, Nemo.data(elem))
    )
end

function load_internal(s::DeserializerState, ::Type{gfp_fmpz_elem}, dict::Dict)
    F = load_type_dispatch(s, Nemo.GaloisFmpzField, dict[:parent])
    return F(load_type_dispatch(s, fmpz, dict[:data]))
end

################################################################################
# SimpleNumField
function save_internal(s::SerializerState, K::SimpleNumField)
    return Dict(
        :def_pol => save_type_dispatch(s, defining_polynomial(K)),
        :var => save_type_dispatch(s, var(K)) 
    )
end

function load_internal(s::DeserializerState, ::Type{<: SimpleNumField}, dict::Dict)
    def_pol = load_unknown_type(s, dict[:def_pol])
    var = load_type_dispatch(s, Symbol, dict[:var])
    K, _ = NumberField(def_pol, var, cached=false)

    return K
end

################################################################################
# FqNmodfinitefield
function save_internal(s::SerializerState, K::FqNmodFiniteField)
    return Dict(
        :def_pol => save_type_dispatch(s, defining_polynomial(K))
    )
end

function load_internal(s::DeserializerState, ::Type{FqNmodFiniteField}, dict::Dict)
    def_pol = load_unknown_type(s, dict[:def_pol])
    K, _ = FiniteField(def_pol, cached=false)

    return K
end

#elements
function save_internal(s::SerializerState, k::Union{nf_elem, fq_nmod, Hecke.NfRelElem})
    K = parent(k)
    polynomial = parent(defining_polynomial(K))(k)
    K_dict = save_type_dispatch(s, K)

    return Dict(
        :parent => K_dict,
        :polynomial => save_type_dispatch(s, polynomial)
    )
end

function load_internal(s::DeserializerState,
                       ::Type{<: Union{nf_elem, fq_nmod, Hecke.NfRelElem}},
                       dict::Dict)
    K = load_unknown_type(s, dict[:parent])
    polynomial = load_unknown_type(s, dict[:polynomial])
    return K(polynomial)
end

################################################################################
# Non Simple Extension
function save_internal(s::SerializerState, K::Union{NfAbsNS, NfRelNS})
    def_pols = defining_polynomials(K)
    return Dict(
        :def_pols => save_type_dispatch(s, def_pols),
        :vars => save_type_dispatch(s, vars(K))
    )
end

function load_internal(s::DeserializerState,
                       ::Type{<: Union{NfAbsNS, NfRelNS}},
                       dict::Dict)
    def_pols = load_unknown_type(s, dict[:def_pols])
    vars = load_type_dispatch(s, Vector{Symbol}, dict[:vars])
    K, _ = NumberField(def_pols, vars, cached=false)
    return K
end

#elements
function save_internal(s::SerializerState, k::Union{NfAbsNSElem, Hecke.NfRelNSElem})
    K = parent(k)
    polynomial = Oscar.Hecke.data(k)
    polynomial_parent = parent(polynomial)
    K_dict = save_type_dispatch(s, K)
    return Dict(
        :parent_field => K_dict,
        :polynomial => save_type_dispatch(s, polynomial),
        :polynomial_parent => save_type_dispatch(s, polynomial_parent)
    )
end

function load_internal(s::DeserializerState,
                       ::Type{<: Union{NfAbsNSElem, Hecke.NfRelNSElem}},
                       dict::Dict)
    K = load_unknown_type(s, dict[:parent_field])
    polynomial = load_unknown_type(s, dict[:polynomial])
    polynomial = evaluate(polynomial, gens(K))

    return K(polynomial)
end


################################################################################
# FracField
function save_internal(s::SerializerState, K::FracField)
    return Dict(
        :base_ring => save_type_dispatch(s, base_ring(K)),
    )
end

function load_internal(s::DeserializerState,
                       ::Type{<: FracField},
                       dict::Dict)
    R, _ = load_unknown_type(s, dict[:base_ring])

    return FractionField(R, cached=false)
end

# elements
function save_internal(s::SerializerState, f::FracElem)
    parent_dict = save_type_dispatch(s, parent(f))
    return Dict(
        :parent => parent_dict,
        :den => save_type_dispatch(s, denominator(f)),
        :num => save_type_dispatch(s, numerator(f))
    )
end

function load_internal(s::DeserializerState, ::Type{<: FracElem}, dict::Dict)
    R = load_unknown_type(s, dict[:parent])
    num = load_unknown_type(s, dict[:num])
    den = load_unknown_type(s, dict[:den])

    return R(num) // R(den)
end
