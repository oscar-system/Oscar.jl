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
        :var => save_type_dispatch(s, String(var(K)))
    )
end
 
function load_internal(s::DeserializerState, ::Type{<: SimpleNumField}, dict::Dict)
    def_pol = load_type_dispatch(s, dict[:def_pol], check_namespace=false)
    var = load_type_dispatch(s, dict[:var], check_namespace=false)
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
    def_pol = load_type_dispatch(s, dict[:def_pol], check_namespace=false)
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
    K = load_type_dispatch(s, dict[:parent], check_namespace=false)
    polynomial = load_type_dispatch(s, dict[:polynomial], check_namespace=false)
    return K(polynomial)
end

################################################################################
# Non Simple Extension
function save_internal(s::SerializerState, K::Union{NfAbsNS, NfRelNS})
    def_pols = defining_polynomials(K)
    return Dict(
        :def_pols => save_type_dispatch(s, def_pols),
    )
end

function load_internal(s::DeserializerState,
                       ::Type{<: Union{NfAbsNS, NfRelNS}},
                       dict::Dict)
    def_pols = load_type_dispatch(s, dict[:def_pols], check_namespace=false)
    K, _ = NumberField(def_pols, cached=false)
    return K
end

#elements
function save_internal(s::SerializerState, k::Union{NfAbsNSElem, Hecke.NfRelNSElem})
  K = parent(k)
  polynomial = Oscar.Hecke.data(k)
  K_dict = save_type_dispatch(s, K)

  return Dict(:parent => K_dict, :polynomial => save_type_dispatch(s, polynomial))
end

function load_internal(s::DeserializerState,
                       ::Type{<: Union{NfAbsNSElem, Hecke.NfRelNSElem}},
                       dict::Dict)
  K = load_type_dispatch(s, dict[:parent]; check_namespace=false)
  polynomial = load_type_dispatch(s, dict[:polynomial]; check_namespace=false)
  polynomial = evaluate(polynomial, gens(parent(K.pol[1])))

  return K(polynomial)
end
