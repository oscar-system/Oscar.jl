################################################################################
# field of rationals (singleton type)
@registerSerializationType(FlintRationalField)


################################################################################
# non-fmpz variant
@registerSerializationType(Nemo.gfp_elem)
@registerSerializationType(Nemo.GaloisField)

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
@registerSerializationType(Nemo.gfp_fmpz_elem)
@registerSerializationType(Nemo.GaloisFmpzField)

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
@registerSerializationType(AnticNumberField)

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
@registerSerializationType(FqNmodFiniteField)

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
@registerSerializationType(fq_nmod)
@registerSerializationType(nf_elem)

function save_internal(s::SerializerState, k::Union{nf_elem, fq_nmod, Hecke.NfRelElem})
    K = parent(k)
    polynomial = parent(defining_polynomial(K))(k)

    return Dict(
        :parent => save_type_dispatch(s, K),
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
@registerSerializationType(NfAbsNS)

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
    return Dict(
        :parent_field => save_type_dispatch(s, K),
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
    return Dict(
        :parent => save_type_dispatch(s, parent(f)),
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

################################################################################
# RealField
function save_internal(s::SerializerState, RR::Nemo.RealField)
    return Dict(
        :precision => save_type_dispatch(s, precision(RR))
    )    
end

function load_internal(s::DeserializerState, ::Type{Nemo.RealField}, dict::Dict)
    prec = load_unknown_type(s, dict[:precision])
    return Nemo.RealField(prec)
end

# elements
function save_internal(s::SerializerState, r::arb)
    c_str = ccall((:arb_dump_str, Nemo.Arb_jll.libarb), Ptr{UInt8}, (Ref{arb},), r)
    arb_unsafe_str = unsafe_string(c_str)

    # free memory
    ccall((:flint_free, Nemo.FLINT_jll.libflint), Nothing, (Ptr{UInt8},), c_str)

    return Dict(
        :parent => save_type_dispatch(s, parent(r)),
        :arb_unsafe_str => save_type_dispatch(s, arb_unsafe_str)
    )    
end


function load_internal(s::DeserializerState, ::Type{arb}, dict::Dict)
    parent = load_type_dispatch(s, Nemo.RealField, dict[:parent])
    arb_unsafe_str = load_type_dispatch(s, String, dict[:arb_unsafe_str])
    r = Nemo.arb()
    ccall((:arb_load_str, Nemo.Arb_jll.libarb),
          Int32, (Ref{arb}, Ptr{UInt8}), r, arb_unsafe_str)
    r.parent = parent
    
    return r
end

################################################################################
# ComplexField

@registerSerializationType(AcbField)

function save_internal(s::SerializerState, CC::AcbField)
    return Dict(
        :precision => save_type_dispatch(s, precision(CC))
    )    
end

function load_internal(s::DeserializerState, ::Type{AcbField}, dict::Dict)
    prec = load_unknown_type(s, dict[:precision])
    return ComplexField(prec)
end

function save_internal(s::SerializerState, c::acb)
    return Dict(
        :parent => save_type_dispatch(s, parent(c)),
        :real => save_type_dispatch(s, real(c)),
        :imag => save_type_dispatch(s, imag(c))
    )
end

function load_internal(s::DeserializerState, ::Type{acb}, dict::Dict)
    CC = load_type_dispatch(s, AcbField, dict[:parent])
    real_part = load_type_dispatch(s, arb, dict[:real])
    imag_part = load_type_dispatch(s, arb, dict[:imag])

    return CC(real_part, imag_part)
end

################################################################################
# Field Embeddings

@registerSerializationType(Hecke.NumFieldEmbNfAbs)

function save_internal(s::SerializerState, E::Hecke.NumFieldEmbNfAbs)
    K = number_field(E)
    g = gen(K)
    g_ball = E(g)
    
    return Dict(
        :num_field => save_type_dispatch(s, K),
        :gen_ball => save_type_dispatch(s, g_ball)
    )
end

function load_internal(s::DeserializerState, ::Type{Hecke.NumFieldEmbNfAbs}, dict::Dict)
    K = load_type_dispatch(s, AnticNumberField, dict[:num_field])
    gen_ball = load_type_dispatch(s, acb, dict[:gen_ball])

    return complex_embedding(K, gen_ball)
end
