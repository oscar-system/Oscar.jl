################################################################################
# field of rationals (singleton type)
@registerSerializationType(QQField)

################################################################################
# non-ZZRingElem variant
@registerSerializationType(Nemo.fpFieldElem)
@registerSerializationType(Nemo.fpField)

function save_internal(s::SerializerState, F::Nemo.fpField)
    return Dict(
        :characteristic => UInt64(characteristic(F))
    )
end

function load_internal(s::DeserializerState, ::Type{Nemo.fpField}, dict::Dict)
    return Nemo.fpField(UInt64(dict[:characteristic]))
end

# elements
function save_internal(s::SerializerState, elem::fpFieldElem)
    return string(elem)
end

function load_internal(s::DeserializerState, z::Type{fpFieldElem}, dict::Dict)
    F = load_type_dispatch(s, Nemo.fpField, dict[:parent])
    return F(UInt64(dict[:data]))
end

function load_internal_with_parent(s::DeserializerState,
                                   z::Type{fpFieldElem},
                                   str::String,
                                   parent::Nemo.fpField)
    return parent(parse(UInt64, str))
end



################################################################################
# ZZRingElem variant
@registerSerializationType(Nemo.FpFieldElem)
@registerSerializationType(Nemo.FpField)

function save_internal(s::SerializerState, F::Nemo.FpField)
    return Dict(
        :characteristic => save_type_dispatch(s, characteristic(F))
    )
end

function load_internal(s::DeserializerState, F::Type{Nemo.FpField}, dict::Dict)
    return F(load_type_dispatch(s, ZZRingElem, dict[:characteristic]))
end

# elements
function save_internal(s::SerializerState, elem::FpFieldElem)
    return string(Nemo.data(elem))
end

function load_internal_with_parent(s::DeserializerState,
                                   ::Type{FpFieldElem},
                                   str::String,
                                   parent::Nemo.FpField)
    return parent(load_type_dispatch(s, ZZRingElem, str))
end

################################################################################
# SimpleNumField

@registerSerializationType(Hecke.NfRel, true)
@registerSerializationType(AnticNumberField, true)

function save_internal(s::SerializerState, K::SimpleNumField)
    return Dict(
        :def_pol => save_type_dispatch(s, defining_polynomial(K)),
        :var => save_type_dispatch(s, var(K)) 
    )
end

function load_internal(s::DeserializerState, ::Type{<: SimpleNumField}, dict::Dict)
    def_pol = load_unknown_type(s, dict[:def_pol])
    var = load_type_dispatch(s, Symbol, dict[:var])
    K, _ = number_field(def_pol, var, cached=false)
    return K
end

################################################################################
# FqNmodfinitefield
@registerSerializationType(fqPolyRepField, true)

function save_internal(s::SerializerState, K::fqPolyRepField)
    return Dict(
        :def_pol => save_type_dispatch(s, defining_polynomial(K))
    )
end

function load_internal(s::DeserializerState,
                       ::Type{<: fqPolyRepField},
                       dict::Dict)
    def_pol = load_unknown_type(s, dict[:def_pol])
    K, _ = FiniteField(def_pol, cached=false)
    return K
end

#elements
@registerSerializationType(fqPolyRepFieldElem)
@registerSerializationType(nf_elem)
@registerSerializationType(Hecke.NfRelElem)

function save_internal(s::SerializerState, k::Union{nf_elem, fqPolyRepFieldElem, Hecke.NfRelElem})
    K = parent(k)
    polynomial = parent(defining_polynomial(K))(k)
    poly_dict = save_internal(s, polynomial)
    parents = poly_dict[:parents]
    push!(parents, save_as_ref(s, K))
    return Dict(
        :parents => parents,
        :terms => poly_dict[:terms]
    )
end

function load_terms(s::DeserializerState, parents::Vector, terms::Vector,
                    parent_ring::Union{fqPolyRepField, SimpleNumField})
    loaded_terms = load_terms(s, parents[1:end - 1], terms, parents[end - 1])
    return parent_ring(loaded_terms)
end

function load_internal(s::DeserializerState,
                       ::Type{<: Union{nf_elem, fqPolyRepFieldElem, Hecke.NfRelElem}},
                       dict::Dict)
    parents = load_parents(s, dict[:parents])
    return load_terms(s, parents[1:end], dict[:terms], parents[end])
end

function load_internal_with_parent(s::DeserializerState,
                                   ::Type{<: Union{nf_elem, fqPolyRepFieldElem, Hecke.NfRelElem}},
                                   dict::Dict,
                                   parent_field::Union{fqPolyRepField, SimpleNumField})
    parents = get_parents(parent(defining_polynomial(parent_field)))
    terms = load_terms(s, parents, dict[:terms], parents[end])
    return parent_field(terms)
end

################################################################################
# FqField
@registerSerializationType(FqField, true)
@registerSerializationType(FqFieldElem)

function save_internal(s::SerializerState, K::FqField)
    if absolute_degree(K) == 1
        return Dict(
            :order => save_type_dispatch(s, order(K))
        )
    end
    return Dict(
        :def_pol => save_type_dispatch(s, defining_polynomial(K))
    )
end

function load_internal(s::DeserializerState,
                       ::Type{<: FqField},
                       dict::Dict)
    if haskey(dict, :order)
        order = load_type_dispatch(s, ZZRingElem, dict[:order])
        return Hecke.Nemo._FiniteField(order)[1]
    end
    def_pol = load_unknown_type(s, dict[:def_pol])
    return Hecke.Nemo._FiniteField(def_pol, cached=false)[1]
end

# elements
function save_internal(s::SerializerState, k::FqFieldElem)
    K = parent(k)
    if absolute_degree(K) == 1
        return string(lift(ZZ, k))
    end
    polynomial = parent(defining_polynomial(K))(lift(ZZ["x"][1], k))
    encoded_polynomial = save_internal(s, polynomial)
    parents = encoded_polynomial[:parents]
    push!(parents, save_as_ref(s, K))
    
    return Dict(
        :parents => parents,
        :terms => encoded_polynomial[:terms]
    )
end

# Field should already be loaded by this point
function load_internal(s::DeserializerState,
                       ::Type{<: FqFieldElem},
                       dict::Dict)
    loaded_parents = load_parents(s, dict[:parents])
    return load_terms(s, loaded_parents, dict[:terms], loaded_parents[end])
end

function load_terms(s::DeserializerState, parents::Vector, terms::Vector,
                    parent_ring::FqField)
    loaded_terms = load_terms(s, parents[1:end - 1], terms, parents[end - 1])
    try
        return parent_ring(loaded_terms)
    catch err
        # hack untill we get updates in nemo
        if err isa ErrorException
            if err.msg == "Polynomial has wrong coefficient ring" && absolute_degree(coefficient_ring(parent(loaded_terms))) == 1
                return parent_ring.forwardmap(loaded_terms)
            end
        end
        throw(err)
    end
end

function load_internal_with_parent(s::DeserializerState,
                                   ::Type{<: FqFieldElem},
                                   str::String,
                                   parent_field::FqField)
    @assert absolute_degree(parent_field) == 1
    return parent_field(ZZ(str))
end


################################################################################
# Non Simple Extension

@registerSerializationType(Hecke.NfRelNS, true)

@registerSerializationType(NfAbsNS, true)
@registerSerializationType(NfAbsNSElem)

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
    K, _ = number_field(def_pols, vars, cached=false)
    return K
end

#elements
@registerSerializationType(Hecke.NfRelNSElem)

function save_internal(s::SerializerState, k::Union{NfAbsNSElem, Hecke.NfRelNSElem})
    K = parent(k)
    polynomial = Oscar.Hecke.data(k)
    polynomial_parent = parent(polynomial)
    poly_dict = save_internal(s, polynomial)
    parents = poly_dict[:parents]
    push!(parents, save_as_ref(s, K))
    return Dict(
        :parents => parents,
        :terms => poly_dict[:terms]
    )
end

function load_terms(s::DeserializerState, parents::Vector, terms::Vector,
                    parent_ring::Union{NfAbsNS, NfRelNS})
    loaded_terms = load_terms(s, parents[1:end - 1], terms, parents[end - 1])
    loaded_terms = evaluate(loaded_terms, gens(parent_ring))
    return parent_ring(loaded_terms)
end

function load_internal(s::DeserializerState,
                       ::Type{<: Union{NfAbsNSElem, Hecke.NfRelNSElem}},
                       dict::Dict)
    K = load_unknown_type(s, dict[:parent_field])
    polynomial = load_unknown_type(s, dict[:polynomial])
    polynomial = evaluate(polynomial, gens(K))

    return K(polynomial)
end

function load_internal_with_parent(s::DeserializerState,
                                   ::Type{<: NfAbsNSElem},
                                   dict::Dict,
                                   parent_field::NfAbsNS)
    polynomial = load_unknown_type(s, dict[:polynomial])
    polynomial = evaluate(polynomial, gens(parent_field))

    return parent_field(polynomial)
end

function load_internal_with_parent(s::DeserializerState,
                                   ::Type{<: Hecke.NfRelNSElem},
                                   dict::Dict,
                                   parent_field::Hecke.NfRelNS)
    ngens = length(gens(parent_field))
    parent_polynomial_ring, _ = polynomial_ring(base_field(parent_field), ngens)
    polynomial = load_unknown_type(s, dict[:polynomial]; parent=parent_polynomial_ring)
    polynomial = evaluate(polynomial, gens(parent_field))

    return parent_field(polynomial)
end


################################################################################
# FracField

@registerSerializationType(FracField, true)

function save_internal(s::SerializerState, K::FracField)
    return Dict(
        :base_ring => save_type_dispatch(s, base_ring(K)),
    )
end

function load_internal(s::DeserializerState,
                       ::Type{<: FracField},
                       dict::Dict)
    R = load_unknown_type(s, dict[:base_ring])

    return fraction_field(R, cached=false)
end

# elements
@registerSerializationType(FracElem)

function save_internal(s::SerializerState, f::FracElem)
    encoded_denominator = save_internal(s, denominator(f))
    encoded_numerator = save_internal(s, numerator(f))
    parents = encoded_numerator[:parents]
    push!(parents, save_as_ref(s, parent(f)))

    return Dict(
        :parents => parents,
        :terms => (encoded_numerator[:terms], encoded_denominator[:terms]),
    )
end

function load_terms(s::DeserializerState, parents::Vector, terms::Vector,
                    parent_ring::FracField)
    num_coeff, den_coeff = terms
    loaded_num = load_terms(s, parents[1:end - 1], num_coeff, parents[end - 1])
    loaded_den = load_terms(s, parents[1:end - 1], den_coeff, parents[end-1])
    return  parent_ring(loaded_num) // parent_ring(loaded_den)
end

function load_internal(s::DeserializerState,
                       ::Type{<: FracElem},
                       dict::Dict)
    parents = load_parents(s, dict[:parents])
    return load_terms(s, parents, dict[:terms], parents[end])
end

function load_internal_with_parent(s::DeserializerState,
                                   ::Type{<: FracElem},
                                   dict::Dict,
                                   parent:: FracField)
    parents = get_parents(parent)
    return load_terms(s, parents, dict[:terms], parents[end])
end

################################################################################
# RationalFunctionField

@registerSerializationType(AbstractAlgebra.Generic.RationalFunctionField,
                           true,
                           "RationalFunctionField")

function save_internal(s::SerializerState,
                       RF::AbstractAlgebra.Generic.RationalFunctionField)
    return Dict(
        :base_ring => save_type_dispatch(s, base_ring(RF)),
        :symbols => save_type_dispatch(s, symbols(RF))
    )
end

function load_internal(s::DeserializerState,
                       ::Type{<: AbstractAlgebra.Generic.RationalFunctionField},
                       dict::Dict)
    R = load_unknown_type(s, dict[:base_ring])
    
    if dict[:symbols] isa String
        symbol = load_type_dispatch(s, Symbol, dict[:symbols])
        return RationalFunctionField(R, symbol, cached=false)[1]
    end

    symbols = load_type_dispatch(s, Vector{Symbol}, dict[:symbols])
    return RationalFunctionField(R, symbols, cached=false)[1]
end

#elements
@registerSerializationType(AbstractAlgebra.Generic.RationalFunctionFieldElem,
                           true,
                           "RationalFunctionFieldElem")

function save_internal(s::SerializerState, f::AbstractAlgebra.Generic.RationalFunctionFieldElem)
    encoded_denominator = save_internal(s, denominator(f))
    encoded_numerator = save_internal(s, numerator(f))
    parents = encoded_numerator[:parents]
    # this removes the unnecessary polynomial ring for the denominator and
    # numerator
    pop!(parents)
    push!(parents, save_as_ref(s, parent(f)))
    return Dict(
        :parents => parents,
        :terms => (encoded_numerator[:terms], encoded_denominator[:terms])
    )
end

function load_terms(s::DeserializerState, parents::Vector, terms::Vector,
                    parent_ring::AbstractAlgebra.Generic.RationalFunctionField)
    num_coeff, den_coeff = terms
    pushfirst!(parents, base_ring(AbstractAlgebra.Generic.fraction_field(parent_ring)))
    loaded_num = load_terms(s, parents[1:end - 1], num_coeff, parents[end - 1])
    loaded_den = load_terms(s, parents[1:end - 1], den_coeff, parents[end - 1])
    return  parent_ring(loaded_num, loaded_den)
end
    
function load_internal(s::DeserializerState,
                       ::Type{<: AbstractAlgebra.Generic.RationalFunctionFieldElem},
                       dict::Dict)
    parents = load_parents(dict[:parents])
    return load_terms(s, parents, dict[:terms], parents[end])
end

function load_internal_with_parent(s::DeserializerState,
                                   ::Type{<: AbstractAlgebra.Generic.RationalFunctionFieldElem},
                                   dict::Dict,
                                   parent:: AbstractAlgebra.Generic.RationalFunctionField)
    forced_parent = base_ring(AbstractAlgebra.Generic.fraction_field(parent))
    num = load_unknown_type(s, dict[:num]; parent=forced_parent)
    den = load_unknown_type(s, dict[:den]; parent=forced_parent)

    return parent(num, den)
end

################################################################################
# ArbField
@registerSerializationType(ArbField, true)
@registerSerializationType(arb)

function save_internal(s::SerializerState, RR::Nemo.ArbField)
    return Dict(
        :precision => save_type_dispatch(s, precision(RR))
    )    
end

function load_internal(s::DeserializerState, ::Type{Nemo.ArbField}, dict::Dict)
    prec = load_type_dispatch(s, Int64, dict[:precision])
    return Nemo.ArbField(prec)
end

# elements
function save_internal(s::SerializerState, r::arb)
    c_str = ccall((:arb_dump_str, Nemo.Arb_jll.libarb), Ptr{UInt8}, (Ref{arb},), r)
    arb_unsafe_str = unsafe_string(c_str)

    # free memory
    ccall((:flint_free, Nemo.libflint), Nothing, (Ptr{UInt8},), c_str)

    return Dict(
        :parent => save_type_dispatch(s, parent(r)),
        :arb_unsafe_str => save_type_dispatch(s, arb_unsafe_str)
    )    
end


function load_internal(s::DeserializerState, ::Type{arb}, dict::Dict)
    parent = load_type_dispatch(s, Nemo.ArbField, dict[:parent])
    arb_unsafe_str = load_type_dispatch(s, String, dict[:arb_unsafe_str])
    r = Nemo.arb()
    ccall((:arb_load_str, Nemo.Arb_jll.libarb),
          Int32, (Ref{arb}, Ptr{UInt8}), r, arb_unsafe_str)
    r.parent = parent
    
    return r
end

function load_internal_with_parent(s::DeserializerState,
                                   ::Type{arb},
                                   dict::Dict,
                                   parent::Nemo.ArbField)
    arb_unsafe_str = load_type_dispatch(s, String, dict[:arb_unsafe_str])
    r = Nemo.arb()
    ccall((:arb_load_str, Nemo.Arb_jll.libarb),
          Int32, (Ref{arb}, Ptr{UInt8}), r, arb_unsafe_str)
    r.parent = parent
    
    return r
end

################################################################################
# AcbField

@registerSerializationType(AcbField, true)
@registerSerializationType(acb)

function save_internal(s::SerializerState, CC::AcbField)
    return Dict(
        :precision => save_type_dispatch(s, precision(CC))
    )    
end

function load_internal(s::DeserializerState, ::Type{AcbField}, dict::Dict)
    prec = load_type_dispatch(s, Int64, dict[:precision])
    return AcbField(prec)
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

function load_internal_with_parent(s::DeserializerState,
                                   ::Type{acb},
                                   dict::Dict,
                                   parent::AcbField)
    real_part = load_type_dispatch(s, arb, dict[:real])
    imag_part = load_type_dispatch(s, arb, dict[:imag])

    return parent(real_part, imag_part)
end

################################################################################
# Field Embeddings

@registerSerializationType(Hecke.NumFieldEmbNfAbs, true)

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

@registerSerializationType(Hecke.NumFieldEmbNfAbsNS, true)

function save_internal(s::SerializerState, E::Hecke.NumFieldEmbNfAbsNS)
    K = number_field(E)
    gen_balls = map(E, gens(K))
    
    return Dict(
        :num_field => save_type_dispatch(s, K),
        :gen_balls => save_type_dispatch(s, gen_balls)
    )
end

function load_internal(s::DeserializerState, ::Type{Hecke.NumFieldEmbNfAbsNS}, dict::Dict)
    K = load_type_dispatch(s, NfAbsNS, dict[:num_field])
    gen_balls = load_type_dispatch(s, Vector{acb}, dict[:gen_balls])

    return complex_embedding(K, gen_balls)
end

################################################################################
# Padic Field

@registerSerializationType(FlintPadicField)
@registerSerializationType(padic)

function save_internal(s::SerializerState, P::FlintPadicField)
    return Dict(
        :prime => save_type_dispatch(s, prime(P)),
        :precision => save_type_dispatch(s, precision(P))
    )
end

function load_internal(s::DeserializerState, ::Type{FlintPadicField}, dict::Dict)
    prime_num = load_type_dispatch(s, ZZRingElem, dict[:prime])
    precision = load_type_dispatch(s, Int64, dict[:precision])

    return PadicField(prime_num, precision)
end

#elements

function save_internal(s::SerializerState, n::padic)
    return string(lift(QQ, n))
end

function load_internal(s::DeserializerState, ::Type{padic}, dict::Dict)
    rational_rep = load_type_dispatch(s, QQFieldElem, dict[:rational_rep])
    parent_field = load_type_dispatch(s, FlintPadicField, dict[:parent])

    return parent_field(rational_rep)
end

function load_internal_with_parent(s::DeserializerState,
                                   ::Type{padic},
                                   str::String,
                                   parent::FlintPadicField)
    rational_rep = load_type_dispatch(s, QQFieldElem, str)
    return parent(rational_rep)
end
