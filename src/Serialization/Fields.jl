################################################################################
# Utility functions for parent tree
function get_parents(parent_ring::Field)
    if has_elem_basic_encoding(parent_ring)
        return Any[]
    end

    if absolute_degree(parent_ring) == 1
        return Any[]
    end
    base = parent(defining_polynomial(parent_ring))
    parents = get_parents(base)
    push!(parents, parent_ring)
    return parents
end

function get_parents(parent_ring::T) where T <: Union{NfAbsNS, NfRelNS}
    n = ngens(parent_ring)
    base = polynomial_ring(base_field(parent_ring), n)[1]
    parents = get_parents(base)
    push!(parents, parent_ring)
    return parents
end

function get_parents(parent_ring::T) where T <: Union{FracField,
                                                      AbstractAlgebra.Generic.RationalFunctionField,
                                                      AbstractAlgebra.Generic.LaurentSeriesField}
    if has_elem_basic_encoding(parent_ring)
        return Any[]
    end
    
    base = base_ring(parent_ring)
    parents = get_parents(base)
    push!(parents, parent_ring)
    return parents
end

################################################################################
# field of rationals (singleton type)
@registerSerializationType(QQField)

################################################################################
# non-ZZRingElem variant
@registerSerializationType(Nemo.fpFieldElem)
@registerSerializationType(Nemo.fpField)
has_elem_basic_encoding(obj::Nemo.fpField) = true

function save_internal(s::SerializerState, F::Nemo.fpField)
    return Dict(
        :characteristic => UInt64(characteristic(F))
    )
end

function load_internal(s::DeserializerState, ::Type{Nemo.fpField}, dict::Dict)
    return Nemo.fpField(UInt64(dict[:characteristic]))
end

# elements
function save_internal(s::SerializerState, elem::fpFieldElem; include_parents::Bool=true)
    if include_parents
        return Dict(
            :parent => save_as_ref(s, parent(elem)),
            :class_rep => string(elem)
        )
    end
    return string(elem)
end

function load_internal(s::DeserializerState, z::Type{fpFieldElem}, dict::Dict)
    F = load_type_dispatch(s, Nemo.fpField, dict[:parent])
    return F(parse(UInt64, dict[:class_rep]))
end

function load_internal_with_parent(s::DeserializerState,
                                   z::Type{fpFieldElem},
                                   int::Int,
                                   parent::Nemo.fpField)
    return parent(UInt64(int))
end

function load_internal_with_parent(s::DeserializerState,
                                   z::Type{fpFieldElem},
                                   str::String,
                                   parent::Nemo.fpField)
    return parent(parse(UInt64,str))
end



################################################################################
# ZZRingElem variant
@registerSerializationType(Nemo.FpFieldElem)
@registerSerializationType(Nemo.FpField)
has_elem_basic_encoding(obj::Nemo.FpField) = true

function save_internal(s::SerializerState, F::Nemo.FpField)
    return Dict(
        :characteristic => save_type_dispatch(s, characteristic(F))
    )
end

function load_internal(s::DeserializerState, F::Type{Nemo.FpField}, dict::Dict)
    return F(load_type_dispatch(s, ZZRingElem, dict[:characteristic]))
end

# elements
function save_internal(s::SerializerState, elem::FpFieldElem; include_parents::Bool=true)
    if include_parents
        return Dict(
            :parent => save_as_ref(s, parent(elem)),
            :class_rep => string(Nemo.data(elem))
        )
    end
    return string(Nemo.data(elem))
end

function load_internal(s::DeserializerState, z::Type{FpFieldElem}, dict::Dict)
    F = load_type_dispatch(s, Nemo.FpField, dict[:parent])
    return F(ZZRingElem(dict[:class_rep]))
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

function save_internal(s::SerializerState, k::Union{nf_elem, fqPolyRepFieldElem, Hecke.NfRelElem};
                       include_parents::Bool=true)
    K = parent(k)
    polynomial = parent(defining_polynomial(K))(k)

    # currently we get parent(defining_polynomial(K)) == parent(defining_polynomial(K)) = false
    # which leads to duplicate refs in the refs section of the file (not in the parent list)
    terms = save_internal(s, polynomial; include_parents=false)
    
    if include_parents
        return Dict(
            :parents => get_parent_refs(s, K),
            :terms => terms
        )
    end
    return terms
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
has_elem_basic_encoding(obj::FqField) = absolute_degree(obj) == 1

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
function save_internal(s::SerializerState, k::FqFieldElem; include_parents::Bool=true)
    K = parent(k)
    
    if absolute_degree(K) == 1
        class_rep = string(lift(ZZ, k))

        if include_parents
            return Dict(
                :parents => get_parent_refs(s, K),
                :class_rep => class_rep
            )
        end
        return class_rep
    end

    poly_parent = parent(defining_polynomial(K))
    parent_base_ring = base_ring(poly_parent)
    # currently this lift won't work for the given types
    # but is necessary for serialization
    polynomial = lift(poly_base_ring, k)
    terms = save_internal(s, polynomial; include_parents=false)
    
    if include_parents
        return Dict(
        :parents => load_parent_refs
        :terms => terms
        )
    end
    return terms
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
    return parent_ring(loaded_terms)
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

function save_internal(s::SerializerState, k::Union{NfAbsNSElem, Hecke.NfRelNSElem};
                       include_parents::Bool=true)
    K = parent(k)
    polynomial = Oscar.Hecke.data(k)
    polynomial_parent = parent(polynomial)
    terms = save_internal(s, polynomial; include_parents=false)

    if include_parents
        return Dict(
            :parents => get_parent_refs(s, K),
            :terms => terms
        )
    end
    return terms
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
    K = load_parents(s, dict[:parents])[end]
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
    n = ngens(parent_field)
    parent_polynomial_ring, _ = polynomial_ring(base_field(parent_field), n)
    polynomial = load_unknown_type(s, dict[:polynomial]; parent=parent_polynomial_ring)
    polynomial = evaluate(polynomial, gens(parent_field))

    return parent_field(polynomial)
end


################################################################################
# FracField

@registerSerializationType(FracField, true)

function save_internal(s::SerializerState, K::FracField)
    return Dict(
        :base_ring => save_as_ref(s, base_ring(K)),
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

function save_internal(s::SerializerState, f::FracElem; include_parents::Bool=true)
    encoded_denominator = save_internal(s, denominator(f); include_parents=false)
    encoded_numerator = save_internal(s, numerator(f); include_parents=false)
    terms = (encoded_numerator, encoded_denominator)
    
    if include_parents
        return Dict(
            :parents => get_parent_refs(s, parent(f)),
            :terms => terms
        )
    end
    return terms
end

function load_terms(s::DeserializerState, parents::Vector, terms::Vector,
                    parent_ring::FracField)
    num_coeff, den_coeff = terms
    loaded_num = load_terms(s, parents[1:end - 1], num_coeff, parents[end - 1])
    loaded_den = load_terms(s, parents[1:end - 1], den_coeff, parents[end-1])
    return  parent_ring(loaded_num) // parent_ring(loaded_den)
end

function load_internal(s::DeserializerState,
                       ::Type{<:FracElem},
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

function save_internal(s::SerializerState, f::AbstractAlgebra.Generic.RationalFunctionFieldElem;
                       include_parents::Bool=true)
    encoded_denominator = save_internal(s, denominator(f); include_parents=false)
    encoded_numerator = save_internal(s, numerator(f); include_parents=false)
    terms = (encoded_numerator, encoded_denominator)

    if include_parents
        return Dict(
            :parents => get_parent_refs(s, parent(f)),
            :terms => terms
        )
    end
    return terms
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
@registerSerializationType(ArbField)
@registerSerializationType(arb)
has_elem_basic_encoding(obj::ArbField) = true

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
function save_internal(s::SerializerState, r::arb; include_parents::Bool=true)
    c_str = ccall((:arb_dump_str, Nemo.Arb_jll.libarb), Ptr{UInt8}, (Ref{arb},), r)
    arb_unsafe_str = unsafe_string(c_str)

    # free memory
    ccall((:flint_free, Nemo.libflint), Nothing, (Ptr{UInt8},), c_str)
    if include_parents
        return Dict(
            :parents => get_parent_refs(s, parent(r)),
            :arb_str => arb_unsafe_str
        )
    end
    return arb_unsafe_str
end

function load_internal(s::DeserializerState, ::Type{arb}, dict::Dict)
    r = Nemo.arb()
    ccall((:arb_load_str, Nemo.Arb_jll.libarb),
          Int32, (Ref{arb}, Ptr{UInt8}), r, dict[:arb_str])

    parent = load_parents(s, dict[:parents])[end]
    r.parent = parent
    return r
end

function load_internal_with_parent(s::DeserializerState,
                                   ::Type{arb},
                                   str::String,
                                   parent::Nemo.ArbField)
    r = Nemo.arb()
    ccall((:arb_load_str, Nemo.Arb_jll.libarb),
          Int32, (Ref{arb}, Ptr{UInt8}), r, str)
    r.parent = parent
    
    return r
end

################################################################################
# AcbField
@registerSerializationType(AcbField)
@registerSerializationType(acb)
has_elem_basic_encoding(obj::AcbField) = true

function save_internal(s::SerializerState, CC::AcbField)
    return Dict(
        :precision => save_type_dispatch(s, precision(CC))
    )    
end

function load_internal(s::DeserializerState, ::Type{AcbField}, dict::Dict)
    prec = load_type_dispatch(s, Int64, dict[:precision])
    return AcbField(prec)
end

# elements
function save_internal(s::SerializerState, c::acb; include_parents::Bool=true)
    encoded_acb = save_internal(s, [real(c), imag(c)]; include_parents=false)
    if include_parents
        return Dict(
            :parent => save_as_ref(s, parent(c)),
            :vector => encoded_acb
        )
    end
    return encoded_acb
end

function load_internal_with_parent(s::DeserializerState,
                                   ::Type{acb},
                                   vec::Vector{Any},
                                   parent::AcbField)
    real_part = load_type_dispatch(s, arb, vec[1], parent=ArbField(precision(parent)))
    imag_part = load_type_dispatch(s, arb, vec[2], parent=ArbField(precision(parent)))
    
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
        :gen_ball => save_internal(s, g_ball)
    )
end

function load_internal(s::DeserializerState, ::Type{Hecke.NumFieldEmbNfAbs}, dict::Dict)
    K = load_type_dispatch(s, AnticNumberField, dict[:num_field])
    parent = load_type_dispatch(s, AcbField, dict[:gen_ball][:parent])
    gen_ball = load_internal_with_parent(s, acb, dict[:gen_ball][:vector], parent)

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
has_elem_basic_encoding(obj::FlintPadicField) = true

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

function save_internal(s::SerializerState, n::padic; include_parents::Bool=true)
    rational_rep = string(lift(QQ, n))
    if include_parents
        return Dict(
            :parent => save_as_ref(s, parent(n)),
            :rational_rep => rational_rep
        )
    end
    return rational_rep
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
