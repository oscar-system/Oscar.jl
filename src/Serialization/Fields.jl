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

function load_internal_with_parent(s::DeserializerState,
                                   z::Type{gfp_elem},
                                   dict::Dict,
                                   parent::Nemo.GaloisField)
    return parent(UInt64(dict[:data]))
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

function load_internal_with_parent(s::DeserializerState,
                                   ::Type{gfp_fmpz_elem},
                                   dict::Dict,
                                   parent::Nemo.GaloisFmpzField)
    return parent(load_type_dispatch(s, fmpz, dict[:data]))
end

################################################################################
# SimpleNumField

@registerSerializationType(Hecke.NfRel)

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

@registerSerializationType(Hecke.NfRelElem)

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

function load_internal_with_parent(s::DeserializerState,
                                   ::Type{<: Union{nf_elem, fq_nmod, Hecke.NfRelElem}},
                                   dict::Dict,
                                   parent_field::Union{FqNmodFiniteField, SimpleNumField})
    polynomial_parent = parent(defining_polynomial(parent_field))
    polynomial = load_unknown_type(s, dict[:polynomial]; parent=polynomial_parent)

    return parent_field(polynomial)
end



################################################################################
# Non Simple Extension

@registerSerializationType(Hecke.NfRelNS)

@registerSerializationType(NfAbsNS)
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
    K, _ = NumberField(def_pols, vars, cached=false)
    return K
end

#elements
@registerSerializationType(Hecke.NfRelNSElem)

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
    parent_polynomial_ring, _ = PolynomialRing(base_field(parent_field), ngens)
    polynomial = load_unknown_type(s, dict[:polynomial]; parent=parent_polynomial_ring)
    polynomial = evaluate(polynomial, gens(parent_field))

    return parent_field(polynomial)
end


################################################################################
# FracField

@registerSerializationType(FracField)

function save_internal(s::SerializerState, K::FracField)
    return Dict(
        :base_ring => save_type_dispatch(s, base_ring(K)),
    )
end

function load_internal(s::DeserializerState,
                       ::Type{<: FracField},
                       dict::Dict)
    R = load_unknown_type(s, dict[:base_ring])

    return FractionField(R, cached=false)
end

# elements
@registerSerializationType(FracElem)

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

function load_internal_with_parent(s::DeserializerState,
                                   ::Type{<: FracElem},
                                   dict::Dict,
                                   parent:: FracField)
    parts_parent = base_ring(parent)
    num = load_unknown_type(s, dict[:num]; parent=parts_parent)
    den = load_unknown_type(s, dict[:den]; parent=parts_parent)
    
    return parent(num, den)
end

################################################################################
# RationalFunctionField

@registerSerializationType(AbstractAlgebra.Generic.RationalFunctionField)

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
    symbols = load_unknown_type(s, dict[:symbols])

    return RationalFunctionField(R, symbols, cached=false)[1]
end

#elements
@registerSerializationType(AbstractAlgebra.Generic.Rat)

function save_internal(s::SerializerState, f::AbstractAlgebra.Generic.Rat)
    frac_elem_parent = save_type_dispatch(s, parent(denominator(f)))
    return Dict(
        :parent => save_type_dispatch(s, parent(f)),
        :frac_elem_parent => frac_elem_parent,
        :den => save_type_dispatch(s, denominator(f)),
        :num => save_type_dispatch(s, numerator(f))
    )
end

function load_internal(s::DeserializerState,
                       ::Type{<: AbstractAlgebra.Generic.Rat},
                       dict::Dict)
    _ = load_unknown_type(s, dict[:frac_elem_parent])
    R = load_type_dispatch(s, AbstractAlgebra.Generic.RationalFunctionField, dict[:parent])
    num = load_unknown_type(s, dict[:num])
    den = load_unknown_type(s, dict[:den])

    if num isa PolyElem
        num = evaluate(num, gen(R))
        den = evaluate(den, gen(R))
    else
        num = evaluate(num, gens(R))
        den = evaluate(den, gens(R))
    end

    return num // den
end

function load_internal_with_parent(s::DeserializerState,
                                   ::Type{<: AbstractAlgebra.Generic.Rat},
                                   dict::Dict,
                                   parent:: AbstractAlgebra.Generic.RationalFunctionField)
    _ = load_unknown_type(s, dict[:frac_elem_parent])
    num = load_unknown_type(s, dict[:num])
    den = load_unknown_type(s, dict[:den])

    if num isa PolyElem
        num = evaluate(num, gen(parent))
        den = evaluate(den, gen(parent))
    else
        num = evaluate(num, gens(parent))
        den = evaluate(den, gens(parent))
    end

    return num // den
end

################################################################################
# RealField
@registerSerializationType(ArbField)
@registerSerializationType(arb)

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

function load_internal_with_parent(s::DeserializerState,
                                   ::Type{arb},
                                   dict::Dict,
                                   parent::Nemo.RealField)
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
@registerSerializationType(acb)

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

@registerSerializationType(Hecke.NumFieldEmbNfAbsNS)

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
    prime_num = load_type_dispatch(s, fmpz, dict[:prime])
    precision = load_type_dispatch(s, Int64, dict[:precision])

    return PadicField(prime_num, precision)
end

#elements

function save_internal(s::SerializerState, n::padic)
    return Dict(
        :rational_rep => save_type_dispatch(s, lift(QQ, n)),
        :parent => save_type_dispatch(s, parent(n))
    )
end

function load_internal(s::DeserializerState, ::Type{padic}, dict::Dict)
    rational_rep = load_type_dispatch(s, fmpq, dict[:rational_rep])
    parent_field = load_type_dispatch(s, FlintPadicField, dict[:parent])

    return parent_field(rational_rep)
end

function load_internal_with_parent(s::DeserializerState,
                                   ::Type{padic},
                                   dict::Dict,
                                   parent::FlintPadicField)
    rational_rep = load_type_dispatch(s, fmpq, dict[:rational_rep])
    parent_field = load_type_dispatch(s, FlintPadicField, dict[:parent])

    # padic num precision is 1 higher than the field it lies in
    if precision(parent_field) > precision(parent)
        @warn("Precision Warning: given parent is less precise than serialized parent",
              maxlog=1)
    end
    
    return parent(rational_rep)
end
