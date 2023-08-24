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
# abstract Field load
function load_object(s::DeserializerState, ::Type{Field}, id::String)
  return load_ref(s, id)
end

function load_object(s:: DeserializerState, ::Type{Field}, dict::Dict{Symbol, Any})
  return load_typed_object(s, dict)
end

################################################################################
# field of rationals (singleton type)
@registerSerializationType(QQField)

################################################################################
# non-ZZRingElem variant
@registerSerializationType(Nemo.fpField)
has_elem_basic_encoding(obj::Nemo.fpField) = true

function save_object(s::SerializerState, F::Nemo.fpField)
  save_object(s, string(characteristic(F)))
end

function load_object(s::DeserializerState, ::Type{Nemo.fpField}, str::String)
  return Nemo.fpField(parse(UInt64, str))
end

# elements
@registerSerializationType(fpFieldElem)
type_needs_params(T::Type{fpFieldElem}) = true

function save_object(s::SerializerState, elem::fpFieldElem)
  data_basic(s, string(elem))
end

function load_object(s::DeserializerState, ::Type{fpFieldElem},
                                 str::String, F::Nemo.fpField)
  return F(parse(UInt64, str))
end

################################################################################
# ZZRingElem variant
@registerSerializationType(Nemo.FpField)
has_elem_basic_encoding(obj::Nemo.FpField) = true

function save_object(s::SerializerState, F::Nemo.FpField)
  save_object(s, string(characteristic(F)))
end

function load_object(s::DeserializerState, ::Type{Nemo.FpField}, str::String)
  return Nemo.FpField(parse(ZZRingElem, str))
end

# elements
@registerSerializationType(FpFieldElem)
type_needs_params(T::Type{FpFieldElem}) = true

function save_object(s::SerializerState, elem::FpFieldElem)
  data_basic(s, string(elem))
end

function load_object(s::DeserializerState, ::Type{FpFieldElem},
                                 str::String, F::Nemo.FpField)
  return F(parse(ZZRingElem, str))
end

################################################################################
# SimpleNumField

@registerSerializationType(Hecke.NfRel, true)
@registerSerializationType(AnticNumberField, true)

function save_object(s::SerializerState, K::SimpleNumField)
  data_dict(s) do 
    save_typed_object(s, defining_polynomial(K), :def_pol)
    save_object(s, var(K), :var)
  end
end

function load_object(s::DeserializerState, ::Type{<: SimpleNumField}, dict::Dict)
  def_pol = load_typed_object(s, dict[:def_pol])
  var = Symbol(dict[:var])
  K, _ = number_field(def_pol, var, cached=false)
  return K
end

################################################################################
# FqNmodfinitefield
@registerSerializationType(fqPolyRepField, true)

function save_object(s::SerializerState, K::fqPolyRepField)
  data_dict(s) do
    save_typed_object(s, defining_polynomial(K), :def_pol)
  end
end

function load_object(s::DeserializerState,
                     ::Type{<: fqPolyRepField},
                     dict::Dict)
  def_pol = load_typed_object(s, dict[:def_pol])
  K, _ = FiniteField(def_pol, cached=false)
  return K
end

#elements
@registerSerializationType(fqPolyRepFieldElem)
@registerSerializationType(nf_elem)
@registerSerializationType(Hecke.NfRelElem)
const NumFieldElemTypeUnion = Union{nf_elem, fqPolyRepFieldElem, Hecke.NfRelElem}
type_needs_params(T::Type{<:NumFieldElemTypeUnion}) = true

function save_object(s::SerializerState, k::NumFieldElemTypeUnion)
  K = parent(k)
  polynomial = parent(defining_polynomial(K))(k)
  save_object(s, polynomial)
end

function load_object(s::DeserializerState, ::Type{<: NumFieldElemTypeUnion},
                                 terms::Vector, parents::Vector)
  K = parents[end]
  polynomial = load_object(s, PolyRingElem, terms, parents[1:end - 1])
  loaded_terms = evaluate(polynomial, gen(K))
  return K(loaded_terms)
end

################################################################################
# FqField
@registerSerializationType(FqField, true)
@registerSerializationType(FqFieldElem)
type_needs_params(::Type{FqFieldElem}) = true
has_elem_basic_encoding(obj::FqField) = absolute_degree(obj) == 1

function save_object(s::SerializerState, K::FqField)
  if absolute_degree(K) == 1
    save_object(s, order(K))
  else
    data_dict(s) do
      save_typed_object(s, defining_polynomial(K))
    end
  end
end

function load_object(s::DeserializerState, ::Type{<: FqField}, str::String)
  order = ZZRingElem(str)
  return Hecke.Nemo._FiniteField(order)[1]
end

function load_object(s::DeserializerState, ::Type{<: FqField}, dict::Dict)
  def_pol = load_typed_object(s, dict)
  return Hecke.Nemo._FiniteField(def_pol, cached=false)[1]
end

# elements
function save_object(s::SerializerState, k::FqFieldElem)
  K = parent(k)
  
  if absolute_degree(K) == 1
    save_object(s, lift(ZZ, k))
  else
    poly_parent = parent(defining_polynomial(K))
    parent_base_ring = base_ring(poly_parent)
    # currently this lift won't work for the given types
    # but is necessary for serialization
    polynomial = lift(poly_parent, k)
    save_object(s, polynomial)
  end
end

function load_object(s::DeserializerState, ::Type{<: FqFieldElem},
                                 terms::Vector, parents::Vector)
  K = parents[end]
  return K(load_object(s, PolyRingElem, terms, parents[1:end - 1]))
end

function load_object(s::DeserializerState, ::Type{<: FqFieldElem},
                                 str::String, parent::FqField)
  return parent(ZZRingElem(str))
end

################################################################################
# Non Simple Extension

@registerSerializationType(Hecke.NfRelNS, true)
@registerSerializationType(NfAbsNS, true)

function save_object(s::SerializerState, K::Union{NfAbsNS, NfRelNS})
  def_pols = defining_polynomials(K)
  data_dict(s) do
    save_typed_object(s, def_pols, :def_pols)
    save_object(s, vars(K), :vars)
  end
end

function load_object(s::DeserializerState,
                     ::Type{<: Union{NfAbsNS, NfRelNS}},
                     dict::Dict{Symbol, Any})
  def_pols = load_typed_object(s, dict[:def_pols])
  vars = map(Symbol, dict[:vars])
  K, _ = number_field(Array(def_pols), vars, cached=false)
  return K
end

#elements
@registerSerializationType(Hecke.NfRelNSElem)
@registerSerializationType(NfAbsNSElem)
type_needs_params(::Type{<:Union{NfAbsNSElem, Hecke.NfRelNSElem}}) = true

function save_object(s::SerializerState, k::Union{NfAbsNSElem, Hecke.NfRelNSElem})
  polynomial = Oscar.Hecke.data(k)
  save_object(s, polynomial)
end

function load_object(s::DeserializerState, ::Type{<: Union{NfAbsNSElem, Hecke.NfRelNSElem}}, terms::Vector, parents::Vector)
  K = parents[end]
  n = ngens(K)
  # forces parent of MPolyElem
  poly_ring = polynomial_ring(base_field(K), n)
  parents[end - 1], _ = poly_ring
  poly_elem_type = elem_type
  polynomial = load_object(s, MPolyRingElem, terms, parents[1:end - 1])
  polynomial = evaluate(polynomial, gens(K))
  return K(polynomial)
end

################################################################################
# FracField

@registerSerializationType(FracField, true)

function save_object(s::SerializerState, K::FracField)
  save_typed_object(s, base_ring(K), :data)
end

function load_object(s::DeserializerState, ::Type{<: FracField}, dict::Dict)
  R = load_typed_object(s, dict)

  return fraction_field(R, cached=false)
end

function load_object(s::DeserializerState, ::Type{<: FracField}, str::String)
  R = load_ref(s, str)

  return fraction_field(R, cached=false)
end

# elements
@registerSerializationType(FracElem)
# might need a better way to block QQFieldElem here so that more FracElems can be serialized
type_needs_params(::Type{<:FracElem{T}}) where T <: Union{MPolyRingElem, PolyRingElem, UniversalPolyRingElem} = true

function save_object(s::SerializerState, f::FracElem)
  data_array(s) do
    save_object(s, numerator(f))
    save_object(s, denominator(f))
  end
end

function load_object(s::DeserializerState, ::Type{<: FracElem},
                                 terms::Vector, parents::Vector)
  parent_ring = parents[end]
  num_coeff, den_coeff = terms
  coeff_type = elem_type(base_ring(parent_ring))
  loaded_num = load_object(s, coeff_type, num_coeff, parents[1:end - 1])
  loaded_den = load_object(s, coeff_type, den_coeff, parents[1:end - 1])
  return  parent_ring(loaded_num, loaded_den)
end

################################################################################
# RationalFunctionField

@registerSerializationType(AbstractAlgebra.Generic.RationalFunctionField,
                           true,
                           "RationalFunctionField")

function save_object(s::SerializerState,
                     RF::AbstractAlgebra.Generic.RationalFunctionField)
  data_dict(s) do
    save_typed_object(s, base_ring(RF), :base_ring)
    syms = symbols(RF)
    save_object(s, syms, :symbols)
  end
end

function load_object(s::DeserializerState,
                     ::Type{<: AbstractAlgebra.Generic.RationalFunctionField},
                     dict::Dict)
  R = load_typed_object(s, dict[:base_ring])
  # ensure proper types of univariate case on load
  if dict[:symbols] isa Vector
    symbols = map(Symbol, dict[:symbols])
  else
    symbols = Symbol(dict[:symbols])
  end
  return RationalFunctionField(R, symbols, cached=false)[1]
end

#elements
@registerSerializationType(AbstractAlgebra.Generic.RationalFunctionFieldElem,
                           true,
                           "RationalFunctionFieldElem")
type_needs_params(::Type{<: AbstractAlgebra.Generic.RationalFunctionFieldElem}) = true

function save_object(s::SerializerState, f::AbstractAlgebra.Generic.RationalFunctionFieldElem)
  data_array(s) do
    save_object(s, numerator(f))
    save_object(s, denominator(f))
  end
end

function load_object(s::DeserializerState, ::Type{<: AbstractAlgebra.Generic.RationalFunctionFieldElem},
                                 terms::Vector, parents::Vector)
  parent_ring = parents[end]
  num_coeff, den_coeff = terms
  base = base_ring(AbstractAlgebra.Generic.fraction_field(parent_ring))
  pushfirst!(parents, base)
  coeff_type = elem_type(base)
  loaded_num = load_object(s, coeff_type, num_coeff, parents[1:end - 1])
  loaded_den = load_object(s, coeff_type, den_coeff, parents[1:end - 1])
  return  parent_ring(loaded_num, loaded_den)
end

################################################################################
# ArbField
@registerSerializationType(ArbField)
@registerSerializationType(arb)
has_elem_basic_encoding(obj::ArbField) = true
type_needs_params(::Type{arb}) = true

function save_object(s::SerializerState, RR::Nemo.ArbField)
  save_object(s, precision(RR))
end

function load_object(s::DeserializerState, ::Type{Nemo.ArbField}, dict::Dict)
  prec = parse(Int64, dict[:precision])
  return Nemo.ArbField(prec)
end

# elements
function save_object(s::SerializerState, r::arb)
  c_str = ccall((:arb_dump_str, Nemo.Arb_jll.libarb), Ptr{UInt8}, (Ref{arb},), r)
  save_object(s, unsafe_string(c_str))
  
  # free memory
  ccall((:flint_free, Nemo.libflint), Nothing, (Ptr{UInt8},), c_str)
end

function load_object(s::DeserializerState, ::Type{arb}, str::String, parent::ArbField)
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
type_needs_params(::Type{acb}) = true

function save_object(s::SerializerState, CC::AcbField)
  save_object(s, precision(CC))
end

function load_object(s::DeserializerState, ::Type{AcbField}, str::String)
  prec = parse(Int, str)
  return AcbField(prec)
end

# elements
function save_object(s::SerializerState, c::acb)
  data_array(s) do
    save_object(s, real(c))
    save_object(s, imag(c))
  end
end

function load_object(s::DeserializerState, ::Type{acb}, vec::Vector{Any}, parent::AcbField)
  real_part = load_object(s, arb, vec[1], ArbField(precision(parent)))
  imag_part = load_object(s, arb, vec[2], ArbField(precision(parent)))
  
  return parent(real_part, imag_part)
end

################################################################################
# Field Embeddings

@registerSerializationType(Hecke.NumFieldEmbNfAbs, true)

function save_object(s::SerializerState, E::Hecke.NumFieldEmbNfAbs)
  K = number_field(E)
  g = gen(K)
  g_ball = E(g)

  data_dict(s) do
    save_typed_object(s, K, :num_field)
    save_typed_object(s, g_ball, :gen_ball)
  end
end

function load_object(s::DeserializerState, ::Type{Hecke.NumFieldEmbNfAbs}, dict::Dict)
  K = load_typed_object(s, dict[:num_field])
  gen_ball = load_typed_object(s, dict[:gen_ball])

  return complex_embedding(K, gen_ball)
end

@registerSerializationType(Hecke.NumFieldEmbNfAbsNS, true)

function save_object(s::SerializerState, E::Hecke.NumFieldEmbNfAbsNS)
  K = number_field(E)
  gen_balls = map(E, gens(K))

  data_dict(s) do
    save_typed_object(s, K, :num_field)
    save_typed_object(s, gen_balls, :gen_balls)
  end
end

function load_object(s::DeserializerState, ::Type{Hecke.NumFieldEmbNfAbsNS}, dict::Dict)
  K = load_typed_object(s, dict[:num_field])
  gen_balls = load_typed_object(s, dict[:gen_balls])

  return complex_embedding(K, gen_balls)
end

################################################################################
# Padic Field
@registerSerializationType(FlintPadicField)

function save_object(s::SerializerState, P::FlintPadicField)
  data_dict(s) do
    save_object(s, prime(P), :prime)
    save_object(s, precision(P), :precision)
  end
end

function load_object(s::DeserializerState, ::Type{FlintPadicField}, dict::Dict)
  prime_num = parse(ZZRingElem, dict[:prime])
  precision = parse(Int64, dict[:precision])

  return PadicField(prime_num, precision)
end

#elements
@registerSerializationType(padic)
type_needs_params(::Type{padic}) = true

function save_object(s::SerializerState, obj::padic)
  # currently it seems padics do not store the underlying polynomial
  save_object(s, lift(QQ, obj))
end

function load_object(s::DeserializerState, ::Type{padic},
                     str::String, parent_field::FlintPadicField)
  rational_rep = load_object(s, QQFieldElem, str)

  return parent_field(rational_rep)
end
