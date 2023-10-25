################################################################################
# Utility functions for parent tree
function get_parents(parent_ring::Field)
  # we have reached the end of the parent references and the current ring
  # can be found as the base_ring of the previous parent without ambiguity
  if !serialize_with_id(parent_ring)
    return RingMatSpaceUnion[]
  end

  if absolute_degree(parent_ring) == 1
    return RingMatSpaceUnion[]
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
  # we have reached the end of the parent references and the current ring
  # can be found as the base_ring of the previous parent without ambiguity
  if !serialize_with_id(parent_ring)
    return RingMatSpaceUnion[]
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
@register_serialization_type QQField

################################################################################
# non-ZZRingElem variant
@register_serialization_type Nemo.fpField

function save_object(s::SerializerState, F::Nemo.fpField)
  save_object(s, string(characteristic(F)))
end

function load_object(s::DeserializerState, ::Type{Nemo.fpField}, str::String)
  return Nemo.fpField(parse(UInt64, str))
end

# elements
@register_serialization_type fpFieldElem uses_params

function save_object(s::SerializerState, elem::fpFieldElem)
  save_data_basic(s, string(elem))
end

function load_object(s::DeserializerState, ::Type{fpFieldElem},
                                 str::String, F::Nemo.fpField)
  return F(parse(UInt64, str))
end

################################################################################
# ZZRingElem variant
@register_serialization_type Nemo.FpField

function save_object(s::SerializerState, F::Nemo.FpField)
  save_object(s, string(characteristic(F)))
end

function load_object(s::DeserializerState, ::Type{Nemo.FpField}, str::String)
  return Nemo.FpField(parse(ZZRingElem, str))
end

# elements
@register_serialization_type FpFieldElem uses_params

function save_object(s::SerializerState, elem::FpFieldElem)
  save_data_basic(s, string(elem))
end

function load_object(s::DeserializerState, ::Type{FpFieldElem},
                                 str::String, F::Nemo.FpField)
  return F(parse(ZZRingElem, str))
end

################################################################################
# SimpleNumField

@register_serialization_type Hecke.NfRel uses_id
@register_serialization_type AnticNumberField uses_id

function save_object(s::SerializerState, K::SimpleNumField)
  save_data_dict(s) do 
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
@register_serialization_type fqPolyRepField uses_id

function save_object(s::SerializerState, K::fqPolyRepField)
  save_data_dict(s) do
    save_typed_object(s, defining_polynomial(K), :def_pol)
  end
end

function load_object(s::DeserializerState,
                     ::Type{<: fqPolyRepField},
                     dict::Dict)
  def_pol = load_typed_object(s, dict[:def_pol])
  K, _ = finite_field(def_pol, cached=false)
  return K
end

#elements
@register_serialization_type fqPolyRepFieldElem uses_params
@register_serialization_type nf_elem uses_params
@register_serialization_type Hecke.NfRelElem uses_params
const NumFieldElemTypeUnion = Union{nf_elem, fqPolyRepFieldElem, Hecke.NfRelElem}

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
@register_serialization_type FqField uses_id
@register_serialization_type FqFieldElem uses_params

function save_object(s::SerializerState, K::FqField)
  if absolute_degree(K) == 1
    save_object(s, order(K))
  else
    save_data_dict(s) do
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

@register_serialization_type Hecke.NfRelNS uses_id
@register_serialization_type NfAbsNS uses_id

function save_object(s::SerializerState, K::Union{NfAbsNS, NfRelNS})
  def_pols = defining_polynomials(K)
  save_data_dict(s) do
    save_typed_object(s, def_pols, :def_pols)
    save_object(s, vars(K), :vars)
  end
end

function load_object(s::DeserializerState,
                     ::Type{<: Union{NfAbsNS, NfRelNS}},
                     dict::Dict{Symbol, Any})
  def_pols = load_typed_object(s, dict[:def_pols])
  vars = map(Symbol, dict[:vars])
  # fix since numberfield doesn't accept PolyRingElem vectors
  array_pols = Array{typeof(def_pols[1]), 1}(def_pols)
  K, _ = number_field(array_pols, vars, cached=false)
  return K
end

#elements
@register_serialization_type Hecke.NfRelNSElem uses_params
@register_serialization_type NfAbsNSElem uses_params

function save_object(s::SerializerState, k::Union{NfAbsNSElem, Hecke.NfRelNSElem})
  polynomial = Oscar.Hecke.data(k)
  save_object(s, polynomial)
end

function load_object(s::DeserializerState, ::Type{<: Union{NfAbsNSElem, Hecke.NfRelNSElem}}, terms::Vector, parents::Vector)
  K = parents[end]
  n = ngens(K)
  # forces parent of MPolyRingElem
  poly_ring = polynomial_ring(base_field(K), n)
  parents[end - 1], _ = poly_ring
  poly_elem_type = elem_type
  polynomial = load_object(s, MPolyRingElem, terms, parents[1:end - 1])
  polynomial = evaluate(polynomial, gens(K))
  return K(polynomial)
end

################################################################################
# FracField

@register_serialization_type FracField uses_id

function save_object(s::SerializerState, K::FracField)
  save_data_dict(s) do 
    save_typed_object(s, base_ring(K), :base_ring)
  end
end

function load_object(s::DeserializerState, ::Type{<: FracField}, dict::Dict)
  R = load_typed_object(s, dict[:base_ring])

  return fraction_field(R, cached=false)
end

function load_object(s::DeserializerState, ::Type{<: FracField}, str::String)
  R = load_ref(s, str)

  return fraction_field(R, cached=false)
end

# elements
# we use the union to prevent QQFieldElem being registered here
@register_serialization_type FracElem{<:Union{MPolyRingElem, PolyRingElem, UniversalPolyRingElem}} "FracElem" uses_params

function save_object(s::SerializerState, f::FracElem)
  save_data_array(s) do
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

@register_serialization_type AbstractAlgebra.Generic.RationalFunctionField "RationalFunctionField" uses_id

function save_object(s::SerializerState,
                     RF::AbstractAlgebra.Generic.RationalFunctionField)
  save_data_dict(s) do
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
  return rational_function_field(R, symbols, cached=false)[1]
end

#elements
@register_serialization_type AbstractAlgebra.Generic.RationalFunctionFieldElem "RationalFunctionFieldElem" uses_params

function save_object(s::SerializerState, f::AbstractAlgebra.Generic.RationalFunctionFieldElem)
  save_data_array(s) do
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
@register_serialization_type ArbField
@register_serialization_type arb uses_params

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
@register_serialization_type AcbField
@register_serialization_type acb uses_params

function save_object(s::SerializerState, CC::AcbField)
  save_object(s, precision(CC))
end

function load_object(s::DeserializerState, ::Type{AcbField}, str::String)
  prec = parse(Int, str)
  return AcbField(prec)
end

# elements
function save_object(s::SerializerState, c::acb)
  save_data_array(s) do
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

@register_serialization_type Hecke.NumFieldEmbNfAbs  uses_id

function save_object(s::SerializerState, E::Hecke.NumFieldEmbNfAbs)
  K = number_field(E)
  g = gen(K)
  g_ball = E(g)

  save_data_dict(s) do
    save_typed_object(s, K, :num_field)
    save_typed_object(s, g_ball, :gen_ball)
  end
end

function load_object(s::DeserializerState, ::Type{Hecke.NumFieldEmbNfAbs}, dict::Dict)
  K = load_typed_object(s, dict[:num_field])
  gen_ball = load_typed_object(s, dict[:gen_ball])

  return complex_embedding(K, gen_ball)
end

@register_serialization_type Hecke.NumFieldEmbNfAbsNS uses_id

function save_object(s::SerializerState, E::Hecke.NumFieldEmbNfAbsNS)
  K = number_field(E)
  gen_balls = map(E, gens(K))

  save_data_dict(s) do
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
@register_serialization_type FlintPadicField

function save_object(s::SerializerState, P::FlintPadicField)
  save_data_dict(s) do
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
@register_serialization_type padic uses_params

function save_object(s::SerializerState, obj::padic)
  # currently it seems padics do not store the underlying polynomial
  save_object(s, lift(QQ, obj))
end

function load_object(s::DeserializerState, ::Type{padic},
                     str::String, parent_field::FlintPadicField)
  rational_rep = load_object(s, QQFieldElem, str)

  return parent_field(rational_rep)
end
