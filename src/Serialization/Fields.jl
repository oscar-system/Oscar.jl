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

function get_parents(e::EmbeddedField)
  base = number_field(e)
  parents = get_parents(base)
  push!(parents, e)
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
function load_object(s:: DeserializerState, ::Type{Field})
  return load_typed_object(s)
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

function load_object(s::DeserializerState, ::Type{Nemo.fpField})
  load_node(s) do str
    return Nemo.fpField(parse(UInt64, str))
  end
end

# elements
@register_serialization_type fpFieldElem uses_params

function save_object(s::SerializerState, elem::fpFieldElem)
  save_data_basic(s, string(elem))
end

function load_object(s::DeserializerState, ::Type{fpFieldElem}, F::Nemo.fpField)
  load_node(s) do str
    return F(parse(UInt64, str))
  end
end

################################################################################
# ZZRingElem variant
@register_serialization_type Nemo.FpField

function save_object(s::SerializerState, F::Nemo.FpField)
  save_object(s, string(characteristic(F)))
end

function load_object(s::DeserializerState, ::Type{Nemo.FpField})
  load_node(s) do str
    Nemo.FpField(parse(ZZRingElem, str))
  end
end

# elements
@register_serialization_type FpFieldElem uses_params

function save_object(s::SerializerState, elem::FpFieldElem)
  save_data_basic(s, string(elem))
end

function load_object(s::DeserializerState, ::Type{FpFieldElem}, F::Nemo.FpField)
  load_node(s) do str
    F(parse(ZZRingElem, str))
  end
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

function load_object(s::DeserializerState, ::Type{<: SimpleNumField})
  def_pol = load_typed_object(s, :def_pol)
  var = load_node(s, :var) do var
    Symbol(var)
  end
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

function load_object(s::DeserializerState, ::Type{<: fqPolyRepField})
  def_pol = load_typed_object(s, :def_pol)
  K, _ = Nemo.Native.finite_field(def_pol, cached=false)
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

function save_object(s::SerializerState, k::Hecke.NfRelElem{NfAbsNSElem})
  K = parent(k)
  polynomial = parent(defining_polynomial(K))(data(k))
  save_object(s, polynomial)
end

function load_object(s::DeserializerState, ::Type{<: NumFieldElemTypeUnion},
                     parents::Vector)
  polynomial = load_node(s) do _
    load_object(s, PolyRingElem, parents[1:end - 1])
  end
  
  K = parents[end]
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

function load_object(s::DeserializerState, ::Type{<: FqField})
  load_node(s) do node
    if node isa String
      order = ZZRingElem(node)
      return finite_field(order)[1]
    else
      def_pol = load_typed_object(s)
      return finite_field(def_pol, cached=false)[1]
    end
  end
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

function load_object(s::DeserializerState, ::Type{<: FqFieldElem}, parents::Vector)
  K = parents[end]
  load_node(s) do _
    K(load_object(s, PolyRingElem, parents[1:end - 1]))
  end
end

function load_object(s::DeserializerState, ::Type{<: FqFieldElem}, parent::FqField)
  load_node(s) do str
    parent(ZZRingElem(str))
  end
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

function load_object(s::DeserializerState, ::Type{<: Union{NfAbsNS, NfRelNS}})
  def_pols = load_typed_object(s, :def_pols)

  vars = load_node(s, :vars) do vars_data
    return map(Symbol, vars_data)
  end
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

function load_object(s::DeserializerState, ::Type{<: Union{NfAbsNSElem, Hecke.NfRelNSElem}},
                     parents::Vector)
  K = parents[end]
  n = ngens(K)
  # forces parent of MPolyRingElem
  poly_ring = polynomial_ring(base_field(K), n)
  parents[end - 1], _ = poly_ring
  poly_elem_type = elem_type
  load_node(s) do _
    polynomial = load_object(s, MPolyRingElem, parents[1:end - 1])
  end
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

function load_object(s::DeserializerState, ::Type{<: FracField})
  R = load_typed_object(s, :base_ring)

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

function load_object(s::DeserializerState, ::Type{<: FracElem}, parents::Vector)
  parent_ring = parents[end]
  load_node(s) do _
    coeff_type = elem_type(base_ring(parent_ring))
    loaded_num = load_node(s, 1) do _
      load_object(s, coeff_type, parents[1:end - 1])
    end
    loaded_den = load_node(s, 2) do _
      load_object(s, coeff_type, parents[1:end - 1])
    end
    return  parent_ring(loaded_num, loaded_den)
  end
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
                     ::Type{<: AbstractAlgebra.Generic.RationalFunctionField})
  R = load_typed_object(s, :base_ring)
  # ensure proper types of univariate case on load
  symbols = load_node(s, :symbols) do symbols_data
    if symbols_data isa Vector
      return Symbol.(symbols_data)
    else
      return Symbol(symbols_data)
    end
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

function load_object(s::DeserializerState,
                     ::Type{<: AbstractAlgebra.Generic.RationalFunctionFieldElem},
                     parents::Vector)
  parent_ring = parents[end]
  base = base_ring(AbstractAlgebra.Generic.fraction_field(parent_ring))
  pushfirst!(parents, base)
  coeff_type = elem_type(base)

  return load_node(s) do _
    loaded_num = load_node(s, 1) do _
      load_object(s, coeff_type, parents[1:end - 1])
    end
    
    loaded_den = load_node(s, 2) do _
      load_object(s, coeff_type, parents[1:end - 1])
    end
    parent_ring(loaded_num, loaded_den)
  end
end

################################################################################
# ArbField
@register_serialization_type ArbField
@register_serialization_type arb uses_params

function save_object(s::SerializerState, RR::Nemo.ArbField)
  save_object(s, precision(RR))
end

function load_object(s::DeserializerState, ::Type{Nemo.ArbField})
  prec = load_object(s, Int64, :precision)
  return Nemo.ArbField(prec)
end

# elements
function save_object(s::SerializerState, r::arb)
  c_str = ccall((:arb_dump_str, Nemo.Arb_jll.libarb), Ptr{UInt8}, (Ref{arb},), r)
  save_object(s, unsafe_string(c_str))
  
  # free memory
  ccall((:flint_free, Nemo.libflint), Nothing, (Ptr{UInt8},), c_str)
end

function load_object(s::DeserializerState, ::Type{arb}, parent::ArbField)
  r = Nemo.arb()
  load_node(s) do str
    ccall((:arb_load_str, Nemo.Arb_jll.libarb),
          Int32, (Ref{arb}, Ptr{UInt8}), r, str)
  end
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

function load_object(s::DeserializerState, ::Type{AcbField})
  prec = load_object(s, Int)
  return AcbField(prec)
end

# elements
function save_object(s::SerializerState, c::acb)
  save_data_array(s) do
    save_object(s, real(c))
    save_object(s, imag(c))
  end
end

function load_object(s::DeserializerState, ::Type{acb}, parent::AcbField)
  (real_part, imag_part) = load_array_node(s) do _
    load_object(s, arb, ArbField(precision(parent)))
  end
  return parent(real_part, imag_part)
end

################################################################################
# Field Embeddings

const FieldEmbeddingTypes = Union{Hecke.NumFieldEmbNfAbs, Hecke.NumFieldEmbNfRel, Hecke.NumFieldEmbNfAbsNS, Hecke.NumFieldEmbNfNS}

@register_serialization_type Hecke.NumFieldEmbNfAbsNS uses_id
@register_serialization_type Hecke.NumFieldEmbNfNS uses_id
@register_serialization_type Hecke.NumFieldEmbNfAbs uses_id
@register_serialization_type Hecke.NumFieldEmbNfRel uses_id

function save_object(s::SerializerState, E::FieldEmbeddingTypes)
  K = number_field(E)
  k = base_field(K)

  save_data_dict(s) do
    save_typed_object(s, K, :num_field)
    if !(base_field(K) isa QQField)
      save_typed_object(s, restrict(E, k), :base_field_emb)
    end
    if is_simple(K)
      a = gen(K)
      data = E(a)
      if any(overlaps(data, e(a)) for e in complex_embeddings(K) if e != E && restrict(E, k) == restrict(e, k))
        error("Internal error in internal serialization.")
      end
    else
      a = gens(K)
      data = E.(a)
      if any(all(overlaps(t[1], t[2]) for t in zip(data, e.(a))) for e in complex_embeddings(K) if e != E && restrict(E, k) == restrict(e, k))
        error("Internal error in internal serialization.")
      end
      data = tuple(data...)
    end
    save_typed_object(s, data, :data)
  end
end

function load_object(s::DeserializerState, ::Type{<:FieldEmbeddingTypes})
  K = load_typed_object(s, :num_field)
  data = load_typed_object(s, :data)
  if data isa Tuple
    data = collect(data)
  end
  if base_field(K) isa QQField
    return complex_embedding(K, data)
  else
    return complex_embedding(K, load_typed_object(s, :base_field_emb), data)
  end
end

@register_serialization_type Hecke.EmbeddedField uses_id

function save_object(s::SerializerState, E::Hecke.EmbeddedField)
  K = number_field(E)
  e = embedding(E)

  save_data_dict(s) do
    save_typed_object(s, K, :num_field)
    save_typed_object(s, e, :embedding)
  end
end

function load_object(s::DeserializerState, ::Type{Hecke.EmbeddedField})
  K = load_typed_object(s, :num_field)
  e = load_typed_object(s, :embedding)

  return Hecke.embedded_field(K, e)[1]
end

@register_serialization_type EmbeddedElem uses_params

function save_object(s::SerializerState, f::EmbeddedElem)
  save_object(s, data(f))
end

function load_object(s::DeserializerState, ::Type{<:EmbeddedElem}, parents::Vector)
  parent_field = parents[end]
  numfield_elem = terms
  coeff_type = elem_type(parents[end - 1])
  loaded_alg_elem = load_object(s, coeff_type, parents[1:end - 1])
  return parent_field(loaded_alg_elem)
end

################################################################################
# Padic Field
@register_serialization_type PadicField

function save_object(s::SerializerState, P::PadicField)
  save_data_dict(s) do
    save_object(s, prime(P), :prime)
    save_object(s, precision(P), :precision)
  end
end

function load_object(s::DeserializerState, ::Type{PadicField})
  prime_num = load_node(s, :prime) do node
    return parse(ZZRingElem, node)
  end
  precision = load_node(s, :precision) do node
    return parse(Int64, node)
  end
  return PadicField(prime_num, precision)
end

#elements
@register_serialization_type PadicFieldElem uses_params

function save_object(s::SerializerState, obj::PadicFieldElem)
  # currently it seems PadicFieldElems do not store the underlying polynomial
  save_object(s, lift(QQ, obj))
end

function load_object(s::DeserializerState, ::Type{PadicFieldElem}, parent_field::PadicField)
  rational_rep = load_object(s, QQFieldElem)
  return parent_field(rational_rep)
end
