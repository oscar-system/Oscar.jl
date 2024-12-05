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

function get_parents(e::EmbeddedNumField)
  base = number_field(e)
  parents = get_parents(base)
  push!(parents, e)
  return parents
end

function get_parents(parent_ring::T) where T <: Union{AbsNonSimpleNumField, RelNonSimpleNumField}
  n = ngens(parent_ring)
  base = polynomial_ring(base_field(parent_ring), n; cached=false)[1]
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
# floats
@register_serialization_type AbstractAlgebra.Floats{Float64} "Floats"

################################################################################
# field of rationals (singleton type)
@register_serialization_type QQField

################################################################################
# non-ZZRingElem variant
@register_serialization_type Nemo.fpField

function save_object(s::SerializerState, F::fpField)
  save_object(s, string(characteristic(F)))
end

function load_object(s::DeserializerState, ::Type{fpField})
  load_node(s) do str
    return fpField(parse(UInt64, str))
  end
end

# elements
@register_serialization_type fpFieldElem uses_params

function save_object(s::SerializerState, elem::fpFieldElem)
  save_data_basic(s, string(elem))
end

function load_object(s::DeserializerState, ::Type{fpFieldElem}, F::fpField)
  load_node(s) do str
    return F(parse(UInt64, str))
  end
end

################################################################################
# ZZRingElem variant
@register_serialization_type Nemo.FpField

function save_object(s::SerializerState, F::FpField)
  save_object(s, string(characteristic(F)))
end

function load_object(s::DeserializerState, ::Type{FpField})
  load_node(s) do str
    FpField(parse(ZZRingElem, str))
  end
end

# elements
@register_serialization_type FpFieldElem uses_params

function save_object(s::SerializerState, elem::FpFieldElem)
  save_data_basic(s, string(elem))
end

function load_object(s::DeserializerState, ::Type{FpFieldElem}, F::FpField)
  load_node(s) do str
    F(parse(ZZRingElem, str))
  end
end

################################################################################
# SimpleNumField

@register_serialization_type Hecke.RelSimpleNumField uses_id
@register_serialization_type AbsSimpleNumField uses_id

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
@register_serialization_type AbsSimpleNumFieldElem uses_params
@register_serialization_type Hecke.RelSimpleNumFieldElem uses_params
const NumFieldElemTypeUnion = Union{AbsSimpleNumFieldElem, fqPolyRepFieldElem, Hecke.RelSimpleNumFieldElem}

function save_object(s::SerializerState, k::NumFieldElemTypeUnion)
  K = parent(k)
  polynomial = parent(defining_polynomial(K))(k)
  save_object(s, polynomial)
end

function save_object(s::SerializerState, k::Hecke.RelSimpleNumFieldElem{AbsNonSimpleNumFieldElem})
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
      save_typed_object(s, defining_polynomial(K), :def_pol)
    end
  end
end

function load_object(s::DeserializerState, ::Type{<: FqField})
  load_node(s) do node
    if node isa String
      order = ZZRingElem(node)
      return finite_field(order)[1]
    else
      def_pol = load_typed_object(s, :def_pol)
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
  if absolute_degree(parent) != 1
    return load_object(s, FqFieldElem, get_parents(parent))
  end
  load_node(s) do str
    return parent(parse(ZZRingElem, str))
  end
end

################################################################################
# Non Simple Extension

@register_serialization_type Hecke.RelNonSimpleNumField uses_id
@register_serialization_type AbsNonSimpleNumField uses_id

function save_object(s::SerializerState, K::Union{AbsNonSimpleNumField, RelNonSimpleNumField})
  def_pols = defining_polynomials(K)
  save_data_dict(s) do
    save_typed_object(s, def_pols, :def_pols)
    save_object(s, vars(K), :vars)
  end
end

function load_object(s::DeserializerState, ::Type{<: Union{AbsNonSimpleNumField, RelNonSimpleNumField}})
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
@register_serialization_type Hecke.RelNonSimpleNumFieldElem uses_params
@register_serialization_type AbsNonSimpleNumFieldElem uses_params

function save_object(s::SerializerState, k::Union{AbsNonSimpleNumFieldElem, Hecke.RelNonSimpleNumFieldElem})
  polynomial = Oscar.Hecke.data(k)
  save_object(s, polynomial)
end

function load_object(s::DeserializerState, ::Type{<: Union{AbsNonSimpleNumFieldElem, Hecke.RelNonSimpleNumFieldElem}},
                     parents::Vector)
  K = parents[end]
  n = ngens(K)
  # forces parent of MPolyRingElem
  poly_ring = polynomial_ring(base_field(K), n; cached=false)
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
@register_serialization_type ArbFieldElem uses_params

function save_object(s::SerializerState, RR::ArbField)
  save_object(s, precision(RR))
end

function load_object(s::DeserializerState, ::Type{ArbField})
  prec = load_object(s, Int64, :precision)
  return ArbField(prec)
end

# elements
function save_object(s::SerializerState, r::ArbFieldElem)
  c_str = ccall((:arb_dump_str, Nemo.libflint), Ptr{UInt8}, (Ref{ArbFieldElem},), r)
  save_object(s, unsafe_string(c_str))

  # free memory
  ccall((:flint_free, Nemo.libflint), Nothing, (Ptr{UInt8},), c_str)
end

function load_object(s::DeserializerState, ::Type{ArbFieldElem}, parent::ArbField)
  r = ArbFieldElem()
  load_node(s) do str
    ccall((:arb_load_str, Nemo.libflint),
          Int32, (Ref{ArbFieldElem}, Ptr{UInt8}), r, str)
  end
  r.parent = parent
  return r
end

################################################################################
# AcbField
@register_serialization_type AcbField
@register_serialization_type AcbFieldElem uses_params

function save_object(s::SerializerState, CC::AcbField)
  save_object(s, precision(CC))
end

function load_object(s::DeserializerState, ::Type{AcbField})
  prec = load_object(s, Int)
  return AcbField(prec)
end

# elements
function save_object(s::SerializerState, c::AcbFieldElem)
  save_data_array(s) do
    save_object(s, real(c))
    save_object(s, imag(c))
  end
end

function load_object(s::DeserializerState, ::Type{AcbFieldElem}, parent::AcbField)
  (real_part, imag_part) = load_array_node(s) do _
    load_object(s, ArbFieldElem, ArbField(precision(parent)))
  end
  return parent(real_part, imag_part)
end

################################################################################
# Field Embeddings

const FieldEmbeddingTypes = Union{Hecke.AbsSimpleNumFieldEmbedding, Hecke.RelSimpleNumFieldEmbedding, Hecke.AbsNonSimpleNumFieldEmbedding, Hecke.RelNonSimpleNumFieldEmbedding}

@register_serialization_type Hecke.AbsNonSimpleNumFieldEmbedding uses_id
@register_serialization_type Hecke.RelNonSimpleNumFieldEmbedding uses_id
@register_serialization_type Hecke.AbsSimpleNumFieldEmbedding uses_id
@register_serialization_type Hecke.RelSimpleNumFieldEmbedding uses_id

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

@register_serialization_type EmbeddedNumField uses_id

function save_object(s::SerializerState, E::EmbeddedNumField)
  K = number_field(E)
  e = embedding(E)

  save_data_dict(s) do
    save_typed_object(s, K, :num_field)
    save_typed_object(s, e, :embedding)
  end
end

function load_object(s::DeserializerState, ::Type{EmbeddedNumField})
  K = load_typed_object(s, :num_field)
  e = load_typed_object(s, :embedding)

  return Hecke.embedded_field(K, e)[1]
end

@register_serialization_type EmbeddedNumFieldElem uses_params

function save_object(s::SerializerState, f::EmbeddedNumFieldElem)
  save_object(s, data(f))
end

function load_object(s::DeserializerState, ::Type{<:EmbeddedNumFieldElem}, parents::Vector)
  parent_field = parents[end]
  numfield_elem = terms
  coeff_type = elem_type(parents[end - 1])
  loaded_alg_elem = load_object(s, coeff_type, parents[1:end - 1])
  return parent_field(loaded_alg_elem)
end

################################################################################
# QQBar

@register_serialization_type QQBarField
@register_serialization_type QQBarFieldElem

function save_object(s::SerializerState, q::QQBarFieldElem)
  is_unique = false
  min_poly_q = minpoly(q)
  roots_min_q = roots(QQBarField(), min_poly_q)
  precision = 30
  approximation = undef

  while(!is_unique)
    CC = AcbField(precision; cached = false)
    approximation = CC(q)
    n_overlaps = length(filter(x -> overlaps(approximation, CC(x)), roots_min_q))
    if n_overlaps == 1
      is_unique = true
    else
      precision *= 2
    end
  end

  save_data_dict(s) do
    save_object(s, min_poly_q, :minpoly)
    save_object(s, approximation, :acb)
    save_object(s, precision, :precision)
  end
end

function load_object(s::DeserializerState, ::Type{QQBarFieldElem})
  Qx, x = polynomial_ring(QQ, :x; cached=false)
  min_poly = load_object(s, PolyRingElem{QQ}, Qx, :minpoly)
  precision = load_object(s, Int, :precision)
  CC = AcbField(precision; cached = false)
  approximation = load_object(s, AcbFieldElem, CC, :acb)
  roots_min_poly = roots(QQBarField(), min_poly)

  try
    only(filter(x -> overlaps(approximation, CC(x)), roots_min_poly))
  catch e
    if e isa ArgumentError
      error("The approximation is not precise enough to determine a unique root")
    end
    rethrow(e)
  end
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
