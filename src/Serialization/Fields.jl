################################################################################
# floats
@register_serialization_type AbstractAlgebra.Floats{Float64} "Floats"

################################################################################
# field of rationals (singleton type)
@register_serialization_type QQField

################################################################################
# type_params for field extension types

################################################################################
# non-ZZRingElem variant
@register_serialization_type fpField "FiniteField"

function save_object(s::SerializerState, F::fpField)
  save_object(s, string(characteristic(F)))
end

function load_object(s::DeserializerState, ::Type{fpField})
  load_node(s) do str
    return fpField(parse(UInt64, str))
  end
end

# elements
@register_serialization_type fpFieldElem 

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
@register_serialization_type FpField "FiniteField"

function save_object(s::SerializerState, F::FpField)
  save_object(s, string(characteristic(F)))
end

function load_object(s::DeserializerState, ::Type{FpField})
  load_node(s) do str
    FpField(parse(ZZRingElem, str))
  end
end

# elements
@register_serialization_type FpFieldElem

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
@register_serialization_type AbsSimpleNumField uses_id [:cyclo]
const SimNumFieldTypeUnion = Union{AbsSimpleNumField, Hecke.RelSimpleNumField}

type_params(obj::T) where T <: SimpleNumField = TypeParams(T, parent(defining_polynomial(obj)))

function save_object(s::SerializerState, K::SimpleNumField)
  save_data_dict(s) do
    save_object(s, defining_polynomial(K), :def_pol)
    save_object(s, var(K), :var)
  end
end

function load_object(s::DeserializerState, ::Type{<: SimpleNumField}, params::PolyRing)
  var = load_object(s, Symbol, :var)
  def_pol = load_object(s, PolyRingElem, params, :def_pol)
  K, _ = number_field(def_pol, var, cached=false)
  return K
end

################################################################################
# FqNmodfinitefield
@register_serialization_type fqPolyRepField uses_id

type_params(K::fqPolyRepField) = TypeParams(fqPolyRepField, parent(defining_polynomial(K)))

function save_object(s::SerializerState, K::fqPolyRepField)
  save_object(s, defining_polynomial(K))
end

function load_object(s::DeserializerState, ::Type{<: fqPolyRepField}, params::PolyRing)
  def_pol = load_object(s, PolyRingElem, params)
  K, _ = Nemo.Native.finite_field(def_pol, cached=false)
  return K
end

#elements
@register_serialization_type fqPolyRepFieldElem
@register_serialization_type AbsSimpleNumFieldElem
@register_serialization_type Hecke.RelSimpleNumFieldElem
const SimNumFieldElemTypeUnion = Union{AbsSimpleNumFieldElem, fqPolyRepFieldElem, Hecke.RelSimpleNumFieldElem}

function save_object(s::SerializerState, k::SimNumFieldElemTypeUnion)
  K = parent(k)
  polynomial = parent(defining_polynomial(K))(k)
  save_object(s, polynomial)
end

function save_object(s::SerializerState, k::Hecke.RelSimpleNumFieldElem{AbsNonSimpleNumFieldElem})
  K = parent(k)
  polynomial = parent(defining_polynomial(K))(data(k))
  save_object(s, polynomial)
end

function load_object(s::DeserializerState, ::Type{<: SimNumFieldElemTypeUnion}, K::Union{SimNumFieldTypeUnion, fqPolyRepField})
  polynomial = load_object(s, PolyRingElem, parent(defining_polynomial(K)))
  loaded_terms = evaluate(polynomial, gen(K))
  return K(loaded_terms)
end

################################################################################
# FqField

@register_serialization_type FqField "FiniteField" uses_id
@register_serialization_type FqFieldElem

function type_params(K::FqField)
  absolute_degree(K) == 1 && return TypeParams(FqField, nothing)
  return TypeParams(FqField, parent(defining_polynomial(K)))
end

function save_object(s::SerializerState, K::FqField)
  if absolute_degree(K) == 1
    save_object(s, order(K))
  else
    save_object(s, defining_polynomial(K))
  end
end

function load_object(s::DeserializerState, ::Type{<: FqField}, params::PolyRing)
  finite_field(load_object(s, PolyRingElem, params), cached=false)[1]
end

function load_object(s::DeserializerState, ::Type{<: FqField})
  finite_field(load_object(s, ZZRingElem, ZZRing()))[1]
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

function load_object(s::DeserializerState, ::Type{<: FqFieldElem}, K::FqField)
  load_node(s) do _
    if absolute_degree(K) != 1
      return K(load_object(s, PolyRingElem, parent(defining_polynomial(K))))
    end
    K(load_object(s, ZZRingElem, ZZRing()))
  end
end

################################################################################
# Non Simple Extension

@register_serialization_type Hecke.RelNonSimpleNumField uses_id
@register_serialization_type AbsNonSimpleNumField uses_id

function type_params(K::T) where T <: Union{AbsNonSimpleNumField, RelNonSimpleNumField}
  return TypeParams(T, parent(defining_polynomials(K)[1]))
end

function save_object(s::SerializerState, K::NonSimpleNumField)
  save_data_dict(s) do
    save_object(s, defining_polynomials(K), :def_pols)
    save_object(s, vars(K), :vars)
  end
end

function load_object(s::DeserializerState,
                     ::Type{<: NonSimpleNumField},
                     params::PolyRing)
  def_pols = load_object(s, Vector{PolyRingElem}, params, :def_pols)
  vars = load_node(s, :vars) do vars_data
    return map(Symbol, vars_data)
  end
  # fix since numberfield doesn't accept PolyRingElem vectors
  array_pols = Array{typeof(def_pols[1]), 1}(def_pols)
  K, _ = number_field(array_pols, vars, cached=false)
  return K
end

#elements
@register_serialization_type Hecke.RelNonSimpleNumFieldElem
@register_serialization_type AbsNonSimpleNumFieldElem

function save_object(s::SerializerState, k::Union{AbsNonSimpleNumFieldElem, Hecke.RelNonSimpleNumFieldElem})
  polynomial = Oscar.Hecke.data(k)
  save_object(s, polynomial)
end

function load_object(s::DeserializerState,
                     ::Type{<: Union{AbsNonSimpleNumFieldElem, Hecke.RelNonSimpleNumFieldElem}},
                     params::Union{AbsNonSimpleNumField, RelNonSimpleNumField})
  K = params
  n = ngens(K)
  # forces parent of MPolyRingElem
  poly_ring, _ = polynomial_ring(base_field(K), n; cached=false)
  poly_elem_type = elem_type
  load_node(s) do _
    polynomial = load_object(s, MPolyRingElem, poly_ring)
  end
  polynomial = evaluate(polynomial, gens(K))
  return K(polynomial)
end

################################################################################
# FracField

@register_serialization_type FracField uses_id

type_params(R::T) where T <: FracField = TypeParams(T, base_ring(R))

const FracUnionTypes = Union{MPolyRingElem, PolyRingElem, UniversalPolyRingElem}
# we use the union to prevent QQField from using these save methods

function save_object(s::SerializerState, K::FracField{T}) where T <: FracUnionTypes
  save_data_dict(s) do
    save_object(s, base_ring(K), :base_ring)
  end
end

load_object(s::DeserializerState, ::Type{<: FracField}, base::Ring) = fraction_field(base, cached=false)

# elements

@register_serialization_type FracElem{<: FracUnionTypes} "FracElem"

function save_object(s::SerializerState, f::FracElem)
  save_data_array(s) do
    save_object(s, numerator(f))
    save_object(s, denominator(f))
  end
end

function load_object(s::DeserializerState, ::Type{<: FracElem}, parent_ring::Ring)
  load_node(s) do _
    base = base_ring(parent_ring)
    coeff_type = elem_type(base)
    return parent_ring(
      load_object(s, coeff_type, base, 1),
      load_object(s, coeff_type, base, 2)
    )
  end
end

################################################################################
# RationalFunctionField

@register_serialization_type AbstractAlgebra.Generic.RationalFunctionField "RationalFunctionField" uses_id

type_params(R::T) where T <: AbstractAlgebra.Generic.RationalFunctionField = TypeParams(T, base_ring(R))

function save_object(s::SerializerState,
                     RF::AbstractAlgebra.Generic.RationalFunctionField{<: FieldElem, <: MPolyRingElem})
  save_data_dict(s) do
    syms = symbols(RF)
    save_object(s, syms, :symbols)
  end
end

function save_object(s::SerializerState,
                     RF::AbstractAlgebra.Generic.RationalFunctionField{<:FieldElem, <: PolyRingElem})
  save_data_dict(s) do
    save_object(s, only(symbols(RF)), :symbol)
  end
end

function load_object(s::DeserializerState,
                     ::Type{<: AbstractAlgebra.Generic.RationalFunctionField}, R::Ring)
  haskey(s, :symbol) && return rational_function_field(R, load_object(s, Symbol, :symbol), cached=false)[1]

  return rational_function_field(R, load_object(s, Vector{Symbol}, :symbols), cached=false)[1]
end

#elements
@register_serialization_type AbstractAlgebra.Generic.RationalFunctionFieldElem "RationalFunctionFieldElem"

function save_object(s::SerializerState, f::AbstractAlgebra.Generic.RationalFunctionFieldElem)
  save_data_array(s) do
    save_object(s, numerator(f))
    save_object(s, denominator(f))
  end
end

function load_object(s::DeserializerState,
                     ::Type{<: AbstractAlgebra.Generic.RationalFunctionFieldElem},
                     parent_ring::AbstractAlgebra.Generic.RationalFunctionField)
  base = base_ring(parent_ring.fraction_field)
  coeff_type = elem_type(base)

  return load_node(s) do _
    loaded_num = load_node(s, 1) do _
      load_object(s, coeff_type, base)
    end

    loaded_den = load_node(s, 2) do _
      load_object(s, coeff_type, base)
    end
    parent_ring(loaded_num, loaded_den)
  end
end

################################################################################
# ArbField
@register_serialization_type ArbField
@register_serialization_type ArbFieldElem

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
@register_serialization_type AcbFieldElem

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

const FieldEmbeddingTypes = Union{
  Hecke.AbsSimpleNumFieldEmbedding,
  Hecke.RelSimpleNumFieldEmbedding,
  Hecke.AbsNonSimpleNumFieldEmbedding,
  Hecke.RelNonSimpleNumFieldEmbedding
}

@register_serialization_type Hecke.AbsNonSimpleNumFieldEmbedding uses_id 
@register_serialization_type Hecke.RelNonSimpleNumFieldEmbedding uses_id
@register_serialization_type Hecke.AbsSimpleNumFieldEmbedding uses_id
@register_serialization_type Hecke.RelSimpleNumFieldEmbedding uses_id

function type_params(E::T) where T <: FieldEmbeddingTypes
  K = number_field(E)
  base_K = base_field(K)
  base_field(K) isa QQField && return TypeParams(T, K)

  base_field_emb = restrict(E, base_K)
  return TypeParams(
    T,
    :num_field => K,
    :base_field_emb => base_field_emb,
  )
end

function save_object(s::SerializerState, E::FieldEmbeddingTypes)
  K = number_field(E)
  base_K = base_field(K)
  if is_simple(K)
    a = gen(K)
    gen_emb_approx = E(a)
    if any(overlaps(gen_emb_approx, e(a)) for e in complex_embeddings(K) if e != E && restrict(E, base_K) == restrict(e, base_K))
      error("Internal error in internal serialization.")
    end
    save_object(s, gen_emb_approx)
  else
    a = gens(K)
    data = E.(a)
    if any(all(overlaps(t[1], t[2]) for t in zip(data, e.(a))) for e in complex_embeddings(K) if e != E && restrict(E, base_field(K)) == restrict(e, base_field(K)))
      error("Internal error in internal serialization.")
    end
    save_object(s, tuple(data...))
  end
end

function load_object(s::DeserializerState, T::Type{<:FieldEmbeddingTypes}, K::Field)
  if !is_simple(K)
    data = load_object(s, Vector{AcbFieldElem}, AcbField())
  else
    data = load_object(s, AcbFieldElem, AcbField())
  end
  return complex_embedding(K, data)
end

function load_object(s::DeserializerState, T::Type{<:FieldEmbeddingTypes}, params::Dict)
  K = params[:num_field]
  if !is_simple(K)
    data = load_object(s, Vector{AcbFieldElem}, AcbField())
  else
    data = load_object(s, AcbFieldElem, AcbField())
  end
  return complex_embedding(K, params[:base_field_emb], data)
end

@register_serialization_type EmbeddedNumField uses_id

type_params(E::T) where T <: EmbeddedNumField = TypeParams(T, embedding(E))

function save_object(s::SerializerState, E::EmbeddedNumField)
  save_data_array(s) do
    # no data required but we leave this function here to generate a valid json
    # probably be a neater to way to do this
  end
end

function load_object(s::DeserializerState, ::Type{<:EmbeddedNumField},
                     E::T) where T <: FieldEmbeddingTypes
  K = number_field(E)
  return Hecke.embedded_field(K, E)[1]
end

@register_serialization_type EmbeddedNumFieldElem

function save_object(s::SerializerState, f::EmbeddedNumFieldElem)
  save_object(s, data(f))
end

function load_object(s::DeserializerState, ::Type{<:EmbeddedNumFieldElem}, E::T) where T <: EmbeddedNumField
  K = number_field(E)
  coeff_type = elem_type(K)
  loaded_alg_elem = load_object(s, coeff_type, K)
  return E(loaded_alg_elem)
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

function load_object(s::DeserializerState, ::Type{QQBarFieldElem}, ::QQBarField)
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
@register_serialization_type PadicFieldElem

function save_object(s::SerializerState, obj::PadicFieldElem)
  # currently it seems PadicFieldElems do not store the underlying polynomial
  save_object(s, lift(QQ, obj))
end

function load_object(s::DeserializerState, ::Type{PadicFieldElem}, parent_field::PadicField)
  rational_rep = load_object(s, QQFieldElem, QQField())
  return parent_field(rational_rep)
end
