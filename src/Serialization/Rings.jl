################################################################################
# Common union types

# this will need a better name at some point
const RingMatElemUnion = Union{RingElem, MatElem, FreeAssAlgElem}

# this union will also need a better name at some point
const RingMatSpaceUnion = Union{Ring, MatSpace, FreeAssAlgebra}

################################################################################
# Utility functions for ring parent tree

# builds parent tree
function get_parents(parent_ring::T) where T <: RingMatSpaceUnion
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
# Handling RingElem MatElem, FieldElem ... Params

function save_type_params(s::SerializerState, x::T) where T <: RingMatElemUnion
  save_data_dict(s) do
    save_object(s, encode_type(T), :name)
    parent_x = parent(x)
    if serialize_with_id(parent_x)
      parent_ref = save_as_ref(s, parent_x)
      save_object(s, parent_ref, :params)
    else
      save_typed_object(s, parent_x, :params)
    end
  end
end

function load_type_params(s::DeserializerState, ::Type{<:RingMatElemUnion}, dict::Dict{Symbol, Any})
  return load_typed_object(s, dict)
end

# fix for polynomial cases
function load_object(s::DeserializerState, T::Type{<:RingMatElemUnion},
                     terms::Vector{Any}, parent_ring::RingMatSpaceUnion) 
  parents = get_parents(parent_ring)
  return load_object(s, T, terms, parents)
end

# fix for series and ideal cases
function load_object(s::DeserializerState, T::Type{<:Union{RingElem, MPolyIdeal, Laurent.LaurentMPolyIdeal}},
                     terms::Dict{Symbol, Any}, parent_ring::S) where S <: Union{Ring, AbstractAlgebra.Generic.LaurentMPolyWrapRing}
  parents = get_parents(parent_ring)
  return load_object(s, T, terms, parents)
end

################################################################################
# ring of integers (singleton type)
@register_serialization_type ZZRing

################################################################################
#  Mod Rings
@register_serialization_type Nemo.zzModRing
@register_serialization_type Nemo.ZZModRing
const ModRingUnion = Union{Nemo.zzModRing, Nemo.ZZModRing}

function save_object(s::SerializerState, R::T) where T <: ModRingUnion
  save_object(s, modulus(R))
end

function load_object(s::DeserializerState, ::Type{Nemo.zzModRing}, str::String)
  modulus = parse(UInt64, str)
  return Nemo.zzModRing(modulus)
end

function load_object(s::DeserializerState, ::Type{Nemo.ZZModRing}, str::String)
  modulus = ZZRingElem(str)
  return Nemo.ZZModRing(modulus)
end

#elements
@register_serialization_type zzModRingElem uses_params
@register_serialization_type ZZModRingElem uses_params
const ModRingElemUnion = Union{zzModRingElem, ZZModRingElem}

function save_object(s::SerializerState, x::ModRingElemUnion)
  save_data_basic(s, string(x))
end

function load_object(s::DeserializerState, ::Type{<:ModRingElemUnion},
                                 str::String, parent_ring::T) where T <: ModRingUnion
  return parent_ring(ZZRingElem(str))
end

################################################################################
#  Polynomial Rings

@register_serialization_type PolyRing uses_id
@register_serialization_type MPolyRing uses_id
@register_serialization_type UniversalPolyRing uses_id
@register_serialization_type AbstractAlgebra.Generic.LaurentMPolyWrapRing uses_id

function save_object(s::SerializerState, R::Union{UniversalPolyRing, MPolyRing, PolyRing, AbstractAlgebra.Generic.LaurentMPolyWrapRing})
  save_data_dict(s) do
    save_typed_object(s, base_ring(R), :base_ring)
    save_object(s, symbols(R), :symbols)
  end
end

function load_object(s::DeserializerState,
                     T::Type{<: Union{UniversalPolyRing, MPolyRing, PolyRing, AbstractAlgebra.Generic.LaurentMPolyWrapRing}},
                     dict::Dict)
  base_ring = load_typed_object(s, dict[:base_ring])
  symbols = map(Symbol, dict[:symbols])
  
  if T <: PolyRing
    return polynomial_ring(base_ring, symbols..., cached=false)[1]
  elseif T <: UniversalPolyRing
    poly_ring = UniversalPolynomialRing(base_ring, cached=false)
    gens(poly_ring, symbols)
    return poly_ring
  elseif T <: AbstractAlgebra.Generic.LaurentMPolyWrapRing
    return LaurentPolynomialRing(base_ring, symbols, cached=false)[1]
  end
  
  return polynomial_ring(base_ring, symbols, cached=false)[1]
end

################################################################################
#  Polynomial Ring Types
@register_serialization_type MPolyRingElem uses_params
@register_serialization_type UniversalPolyRingElem uses_params
@register_serialization_type AbstractAlgebra.Generic.LaurentMPolyWrap uses_params
const PolyElemUniontype = Union{MPolyRingElem, UniversalPolyRingElem, AbstractAlgebra.Generic.LaurentMPolyWrap}

# elements
function save_object(s::SerializerState, p::Union{UniversalPolyRingElem, MPolyRingElem})
  # we use this line instead of typeof(coeff(p, 1)) to catch the 0 polynomial
  coeff_type = elem_type(base_ring(parent(p)))
  save_data_array(s) do
    for i in 1:length(p)
      save_data_array(s) do 
        save_object(s, map(string, exponent_vector(p, i)))
        save_object(s, coeff(p, i))
      end
    end
  end
end

function save_object(s::SerializerState, p::AbstractAlgebra.Generic.LaurentMPolyWrap)
  exponent_vectors_gen = AbstractAlgebra.exponent_vectors(p)
  index = 0
  save_data_array(s) do
    for c in coefficients(p)
      exponent_vector, index = iterate(exponent_vectors_gen, index)
      save_data_array(s) do
        save_object(s, map(string, exponent_vector))
        save_object(s, c)
      end
    end
  end
end

################################################################################
# Univariate Polynomials

@register_serialization_type PolyRingElem uses_params

function save_object(s::SerializerState, p::PolyRingElem)
  coeffs = coefficients(p)
  exponent = 0
  save_data_array(s) do
    for coeff in coeffs
      # collect only non trivial terms
      if is_zero(coeff)
        exponent += 1
        continue
      end
      save_data_array(s) do
        save_object(s, string(exponent))
        save_object(s, coeff)
      end
      exponent += 1
    end
  end
end

function load_object(s::DeserializerState, ::Type{<: PolyRingElem},
                     terms::Vector, parents::Vector)
  parent_ring = parents[end]
  if isempty(terms)
    return parent_ring(0)
  end
  # load exponent
  terms = map(x->(parse(Int, x[1]), x[2]), terms)
  # shift so constant starts at 1
  degree = max(map(x->x[1] + 1, terms)...)
  base = base_ring(parent_ring)
  loaded_terms = zeros(base, degree)
  coeff_type = elem_type(base)

  for term in terms
    exponent, coeff = term
    # account for shift
    exponent += 1
    if serialize_with_params(coeff_type)
      if length(parents) == 1
        params = coefficient_ring(parent_ring)
      else
        params = parents[1:end - 1]
      end
      loaded_terms[exponent] = load_object(s, coeff_type, coeff, params)
    else
      loaded_terms[exponent] = load_object(s, coeff_type, coeff)
    end
  end
  return parent_ring(loaded_terms)
end


function load_object(s::DeserializerState,
                                 ::Type{<:Union{MPolyRingElem, UniversalPolyRingElem, AbstractAlgebra.Generic.LaurentMPolyWrap}},
                                 terms::Vector, parents::Vector)
  parent_ring = parents[end]
  base = base_ring(parent_ring)
  polynomial = MPolyBuildCtx(parent_ring)
  coeff_type = elem_type(base)
  for (e, coeff) in terms
    if serialize_with_params(coeff_type)
      if length(parents) == 1
        params = coefficient_ring(parent_ring)
      else
        params = parents[1:end - 1]
      end
      c = load_object(s, coeff_type, coeff, params)
    else
      c = load_object(s, coeff_type, coeff)
    end
    e_int = [parse(Int, x) for x in e]
    push_term!(polynomial, c, e_int)
  end
  return finish(polynomial)
end

################################################################################
# Polynomial Ideals

@register_serialization_type MPolyIdeal uses_params
@register_serialization_type Laurent.LaurentMPolyIdeal uses_params
const IdealUnionType = Union{MPolyIdeal, Laurent.LaurentMPolyIdeal, FreeAssAlgIdeal}

function save_type_params(s::SerializerState, x::T) where T <: IdealUnionType
  save_data_dict(s) do
    save_object(s, encode_type(T), :name)
    ref = save_as_ref(s, parent(gens(x)[1]))
    save_object(s, ref, :params)
  end
end

function load_type_params(s::DeserializerState, ::Type{<: IdealUnionType},
                          params::T) where T <: RingMatSpaceUnion
  return load_type_params(s, RingElem, params)
end

function save_object(s::SerializerState, I::T) where T <: IdealUnionType
  save_object(s, gens(I))
end

function load_object(s::DeserializerState, T::Type{<: IdealUnionType},
                     gen_terms::Vector, params::Vector)
  return load_object(s, T, gen_terms, params[end])
end

function load_object(s::DeserializerState, ::Type{<: IdealUnionType},
                     gen_terms::Vector, parent_ring::RingMatSpaceUnion)
  gens = [
    load_object(s, elem_type(parent_ring), g, parent_ring) for g in gen_terms
      ]
  return ideal(parent_ring, gens)
end

################################################################################
# Matrices
@register_serialization_type MatSpace uses_id
@register_serialization_type MatElem uses_params

function save_object(s::SerializerState, obj::MatSpace)
  save_data_dict(s) do
    save_typed_object(s, base_ring(obj), :base_ring)
    save_object(s, ncols(obj), :ncols)
    save_object(s, nrows(obj), :nrows)
  end
end

function load_object(s::DeserializerState, ::Type{<:MatSpace}, dict::Dict)
  base_ring = load_typed_object(s, dict[:base_ring])
  ncols = parse(Int, dict[:ncols])
  nrows = parse(Int, dict[:nrows])
  return matrix_space(base_ring, nrows, ncols)
end

function save_object(s::SerializerState, obj::MatElem)
  save_object(s, Array(obj))
end

function load_object(s::DeserializerState, ::Type{<:MatElem},
                     entries::Vector, parents::Vector)
  parent = parents[end]
  T = elem_type(base_ring(parent))
  if serialize_with_params(T)
    if length(parents) == 1
      params = base_ring(parent)
    else
      params = parents[1:end - 1]
    end
    m = load_object(s, Matrix, entries, (T, params))
  else
    m = load_object(s, Matrix, entries, T)
  end
  return parent(m)
end

################################################################################
# Power Series
@register_serialization_type SeriesRing uses_id

function save_object(s::SerializerState, R::Union{
  Generic.RelPowerSeriesRing,
  QQRelPowerSeriesRing,
  ZZRelPowerSeriesRing,
  fqPolyRepRelPowerSeriesRing,
  FqRelPowerSeriesRing,
  zzModRelPowerSeriesRing})
  save_data_dict(s) do
    save_typed_object(s, base_ring(R), :base_ring)
    save_object(s, var(R), :var)
    save_object(s, max_precision(R), :max_precision)
    save_object(s, :capped_relative, :model)
  end
end

function save_object(s::SerializerState, R::Union{
  Generic.AbsPowerSeriesRing,
  QQAbsPowerSeriesRing,
  ZZAbsPowerSeriesRing,
  FqAbsPowerSeriesRing,
  fqPolyRepAbsPowerSeriesRing,
  zzModAbsPowerSeriesRing})

  save_data_dict(s) do
    save_typed_object(s, base_ring(R), :base_ring)
    save_object(s, var(R), :var)
    save_object(s, max_precision(R), :max_precision)
    save_object(s, :capped_absolute, :model)
  end
end

function load_object(s::DeserializerState, ::Type{<: SeriesRing}, dict::Dict)
  base_ring = load_typed_object(s, dict[:base_ring])
  var = load_object(s, Symbol, dict[:var])
  max_precision = load_object(s, Int, dict[:max_precision])
  model = load_object(s, Symbol, dict[:model])
  
  return power_series_ring(base_ring, max_precision, var; cached=false, model=model)[1]
end

# elements
@register_serialization_type RelPowerSeriesRingElem uses_params
@register_serialization_type AbsPowerSeriesRingElem uses_params

function save_object(s::SerializerState, r::RelPowerSeriesRingElem)
  v = valuation(r)
  pl = pol_length(r)
  encoded_terms = []
  save_data_dict(s) do
    save_data_array(s, :terms) do
      for exponent in v: v + pl
        coefficient = coeff(r, exponent)
        #collect only non trivial values
        if is_zero(coefficient)
          continue
        end

        save_data_array(s) do
          save_object(s, exponent)
          save_object(s, coefficient)
        end
      end
    end
    save_object(s, pl, :pol_length)
    save_object(s, precision(r), :precision)
    save_object(s, v, :valuation)
  end
end

function save_object(s::SerializerState, r::AbsPowerSeriesRingElem)
  pl = pol_length(r)
  encoded_terms = []
  parents = []
  parent_ring = parent(r)
  save_data_dict(s) do
    save_data_array(s, :terms) do
      for exponent in 0:pl
        coefficient = coeff(r, exponent)
        #collect only non trivial values
        if is_zero(coefficient)
          continue
        end
        save_data_array(s) do
          save_object(s, exponent)
          save_object(s, coefficient)
        end
      end
    end
    save_object(s, pl, :pol_length)
    save_object(s, precision(r),:precision)
  end
end

function load_object(s::DeserializerState, ::Type{<: RelPowerSeriesRingElem},
                     dict::Dict, parents::Vector)
  parent_ring = parents[end]
  valuation = parse(Int, dict[:valuation])
  pol_length = parse(Int, dict[:pol_length])
  precision = parse(Int, dict[:precision])
  base = base_ring(parent_ring)
  loaded_terms = zeros(base, pol_length)
  coeff_type = elem_type(base)
  for (exponent, coeff) in dict[:terms]
    if serialize_with_params(coeff_type)
      if length(parents) == 1
        params = base
      else
        params = parents[1:end - 1]
      end
      c = load_object(s, coeff_type, coeff, params)
    else
      c = load_object(s, coeff_type, coeff)
    end
    e = parse(Int, exponent)
    loaded_terms[e] = c
  end
  
  return parent_ring(loaded_terms, pol_length, precision, valuation)
end

function load_object(s::DeserializerState, ::Type{<: AbsPowerSeriesRingElem},
                     dict::Dict, parents::Vector)
  parent_ring = parents[end]
  pol_length = parse(Int, dict[:pol_length])
  precision = parse(Int, dict[:precision])
  base = base_ring(parent_ring)
  loaded_terms = zeros(base, pol_length)
  coeff_type = elem_type(base)
  for (exponent, coeff) in dict[:terms]
    if serialize_with_params(coeff_type)
      if length(parents) == 1
        params = base
      else
        params = parents[1:end - 1]
      end
      c = load_object(s, coeff_type, coeff, params)
    else
      c = load_object(s, coeff_type, coeff)
    end
    e = parse(Int, exponent)
    e += 1
    loaded_terms[e] = c
  end
  
  return parent_ring(loaded_terms, pol_length, precision)
end

################################################################################
# Laurent Series
@register_serialization_type Generic.LaurentSeriesRing  "LaurentSeriesRing" uses_id
@register_serialization_type Generic.LaurentSeriesField "LaurentSeriesField" uses_id
@register_serialization_type ZZLaurentSeriesRing uses_id

function save_object(s::SerializerState, R::Union{
  Generic.LaurentSeriesRing,
  Generic.LaurentSeriesField,
  ZZLaurentSeriesRing})
  save_data_dict(s) do
    save_typed_object(s, base_ring(R), :base_ring)
    save_object(s, var(R), :var)
    save_object(s, max_precision(R), :max_precision)
  end
end

function load_object(s::DeserializerState,
                     ::Type{<: Union{
                       Generic.LaurentSeriesRing,
                       Generic.LaurentSeriesField,
                       ZZLaurentSeriesRing}},
                     dict::Dict)
  base_ring = load_typed_object(s, dict[:base_ring])
  var = Symbol(dict[:var])
  max_precision = parse(Int, dict[:max_precision])

  return laurent_series_ring(base_ring, max_precision, var; cached=false)[1]
end

# elements
@register_serialization_type Generic.LaurentSeriesFieldElem "LaurentSeriesFieldElem" uses_params
@register_serialization_type Generic.LaurentSeriesRingElem "LaurentSeriesRingElem" uses_params
@register_serialization_type ZZLaurentSeriesRingElem uses_params

function save_object(s::SerializerState, r:: Union{Generic.LaurentSeriesElem, ZZLaurentSeriesRingElem})
  v = valuation(r)
  pl = pol_length(r)
  encoded_terms = []
  save_data_dict(s) do
    save_data_array(s, :terms) do
      for exponent in v: v + pl
        coefficient = coeff(r, exponent)
        #collect only non trivial values
        if is_zero(coefficient)
          continue
        end

        save_data_array(s) do
          save_object(s, exponent)
          save_object(s, coefficient)
        end
      end
    end
    save_object(s, pl, :pol_length)
    save_object(s, precision(r), :precision)
    save_object(s, v, :valuation)
    save_object(s, Generic.scale(r), :scale)
  end
end

function load_object(s::DeserializerState, ::Type{<: Union{Generic.LaurentSeriesElem, ZZLaurentSeriesRingElem}},
                     dict::Dict, parents::Vector)
  parent_ring = parents[end]
  terms = dict[:terms]
  highest_degree = max(map(x->parse(Int, x[1]), terms)...)
  lowest_degree = min(map(x->parse(Int, x[1]), terms)...)
  base = base_ring(parent_ring)
  coeff_type = elem_type(base)
  # account for index shift
  loaded_terms = zeros(base, highest_degree - lowest_degree + 1)
  for (exponent, coeff) in terms
    e = parse(Int, exponent)
    e -= lowest_degree - 1
    if serialize_with_params(coeff_type)
      if length(parents) == 1
        params = base
      else
        params = parents[1:end - 1]
      end
      c = load_object(s, coeff_type, coeff, params)
    else
      c = load_object(s, coeff_type, coeff)
    end
    loaded_terms[e] = c
  end
  valuation = parse(Int, dict[:valuation])
  pol_length = parse(Int, dict[:pol_length])
  precision = parse(Int, dict[:precision])
  scale = parse(Int, dict[:scale])
  return parent_ring(loaded_terms, pol_length, precision, valuation, scale)
end
