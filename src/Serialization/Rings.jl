################################################################################
# Utility functions for ring parent tree
# this union will also need a better name at some point
RingMatSpaceUnion = Union{Ring, MatSpace, FreeAssAlgebra}

# builds parent tree
function get_parents(parent_ring::T) where T <: RingMatSpaceUnion
  # with new structure it seems like we may be able to remove the
  # array of parents encoded, or at least it seems like it may
  # make the code a bit cleaner, which will probably allow for
  # this line, and the function it uses below to be removed.
  # However i will leave this to a later PR.
  if has_elem_basic_encoding(parent_ring)
    return Any[]
  end
  base = base_ring(parent_ring)

  parents = get_parents(base)
  push!(parents, parent_ring)
  return parents
end

function save_parents(s::SerializerState, parent_ring::T) where T <: RingMatSpaceUnion
  parents = get_parents(parent_ring)
  refs = []
  for p in parents
    push!(refs, save_as_ref(s, p))
  end
  return refs
end

################################################################################
# Handling RingElem Params
# this will need a better name at some point
RingMatElemUnion = Union{RingElem, MatElem, FreeAssAlgElem}

function save_type_params(s::SerializerState, x::T) where T <: RingMatElemUnion
  data_dict(s) do
    save_object(s, encode_type(T), :name)
    parent_x = parent(x)
    if serialize_with_id(parent_x)
      parent_refs = save_parents(s, parent_x)
      save_object(s, parent_refs, :params)
    else
      save_typed_object(s, parent_x, :params)
    end
  end
end

function load_type_params(s::DeserializerState, ::Type{<:RingMatElemUnion}, dict::Dict{Symbol, Any})
  return load_typed_object(s, dict)
end

function load_type_params(s::DeserializerState, ::Type{<:RingMatElemUnion}, refs::Vector{Any})
  return load_parents(s, refs)
end

function load_type_params(s::DeserializerState, ::Type{<:RingMatElemUnion}, parent_ring::T)  where T <: RingMatSpaceUnion
  return get_parents(parent_ring)
end

# this should be properly dealt with later and is only an intermediate solution
# to getting the tests to pass when forcing a type and passing params
# ideally all load should look like this without passing a vector of parents

# fix for polynomial cases
function load_object(s::DeserializerState, T::Type{<:RingMatElemUnion}, terms::Vector{Any}, parent_ring::RingMatSpaceUnion) 
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
@registerSerializationType(ZZRing)

################################################################################
#  non simpleton base rings
@registerSerializationType(Nemo.zzModRing, "Nemo.zzModRing")
has_elem_basic_encoding(obj::Nemo.zzModRing) = true

function save_object(s::SerializerState, R::Nemo.zzModRing)
  save_object(s, string(modulus(R)))
end

function load_object(s::DeserializerState, ::Type{Nemo.zzModRing}, str::String)
  modulus = parse(UInt64, str)
  return Nemo.zzModRing(modulus)
end

#elements
@registerSerializationType(zzModRingElem)
type_needs_params(T::Type{zzModRingElem}) = true

function save_object(s::SerializerState, x::zzModRingElem)
  data_basic(s, string(x))
end

function load_object(s::DeserializerState, ::Type{zzModRingElem},
                                 str::String, parent_ring::Nemo.zzModRing)
  return parent_ring(ZZRingElem(str))
end

################################################################################
#  Polynomial Rings

@registerSerializationType(PolyRing, true)
@registerSerializationType(MPolyRing, true)
@registerSerializationType(UniversalPolyRing, true)
@registerSerializationType(AbstractAlgebra.Generic.LaurentMPolyWrapRing, true)

function save_object(s::SerializerState, R::Union{UniversalPolyRing, MPolyRing, PolyRing, AbstractAlgebra.Generic.LaurentMPolyWrapRing})
  data_dict(s) do
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
@registerSerializationType(MPolyRingElem)
@registerSerializationType(UniversalPolyRingElem)
@registerSerializationType(AbstractAlgebra.Generic.LaurentMPolyWrap)
PolyElemUniontype = Union{MPolyRingElem, UniversalPolyRingElem, AbstractAlgebra.Generic.LaurentMPolyWrap}
type_needs_params(::Type{<:PolyElemUniontype}) = true

# elements
function save_object(s::SerializerState, p::Union{UniversalPolyRingElem, MPolyRingElem})
  # we use this line instead of typeof(coeff(p, 1)) to catch the 0 polynomial
  coeff_type = elem_type(base_ring(parent(p)))
  data_array(s) do
    for i in 1:length(p)
      data_array(s) do 
        save_object(s, map(string, exponent_vector(p, i)))
        save_object(s, coeff(p, i))
      end
    end
  end
end

function save_object(s::SerializerState, p::AbstractAlgebra.Generic.LaurentMPolyWrap)
  exponent_vectors_gen = AbstractAlgebra.exponent_vectors(p)
  index = 0
  data_array(s) do
    for c in coefficients(p)
      exponent_vector, index = iterate(exponent_vectors_gen, index)
      data_array(s) do
        save_object(s, map(string, exponent_vector))
        save_object(s, c)
      end
    end
  end
end

################################################################################
# Univariate Polynomials

@registerSerializationType(PolyRingElem)
type_needs_params(::Type{<:PolyRingElem}) = true

function save_object(s::SerializerState, p::PolyRingElem)
  coeffs = coefficients(p)
  exponent = 0
  data_array(s) do
    for coeff in coeffs
      # collect only non trivial terms
      if is_zero(coeff)
        exponent += 1
        continue
      end
      data_array(s) do
        save_object(s, string(exponent))
        save_object(s, coeff)
      end
      exponent += 1
    end
  end
end

function load_object(s::DeserializerState,
                                 ::Type{<: PolyRingElem},
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
    if type_needs_params(coeff_type)
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
    if type_needs_params(coeff_type)
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

@registerSerializationType(MPolyIdeal)
@registerSerializationType(Laurent.LaurentMPolyIdeal)
const IdealUnionType = Union{MPolyIdeal, Laurent.LaurentMPolyIdeal, FreeAssAlgIdeal}
type_needs_params(::Type{<: IdealUnionType}) = true

function save_type_params(s::SerializerState, x::T) where T <: IdealUnionType
  data_dict(s) do
    save_object(s, encode_type(T), :name)
    refs = save_parents(s, parent(gens(x)[1]))
    save_object(s, refs, :params)
  end
end

function load_type_params(s::DeserializerState, ::Type{<: IdealUnionType}, params::Any)
  return load_type_params(s, RingElem, params)
end

function save_object(s::SerializerState, I::T) where T <: IdealUnionType
  data_dict(s) do
    save_object(s, gens(I), :gens)
  end
end

function load_object(s::DeserializerState, ::Type{<: IdealUnionType},
                                 dict::Dict{Symbol, Any}, params::Vector)
  parent_ring = params[end]
  gens = [
    load_object(s, elem_type(parent_ring), g, params) for g in dict[:gens]
      ]
  return ideal(parent_ring, gens)
end

################################################################################
# Matrices
@registerSerializationType(MatSpace, true)
@registerSerializationType(MatElem)
type_needs_params(::Type{<:MatElem}) = true

# not all 
function save_object(s::SerializerState, obj::MatSpace)
  data_dict(s) do
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
  if type_needs_params(T)
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
@registerSerializationType(SeriesRing, true)

function save_object(s::SerializerState, R::Union{
  Generic.RelPowerSeriesRing,
  QQRelPowerSeriesRing,
  ZZRelPowerSeriesRing,
  fqPolyRepRelPowerSeriesRing,
  FqRelPowerSeriesRing,
  zzModRelPowerSeriesRing})
  data_dict(s) do
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

  data_dict(s) do
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
@registerSerializationType(RelPowerSeriesRingElem)
@registerSerializationType(AbsPowerSeriesRingElem)
type_needs_params(::Type{<: Union{RelPowerSeriesRingElem, AbsPowerSeriesRingElem}}) = true

function save_object(s::SerializerState, r::RelPowerSeriesRingElem)
  v = valuation(r)
  pl = pol_length(r)
  encoded_terms = []
  data_dict(s) do
    s.key = :terms
    data_array(s) do
      for exponent in v: v + pl
        coefficient = coeff(r, exponent)
        #collect only non trivial values
        if is_zero(coefficient)
          continue
        end

        data_array(s) do
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
  data_dict(s) do
    s.key = :terms
    data_array(s) do
      for exponent in 0:pl
        coefficient = coeff(r, exponent)
        #collect only non trivial values
        if is_zero(coefficient)
          continue
        end
        data_array(s) do
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
    if type_needs_params(coeff_type)
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
    if type_needs_params(coeff_type)
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
@registerSerializationType(Generic.LaurentSeriesRing, true, "LaurentSeriesRing")
@registerSerializationType(Generic.LaurentSeriesField, true, "LaurentSeriesField")
@registerSerializationType(ZZLaurentSeriesRing)

function save_object(s::SerializerState, R::Union{
  Generic.LaurentSeriesRing,
  Generic.LaurentSeriesField,
  ZZLaurentSeriesRing})
  data_dict(s) do
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
@registerSerializationType(Generic.LaurentSeriesFieldElem, "LaurentSeriesFieldElem")
@registerSerializationType(Generic.LaurentSeriesRingElem, "LaurentSeriesRingElem")
@registerSerializationType(ZZLaurentSeriesRingElem)
type_needs_params(::Type{<: Union{ZZLaurentSeriesRingElem,
                                  Generic.LaurentSeriesFieldElem,
                                  Generic.LaurentSeriesRingElem}}) = true

function save_object(s::SerializerState, r:: Union{Generic.LaurentSeriesElem, ZZLaurentSeriesRingElem})
  v = valuation(r)
  pl = pol_length(r)
  encoded_terms = []
  data_dict(s) do
    s.key = :terms
    data_array(s) do
      for exponent in v: v + pl
        coefficient = coeff(r, exponent)
        #collect only non trivial values
        if is_zero(coefficient)
          continue
        end

        data_array(s) do
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

function load_object(s::DeserializerState,
                                 ::Type{<: Union{Generic.LaurentSeriesElem, ZZLaurentSeriesRingElem}},
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
    if type_needs_params(coeff_type)
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
