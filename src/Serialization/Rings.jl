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

function load_type_params(s::DeserializerState, ::Type{<:RingMatElemUnion})
  return load_typed_object(s)
end

# fix for polynomial cases
function load_object(s::DeserializerState, T::Type{<:RingMatElemUnion}, parent_ring::RingMatSpaceUnion) 
  parents = get_parents(parent_ring)
  return load_object(s, T, parents)
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

function load_object(s::DeserializerState, ::Type{Nemo.zzModRing})
  modulus = load_object(s, UInt64)
  return Nemo.zzModRing(modulus)
end

function load_object(s::DeserializerState, ::Type{Nemo.ZZModRing})
  modulus = load_object(s, ZZRingElem)
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
                     parent_ring::T) where T <: ModRingUnion
  return parent_ring(load_object(s, ZZRingElem))
end

################################################################################
#  Polynomial Rings

@register_serialization_type PolyRing uses_id
@register_serialization_type MPolyRing uses_id
@register_serialization_type UniversalPolyRing uses_id
@register_serialization_type MPolyDecRing uses_id
@register_serialization_type AbstractAlgebra.Generic.LaurentMPolyWrapRing uses_id

function save_object(s::SerializerState, R::Union{UniversalPolyRing, MPolyRing, PolyRing, AbstractAlgebra.Generic.LaurentMPolyWrapRing})
  save_data_dict(s) do
    save_typed_object(s, base_ring(R), :base_ring)
    save_object(s, symbols(R), :symbols)
  end
end

function load_object(s::DeserializerState,
                     T::Type{<: Union{UniversalPolyRing, MPolyRing, PolyRing, AbstractAlgebra.Generic.LaurentMPolyWrapRing}})
  base_ring = load_typed_object(s, :base_ring)
  symbols = load_object(s, Vector, Symbol, :symbols)

  if T <: PolyRing
    return polynomial_ring(base_ring, symbols..., cached=false)[1]
  elseif T <: UniversalPolyRing
    poly_ring = UniversalPolynomialRing(base_ring, cached=false)
    gens(poly_ring, symbols)
    return poly_ring
  elseif T <: AbstractAlgebra.Generic.LaurentMPolyWrapRing
    return laurent_polynomial_ring(base_ring, symbols, cached=false)[1]
  end
  return polynomial_ring(base_ring, symbols, cached=false)[1]
end

# with grading

function save_object(s::SerializerState, R::MPolyDecRing)
  save_data_dict(s) do
    save_typed_object(s, _grading(R), :grading)
    save_typed_object(s, forget_grading(R), :ring)
  end
end

function load_object(s::DeserializerState, ::Type{<:MPolyDecRing})
  ring = load_typed_object(s, :ring)
  grading = load_typed_object(s, :grading)
  return grade(ring, grading)[1]
end

################################################################################
#  Polynomial Ring Elem Types
@register_serialization_type MPolyRingElem uses_params
@register_serialization_type MPolyDecRingElem uses_params
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

function load_object(s::DeserializerState, ::Type{<: PolyRingElem}, parents::Vector)
  parent_ring = parents[end]
  load_node(s) do terms
    if isempty(terms)
      return parent_ring(0)
    end
    # load exponents and account for shift
    exponents = []
    for i in 1:length(terms)
      e = load_node(s, i) do _
        load_object(s, Int, 1) + 1
      end
      push!(exponents, e)
    end
    degree = max(exponents...)
    base = base_ring(parent_ring)
    loaded_terms = zeros(base, degree)
    coeff_type = elem_type(base)

    for (i, exponent) in enumerate(exponents)
      load_node(s, i) do term
        if serialize_with_params(coeff_type)
          if length(parents) == 1
            params = coefficient_ring(parent_ring)
          else
            params = parents[1:end - 1]
          end
          # place coefficient at s.obj
          load_node(s, 2) do _
            loaded_terms[exponent] = load_object(s, coeff_type, params)
          end
        else
          load_node(s, 2) do _
            loaded_terms[exponent] = load_object(s, coeff_type)
          end
        end
      end
    end
    return parent_ring(loaded_terms)
  end
end


function load_object(s::DeserializerState,
                     ::Type{<:Union{MPolyRingElem, UniversalPolyRingElem, AbstractAlgebra.Generic.LaurentMPolyWrap}},
                     parents::Vector)
  load_node(s) do terms
    exponents = [term[1] for term in terms]
    parent_ring = parents[end]
    base = base_ring(parent_ring)
    polynomial = MPolyBuildCtx(parent_ring)
    coeff_type = elem_type(base)

    for (i, e) in enumerate(exponents)
      load_node(s, i) do _
        c = nothing
        if serialize_with_params(coeff_type)
          if length(parents) == 1
            params = coefficient_ring(parent_ring)
          else
            params = parents[1:end - 1]
          end
          load_node(s, 2) do _
            c = load_object(s, coeff_type, params)
          end
        else
          load_node(s, 2) do _
            c = load_object(s, coeff_type)
          end
        end
        e_int = [parse(Int, x) for x in e]
        push_term!(polynomial, c, e_int)
      end
    end
    return finish(polynomial)
  end
end

function load_object(s::DeserializerState, ::Type{<:MPolyDecRingElem}, parents::Vector)
  parent_ring = parents[end]
  new_parents = push!(parents[1:end - 1], forget_grading(parent_ring))
  poly = load_object(s, MPolyRingElem, new_parents)
  return parent_ring(poly)
end


################################################################################
# Polynomial Ideals

@register_serialization_type MPolyIdeal uses_params
@register_serialization_type LaurentMPolyIdeal uses_params

# we should avoid this list getting too long and find a
# way to abstract saving params soon
const IdealOrdUnionType = Union{MPolyIdeal,
                                LaurentMPolyIdeal,
                                FreeAssAlgIdeal,
                                IdealGens,
                                MonomialOrdering}

function save_type_params(s::SerializerState, x::T) where T <: IdealOrdUnionType
  save_data_dict(s) do
    save_object(s, encode_type(T), :name)
    ref = save_as_ref(s, base_ring(x))
    save_object(s, ref, :params)
  end
end

function load_type_params(s::DeserializerState, ::Type{<: IdealOrdUnionType})
  return load_type_params(s, RingElem)
end

function save_object(s::SerializerState, I::T) where T <: IdealOrdUnionType
  save_object(s, gens(I))
end

function load_object(s::DeserializerState, ::Type{<: IdealOrdUnionType}, parent_ring::RingMatSpaceUnion)
  gens = elem_type(parent_ring)[]
  load_node(s) do gens_data
    for i in 1:length(gens_data)
      gen = load_node(s, i) do _
        load_object(s, elem_type(parent_ring), parent_ring)
      end
      push!(gens, gen)
    end
  end
  return ideal(parent_ring, gens)
end

################################################################################
# IdealGens

# this will need adjustments to cover the NCRing case

@register_serialization_type IdealGens uses_params

function save_object(s::SerializerState, obj::IdealGens)
  save_data_dict(s) do
    save_object(s, ordering(obj), :ordering)
    save_object(s, gens(obj), :gens)
    save_object(s, is_groebner_basis(obj), :is_gb)
    save_object(s, obj.isReduced, :is_reduced)
    save_object(s, obj.keep_ordering, :keep_ordering)
  end
end

function load_object(s::DeserializerState, ::Type{<:IdealGens}, base_ring::MPolyRing)
  ord = load_object(s, MonomialOrdering, base_ring, :ordering)
  generators = load_object(s, Vector{MPolyRingElem}, base_ring, :gens)
  is_gb = load_object(s, Bool, :is_gb)
  is_reduced = load_object(s, Bool, :is_reduced)
  keep_ordering = load_object(s, Bool, :keep_ordering)
  return IdealGens(base_ring, generators, ord;
                   keep_ordering=keep_ordering,
                   isReduced=is_reduced,
                   isGB=is_gb)
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

function load_object(s::DeserializerState, ::Type{<:MatSpace})
  base_ring = load_typed_object(s, :base_ring)
  ncols = load_object(s, Int, :ncols)
  nrows = load_object(s, Int, :nrows)
  return matrix_space(base_ring, nrows, ncols)
end

function save_object(s::SerializerState, obj::MatElem)
  save_object(s, Array(obj))
end

function load_object(s::DeserializerState, ::Type{<:MatElem}, parents::Vector)
  parent = parents[end]
  T = elem_type(base_ring(parent))
  if serialize_with_params(T)
    if length(parents) == 1
      params = base_ring(parent)
    else
      params = parents[1:end - 1]
    end
    m = load_object(s, Matrix, (T, params))
  else
    m = load_object(s, Matrix, T)
  end
  if isempty(m)
    return parent()
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

function load_object(s::DeserializerState, ::Type{<: SeriesRing})
  base_ring = load_typed_object(s, :base_ring)
  var = load_object(s, Symbol, :var)
  max_precision = load_object(s, Int, :max_precision)
  model = load_object(s, Symbol, :model)
  
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

function load_object(s::DeserializerState, ::Type{<:RelPowerSeriesRingElem}, parents::Vector)
  parent_ring = parents[end]
  valuation = load_object(s, Int, :valuation)
  pol_length = load_object(s, Int, :pol_length)
  precision = load_object(s, Int, :precision)
  base = base_ring(parent_ring)
  loaded_terms = zeros(base, pol_length)
  coeff_type = elem_type(base)
  
  load_node(s, :terms) do terms
    for i in 1:length(terms)
      load_node(s, i) do (exponent, _)
        if serialize_with_params(coeff_type)
          if length(parents) == 1
            params = base
          else
            params = parents[1:end - 1]
          end
          c = load_object(s, coeff_type, params, 2)
        else
          c = load_object(s, coeff_type, 2)
        end
        e = parse(Int, exponent)
        loaded_terms[e] = c
      end
    end
  end
  return parent_ring(loaded_terms, pol_length, precision, valuation)
end

function load_object(s::DeserializerState, ::Type{<:AbsPowerSeriesRingElem}, parents::Vector)
  parent_ring = parents[end]
  pol_length = load_object(s, Int, :pol_length)
  precision = load_object(s, Int, :precision)
  base = base_ring(parent_ring)
  loaded_terms = zeros(base, pol_length)
  coeff_type = elem_type(base)

  load_node(s, :terms) do terms
    for i in 1:length(terms)
      load_node(s, i) do (exponent, _)
        if serialize_with_params(coeff_type)
          if length(parents) == 1
            params = base
          else
            params = parents[1:end - 1]
          end
          c = load_object(s, coeff_type, params, 2)
        else
          c = load_object(s, coeff_type, 2)
        end
        e = parse(Int, exponent)
        e += 1
        loaded_terms[e] = c
      end
    end
  end
  return parent_ring(loaded_terms, pol_length, precision)
end

################################################################################
# Laurent Series
@register_serialization_type Generic.LaurentSeriesRing "LaurentSeriesRing" uses_id
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
                       ZZLaurentSeriesRing}})
  base_ring = load_typed_object(s, :base_ring)
  var = load_object(s, Symbol, :var)
  max_precision = load_object(s, Int, :max_precision)

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

function load_object(s::DeserializerState,
                     ::Type{<: Union{Generic.LaurentSeriesElem, ZZLaurentSeriesRingElem}},
                     parents::Vector)
  parent_ring = parents[end]

  terms = load_node(s, :terms) do terms_data
    exponents = []
    for i in 1:length(terms_data)
      load_node(s, i) do _
        push!(exponents, load_object(s, Int, 1))
      end
    end
    
    highest_degree = max(exponents...)
    lowest_degree = min(exponents...)
    base = base_ring(parent_ring)
    coeff_type = elem_type(base)
    # account for index shift
    loaded_terms = zeros(base, highest_degree - lowest_degree + 1)
    for (i, e) in enumerate(exponents)
      load_node(s, i) do _
        e -= lowest_degree - 1
        if serialize_with_params(coeff_type)
          if length(parents) == 1
            params = base
          else
            params = parents[1:end - 1]
          end
          c = load_object(s, coeff_type, params, 2)
        else
          c = load_object(s, coeff_type, 2)
        end
        loaded_terms[e] = c
      end
    end
    return loaded_terms
  end
  valuation = load_object(s, Int, :valuation)
  pol_length = load_object(s, Int, :pol_length)
  precision = load_object(s, Int, :precision)
  scale = load_object(s, Int, :scale)
  return parent_ring(terms, pol_length, precision, valuation, scale)
end

### Affine algebras
@register_serialization_type MPolyQuoRing uses_id

function save_object(s::SerializerState, A::MPolyQuoRing)
  save_data_dict(s) do # Saves stuff in a JSON dictionary. This opens a `{`, puts stuff 
                       # inside there for the various keys and then closes it with `}`.
                       # It's not using Julia Dicts.
    save_typed_object(s, modulus(A), :modulus)
    save_typed_object(s, ordering(A), :ordering) # Does this already serialize???
  end
end

function load_object(s::DeserializerState, ::Type{MPolyQuoRing})
  I = load_typed_object(s, :modulus) 
  R = base_ring(I)
  o = load_typed_object(s, :ordering)
  return MPolyQuoRing(R, I, o)
end

### Serialization of Monomial orderings
@register_serialization_type MonomialOrdering uses_params

function save_object(s::SerializerState, o::MonomialOrdering)
  save_data_dict(s) do
    save_typed_object(s, o.o, :internal_ordering) # TODO: Is there a getter for this?
    if isdefined(o, :is_total)
      save_object(s, o.is_total, :is_total)
    end
  end
end

function load_object(s::DeserializerState, ::Type{MonomialOrdering}, ring::MPolyRing)
  ord = load_typed_object(s, :internal_ordering)
  result = MonomialOrdering(ring, ord)

  if haskey(s, :is_total)
    result.is_total = load_object(s, Bool, :is_total)
  end
  return result
end

# we will need to extend this to more orderings at some point
@register_serialization_type Orderings.SymbOrdering

function save_object(s::SerializerState, o::Orderings.SymbOrdering{S}) where {S}
  save_data_dict(s) do
    save_typed_object(s, S, :ordering_symbol_as_type)
    save_typed_object(s, o.vars, :vars) # TODO: Is there a getter?
  end
end

function load_object(s::DeserializerState, ::Type{Orderings.SymbOrdering})
  S = load_typed_object(s, :ordering_symbol_as_type)
  vars = load_typed_object(s, :vars)
  return Orderings.SymbOrdering(S, vars)
end

