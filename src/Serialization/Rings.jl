################################################################################
# Common union types

const RingMatElemUnion = Union{RingElem, MatElem, FreeAssociativeAlgebraElem,
                               SMat, TropicalSemiringElem}
const RingMatSpaceUnion = Union{Ring, MatSpace, SMatSpace,
                                FreeAssociativeAlgebra, TropicalSemiring}
const ModRingUnion = Union{zzModRing, ZZModRing}
const ModRingElemUnion = Union{zzModRingElem, ZZModRingElem}

const PolyRingUnionType = Union{UniversalPolyRing,
                            MPolyRing,
                            PolyRing,
                            AbstractAlgebra.Generic.LaurentMPolyWrapRing}

const IdealUnionType = Union{MPolyIdeal,
                                MPolyQuoIdeal,
                                MPolyLocalizedIdeal,
                                MPolyQuoLocalizedIdeal,
                                LaurentMPolyIdeal,
                                FreeAssociativeAlgebraIdeal,
                                IdealGens
                            }

const RelPowerSeriesUnionType = Union{Generic.RelPowerSeriesRing,
                                      QQRelPowerSeriesRing,
                                      ZZRelPowerSeriesRing,
                                      fqPolyRepRelPowerSeriesRing,
                                      FqRelPowerSeriesRing,
                                      zzModRelPowerSeriesRing}
const AbsPowerSeriesUnionType = Union{Generic.AbsPowerSeriesRing,
                                      QQAbsPowerSeriesRing,
                                      ZZAbsPowerSeriesRing,
                                      FqAbsPowerSeriesRing,
                                      fqPolyRepAbsPowerSeriesRing,
                                      zzModAbsPowerSeriesRing}

const LaurentUnionType = Union{Generic.LaurentSeriesRing,
                               Generic.LaurentSeriesField,
                               ZZLaurentSeriesRing}

################################################################################
# type_params functions

type_params(x::T) where T <: RingMatElemUnion = TypeParams(T, parent(x))
type_params(R::T) where T <: RingMatSpaceUnion = TypeParams(T, base_ring(R))
type_params(x::T) where T <: IdealUnionType = TypeParams(T, base_ring(x))
# exclude from ring union
type_params(::ZZRing) = TypeParams(ZZRing, nothing)
type_params(::ZZRingElem) = TypeParams(ZZRingElem, nothing)
type_params(R::T) where T <: ModRingUnion = TypeParams(T, nothing)

################################################################################
# ring of integers (singleton type)
@register_serialization_type ZZRing

################################################################################
#  Mod Rings
@register_serialization_type Nemo.zzModRing
@register_serialization_type Nemo.ZZModRing

function save_object(s::SerializerState, R::T) where T <: ModRingUnion
  save_object(s, modulus(R))
end

function load_object(s::DeserializerState, ::Type{zzModRing})
  modulus = load_object(s, UInt64)
  return zzModRing(modulus)
end

function load_object(s::DeserializerState, ::Type{ZZModRing})
  modulus = load_object(s, ZZRingElem)
  return ZZModRing(modulus)
end

#elements
@register_serialization_type zzModRingElem
@register_serialization_type ZZModRingElem

function save_object(s::SerializerState, x::ModRingElemUnion)
  save_data_basic(s, string(x))
end

function load_object(s::DeserializerState, ::Type{<:ModRingElemUnion},
                     parent_ring::T) where T <: ModRingUnion
  return parent_ring(load_object(s, ZZRingElem, ZZRing()))
end

################################################################################
#  Polynomial Rings

@register_serialization_type PolyRing uses_id 
@register_serialization_type MPolyRing uses_id
@register_serialization_type UniversalPolyRing uses_id
@register_serialization_type MPolyDecRing uses_id
@register_serialization_type AbstractAlgebra.Generic.LaurentMPolyWrapRing uses_id

function save_object(s::SerializerState, R::PolyRingUnionType)
  base = base_ring(R)
  save_data_dict(s) do
    save_object(s, symbols(R), :symbols)
  end
end

function load_object(s::DeserializerState,
                     T::Type{<: PolyRingUnionType},
                     params::Ring)
  symbols = load_object(s, Vector{Symbol}, :symbols)
  if T <: PolyRing
    return polynomial_ring(params, symbols..., cached=false)[1]
  elseif T <: UniversalPolyRing
    poly_ring = universal_polynomial_ring(params, cached=false)
    gens(poly_ring, symbols)
    return poly_ring
  elseif T <: AbstractAlgebra.Generic.LaurentMPolyWrapRing
    return laurent_polynomial_ring(params, symbols, cached=false)[1]
  end
  return polynomial_ring(params, symbols, cached=false)[1]
end

# with grading
type_params(R::MPolyDecRing) = TypeParams(
  MPolyDecRing,
  :grading_group => grading_group(R),
  :ring => forget_grading(R),
)

function save_object(s::SerializerState, R::MPolyDecRing)
  save_object(s, _grading(R))
end

function load_object(s::DeserializerState, ::Type{<:MPolyDecRing}, d::Dict)
  ring = d[:ring]
  grading = load_object(s, Vector{elem_type(d[:grading_group])}, d[:grading_group])
  return grade(ring, grading)[1]
end

################################################################################
#  Polynomial Ring Elem Types
@register_serialization_type MPolyRingElem
@register_serialization_type MPolyDecRingElem
@register_serialization_type UniversalPolyRingElem
@register_serialization_type AbstractAlgebra.Generic.LaurentMPolyWrap

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

@register_serialization_type PolyRingElem

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
                     parent_ring::PolyRing)
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
      load_node(s, i) do _
        load_node(s, 2) do _
          loaded_terms[exponent] = load_object(s, coeff_type, base)
        end
      end
    end
    return parent_ring(loaded_terms)
  end
end


function load_object(s::DeserializerState,
                     ::Type{<:Union{MPolyRingElem, UniversalPolyRingElem, AbstractAlgebra.Generic.LaurentMPolyWrap}},
                     parent_ring::PolyRingUnionType)
  load_node(s) do terms
    exponents = [term[1] for term in terms]
    base = base_ring(parent_ring)
    polynomial = MPolyBuildCtx(parent_ring)
    coeff_type = elem_type(base)
    for (i, e) in enumerate(exponents)
      load_node(s, i) do _
        c = load_object(s, coeff_type, base, 2)
        e_int = load_array_node(s, 1) do _
          load_object(s, Int)
        end
        push_term!(polynomial, c, e_int)
      end
    end
    return finish(polynomial)
  end
end

function load_object(s::DeserializerState, ::Type{<:MPolyDecRingElem}, parent_ring::MPolyDecRingElem)
  poly = load_object(s, MPolyRingElem, forget_grading(parent_ring))
  return parent_ring(poly)
end

################################################################################
# Polynomial Ideals
@register_serialization_type MPolyIdeal
@register_serialization_type LaurentMPolyIdeal
@register_serialization_type MPolyLocalizedIdeal
@register_serialization_type MPolyQuoLocalizedIdeal
@register_serialization_type MPolyQuoIdeal

function save_object(s::SerializerState, I::T) where T <: IdealUnionType
  # we might want to serialize generating_system(I) and I.gb
  # in the future
  save_object(s, gens(I))
end

function load_object(s::DeserializerState, ::Type{<: IdealUnionType}, parent_ring::RingMatSpaceUnion)
  gens = elem_type(parent_ring)[]
  load_array_node(s) do _
    push!(gens, load_object(s, elem_type(parent_ring), parent_ring))
  end
  return ideal(parent_ring, gens)
end

################################################################################
# IdealGens

# this will need adjustments to cover the NCRing case

@register_serialization_type IdealGens

type_params(ig::IdealGens) = TypeParams(
  IdealGens,
  :base_ring => base_ring(ig),
  :ordering_type => TypeParams(typeof(ordering(ig)), nothing)
)

function save_object(s::SerializerState, obj::IdealGens)
  save_data_dict(s) do
    save_object(s, ordering(obj), :ordering)
    save_object(s, gens(obj), :gens)
    save_object(s, is_groebner_basis(obj), :is_gb)
    save_object(s, obj.isReduced, :is_reduced)
    save_object(s, obj.keep_ordering, :keep_ordering)
  end
end

function load_object(s::DeserializerState, ::Type{<:IdealGens}, params::Dict)
  base_ring = params[:base_ring]
  ordering_type = params[:ordering_type]

  if ordering_type <: MonomialOrdering
    ord = load_object(s, ordering_type, base_ring, :ordering)
  else
    ord = load_node(s, :ordering) do _
      MonomialOrdering(base_ring, load_object(s, ordering_type, :internal_ordering))
    end
  end
  generators = load_object(s, Vector{elem_type(base_ring)}, base_ring, :gens)
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
@register_serialization_type MatElem
@register_serialization_type SMatSpace uses_id
@register_serialization_type SMat

function save_object(s::SerializerState, obj::MatSpace{T}) where T
  save_data_dict(s) do
    save_object(s, ncols(obj), :ncols)
    save_object(s, nrows(obj), :nrows)
  end
end

function save_object(s::SerializerState, obj::SMatSpace)
  save_data_dict(s) do
    # getters currently do not seem to exist
    save_object(s, obj.cols, :ncols)
    save_object(s, obj.rows, :nrows)
  end
end

function load_object(s::DeserializerState, ::Type{MatSpace}, base_ring::Ring)
  ncols = load_object(s, Int, :ncols)
  nrows = load_object(s, Int, :nrows)
  return matrix_space(base_ring, nrows, ncols)
end

function load_object(s::DeserializerState, ::Type{SMatSpace}, base_ring::Ring)
  ncols = load_object(s, Int, :ncols)
  nrows = load_object(s, Int, :nrows)
  return SMatSpace(base_ring, nrows, ncols)
end

# elems
function save_object(s::SerializerState, obj::MatElem)
  save_object(s, Array(obj))
end

function save_object(s::SerializerState, obj::SMat)
  save_data_array(s) do
    for r in obj
      save_object(s, collect(r))
    end
  end
end

function load_object(s::DeserializerState, ::Type{<:MatElem}, parent::MatSpace{T}) where T
  m = load_object(s, Matrix{T}, base_ring(parent))
  if isempty(m)
    return parent()
  end
  return parent(m)
end

function load_object(s::DeserializerState, ::Type{<:SMat}, parent::SMatSpace{T}) where T
  base = base_ring(parent)
  M = sparse_matrix(base)

  load_array_node(s) do _
    row_entries = Tuple{Int, T}[]
    load_array_node(s) do _
      push!(row_entries, load_object(s, Tuple{Int, T}, (nothing, base)))
    end
    push!(M, sparse_row(base, row_entries))
  end
  return M
end

################################################################################
# Power Series
@register_serialization_type SeriesRing uses_id


function save_object(s::SerializerState, R::RelPowerSeriesUnionType)
  save_data_dict(s) do
    save_object(s, var(R), :var)
    save_object(s, max_precision(R), :max_precision)
    save_object(s, :capped_relative, :model)
  end
end

function save_object(s::SerializerState, R::AbsPowerSeriesUnionType)
  save_data_dict(s) do
    save_object(s, var(R), :var)
    save_object(s, max_precision(R), :max_precision)
    save_object(s, :capped_absolute, :model)
  end
end

function load_object(s::DeserializerState, ::Type{<: SeriesRing}, base_ring::Ring)
  var = load_object(s, Symbol, :var)
  max_precision = load_object(s, Int, :max_precision)
  model = load_object(s, Symbol, :model)
  
  return power_series_ring(base_ring, max_precision, var; cached=false, model=model)[1]
end

# elements
@register_serialization_type RelPowerSeriesRingElem
@register_serialization_type AbsPowerSeriesRingElem

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

function load_object(s::DeserializerState, ::Type{<:RelPowerSeriesRingElem},
                     parent_ring::RelPowerSeriesUnionType)
  valuation = load_object(s, Int, :valuation)
  pol_length = load_object(s, Int, :pol_length)
  precision = load_object(s, Int, :precision)
  base = base_ring(parent_ring)
  loaded_terms = zeros(base, pol_length)
  coeff_type = elem_type(base)
  
  load_node(s, :terms) do _
    load_array_node(s) do _
      e = load_object(s, Int, 1)
      loaded_terms[e] = load_object(s, coeff_type, base, 2)
    end
  end
  return parent_ring(loaded_terms, pol_length, precision, valuation)
end

function load_object(s::DeserializerState, ::Type{<:AbsPowerSeriesRingElem},
                     parent_ring::AbsPowerSeriesUnionType)
  pol_length = load_object(s, Int, :pol_length)
  precision = load_object(s, Int, :precision)
  base = base_ring(parent_ring)
  loaded_terms = zeros(base, pol_length)
  coeff_type = elem_type(base)

  load_node(s, :terms) do _
    load_array_node(s) do _
      e = load_object(s, Int, 1)
      loaded_terms[e + 1] = load_object(s, coeff_type, base, 2)
    end
  end
  return parent_ring(loaded_terms, pol_length, precision)
end

################################################################################
# Laurent Series
@register_serialization_type Generic.LaurentSeriesRing "LaurentSeriesRing" uses_id
@register_serialization_type Generic.LaurentSeriesField "LaurentSeriesField" uses_id
@register_serialization_type ZZLaurentSeriesRing uses_id

function save_object(s::SerializerState, R::LaurentUnionType)
  save_data_dict(s) do
    save_object(s, var(R), :var)
    save_object(s, max_precision(R), :max_precision)
  end
end

function load_object(s::DeserializerState, ::Type{<: LaurentUnionType}, base_ring::Ring)
  var = load_object(s, Symbol, :var)
  max_precision = load_object(s, Int, :max_precision)

  return laurent_series_ring(base_ring, max_precision, var; cached=false)[1]
end

# elements
@register_serialization_type Generic.LaurentSeriesFieldElem "LaurentSeriesFieldElem"
@register_serialization_type Generic.LaurentSeriesRingElem "LaurentSeriesRingElem"
@register_serialization_type ZZLaurentSeriesRingElem

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
                     parent_ring::LaurentUnionType)
  terms = load_node(s, :terms) do terms_data
    # reading all exponents before ...
    # might be more efficient way ...
    exponents = Int[]
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
      e -= lowest_degree - 1
      load_node(s, i) do _
        loaded_terms[e] = load_object(s, coeff_type, base, 2)
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

type_params(A::MPolyQuoRing) = TypeParams(
  MPolyQuoRing,
  :base_ring => base_ring(A),
  :ordering => typeof(ordering(A))
)

function save_object(s::SerializerState, A::MPolyQuoRing)
  save_data_dict(s) do # Saves stuff in a JSON dictionary. This opens a `{`, puts stuff 
                       # inside there for the various keys and then closes it with `}`.
                       # It's not using Julia Dicts.
    save_object(s, modulus(A), :modulus)
    save_object(s, ordering(A), :ordering)
  end
end

function load_object(s::DeserializerState, ::Type{MPolyQuoRing}, params::Dict)
  R = params[:base_ring]
  ordering_type = params[:ordering]
  o = load_object(s, ordering_type, R, :ordering)
  I = load_object(s, ideal_type(R), R, :modulus)

  return MPolyQuoRing(R, I, o)
end

@register_serialization_type MPolyQuoRingElem

function save_object(s::SerializerState, a::MPolyQuoRingElem)
  save_object(s, lift(a))
end

function load_object(s::DeserializerState, ::Type{<:MPolyQuoRingElem}, Q::MPolyQuoRing)
  R = base_ring(Q)
  rep = load_object(s, elem_type(R), R)
  return Q(rep)
end

### Serialization of Monomial orderings
@register_serialization_type MonomialOrdering

function save_object(s::SerializerState, o::MonomialOrdering)
  save_data_dict(s) do
    save_object(s, o.o, :internal_ordering) # TODO: Is there a getter for this?
    if isdefined(o, :is_total)
      save_object(s, o.is_total, :is_total)
    end
  end
end

function load_object(s::DeserializerState, ::Type{<:MonomialOrdering}, ring::MPolyRing)
  # this will need to be changed to include other orderings, see below
  ord = load_object(s, Orderings.SymbOrdering, :internal_ordering)
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
    save_object(s, S, :ordering_symbol_as_type)
    save_object(s, o.vars, :vars) # TODO: Is there a getter?
  end
end

function load_object(s::DeserializerState, ::Type{Orderings.SymbOrdering})
  S = load_object(s, Symbol, :ordering_symbol_as_type)
  vars = load_object(s, Vector{Int}, :vars) # are these always Vector{Int} ?
  return Orderings.SymbOrdering(S, vars)
end


# localizations of polynomial rings
@register_serialization_type MPolyPowersOfElement uses_id

type_params(U::T) where T <: MPolyPowersOfElement = TypeParams(T, ring(U))

function save_object(s::SerializerState, U::MPolyPowersOfElement)
  save_data_dict(s) do
    save_object(s, denominators(U), :dens)
  end
end

function load_object(s::DeserializerState, ::Type{<:MPolyPowersOfElement}, R::MPolyRing)
  dens = Vector{elem_type(R)}(load_object(s, Vector{elem_type(R)}, R, :dens)) # casting is necessary for empty arrays
  return MPolyPowersOfElement(R, dens)
end


@register_serialization_type MPolyComplementOfPrimeIdeal uses_id

type_params(U::MPolyComplementOfPrimeIdeal) = TypeParams(typeof(U), ring(U))

function save_object(s::SerializerState, U::MPolyComplementOfPrimeIdeal)
  save_data_dict(s) do
    save_object(s, prime_ideal(U), :ideal)
  end
end

function load_object(s::DeserializerState, ::Type{<:MPolyComplementOfPrimeIdeal}, R::MPolyRing)
  id = load_object(s, ideal_type(R), R, :ideal)
  return MPolyComplementOfPrimeIdeal(id)
end

@register_serialization_type MPolyLocRing uses_id

type_params(W::T) where {T <: MPolyLocRing} = TypeParams(T, :base_ring => base_ring(W), :mult_set_type => TypeParams(typeof(inverted_set(W)), nothing)) # TODO: This seems to cause the trouble!

function save_object(s::SerializerState, L::MPolyLocRing)
  save_object(s, inverted_set(L))
end

function load_object(
    s::DeserializerState, 
    ::Type{<:MPolyLocRing}, params::Dict
  ) 
  U = params[:mult_set_type]
  R = params[:base_ring]
  mult_set = load_object(s, U, R)
  return MPolyLocRing(R, mult_set)
end

@register_serialization_type MPolyLocRingElem uses_params

type_params(a::MPolyLocRingElem) = TypeParams(MPolyLocRingElem, parent(a))

function save_object(s::SerializerState, a::MPolyLocRingElem)
  # Because the `parent` of `a` is a `Ring` the generic implementation
  # for `uses_params` above calls `save_type_params` and that stores 
  # the ring. Hopefully. 
  save_data_array(s) do
    save_object(s, numerator(a))
    save_object(s, denominator(a))
  end
end

function load_object(s::DeserializerState, ::Type{<:MPolyLocRingElem}, parent::MPolyLocRing)
  P = base_ring(parent)
  RET = elem_type(P)
  num = load_object(s, RET, P, 1)
  den = load_object(s, RET, P, 2)
  return parent(num, den; check=false)
end

@register_serialization_type MPolyQuoLocRing uses_id

type_params(L::T) where {T <: MPolyQuoLocRing} = TypeParams(T, :base_ring=>base_ring(L), :loc_ring=>localized_ring(L), :quo_ring=>underlying_quotient(L))

function save_object(s::SerializerState, L::MPolyQuoLocRing)
  save_data_dict(s) do
    # Everything happens in the type_params.
    # We still need to do something here, because otherwise 
    # we get an error. TODO: Make this better!
  end
end

function load_object(s::DeserializerState, ::Type{<:MPolyQuoLocRing}, params::Dict)
  R = params[:base_ring]::MPolyRing
  L = params[:loc_ring]::MPolyLocRing
  Q = params[:quo_ring]::MPolyQuoRing
  return MPolyQuoLocRing(R, modulus(Q), inverted_set(L), Q, L)
end

@register_serialization_type MPolyQuoLocRingElem uses_params

type_params(a::T) where {T<:MPolyQuoLocRingElem} = TypeParams(T, parent(a))

function save_object(s::SerializerState, a::MPolyQuoLocRingElem)
  save_data_array(s) do
    save_object(s, lifted_numerator(a))
    save_object(s, lifted_denominator(a))
  end
end

function load_object(s::DeserializerState, ::Type{<:MPolyQuoLocRingElem}, parent::MPolyQuoLocRing)
  P = base_ring(parent)
  RET = elem_type(P)
  num = load_object(s, RET, P, 1)
  den = load_object(s, RET, P, 2)
  return parent(num, den; check=false)
end

@register_serialization_type MPolyComplementOfKPointIdeal uses_id

type_params(U::T) where {T<:MPolyComplementOfKPointIdeal} = TypeParams(T, ring(U))

function save_object(s::SerializerState, U::MPolyComplementOfKPointIdeal)
  save_object(s, point_coordinates(U))
end

function load_object(s::DeserializerState, ::Type{<:MPolyComplementOfKPointIdeal}, R::Ring)
  kk = coefficient_ring(R)
  T = elem_type(kk)
  a = load_object(s, Vector{T}, kk)
  return MPolyComplementOfKPointIdeal(R, a)
end

### Morphisms of the four types of rings

@register_serialization_type MPolyLocalizedRingHom uses_id

type_params(phi::T) where {T<:MPolyLocalizedRingHom} = TypeParams(T, :domain=>domain(phi), :codomain=>codomain(phi), :restricted_map_params=>type_params(restricted_map(phi)))

function save_object(s::SerializerState, phi::MPolyLocalizedRingHom)
  save_object(s, restricted_map(phi))
end

function load_object(s::DeserializerState, ::Type{T}, params::Dict) where {T<:MPolyLocalizedRingHom} # RT is the type of the `restricted_map`
  dom = params[:domain]
  cod = params[:codomain]
  rm_tp = params[:restricted_map_params]
  res = load_object(s, MPolyAnyMap, rm_tp)
  return MPolyLocalizedRingHom(dom, cod, res; check=false)
end

@register_serialization_type MPolyQuoLocalizedRingHom uses_id

type_params(phi::T) where {T<:MPolyQuoLocalizedRingHom} = TypeParams(T, :domain=>domain(phi), :codomain=>codomain(phi), :restricted_map_params=>type_params(restricted_map(phi)))

function save_object(s::SerializerState, phi::MPolyQuoLocalizedRingHom)
  save_object(s, restricted_map(phi))
end

function load_object(s::DeserializerState, ::Type{T}, params::Dict) where {T<:MPolyQuoLocalizedRingHom} # RT is the type of the `restricted_map`
  dom = params[:domain]
  cod = params[:codomain]
  rm_tp = params[:restricted_map_params]
  res = load_object(s, MPolyAnyMap, rm_tp)
  return MPolyQuoLocalizedRingHom(dom, cod, res; check=false)
end

