################################################################################
# Common union types

const ModRingUnion = Union{zzModRing, ZZModRing}
const ModRingElemUnion = Union{zzModRingElem, ZZModRingElem}

const PolyRingUnionType = Union{UniversalPolyRing,
                            MPolyRing,
                            PolyRing,
                            AbstractAlgebra.Generic.LaurentMPolyWrapRing}

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

# element types by default use their parent as reference object
type_params(x::T) where T <: SetElem = TypeParams(T, parent(x))
type_params(x::T) where T <: SMat = TypeParams(T, parent(x))

# rings, groups etc. default have no reference object
type_params(R::T) where T <: AbstractAlgebra.Set = TypeParams(T, nothing)

# ideals and matrix spaces have their base ring as reference object
type_params(x::T) where T <: Ideal = TypeParams(T, base_ring(x))
type_params(x::T) where T <: MatSpace = TypeParams(T, base_ring(x))
type_params(x::T) where T <: SMatSpace = TypeParams(T, base_ring(x))


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

function load_object(s::DeserializerState, tp::TypeParams{<:ModRingElemUnion, <:ModRingUnion})
  parent_ring = Oscar.params(tp)
  return parent_ring(load_object(s, ZZRingElem))
end

################################################################################
#  Polynomial Rings

@register_serialization_type PolyRing uses_id 
@register_serialization_type MPolyRing uses_id
@register_serialization_type UniversalPolyRing uses_id
@register_serialization_type MPolyDecRing uses_id
@register_serialization_type AbstractAlgebra.Generic.LaurentMPolyWrapRing uses_id

# polynomial-like rings use their coefficient ring as reference object
type_params(R::T) where T <: PolyRingUnionType = TypeParams(T, coefficient_ring(R))

function save_object(s::SerializerState, R::PolyRingUnionType)
  save_data_dict(s) do
    save_object(s, symbols(R), :symbols)
  end
end

function load_object(s::DeserializerState, tp::TypeParams{<:PolyRing, <:Ring})
  params = Oscar.params(tp)
  symbols = load_object(s, Vector{Symbol}, :symbols)
  return polynomial_ring(params, only(symbols); cached=false)[1]
end

function load_object(s::DeserializerState, tp::TypeParams{<:MPolyRing, <:Ring})
  params = Oscar.params(tp)
  symbols = load_object(s, Vector{Symbol}, :symbols)
  return polynomial_ring(params, symbols; cached=false)[1]
end

function load_object(s::DeserializerState, tp::TypeParams{<:UniversalPolyRing, <:Ring})
  params = Oscar.params(tp)
  symbols = load_object(s, Vector{Symbol}, :symbols)
  return universal_polynomial_ring(params, symbols; cached=false)[1]
end

function load_object(s::DeserializerState, tp::TypeParams{<:AbstractAlgebra.Generic.LaurentMPolyWrapRing, <:Ring})
  params = Oscar.params(tp)
  symbols = load_object(s, Vector{Symbol}, :symbols)
  return laurent_polynomial_ring(params, symbols; cached=false)[1]
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

function load_object(s::DeserializerState, tp::TypeParams{<:MPolyDecRing, <:Tuple{Vararg{Pair}}})
  ring = tp[:ring]
  grading_group = tp[:grading_group]
  grading = load_object(s, TypeParams(Vector{elem_type(grading_group)}, grading_group))
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
  save_object(s, [(exponent_vector(p, i), coeff(p, i)) for i in 1:length(p)])
end

function save_object(s::SerializerState, p::AbstractAlgebra.Generic.LaurentMPolyWrap)
  save_object(s, [(e, c) for (c, e) in zip(coefficients(p), AbstractAlgebra.exponent_vectors(p))])
end

function load_object(s::DeserializerState,
                     tp::TypeParams{<:Union{MPolyRingElem, UniversalPolyRingElem, AbstractAlgebra.Generic.LaurentMPolyWrap}, <:PolyRingUnionType})
  parent_ring = Oscar.params(tp)
  coeff_ring = coefficient_ring(parent_ring)
  polynomial = MPolyBuildCtx(parent_ring)
  coeff_type = elem_type(coeff_ring)
  exps_coeffs = load_object(s, TypeParams(Vector{Tuple{Vector{Int}, coeff_type}}, (nothing, coeff_ring)))

  for (e, c) in exps_coeffs
    push_term!(polynomial, c, e)
  end
  return finish(polynomial)
end

function load_object(s::DeserializerState, tp::TypeParams{<:MPolyDecRingElem, <:MPolyDecRing})
  parent_ring = Oscar.params(tp)
  poly = load_object(s, TypeParams(MPolyRingElem, forget_grading(parent_ring)))
  return parent_ring(poly)
end

################################################################################
# Univariate Polynomials

@register_serialization_type PolyRingElem

function save_object(s::SerializerState, p::PolyRingElem)
  save_object(s, [
    (i - 1, c) for (i, c) in enumerate(coefficients(p))
      if !is_zero(c) ])
end

function save_object(s::SerializerState{IPCSerializer},
                     p::PolyRingElem)
  save_object(s, collect(coefficients(p)))
end

function load_object(s::DeserializerState{IPCSerializer},
                     tp::TypeParams{<:PolyRingElem, <:PolyRing})
  parent_ring = Oscar.params(tp)
  CR = coefficient_ring(parent_ring)
  parent_ring(load_object(s, TypeParams(Vector{elem_type(CR)}, CR)))
end

function load_object(s::DeserializerState,
                     tp::TypeParams{<:PolyRingElem, <:PolyRing})
  parent_ring = Oscar.params(tp)
  coeff_ring = coefficient_ring(parent_ring)
  coeff_type = elem_type(coeff_ring)
  exps_coeffs = load_object(s, TypeParams(Vector{Tuple{Int, coeff_type}}, (nothing, coeff_ring)))

  isempty(exps_coeffs) && return parent_ring(0)

  degree = max([e for (e, _) in exps_coeffs]...)
  loaded_terms = Hecke.zeros_array(coeff_ring, degree + 1)
  for (e, c) in exps_coeffs
    loaded_terms[e + 1] = c
  end
  return parent_ring(loaded_terms)
end


################################################################################
# Polynomial Ideals
@register_serialization_type MPolyIdeal
@register_serialization_type LaurentMPolyIdeal
@register_serialization_type MPolyLocalizedIdeal
@register_serialization_type MPolyQuoLocalizedIdeal
@register_serialization_type MPolyQuoIdeal
@register_serialization_type Hecke.PIDIdeal

function save_object(s::SerializerState, I::Ideal)
  # we might want to serialize generating_system(I) and I.gb
  # in the future
  save_object(s, gens(I))
end

function load_object(s::DeserializerState, tp::TypeParams{<:Ideal, <:NCRing})
  parent_ring = Oscar.params(tp)
  gens = elem_type(parent_ring)[]
  load_array_node(s) do _
    push!(gens, load_object(s, TypeParams(elem_type(parent_ring), parent_ring)))
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

function load_object(s::DeserializerState, tp::TypeParams{<:IdealGens, <:Tuple{Vararg{Pair}}})
  base_ring = tp[:base_ring]
  ordering_type = tp[:ordering_type]

  if ordering_type <: MonomialOrdering
    ord = load_object(s, TypeParams(ordering_type, base_ring), :ordering)
  else
    ord = load_node(s, :ordering) do _
      MonomialOrdering(base_ring, load_object(s, ordering_type, :internal_ordering))
    end
  end
  generators = load_object(s, TypeParams(Vector{elem_type(base_ring)}, base_ring), :gens)
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

function load_object(s::DeserializerState, tp::TypeParams{MatSpace, <:NCRing})
  base_ring = Oscar.params(tp)
  ncols = load_object(s, Int, :ncols)
  nrows = load_object(s, Int, :nrows)
  return matrix_space(base_ring, nrows, ncols)
end

function load_object(s::DeserializerState, tp::TypeParams{SMatSpace, <:Ring})
  base_ring = Oscar.params(tp)
  ncols = load_object(s, Int, :ncols)
  nrows = load_object(s, Int, :nrows)
  return SMatSpace(base_ring, nrows, ncols)
end

# elems
function save_object(s::SerializerState, obj::MatElem)
  save_object(s, Matrix(obj))
end

function save_object(s::SerializerState, obj::SMat)
  save_data_array(s) do
    for r in obj
      save_object(s, collect(r))
    end
  end
end

function load_object(s::DeserializerState, tp::TypeParams{<:MatElem, <:MatSpace{T}}) where T
  parent = Oscar.params(tp)
  m = load_object(s, TypeParams(Matrix{T}, base_ring(parent)))
  if isempty(m)
    return parent()
  end
  return parent(m)
end

function load_object(s::DeserializerState, tp::TypeParams{<:SMat, <:SMatSpace{T}}) where T
  parent = Oscar.params(tp)
  base = base_ring(parent)
  M = sparse_matrix(base)

  load_array_node(s) do _
    row_entries = Tuple{Int, T}[]
    load_array_node(s) do _
      push!(row_entries, load_object(s, TypeParams(Tuple{Int, T}, (nothing, base))))
    end
    push!(M, sparse_row(base, row_entries))
  end
  M.c = parent.cols
  @assert nrows(M) == parent.rows
  return M
end

################################################################################
# Power Series
@register_serialization_type SeriesRing uses_id

type_params(R::T) where T <: SeriesRing = TypeParams(T, base_ring(R))

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

function load_object(s::DeserializerState, tp::TypeParams{<:SeriesRing, <:Ring})
  base_ring = Oscar.params(tp)
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

function load_object(s::DeserializerState, tp::TypeParams{<:RelPowerSeriesRingElem, <:RelPowerSeriesUnionType})
  parent_ring = Oscar.params(tp)
  valuation = load_object(s, Int, :valuation)
  pol_length = load_object(s, Int, :pol_length)
  precision = load_object(s, Int, :precision)
  base = base_ring(parent_ring)
  loaded_terms = Hecke.zeros_array(base, pol_length)
  coeff_type = elem_type(base)

  load_node(s, :terms) do _
    load_array_node(s) do _
      e = load_object(s, Int, 1)
      loaded_terms[e] = load_object(s, TypeParams(coeff_type, base), 2)
    end
  end
  return parent_ring(loaded_terms, pol_length, precision, valuation)
end

function load_object(s::DeserializerState, tp::TypeParams{<:AbsPowerSeriesRingElem, <:AbsPowerSeriesUnionType})
  parent_ring = Oscar.params(tp)
  pol_length = load_object(s, Int, :pol_length)
  precision = load_object(s, Int, :precision)
  base = base_ring(parent_ring)
  loaded_terms = Hecke.zeros_array(base, pol_length)
  coeff_type = elem_type(base)

  load_node(s, :terms) do _
    load_array_node(s) do _
      e = load_object(s, Int, 1)
      loaded_terms[e + 1] = load_object(s, TypeParams(coeff_type, base), 2)
    end
  end
  return parent_ring(loaded_terms, pol_length, precision)
end

################################################################################
# Laurent Series
@register_serialization_type Generic.LaurentSeriesRing "LaurentSeriesRing" uses_id
@register_serialization_type Generic.LaurentSeriesField "LaurentSeriesField" uses_id
@register_serialization_type ZZLaurentSeriesRing uses_id

type_params(R::T) where T <: LaurentUnionType = TypeParams(T, base_ring(R))

function save_object(s::SerializerState, R::LaurentUnionType)
  save_data_dict(s) do
    save_object(s, var(R), :var)
    save_object(s, max_precision(R), :max_precision)
  end
end

function load_object(s::DeserializerState, tp::TypeParams{<:LaurentUnionType, <:Ring})
  base_ring = Oscar.params(tp)
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
                     tp::TypeParams{<:Union{Generic.LaurentSeriesElem, ZZLaurentSeriesRingElem}, <:LaurentUnionType})
  parent_ring = Oscar.params(tp)
  terms = load_node(s, :terms) do _
    exponents = Int[]
    for i in 1:length(s.obj)
      load_node(s, i) do _
        push!(exponents, load_object(s, Int, 1))
      end
    end

    highest_degree = max(exponents...)
    lowest_degree = min(exponents...)
    base = base_ring(parent_ring)
    coeff_type = elem_type(base)
    loaded_terms = Hecke.zeros_array(base, highest_degree - lowest_degree + 1)
    for (i, e) in enumerate(exponents)
      e -= lowest_degree - 1
      load_node(s, i) do _
        loaded_terms[e] = load_object(s, TypeParams(coeff_type, base), 2)
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

function load_object(s::DeserializerState, tp::TypeParams{<:MPolyQuoRing, <:Tuple{Vararg{Pair}}})
  R = tp[:base_ring]
  ordering_type = tp[:ordering]
  o = load_object(s, TypeParams(ordering_type, R), :ordering)
  I = load_object(s, TypeParams(ideal_type(R), R), :modulus)

  return MPolyQuoRing(R, I, o)
end

@register_serialization_type MPolyQuoRingElem

function save_object(s::SerializerState, a::MPolyQuoRingElem)
  save_object(s, lift(a))
end

function load_object(s::DeserializerState, tp::TypeParams{<:MPolyQuoRingElem, <:MPolyQuoRing})
  Q = Oscar.params(tp)
  R = base_ring(Q)
  rep = load_object(s, TypeParams(elem_type(R), R))
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

function load_object(s::DeserializerState, tp::TypeParams{<:MonomialOrdering, <:MPolyRing})
  ring = Oscar.params(tp)
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
  save_object(s, denominators(U))
end

function load_object(s::DeserializerState, tp::TypeParams{<:MPolyPowersOfElement, <:MPolyRing})
  R = Oscar.params(tp)
  dens = load_object(s, TypeParams(Vector{elem_type(R)}, R))
  return MPolyPowersOfElement(R, dens)
end


@register_serialization_type MPolyComplementOfPrimeIdeal uses_id

type_params(U::MPolyComplementOfPrimeIdeal) = TypeParams(typeof(U), ring(U))

function save_object(s::SerializerState, U::MPolyComplementOfPrimeIdeal)
  save_object(s, prime_ideal(U))
end

function load_object(s::DeserializerState, tp::TypeParams{<:MPolyComplementOfPrimeIdeal, <:MPolyRing})
  R = Oscar.params(tp)
  id = load_object(s, TypeParams(ideal_type(R), R))
  return MPolyComplementOfPrimeIdeal(id)
end

@register_serialization_type MPolyLocRing uses_id

type_params(W::T) where {T <: MPolyLocRing} = TypeParams(T, :base_ring => base_ring(W), :mult_set_type => TypeParams(typeof(inverted_set(W)), nothing)) # TODO: make this neater!

function save_object(s::SerializerState, L::MPolyLocRing)
  save_object(s, inverted_set(L))
end

function load_object(s::DeserializerState, tp::TypeParams{<:MPolyLocRing, <:Tuple{Vararg{Pair}}})
  U = tp[:mult_set_type]
  R = tp[:base_ring]
  mult_set = load_object(s, TypeParams(U, R))
  return MPolyLocRing(R, mult_set)
end

@register_serialization_type MPolyLocRingElem

function save_object(s::SerializerState, a::MPolyLocRingElem)
  # `save_type_params` will store the output of type_params
  # in this case the parent ring
  save_data_array(s) do
    save_object(s, numerator(a))
    save_object(s, denominator(a))
  end
end

function load_object(s::DeserializerState, tp::TypeParams{<:MPolyLocRingElem, <:MPolyLocRing})
  parent = Oscar.params(tp)
  P = base_ring(parent)
  RET = elem_type(P)
  num = load_object(s, TypeParams(RET, P), 1)
  den = load_object(s, TypeParams(RET, P), 2)
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

function load_object(s::DeserializerState, tp::TypeParams{<:MPolyQuoLocRing, <:Tuple{Vararg{Pair}}})
  R = tp[:base_ring]::MPolyRing
  L = tp[:loc_ring]::MPolyLocRing
  Q = tp[:quo_ring]::MPolyQuoRing
  return MPolyQuoLocRing(R, modulus(Q), inverted_set(L), Q, L)
end

@register_serialization_type MPolyQuoLocRingElem

function save_object(s::SerializerState, a::MPolyQuoLocRingElem)
 save_object(s, [lifted_numerator(a), lifted_denominator(a)])
end

function load_object(s::DeserializerState, tp::TypeParams{<:MPolyQuoLocRingElem, <:MPolyQuoLocRing})
  parent = Oscar.params(tp)
  P = base_ring(parent)
  RET = elem_type(P)
  num = load_object(s, TypeParams(RET, P), 1)
  den = load_object(s, TypeParams(RET, P), 2)
  return parent(num, den; check=false)
end

@register_serialization_type MPolyComplementOfKPointIdeal uses_id

type_params(U::T) where {T<:MPolyComplementOfKPointIdeal} = TypeParams(T, ring(U))

function save_object(s::SerializerState, U::MPolyComplementOfKPointIdeal)
  save_object(s, point_coordinates(U))
end

function load_object(s::DeserializerState, tp::TypeParams{<:MPolyComplementOfKPointIdeal, <:Ring})
  R = Oscar.params(tp)
  kk = coefficient_ring(R)
  T = elem_type(kk)
  a = load_object(s, TypeParams(Vector{T}, kk))
  return MPolyComplementOfKPointIdeal(R, a)
end

### Morphisms of the four types of rings

@register_serialization_type MPolyLocalizedRingHom uses_id

type_params(phi::T) where {T<:MPolyLocalizedRingHom} = TypeParams(T, :domain=>domain(phi), :codomain=>codomain(phi), :restricted_map_params=>type_params(restricted_map(phi)))

function save_object(s::SerializerState, phi::MPolyLocalizedRingHom)
  save_object(s, restricted_map(phi))
end

function load_object(s::DeserializerState, tp::TypeParams{T, <:Tuple{Vararg{Pair}}}) where {T<:MPolyLocalizedRingHom}
  dom = tp[:domain]
  cod = tp[:codomain]
  rm_tp = tp[:restricted_map_params]
  res = load_object(s, TypeParams(MPolyAnyMap, rm_tp))
  return MPolyLocalizedRingHom(dom, cod, res; check=false)
end

@register_serialization_type MPolyQuoLocalizedRingHom uses_id

type_params(phi::T) where {T<:MPolyQuoLocalizedRingHom} = TypeParams(T, :domain=>domain(phi), :codomain=>codomain(phi), :restricted_map_params=>type_params(restricted_map(phi)))

function save_object(s::SerializerState, phi::MPolyQuoLocalizedRingHom)
  save_object(s, restricted_map(phi))
end

function load_object(s::DeserializerState, tp::TypeParams{T, <:Tuple{Vararg{Pair}}}) where {T<:MPolyQuoLocalizedRingHom}
  dom = tp[:domain]
  cod = tp[:codomain]
  rm_tp = tp[:restricted_map_params]
  res = load_object(s, TypeParams(MPolyAnyMap, rm_tp))
  return MPolyQuoLocalizedRingHom(dom, cod, res; check=false)
end

