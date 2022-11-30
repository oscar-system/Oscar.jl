################################################################################
# ring of integers (singleton type)
@registerSerializationType(FlintIntegerRing)


################################################################################
#  non simpleton base rings
@registerSerializationType(Nemo.NmodRing, "Nemo.NmodRing")

function save_internal(s::SerializerState, R::Nemo.NmodRing)
    return Dict(
        :modulus => save_type_dispatch(s, modulus(R))
    )
end

function load_internal(s::DeserializerState, ::Type{Nemo.NmodRing}, dict::Dict)
    modulus = load_type_dispatch(s, UInt64, dict[:modulus])
    return Nemo.NmodRing(modulus)
end

#elements
@registerSerializationType(nmod)

function save_internal(s::SerializerState, r::nmod)
    return Dict(
        :parent => save_type_dispatch(s, parent(r)),
        :class_val => save_type_dispatch(s, fmpz(r))
    )
end

function load_internal(s::DeserializerState, ::Type{nmod}, dict::Dict)
    parent_ring = load_type_dispatch(s, Nemo.NmodRing, dict[:parent])
    class_val = load_type_dispatch(s, fmpz, dict[:class_val])
    return parent_ring(class_val)
end

function load_internal_with_parent(s::DeserializerState,
                                   ::Type{nmod},
                                   dict::Dict,
                                   parent_ring::Nemo.NmodRing)
    class_val = load_type_dispatch(s, fmpz, dict[:class_val])
    return parent_ring(class_val)
end


################################################################################
#  Polynomial Rings

@registerSerializationType(PolyRing)

@registerSerializationType(MPolyRing)

@registerSerializationType(FmpqMPolyRing)
@registerSerializationType(FmpqPolyRing)
@registerSerializationType(FmpzMPolyRing)
@registerSerializationType(FmpzPolyRing)
@registerSerializationType(FqNmodMPolyRing)
@registerSerializationType(FqNmodPolyRing)
@registerSerializationType(GFPPolyRing)
@registerSerializationType(NmodMPolyRing)
@registerSerializationType(NmodPolyRing)

function save_internal(s::SerializerState, R::Union{MPolyRing, PolyRing})
    return Dict(
        :symbols => save_type_dispatch(s, symbols(R)),
        :base_ring => save_type_dispatch(s, base_ring(R)),
    )
end

function load_internal(s::DeserializerState,
                       T::Type{<: Union{MPolyRing, PolyRing}},
                       dict::Dict)
    base_ring = load_unknown_type(s, dict[:base_ring])
    symbols = load_type_dispatch(s, Vector{Symbol}, dict[:symbols])

    if T <: PolyRing
        return PolynomialRing(base_ring, symbols..., cached=false)[1]
    end

    return PolynomialRing(base_ring, symbols, cached=false)[1]
end

################################################################################
# Multivariate Polynomials
@registerSerializationType(fmpq_mpoly)
@registerSerializationType(fmpz_mpoly)
@registerSerializationType(fq_nmod_mpoly)
@registerSerializationType(nmod_mpoly)
@registerSerializationType(MPolyElem)

function save_internal(s::SerializerState, p::MPolyElem)
    parent_ring = parent(p)
    parent_ring = save_type_dispatch(s, parent_ring)
    terms = []

    for i in 1:length(p)
        term = Dict(
            :exponent => save_type_dispatch(s, exponent_vector(p, i)),
            :coeff => save_type_dispatch(s, coeff(p, i))
        )
        push!(terms, term)
    end

    return Dict(
        :terms => terms,
        :parent => parent_ring,
    )
end

function load_internal(s::DeserializerState, ::Type{<: MPolyElem}, dict::Dict)
    R = load_unknown_type(s, dict[:parent])
    coeff_ring = coefficient_ring(R)
    coeff_type = elem_type(coeff_ring)
    polynomial = MPolyBuildCtx(R)

    for term in dict[:terms]
        c = load_type_dispatch(s, coeff_type, term[:coeff])
        e = load_type_dispatch(s, Vector{Int}, term[:exponent])
        push_term!(polynomial, c, e)
    end
    return finish(polynomial)
end

function load_internal_with_parent(s::DeserializerState,
                                   ::Type{<: MPolyElem},
                                   dict::Dict,
                                   parent_ring::MPolyRing)
    coeff_ring = coefficient_ring(parent_ring)
    coeff_type = elem_type(coeff_ring)
    polynomial = MPolyBuildCtx(parent_ring)

    for term in dict[:terms]
        c = load_type_dispatch(s, coeff_type, term[:coeff]; parent=coeff_ring)
        e = load_type_dispatch(s, Vector{Int}, term[:exponent])
        push_term!(polynomial, c, e)
    end
    return finish(polynomial)
end

################################################################################
# Univariate Polynomials

@registerSerializationType(PolyElem)

@registerSerializationType(fmpq_poly)
@registerSerializationType(fmpz_poly)
@registerSerializationType(fq_nmod_poly)
@registerSerializationType(gfp_poly)
@registerSerializationType(nmod_poly)

function save_internal(s::SerializerState, p::PolyElem)
    parent_ring = parent(p)
    parent_ring = save_type_dispatch(s, parent_ring)

    return Dict(
        :parent => parent_ring,
        :coeffs => save_type_dispatch(s, collect(coefficients(p)))
    )
end

function load_internal(s::DeserializerState, ::Type{<: PolyElem}, dict::Dict)
    R = load_unknown_type(s, dict[:parent])
    coeff_ring = coefficient_ring(R)
    coeff_type = elem_type(coeff_ring)
    coeffs = load_type_dispatch(s, Vector{coeff_type}, dict[:coeffs])

    return R(coeffs)
end

function load_internal_with_parent(s::DeserializerState,
                                   ::Type{<: PolyElem},
                                   dict::Dict,
                                   parent_ring::PolyRing)
    coeff_ring = coefficient_ring(parent_ring)
    coeff_type = elem_type(coeff_ring)
    coeffs = load_type_dispatch(s, Vector{coeff_type}, dict[:coeffs]; parent=coeff_ring)

    return parent_ring(coeffs)
end

################################################################################
# Polynomial Ideals

@registerSerializationType(MPolyIdeal)

function save_internal(s::SerializerState, i::MPolyIdeal)
    generators = gens(i)
    parent_ring = save_type_dispatch(s, parent(generators[1]))

    return Dict(
        :parent => parent_ring,
        :gens => save_type_dispatch(s, generators),
    )
end

function load_internal(s::DeserializerState, ::Type{<: MPolyIdeal}, dict::Dict)
    parent_ring = load_unknown_type(s, dict[:parent])
    gens = load_type_dispatch(s, Vector{elem_type(parent_ring)}, dict[:gens])

    return ideal(parent_ring, gens)
end

function load_internal_with_parent(s::DeserializerState,
                                   ::Type{<: MPolyIdeal},
                                   dict::Dict,
                                   parent_ring::MPolyRing)
    gens = load_type_dispatch(s, Vector{elem_type(parent_ring)}, dict[:gens])
    return ideal(parent_ring, gens)
end

################################################################################
# Matrices

@registerSerializationType(MatElem)

@registerSerializationType(fmpz_mat)
@registerSerializationType(fmpq_mat)
@registerSerializationType(fq_nmod_mat)
@registerSerializationType(nmod_mat)

function save_internal(s::SerializerState, m::MatrixElem)
    m_base_ring = save_type_dispatch(s, base_ring(parent(m)))
    return Dict(
        :matrix => save_type_dispatch(s, Array(m)),
        :base_ring => m_base_ring
    )
end

function load_internal(s::DeserializerState,
                       ::Type{<: MatElem},
                       dict::Dict)
    entries_ring = load_unknown_type(s, dict[:base_ring])
    mat = load_type_dispatch(s, Matrix, dict[:matrix])

    return matrix(entries_ring, mat)
end

function load_internal_with_parent(s::DeserializerState,
                                   ::Type{<: MatElem},
                                   dict::Dict,
                                   parent_ring::T) where T <: MatSpace
    mat = load_type_dispatch(s, Matrix, dict[:matrix], parent=base_ring(parent_ring))
    return matrix(base_ring(parent_ring), mat)
end

################################################################################
# Power Series
@registerSerializationType(SeriesRing)

function save_internal(s::SerializerState, R::Union{
    Generic.RelSeriesRing,
    FmpqRelSeriesRing,
    FmpzRelSeriesRing,
    FqNmodRelSeriesRing,
    NmodRelSeriesRing})
    return Dict(
        :base_ring => save_type_dispatch(s, base_ring(R)),
        :var => save_type_dispatch(s, var(R)),
        :max_precision => save_type_dispatch(s, max_precision(R)),
        :model => save_type_dispatch(s, :capped_relative)
    )
end

function save_internal(s::SerializerState, R::Union{
    Generic.AbsSeriesRing,
    FmpqAbsSeriesRing,
    FmpzAbsSeriesRing,
    FqNmodAbsSeriesRing,
    NmodAbsSeriesRing})
    return Dict(
        :base_ring => save_type_dispatch(s, base_ring(R)),
        :var => save_type_dispatch(s, var(R)),
        :max_precision => save_type_dispatch(s, max_precision(R)),
        :model => save_type_dispatch(s, :capped_absolute)
    )
end

function load_internal(s::DeserializerState, ::Type{<: SeriesRing}, dict::Dict)
    base_ring = load_unknown_type(s, dict[:base_ring])
    var = load_type_dispatch(s, Symbol, dict[:var])
    max_precision = load_type_dispatch(s, Int, dict[:max_precision])
    model = load_type_dispatch(s, Symbol, dict[:model])
    
    return PowerSeriesRing(base_ring, max_precision, var; cached=false, model=model)[1]
end

# elements
@registerSerializationType(RelSeriesElem)

@registerSerializationType(AbsSeriesElem)

function save_internal(s::SerializerState, r::RelSeriesElem)
    v = valuation(r)
    pl = pol_length(r)
    coeffs = map(x -> coeff(r, x), v:v + pl)
    return Dict(
        :parent => save_type_dispatch(s, parent(r)),
        :coeffs => save_type_dispatch(s, coeffs),
        :valuation => save_type_dispatch(s, v),
        :pol_length => save_type_dispatch(s, pl),
        :precision => save_type_dispatch(s, precision(r))
    )
end

function save_internal(s::SerializerState, r::AbsSeriesElem)
    coeffs = map(x -> coeff(r, x), 0:pol_length(r))
    return Dict(
        :parent => save_type_dispatch(s, parent(r)),
        :coeffs => save_type_dispatch(s, coeffs),
        :pol_length => save_type_dispatch(s, pol_length(r)),
        :precision => save_type_dispatch(s, precision(r))
    )
end

function load_internal(s::DeserializerState, ::Type{<: RelSeriesElem}, dict::Dict)
    parent = load_type_dispatch(s, SeriesRing, dict[:parent])
    coeffs = load_type_dispatch(s, Vector, dict[:coeffs])
    valuation = load_type_dispatch(s, Int, dict[:valuation])
    pol_length = load_type_dispatch(s, Int, dict[:pol_length])
    precision = load_type_dispatch(s, Int, dict[:precision])
    
    return parent(coeffs, pol_length, precision, valuation)
end

function load_internal_with_parent(s::DeserializerState,
                                   ::Type{<: RelSeriesElem},
                                   dict::Dict,
                                   parent_ring::SeriesRing)
    # cache parent inside serializer state in case the coefficient ring
    # needs to be checked against the coefficient ring of the passed parent
    # to ensure neccessary backrefs are loaded
    _ = load_unknown_type(s, dict[:parent])

    coeff_ring = base_ring(parent_ring)
    coeff_type = elem_type(coeff_ring)
    coeffs = load_type_dispatch(s, Vector{coeff_type}, dict[:coeffs]; parent=coeff_ring)
    valuation = load_type_dispatch(s, Int, dict[:valuation])
    pol_length = load_type_dispatch(s, Int, dict[:pol_length])
    precision = load_type_dispatch(s, Int, dict[:precision])

    if precision > max_precision(parent_ring)
        @warn("Precision Warning: given parent is less precise than serialized elem",
              maxlog=1)
    end

    return parent_ring(coeffs, pol_length, precision, valuation)
end

function load_internal(s::DeserializerState, ::Type{<: AbsSeriesElem}, dict::Dict)
    parent = load_type_dispatch(s, SeriesRing, dict[:parent])
    coeffs = load_type_dispatch(s, Vector, dict[:coeffs])
    pol_length = load_type_dispatch(s, Int, dict[:pol_length])
    precision = load_type_dispatch(s, Int, dict[:precision])
    
    return parent(coeffs, pol_length, precision)
end

function load_internal_with_parent(s::DeserializerState,
                                   ::Type{<: AbsSeriesElem},
                                   dict::Dict,
                                   parent_ring::SeriesRing)
    # cache parent inside serializer state in case parent needs
    # to be checked against the passed parent
    # to ensure neccessary backrefs are loaded
    _ = load_unknown_type(s, dict[:parent])

    coeff_ring = base_ring(parent_ring)
    coeff_type = elem_type(coeff_ring)
    coeffs = load_type_dispatch(s, Vector{coeff_type}, dict[:coeffs]; parent=coeff_ring)
    pol_length = load_type_dispatch(s, Int, dict[:pol_length])
    precision = load_type_dispatch(s, Int, dict[:precision])

    if precision > max_precision(parent_ring)
        @warn("Precision Warning: given parent is less precise than serialized elem",
              maxlog=1)
    end

    return parent_ring(coeffs, pol_length, precision)
end

################################################################################
# Laurent Series
@registerSerializationType(Generic.LaurentSeriesRing, "LaurentSeriesRing")

@registerSerializationType(Generic.LaurentSeriesField, "LaurentSeriesField")

@registerSerializationType(FmpzLaurentSeriesRing)

function save_internal(s::SerializerState, R::Union{
    Generic.LaurentSeriesRing,
    Generic.LaurentSeriesField,
    FmpzLaurentSeriesRing})
    return Dict(
        :base_ring => save_type_dispatch(s, base_ring(R)),
        :var => save_type_dispatch(s, var(R)),
        :max_precision => save_type_dispatch(s, max_precision(R)),
    )
end

function load_internal(s::DeserializerState,
                       ::Type{<: Union{
                           Generic.LaurentSeriesRing,
                           Generic.LaurentSeriesField,
                           FmpzLaurentSeriesRing}},
                       dict::Dict)
    base_ring = load_unknown_type(s, dict[:base_ring])
    var = load_type_dispatch(s, Symbol, dict[:var])
    max_precision = load_type_dispatch(s, Int, dict[:max_precision])

    return LaurentSeriesRing(base_ring, max_precision, var; cached=false)[1]
end

# elements
@registerSerializationType(Generic.LaurentSeriesFieldElem, "LaurentSeriesFieldElem")

@registerSerializationType(Generic.LaurentSeriesRingElem, "LaurentSeriesRingElem")

@registerSerializationType(fmpz_laurent_series)

function save_internal(s::SerializerState, r:: fmpz_laurent_series)
    v = valuation(r)
    l = pol_length(r)
    coeffs = map(x -> coeff(r, x), v:l + v)
    return Dict(
        :parent => save_type_dispatch(s, parent(r)),
        :coeffs => save_type_dispatch(s, coeffs),
        :valuation => save_type_dispatch(s, v),
        :pol_length => save_type_dispatch(s, l),
        :precision => save_type_dispatch(s, precision(r)),
        :scale => save_type_dispatch(s, Nemo.scale(r))
    )
end

function save_internal(s::SerializerState, r:: Generic.LaurentSeriesElem)
    v = valuation(r)
    l = pol_length(r)
    coeffs = map(x -> coeff(r, x), v:l + v)
    return Dict(
        :parent => save_type_dispatch(s, parent(r)),
        :coeffs => save_type_dispatch(s, coeffs),
        :valuation => save_type_dispatch(s, v),
        :pol_length => save_type_dispatch(s, l),
        :precision => save_type_dispatch(s, precision(r)),
        :scale => save_type_dispatch(s, Generic.scale(r))
    )
end

function load_internal(s::DeserializerState,
                       T::Type{<: Union{
                           Generic.LaurentSeriesElem,                           
                           fmpz_laurent_series}},
                       dict::Dict)
    parent = load_unknown_type(s, dict[:parent])
    coeffs = load_type_dispatch(s, Vector, dict[:coeffs])
    pol_length = load_type_dispatch(s, Int, dict[:pol_length])
    precision = load_type_dispatch(s, Int, dict[:precision])
    valuation = load_type_dispatch(s, Int, dict[:valuation])
    scale = load_type_dispatch(s, Int, dict[:scale])
    
    return parent(coeffs, pol_length, precision, valuation, scale)
end

function load_internal_with_parent(s::DeserializerState,
                                   T::Type{<: Union{
                                       Generic.LaurentSeriesElem, fmpz_laurent_series}},
                                   dict::Dict,
                                   parent_ring)
    # cache parent inside serializer state in case parent needs
    # to be checked against the passed parent
    # to ensure neccessary backrefs are loaded
    _ = load_unknown_type(s, dict[:parent])

    coeff_ring = base_ring(parent_ring)
    coeff_type = elem_type(coeff_ring)
    coeffs = load_type_dispatch(s, Vector{coeff_type}, dict[:coeffs]; parent=coeff_ring)
    pol_length = load_type_dispatch(s, Int, dict[:pol_length])
    precision = load_type_dispatch(s, Int, dict[:precision])
    valuation = load_type_dispatch(s, Int, dict[:valuation])
    scale = load_type_dispatch(s, Int, dict[:scale])

    if precision > max_precision(parent_ring)
        @warn("Precision Warning: given parent is less precise than serialized elem",
              maxlog=1)
    end

    return parent_ring(coeffs, pol_length, precision, valuation, scale)
end

