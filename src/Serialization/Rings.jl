################################################################################
# Utility functions for ring parent tree

# builds parent tree
function get_parents(parent_ring::Ring)
    if has_elem_basic_encoding(parent_ring)
        return Any[]
    end
    base = base_ring(parent_ring)

    parents = get_parents(base)
    push!(parents, parent_ring)
    return parents
end

function get_parent_refs(s, parent_ring::Ring)
    parents = get_parents(parent_ring)
    return [save_as_ref(s, parent) for parent in parents]
end


################################################################################
# ring of integers (singleton type)
@registerSerializationType(ZZRing)


################################################################################
#  non simpleton base rings
@registerSerializationType(Nemo.zzModRing, "Nemo.zzModRing")
has_elem_basic_encoding(obj::Nemo.zzModRing) = true

function save_internal(s::SerializerState, R::Nemo.zzModRing)
    return Dict(
        :modulus => save_type_dispatch(s, modulus(R))
    )
end

function load_internal(s::DeserializerState, ::Type{Nemo.zzModRing}, dict::Dict)
    modulus = load_type_dispatch(s, UInt64, dict[:modulus])
    return Nemo.zzModRing(modulus)
end

#elements
@registerSerializationType(zzModRingElem)

function save_internal(s::SerializerState, r::zzModRingElem; include_parents::Bool=true)
    class_rep = string(r)
    if include_parents
        return Dict(
            :parent => save_as_ref(s, parent(r)),
            :class_rep => class_rep
        )
    end
    return string(r)
end

function load_internal(s::DeserializerState, ::Type{zzModRingElem}, dict::Dict)
    parent_ring = load_unknown_type(s, dict[:parent])
    class_rep = load_type_dispatch(s, ZZRingElem, dict[:class_rep])
    return parent_ring(class_rep)
end

function load_internal_with_parent(s::DeserializerState,
                                   ::Type{zzModRingElem},
                                   str::String,
                                   parent_ring::Nemo.zzModRing)
    return parent_ring(ZZRingElem(str))
end


################################################################################
#  Polynomial Rings

@registerSerializationType(PolyRing, true)
@registerSerializationType(MPolyRing, true)
@registerSerializationType(UniversalPolyRing, true)
@registerSerializationType(AbstractAlgebra.Generic.LaurentMPolyWrapRing, true)

function save_internal(s::SerializerState, R::Union{UniversalPolyRing, MPolyRing, PolyRing, AbstractAlgebra.Generic.LaurentMPolyWrapRing})
    base = base_ring(R)

    return Dict(
        :symbols => save_type_dispatch(s, symbols(R)),
        :base_ring => save_as_ref(s, base),
    )
end

function load_internal(s::DeserializerState,
                       T::Type{<: Union{UniversalPolyRing, MPolyRing, PolyRing, AbstractAlgebra.Generic.LaurentMPolyWrapRing}},
                       dict::Dict)
    base_ring = load_unknown_type(s, dict[:base_ring])
    symbols = load_type_dispatch(s, Vector{Symbol}, dict[:symbols])
    
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
# Multivariate and Universal Polynomials
@registerSerializationType(MPolyRingElem)
@registerSerializationType(UniversalPolyRingElem)

function save_internal(s::SerializerState, p::Union{UniversalPolyRingElem, MPolyRingElem};
                       include_parents::Bool=true)
    parent_ring = parent(p)
    base = base_ring(parent_ring)
    encoded_terms = []

    for i in 1:length(p)
        encoded_coeff = save_internal(s, coeff(p, i); include_parents=false)
        push!(encoded_terms,  (exponent_vector(p, i), encoded_coeff))
    end

    if include_parents
        return Dict(
            :terms => encoded_terms,
            :parents => get_parent_refs(s, parent_ring),
        )
    end
    return encoded_terms
end

@registerSerializationType(AbstractAlgebra.Generic.LaurentMPolyWrap)
function save_internal(s::SerializerState, p::AbstractAlgebra.Generic.LaurentMPolyWrap;
                       include_parents::Bool=true)
    parent_ring = parent(p)
    base = base_ring(parent_ring)
    encoded_terms = []

    exponent_vectors_gen = AbstractAlgebra.exponent_vectors(p)
    index = 0
    for c in coefficients(p)
        exponent_vector, index = iterate(exponent_vectors_gen, index)
        encoded_coeff = save_internal(s, c; include_parents=false)
        push!(encoded_terms,  (exponent_vector, encoded_coeff))
    end

    if include_parents
        return Dict(
            :terms => encoded_terms,
            :parents => get_parent_refs(s, parent_ring),
        )
    end
    return encoded_terms
end

################################################################################
# Univariate Polynomials

@registerSerializationType(PolyRingElem)

function save_internal(s::SerializerState, p::PolyRingElem; include_parents::Bool=true)
    # store a polynomial over a ring provided we can store elements in that ring
    parent_ring = parent(p)
    base = base_ring(parent_ring)
    coeffs = coefficients(p)
    encoded_terms = []
    exponent = 0
    for coeff in coeffs
        # collect only non trivial terms
        if is_zero(coeff)
            exponent += 1
            continue
        end
        
        encoded_coeff = save_internal(s, coeff; include_parents=false)
        push!(encoded_terms,  (exponent, encoded_coeff))
        exponent += 1
    end

    if include_parents
        return Dict(
            :parents => get_parent_refs(s, parent_ring),
            :terms => encoded_terms
        )
    end
    return encoded_terms
end

function load_terms(s::DeserializerState, parents::Vector, terms::Vector,
                    parent_ring::PolyRing)
    if isempty(terms)
        return parent_ring(0)
    end
    
    # shift so constant starts at 1
    degree = max(map(x->x[1] + 1, terms)...)
    base = base_ring(parent_ring)
    loaded_terms = zeros(base, degree)

    for term in terms
        exponent, coeff = term
        # account for shift
        exponent += 1

        if length(parents) == 1
            @assert has_elem_basic_encoding(base)
            coeff_type = elem_type(base)
            loaded_terms[exponent] = load_internal_with_parent(s, coeff_type, coeff, base)
        else
            loaded_terms[exponent] = load_terms(s, parents[1:end - 1], coeff, parents[end - 1])
        end
    end
    return parent_ring(loaded_terms)
end

function load_terms(s::DeserializerState, parents::Vector, terms::Vector,
                    parent_ring::Union{MPolyRing, UniversalPolyRing, AbstractAlgebra.Generic.LaurentMPolyWrapRing})
    base = base_ring(parent_ring)
    polynomial = MPolyBuildCtx(parent_ring)
    for (e, coeff) in terms
        if length(parents) == 1
            coeff_type = elem_type(base)
            @assert has_elem_basic_encoding(base)
            c = load_type_dispatch(s, coeff_type, coeff, parent=base)
        else
            c = load_terms(s, parents[1:end - 1], coeff, parents[end - 1])
        end
        
        push_term!(polynomial, c, Vector{Int}(e))
    end
    return finish(polynomial)
end


function load_internal(s::DeserializerState, ::Type{<: Union{
    PolyRingElem, UniversalPolyRingElem, MPolyRingElem, AbstractAlgebra.Generic.LaurentMPolyWrap}}, dict::Dict)
    loaded_parents = load_parents(s, dict[:parents])
    return load_terms(s, loaded_parents, dict[:terms], loaded_parents[end])
end

function load_internal_with_parent(s::DeserializerState, ::Type{<: Union{
    PolyRingElem, UniversalPolyRingElem, MPolyRingElem, AbstractAlgebra.Generic.LaurentMPolyWrap}}, dict::Dict,
                                   parent_ring::Union{PolyRing, MPolyRing, UniversalPolyRing, AbstractAlgebra.Generic.LaurentMPolyWrapRing})
    parents = get_parents(parent_ring)
    return load_terms(s, parents, dict[:terms], parents[end])
end

################################################################################
# Polynomial Ideals

@registerSerializationType(MPolyIdeal)
@registerSerializationType(Laurent.LaurentMPolyIdeal)

function save_internal(s::SerializerState,
                       I::Union{MPolyIdeal, Laurent.LaurentMPolyIdeal})
    generators = gens(I)

    return Dict(
        :gens => save_type_dispatch(s, generators),
    )
end

function load_internal(s::DeserializerState, ::Type{<: Union{
    MPolyIdeal, Laurent.LaurentMPolyIdeal}}, dict::Dict)
    gens = load_type_dispatch(s, Vector, dict[:gens])

    return ideal(parent(gens[1]), gens)
end

function load_internal_with_parent(s::DeserializerState,
                                   ::Type{<: Union{MPolyIdeal, Laurent.LaurentMPolyIdeal}},
                                   dict::Dict,
                                   parent_ring::Union{MPolyRing, AbstractAlgebra.Generic.LaurentMPolyWrapRing})
    gens = load_type_dispatch(s, Vector{elem_type(parent_ring)},
                              dict[:gens], parent=parent_ring)
    return ideal(parent_ring, gens)
end

################################################################################
# Matrices
@registerSerializationType(MatElem)

function save_internal(s::SerializerState, m::MatrixElem)
    return Dict(
        :base_ring => save_as_ref(s, base_ring(parent(m))),
        :matrix => save_type_dispatch(s, Array(m)),
    )
end

function load_internal(s::DeserializerState,
                       ::Type{<: MatElem},
                       dict::Dict)
    entries_ring = load_unknown_type(s, dict[:base_ring])
    mat = load_type_dispatch(s, Matrix, dict[:matrix]; parent=entries_ring)

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
@registerSerializationType(SeriesRing, true)

function save_internal(s::SerializerState, R::Union{
    Generic.RelPowerSeriesRing,
    QQRelPowerSeriesRing,
    ZZRelPowerSeriesRing,
    fqPolyRepRelPowerSeriesRing,
    zzModRelPowerSeriesRing})
    return Dict(
        :base_ring => save_as_ref(s, base_ring(R)),
        :var => save_type_dispatch(s, var(R)),
        :max_precision => save_type_dispatch(s, max_precision(R)),
        :model => save_type_dispatch(s, :capped_relative)
    )
end

function save_internal(s::SerializerState, R::Union{
    Generic.AbsPowerSeriesRing,
    QQAbsPowerSeriesRing,
    ZZAbsPowerSeriesRing,
    fqPolyRepAbsPowerSeriesRing,
    zzModAbsPowerSeriesRing})

    return Dict(
        :base_ring => save_as_ref(s, base_ring(R)),
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
    
    return power_series_ring(base_ring, max_precision, var; cached=false, model=model)[1]
end

# elements
@registerSerializationType(RelPowerSeriesRingElem)
@registerSerializationType(AbsPowerSeriesRingElem)

function save_internal(s::SerializerState, r::RelPowerSeriesRingElem;
                       include_parents::Bool=true)
    v = valuation(r)
    pl = pol_length(r)
    encoded_terms = []
    parent_ring = parent(r)
    base = base_ring(parent_ring)
    for exponent in v: v + pl
        coefficient = coeff(r, exponent)
        #collect only non trivial values
        if is_zero(coefficient)
            continue
        end

        encoded_coeff = save_internal(s, coefficient; include_parents=false)
        push!(encoded_terms,  (exponent, encoded_coeff))
    end

    if include_parents
        return Dict(
            :parents => get_parent_refs(s, parent_ring),
            :terms => encoded_terms,
            :valuation => save_type_dispatch(s, v),
            :pol_length => save_type_dispatch(s, pl),
            :precision => save_type_dispatch(s, precision(r))
        )
    end
    return encoded_terms
end

function save_internal(s::SerializerState, r::AbsPowerSeriesRingElem;
                       include_parents::Bool=true)
    pl = pol_length(r)
    encoded_terms = []
    parents = []
    parent_ring = parent(r)
    base = base_ring(parent_ring)
    for exponent in 0:pl
        coefficient = coeff(r, exponent)
        #collect only non trivial values
        if is_zero(coefficient)
            continue
        end

        encoded_coeff = save_internal(s, coefficient; include_parents=false)
        push!(encoded_terms,  (exponent, encoded_coeff))
    end

    if include_parents
        return Dict(
            :parents => get_parent_refs(s, parent_ring),
            :terms => encoded_terms,
            :pol_length => save_type_dispatch(s, pl),
            :precision => save_type_dispatch(s, precision(r))
        )
    end
    return encoded_terms
end

function load_internal(s::DeserializerState, ::Type{<: RelPowerSeriesRingElem}, dict::Dict)
    parents = load_parents(s, dict[:parents])
    parent_ring = parents[end]
    valuation = load_type_dispatch(s, Int, dict[:valuation])
    pol_length = load_type_dispatch(s, Int, dict[:pol_length])
    precision = load_type_dispatch(s, Int, dict[:precision])
    base = base_ring(parent_ring)
    loaded_terms = zeros(base, pol_length)

    for term in dict[:terms]
        exponent, coeff = term
        
        if length(parents) == 1
            @assert has_elem_basic_encoding(base)
            coeff_type = elem_type(base)
            loaded_terms[exponent] = load_type_dispatch(s, coeff_type, coeff,
                                                            parent=base)
        else
            loaded_terms[exponent] = load_terms(s, parents[1:end - 1], coeff, parents[end - 1])
        end
    end
    
    return parent_ring(loaded_terms, pol_length, precision, valuation)
end

function load_internal_with_parent(s::DeserializerState,
                                   ::Type{<: RelPowerSeriesRingElem},
                                   dict::Dict,
                                   parent_ring::SeriesRing)
    parents = get_parents(parent_ring)
    terms = load_terms(s, parents, dict[:terms], parents[end])
    valuation = load_type_dispatch(s, Int, dict[:valuation])
    pol_length = load_type_dispatch(s, Int, dict[:pol_length])
    precision = load_type_dispatch(s, Int, dict[:precision])

    return parent_ring(terms[valuation + 1: end], pol_length, precision, valuation)
end

function load_terms(s::DeserializerState, parents::Vector, terms::Vector, parent_ring::SeriesRing)
    highest_degree = max(map(x->x[1] + 1, terms)...)
    base = base_ring(parent_ring)
    # account for index shift
    loaded_terms = zeros(base, highest_degree + 1)
    for term in terms
        exponent, coeff = term
        exponent += 1
        if length(parents) == 1
            @assert has_elem_basic_encoding(base)
            coeff_type = elem_type(base)
            loaded_terms[exponent] = load_type_dispatch(s, coeff_type, coeff,
                                                            parent=base)
        else
            loaded_terms[exponent] = load_terms(s, parents[1:end - 1], coeff, parents[end - 1])
        end
    end
    
    return loaded_terms
end

function load_internal(s::DeserializerState, ::Type{<: AbsPowerSeriesRingElem}, dict::Dict)
    parents = load_parents(s, dict[:parents])
    parent_ring = parents[end]
    pol_length = load_type_dispatch(s, Int, dict[:pol_length])
    precision = load_type_dispatch(s, Int, dict[:precision])
    terms = load_terms(s, parents, dict[:terms], parents[end])

    parent_ring(terms, pol_length, precision)
end

function load_internal_with_parent(s::DeserializerState,
                                   ::Type{<: AbsPowerSeriesRingElem},
                                   dict::Dict,
                                   parent_ring::SeriesRing)
    parents = get_parents(parent_ring)
    terms = load_terms(s, parents, dict[:terms], parents[end])
    pol_length = load_type_dispatch(s, Int, dict[:pol_length])
    precision = load_type_dispatch(s, Int, dict[:precision])

    return parent_ring(terms, pol_length, precision)
end

################################################################################
# Laurent Series
@registerSerializationType(Generic.LaurentSeriesRing, true, "LaurentSeriesRing")
@registerSerializationType(Generic.LaurentSeriesField, true, "LaurentSeriesField")
@registerSerializationType(ZZLaurentSeriesRing)

function save_internal(s::SerializerState, R::Union{
    Generic.LaurentSeriesRing,
    Generic.LaurentSeriesField,
    ZZLaurentSeriesRing})
    return Dict(
        :base_ring => save_as_ref(s, base_ring(R)),
        :var => save_type_dispatch(s, var(R)),
        :max_precision => save_type_dispatch(s, max_precision(R)),
    )
end

function load_internal(s::DeserializerState,
                       ::Type{<: Union{
                           Generic.LaurentSeriesRing,
                           Generic.LaurentSeriesField,
                           ZZLaurentSeriesRing}},
                       dict::Dict)
    base_ring = load_unknown_type(s, dict[:base_ring])
    var = load_type_dispatch(s, Symbol, dict[:var])
    max_precision = load_type_dispatch(s, Int, dict[:max_precision])

    return laurent_series_ring(base_ring, max_precision, var; cached=false)[1]
end

# elements
@registerSerializationType(Generic.LaurentSeriesFieldElem, "LaurentSeriesFieldElem")
@registerSerializationType(Generic.LaurentSeriesRingElem, "LaurentSeriesRingElem")
@registerSerializationType(ZZLaurentSeriesRingElem)

function save_internal(s::SerializerState, r:: Union{Generic.LaurentSeriesElem, ZZLaurentSeriesRingElem};
                       include_parents::Bool=true)
    v = valuation(r)
    pl = pol_length(r)
    encoded_terms = []
    for exponent in v: v + pl
        coefficient = coeff(r, exponent)
        #collect only non trivial values
        if is_zero(coefficient)
            continue
        end

        encoded_coeff = save_internal(s, coefficient; include_parents=false)
        push!(encoded_terms,  (exponent, encoded_coeff))
    end

    return Dict(
        :parents => get_parent_refs(s, parent(r)),
        :terms => encoded_terms,
        :valuation => save_type_dispatch(s, v),
        :pol_length => save_type_dispatch(s, pl),
        :precision => save_type_dispatch(s, precision(r)),
        :scale => save_type_dispatch(s, Generic.scale(r))
    )
end

function load_terms(s::DeserializerState,
                    parents::Vector,
                    terms::Vector,
                    parent_ring::Union{
                           Generic.LaurentSeriesRing,
                           Generic.LaurentSeriesField,
                           ZZLaurentSeriesRing})
    highest_degree = max(map(x->x[1], terms)...)
    lowest_degree = min(map(x->x[1], terms)...)
    base = base_ring(parent_ring)
    # account for index shift
    loaded_terms = zeros(base, highest_degree - lowest_degree + 1)
    for term in terms
        exponent, coeff = term
        exponent -= lowest_degree - 1
        if length(parents) == 1
            @assert has_elem_basic_encoding(base)
            coeff_type = elem_type(base)
            loaded_terms[exponent] = load_type_dispatch(s, coeff_type, coeff,
                                                            parent=base)
        else
            loaded_terms[exponent] = load_terms(s, parents[1:end - 1], coeff, parents[end - 1])
        end
    end
    
    return loaded_terms
end

function load_internal(s::DeserializerState,
                       T::Type{<: Union{
                           Generic.LaurentSeriesElem,                           
                           ZZLaurentSeriesRingElem}},
                       dict::Dict)
    parents = load_parents(s, dict[:parents])
    pol_length = load_type_dispatch(s, Int, dict[:pol_length])
    precision = load_type_dispatch(s, Int, dict[:precision])
    valuation = load_type_dispatch(s, Int, dict[:valuation])
    scale = load_type_dispatch(s, Int, dict[:scale])
    terms = load_terms(s, parents, dict[:terms], parents[end])
    parent_ring = parents[end]
    return parent_ring(terms, pol_length, precision, valuation, scale)
end

function load_internal_with_parent(s::DeserializerState,
                                   T::Type{<: Union{Generic.LaurentSeriesElem, ZZLaurentSeriesRingElem}},
                                   dict::Dict,
                                   parent_ring::Union{Generic.LaurentSeriesRing,
                                                      Generic.LaurentSeriesField,
                                                      ZZLaurentSeriesRing})
    parents = get_parents(parent_ring)
    terms = load_terms(s, parents, dict[:terms], parents[end])
    pol_length = load_type_dispatch(s, Int, dict[:pol_length])
    precision = load_type_dispatch(s, Int, dict[:precision])
    valuation = load_type_dispatch(s, Int, dict[:valuation])
    scale = load_type_dispatch(s, Int, dict[:scale])

    if precision > max_precision(parent_ring)
        @warn("Precision Warning: given parent is less precise than serialized elem",
              maxlog=1)
    end

    return parent_ring(terms, pol_length, precision, valuation, scale)
end

