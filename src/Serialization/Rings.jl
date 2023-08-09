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

function save_parents(s::SerializerState, parent_ring::Ring)
    parents = get_parents(parent_ring)
    refs = []
    for p in parents
        push!(refs, save_as_ref(s, p))
    end
    return refs
end

################################################################################
# Handling RingElem Params

function save_type_params(s::SerializerState, x::T, key::Symbol) where T <: RingElem
    s.key = key
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

function load_type_params(s::DeserializerState, ::Type{<:RingElem}, dict::Dict{Symbol, Any})
    return load_typed_object(s, dict)
end

function load_type_params(s::DeserializerState, ::Type{<:RingElem}, refs::Vector{Any})
    return load_parents(s, refs)
end

function load_type_params(s::DeserializerState, ::Type{<:RingElem}, parent_ring::Ring) 
    return get_parents(parent_ring)
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

function load_object_with_params(s::DeserializerState, ::Type{zzModRingElem},
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
        save_type_object(s, base_ring(R), :base_ring)
        save_object(s, symbols(x), :symbols)
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
@registerSerializationType(PolyRingElem)

PolyElemUniontype = Union{UniversalPolyRingElem, MPolyRingElem, PolyRingElem}
type_needs_params(::Type{<:PolyElemUniontype}) = true

# elements
function save_object(s::SerializerState, p::Union{UniversalPolyRingElem, MPolyRingElem})
    coeff_type = typeof(coeff(p, 1))
    terms = Tuple{Vector{UInt}, coeff_type}[]
    for i in 1:length(p)
        push!(terms, (exponent_vector(p, i), coeff(p, i)))
    end
    save_object(s, terms)
end

@registerSerializationType(AbstractAlgebra.Generic.LaurentMPolyWrap)
function save_object(s::SerializerState, p::AbstractAlgebra.Generic.LaurentMPolyWrap)
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

function load_object_with_params(s::DeserializerState,
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
            loaded_terms[exponent] = load_object_with_params(s, coeff_type, coeff, params)
        else
            loaded_terms[exponent] = load_object(s, coeff_type, coeff)
        end
    end
    return parent_ring(loaded_terms)
end


function load_object_with_params(s::DeserializerState,
                                 ::Type{<:Union{MPolyRing, UniversalPolyRing, AbstractAlgebra.Generic.LaurentMPolyWrapRing}},
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
            c = load_object_with_params(s, coeff_type, coeff, params)
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
IdealUnionType = Union{MPolyIdeal, Laurent.LaurentMPolyIdeal}
type_needs_params(::Type{<: IdealUnionType) = true

function save_type_params(s::SerializerState, x::T, key::Symbol) where T <: IdealUnionType
    s.key = key
    data_dict(s) do
        save_object(s, encode_type(T), :name)
        save_typed_object(s, parent(gens(x)[1]), :params)
    end
end

function save_object(s::SerializerState, I::MPolyIdeal)
    data_dict(s) do
        save_typed_object(s, gens(I), :gens)
    end
end

function load_object(s::DeserializerState, ::Type{<: IdealUnionType}, dict::Dict)
    gens = load_typed_object(s, dict[:gens])
    return ideal(parent(gens[1]), gens)
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

    
