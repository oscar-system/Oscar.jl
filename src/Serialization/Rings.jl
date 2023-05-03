################################################################################
# ring of integers (singleton type)
@registerSerializationType(ZZRing)


################################################################################
#  non simpleton base rings
@registerSerializationType(Nemo.zzModRing,
                           "Nemo.zzModRing")

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

function save_internal(s::SerializerState, r::zzModRingElem)
    return string(r)
end

function load_internal(s::DeserializerState, ::Type{zzModRingElem}, dict::Dict)
    parent_ring = load_type_dispatch(s, Nemo.zzModRing, dict[:parent])
    class_val = load_type_dispatch(s, ZZRingElem, dict[:class_val])
    return parent_ring(class_val)
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

function save_internal(s::SerializerState, R::Union{UniversalPolyRing, MPolyRing, PolyRing})
    base = base_ring(R)

    if has_elem_basic_encoding(base)
        base_dict = save_type_dispatch(s, base)
    else
        base_dict = save_as_ref(s, base)
    end
    return Dict(
        :symbols => save_type_dispatch(s, symbols(R)),
        :base_ring => base_dict,
    )
end

function load_internal(s::DeserializerState,
                       T::Type{<: Union{UniversalPolyRing, MPolyRing, PolyRing}},
                       dict::Dict)
    base_ring = load_unknown_type(s, dict[:base_ring])
    symbols = load_type_dispatch(s, Vector{Symbol}, dict[:symbols])
    
    if T <: PolyRing
        return polynomial_ring(base_ring, symbols..., cached=false)[1]
    elseif T <: UniversalPolyRing
        poly_ring = UniversalPolynomialRing(base_ring, cached=false)
        gens(poly_ring, symbols)
        return poly_ring
    end

    return polynomial_ring(base_ring, symbols, cached=false)[1]
end

################################################################################
# Multivariate and Universal Polynomials
@registerSerializationType(MPolyRingElem)
@registerSerializationType(UniversalPolyRingElem)

function save_internal(s::SerializerState, p::Union{UniversalPolyRingElem, MPolyRingElem})
    parent_ring = parent(p)
    base = base_ring(parent_ring)
    parents = []
    encoded_terms = []

    for i in 1:length(p)
        if coeff == parent_ring(0)
            continue
        end
        encoded_coeff = save_internal(s, coeff(p, i))
        if has_elem_basic_encoding(base)
            push!(encoded_terms,  (exponent_vector(p, i), encoded_coeff))
        else
            parents = encoded_coeff[:parents]
            push!(encoded_terms, (exponent_vector(p, i), encoded_coeff[:terms]))
        end
    end
    parent_ring = save_as_ref(s, parent_ring)
    # end of list should be loaded first
    push!(parents, parent_ring)

    return Dict(
        :terms => encoded_terms,
        :parents => parents,
    )
end


################################################################################
# Univariate Polynomials

@registerSerializationType(PolyRingElem)

function save_internal(s::SerializerState, p::PolyRingElem)
    # store a polynomial over a ring provided we can store elements in that ring
    parent_ring = parent(p)
    base = base_ring(parent_ring)
    coeffs = coefficients(p)
    encoded_terms = []
    # flattened tree required for coefficients having parents
    parents = []
    exponent = 0
    for coeff in coeffs
        # collect not trivial terms
        if is_zero(coeff)
            continue
        end

        encoded_coeff = save_internal(s, coeff)
        if has_elem_basic_encoding(base)
            push!(encoded_terms,  (exponent, encoded_coeff))
        else
            parents = encoded_coeff[:parents]
            # there should not be a type check within such a function
            if typeof(base) <: Union{AbstractAlgebra.Generic.RationalFunctionField,
                             FracField}
                push!(encoded_terms,  (exponent, (encoded_coeff[:num_terms],
                                                  encoded_coeff[:den_terms])))
            else
                push!(encoded_terms,  (exponent, encoded_coeff[:terms]))
            end
        end
        exponent += 1
    end
    parent_ring = save_as_ref(s, parent_ring)
    # end of list should be loaded first
    push!(parents, parent_ring)
    
    return Dict(
        :parents => parents,
        :terms =>  encoded_terms
    )
end

function load_terms(s::DeserializerState, parents::Vector, terms::Vector,
                    parent_ring::PolyRing)
    # shift so constant starts at 1
    degree = max(map(x->x[1] + 1, terms)...)
    base = base_ring(parent_ring)
    loaded_terms = zeros(base, degree)
    
    for term in terms
        exponent, coeff = term
        # account for shift
        exponent += 1

        if length(parents) == 1
            coeff_type = elem_type(base)
            if has_elem_basic_encoding(base)
                loaded_terms[exponent] = load_type_dispatch(s, coeff_type, coeff,
                                                            parent=base)
            else
                loaded_terms[exponent] = base(coeff)
            end
        else
            loaded_terms[exponent] = load_terms(s, parents[1:end - 1], coeff, parents[end - 1])
        end
    end
    return parent_ring(loaded_terms)
end

function load_terms(s::DeserializerState, parents::Vector, terms::Vector,
                    parent_ring::MPolyRing)
    base = base_ring(parent_ring)
    polynomial = MPolyBuildCtx(parent_ring)
    for term in terms
        e, coeff = term
        if length(parents) == 1
            coeff_type = elem_type(base)
            if has_elem_basic_encoding(base)
                c = load_type_dispatch(s, coeff_type, coeff, parent=base)
            else
                c = base(c)
            end
        else
            c = load_terms(s, parents[1:end - 1], coeff, parents[end - 1])
        end
        
        push_term!(polynomial, c, Vector{Int}(e))
    end
    return finish(polynomial)
end

function load_parents(s::DeserializerState, parents::Vector)
    parent_ids = [parent[:id] for parent in parents]
    loaded_parents = []
    
    for id in parent_ids
        if haskey(s.objs, UUID(id))
            loaded_parent = s.objs[UUID(id)]
        else
            parent_dict = s.refs[Symbol(id)]
            parent_dict[:id] = id
            loaded_parent = load_unknown_type(s, parent_dict)
        end            
        push!(loaded_parents, loaded_parent)
    end
    return loaded_parents
end
function load_internal(s::DeserializerState, ::Type{<: Union{
    PolyRingElem, UniversalPolyRingElem, MPolyRingElem}}, dict::Dict)
    loaded_parents = load_parents(s, dict[:parents])
    return load_terms(s, loaded_parents, dict[:terms], loaded_parents[end])
end

function load_internal_with_parent(s::DeserializerState, ::Type{<: Union{
    PolyRingElem, UniversalPolyRingElem, MPolyRingElem}}, dict::Dict,
                                   parent_ring::Union{PolyRing, MPolyRing, UniversalPolyRing})
    parents = Any[parent_ring]
    base = base_ring(parent_ring)
    
    while (!has_elem_basic_encoding(base))
        if base isa Field && !(base isa FracField)
            if absolute_degree(base) == 1
                break
            end
            if typeof(base) <: Union{NfAbsNS, NfRelNS}
                pushfirst!(parents, base)
                ngens = length(gens(base))
                base = polynomial_ring(base_field(base), ngens)[1]
            else
                pushfirst!(parents, base)
                base = parent(defining_polynomial(base))
            end
            continue
        end
        pushfirst!(parents, base)
        base = base_ring(base)
    end
    return load_terms(s, parents, dict[:terms], parents[end])
end

################################################################################
# Polynomial Ideals

@registerSerializationType(MPolyIdeal)

function save_internal(s::SerializerState, I::MPolyIdeal)
    generators = gens(I)
    parent_ring = save_type_dispatch(s, base_ring(I))

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
    gens = load_type_dispatch(s, Vector{elem_type(parent_ring)},
                              dict[:gens], parent=parent_ring)
    return ideal(parent_ring, gens)
end

################################################################################
# Matrices
@registerSerializationType(MatElem)

function save_internal(s::SerializerState, m::MatrixElem)
    return Dict(
        :base_ring => save_type_dispatch(s, base_ring(parent(m))),
        :matrix => save_type_dispatch(s, Array(m)),
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
@registerSerializationType(SeriesRing, true)

function save_internal(s::SerializerState, R::Union{
    Generic.RelPowerSeriesRing,
    QQRelPowerSeriesRing,
    ZZRelPowerSeriesRing,
    fqPolyRepRelPowerSeriesRing,
    zzModRelPowerSeriesRing})
    return Dict(
        :base_ring => save_type_dispatch(s, base_ring(R)),
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
    
    return power_series_ring(base_ring, max_precision, var; cached=false, model=model)[1]
end

# elements
@registerSerializationType(RelPowerSeriesRingElem)
@registerSerializationType(AbsPowerSeriesRingElem)

function save_internal(s::SerializerState, r::RelPowerSeriesRingElem)
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

function save_internal(s::SerializerState, r::AbsPowerSeriesRingElem)
    coeffs = map(x -> coeff(r, x), 0:pol_length(r))
    return Dict(
        :parent => save_type_dispatch(s, parent(r)),
        :coeffs => save_type_dispatch(s, coeffs),
        :pol_length => save_type_dispatch(s, pol_length(r)),
        :precision => save_type_dispatch(s, precision(r))
    )
end

function load_internal(s::DeserializerState, ::Type{<: RelPowerSeriesRingElem}, dict::Dict)
    parent = load_type_dispatch(s, SeriesRing, dict[:parent])
    coeffs = load_type_dispatch(s, Vector, dict[:coeffs])
    valuation = load_type_dispatch(s, Int, dict[:valuation])
    pol_length = load_type_dispatch(s, Int, dict[:pol_length])
    precision = load_type_dispatch(s, Int, dict[:precision])
    
    return parent(coeffs, pol_length, precision, valuation)
end

function load_internal_with_parent(s::DeserializerState,
                                   ::Type{<: RelPowerSeriesRingElem},
                                   dict::Dict,
                                   parent_ring::SeriesRing)
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

function load_internal(s::DeserializerState, ::Type{<: AbsPowerSeriesRingElem}, dict::Dict)
    parent = load_type_dispatch(s, SeriesRing, dict[:parent])
    coeffs = load_type_dispatch(s, Vector, dict[:coeffs])
    pol_length = load_type_dispatch(s, Int, dict[:pol_length])
    precision = load_type_dispatch(s, Int, dict[:precision])
    
    return parent(coeffs, pol_length, precision)
end

function load_internal_with_parent(s::DeserializerState,
                                   ::Type{<: AbsPowerSeriesRingElem},
                                   dict::Dict,
                                   parent_ring::SeriesRing)
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
@registerSerializationType(Generic.LaurentSeriesRing, true, "LaurentSeriesRing")
@registerSerializationType(Generic.LaurentSeriesField, true, "LaurentSeriesField")
@registerSerializationType(ZZLaurentSeriesRing)

function save_internal(s::SerializerState, R::Union{
    Generic.LaurentSeriesRing,
    Generic.LaurentSeriesField,
    ZZLaurentSeriesRing})
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

function save_internal(s::SerializerState, r:: ZZLaurentSeriesRingElem)
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
                           ZZLaurentSeriesRingElem}},
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
                                       Generic.LaurentSeriesElem, ZZLaurentSeriesRingElem}},
                                   dict::Dict,
                                   parent_ring)
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

