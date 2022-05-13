################################################################################
#  non simpleton base rings

function load_internal(s::DeserializerState, t::Type{Nemo.NmodRing}, dict::Dict)
    modulus = load_type_dispatch(s, UInt64, dict[:modulus])
    
    return t(modulus)
end

function save_internal(s::SerializerState, R::Nemo.NmodRing)
    return Dict(
        :modulus => save_type_dispatch(s, modulus(R))
    )
end

################################################################################
#  Polynomial Rings
function save_internal(s::SerializerState, R::Union{MPolyRing, PolyRing})
    return Dict(
        :symbols => save_type_dispatch(s, symbols(R)),
        :base_ring => save_type_dispatch(s, base_ring(R)),
    )
end

function load_internal(s::DeserializerState,
                       T::Type{<: Union{MPolyRing, PolyRing}},
                       dict::Dict)
    base_ring = load_type_dispatch(s, dict[:base_ring], check_namespace=false)
    symbols = load_type_dispatch(s, Vector{Symbol}, dict[:symbols]) 

    if T <: PolyRing
        return PolynomialRing(base_ring, symbols...)
    end

    return PolynomialRing(base_ring, symbols)
end

################################################################################
# Multivariate Polynomials
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
    R, symbols = load_type_dispatch(s, dict[:parent], check_namespace=false)
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

################################################################################
# Univariate Polynomials
function save_internal(s::SerializerState, p::PolyElem)
    parent_ring = parent(p)
    parent_ring = save_type_dispatch(s, parent_ring)

    return Dict(
        :parent => parent_ring,
        :coeffs => save_type_dispatch(s, collect(coefficients(p)))
    )
end

function load_internal(s::DeserializerState, ::Type{<: PolyElem}, dict::Dict)
    R, y = load_type_dispatch(s, dict[:parent], check_namespace=false)
    coeff_ring = coefficient_ring(R)
    coeff_type = elem_type(coeff_ring)
    coeffs = load_type_dispatch(s, Vector{coeff_type}, dict[:coeffs])

    return polynomial(coeff_ring, coeffs, String(symbols(R)[1]))
end

################################################################################
# Polynomial Ideals
function save_internal(s::SerializerState, i::MPolyIdeal)
    generators = gens(i)
    parent_ring = save_type_dispatch(s, parent(generators[1]))
    
    return Dict(
        :parent => parent_ring,
        :gens => save_type_dispatch(s, generators),
    )
end
                       
function load_internal(s::DeserializerState, ::Type{<: MPolyIdeal}, dict::Dict)
    parent_ring, _ = load_type_dispatch(s, dict[:parent], check_namespace=false)
    gens = load_type_dispatch(s, Vector{MPolyElem}, dict[:gens])

    return ideal(parent_ring, gens)
end
