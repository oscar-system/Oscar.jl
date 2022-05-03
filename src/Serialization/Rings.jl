################################################################################
#  base ring
SimpleRingType = Union{FlintRationalField, FlintIntegerRing}

function load_internal(s::SerializerState, t::Type{SimpleRingType}, dict::Dict)
    return t()
end

function load_internal(s::SerializerState, t::Type{Nemo.NmodRing}, dict::Dict)
    return t(dict[:data][:n][:data])
end

################################################################################
#  polynomial ring
function save_internal(s::SerializerState, R::MPolyRing)
    base_ring = R.base_ring
    return Dict(
        :symbols => save_type_dispatch(s, R.S),
        :base_ring => save_type_dispatch(s, base_ring),
    )
end

function load_internal(s::DeserializerState, ::Type{<: MPolyRing}, dict::Dict)
    ring_type = decodeType(dict[:base_ring][:type])
    base_ring = load_type_dispatch(s, ring_type, dict[:base_ring])
    symbols = load_type_dispatch(s, Vector{Symbol}, dict[:symbols]) 
    
    return PolynomialRing(base_ring, symbols)
end

################################################################################
# Polynomials
function save_internal(s::SerializerState, p::MPolyElem)
    parent_ring = parent(p)
    parent_ring = save_type_dispatch(s, parent_ring)
    terms = []
    
    for (i, c) in enumerate(coeffs(p))
        term = Dict(
            :exponent => save_type_dispatch(s, exponent_vector(p, i)),
            :coeff => save_type_dispatch(s, c)
        )
        push!(terms, term)
    end

    return Dict(
        :terms => terms,
        :parent => parent_ring,
    )
end

function load_internal(s::DeserializerState, ::Type{<: MPolyElem}, dict::Dict)
    ring_dict = dict[:parent]
    ring_type = decodeType(ring_dict[:type])
    R, symbols = load_type_dispatch(s, ring_type, ring_dict)
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


