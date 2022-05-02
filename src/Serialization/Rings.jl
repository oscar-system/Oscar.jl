################################################################################
#  polynomial ring
function save_internal(s::SerializerState, R::MPolyRing)
    base_ring = R.base_ring
    return Dict(
        :symbols => save_type_dispatch(s, R.S),
        :base_ring => save_type_dispatch(s, base_ring),
        :nvars => save_type_dispatch(s, R.nvars),
        :nfields => save_type_dispatch(s, R.nfields),
    )
end

function load_internal(s::DeserializerState, ::Type{MPolyRing}, dict::Dict)
    ring_type = decodeType(dict[:base_ring][:type])
    base_ring = ring_type()
    symbols = load_type_dispatch(s, Vector{Symbol}, dict[:symbols]) 
    
    return PolynomialRing(base_ring, symbols)
end

################################################################################
# Polynomials
function save_internal(s::SerializerState, p::MPolyElem)
    ring = parent(p)
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
        :parent => save_type_dispatch(s, ring),
        :content_den => save_type_dispatch(s, p.content_den),
        :content_num => save_type_dispatch(s, p.content_num),
    )
end

function load_internal(s::DeserializerState, ::Type{MPolyElem}, dict::Dict)
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


