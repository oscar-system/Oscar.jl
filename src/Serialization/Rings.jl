################################################################################
# fmpq polynomial ring
function save_internal(s::SerializerState, R::FmpqMPolyRing)
    base_ring = R.base_ring
    return Dict(
        :symbols => save_type_dispatch(s, [String(x) for x in R.S]),
        :base_ring => save_type_dispatch(s, base_ring),
        :nvars => save_type_dispatch(s, R.nvars),
        :nfields => save_type_dispatch(s, R.nfields),
    )
end

function load_internal(s::DeserializerState, ::Type{FmpqMPolyRing}, dict::Dict)
    base_ring = eval(Meta.parse(dict[:base_ring][:type] * "()"))
    symbols = [x[:data] for x in dict[:symbols][:data][:vector]]
    
    return PolynomialRing(base_ring, symbols)
end

################################################################################
# Polynomials
function save_internal(s::SerializerState, p::fmpq_mpoly)
    ring = parent(p)
    monomials = []
    
    for (i, c) in enumerate(coeffs(p))
        monomial = Dict(
            :exponent => exponent_vector(p, i),
            :coeff => c
        )
        push!(monomials, monomial)
    end

    return Dict(
        :monomials => monomials,
        :parent => save_type_dispatch(s, ring)
    )
end

function load_internal(s::DeserializerState, ::Type{fmpq_mpoly}, dict::Dict)
    ring_dict = dict[:parent]
    ring_type = eval(Meta.parse(ring_dict[:type]))
    R, symbols = load_type_dispatch(s, ring_type, ring_dict)
    polynomial = R(0)

    for monomial in dict[:monomials]
        c = monomial[:coeff][:num] // monomial[:coeff][:den]
        term = reduce((a, (s, e)) -> a * s^e,
                      zip(symbols, monomial[:exponent]),
                      init=c)
        polynomial += term
    end
    return polynomial
end


