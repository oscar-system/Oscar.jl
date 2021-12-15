@doc Markdown.doc"""
    tropical_polynomial{M}(f::AbstractAlgebra.Generic.MPoly{Any})

Returns the tropicalization of a polynomial in form of a polynomial over the tropical numbers.
If M=min, returns a polynomial over the min-plus semiring.
If M=max, returns a polynomial over the max-plus semiring.

# Examples
```jldoctest
julia> K = PadicField(7, 2);

julia> Kxy, (x,y) = K["x", "y"];

julia> f = 7*x+y+49;

julia> tropical_polynomial(f,min)
```
"""
function tropical_polynomial(f::AbstractAlgebra.Generic.MPoly{<:RingElement}, M::Union{typeof(min),typeof(max)})
    T = tropical_numbers(M)
    Tx = PolynomialRing(T,[repr(x) for x in gens(parent(f))])
    tropf = inf(T)

    if base_ring(parent(f)) isa NonArchLocalField
        for expv in exponent_vectors(f)
            tropf = tropf + T(valuation(coeff(f,expv)))*monomial(Tx[1],expv)
        end
    else
        for expv in exponent_vectors(f)
            tropf = tropf + T(0)*monomial(Tx[1],expv)
        end
    end

    return tropf
end
export tropical_polynomial
