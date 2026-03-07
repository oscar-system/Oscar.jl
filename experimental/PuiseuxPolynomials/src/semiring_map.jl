#################################################################################
#
# tropical semiring maps for puiseux polynomials
#
#################################################################################


@doc raw"""
    tropical_semiring_map(K::Field, minOrMax::Union{typeof(min),typeof(max)}=min)

Return a map `nu` from `K` to the min (default) or max tropical semiring `T` such that `nu(0)=zero(T)` and `nu(c)=one(T)` for `c` non-zero.  In other words, `nu` extends the trivial valuation on `K`.

# Examples
```jldoctest
julia> nu = tropical_semiring_map(QQ) # arbitrary rings possible
Map into Min tropical semiring encoding the trivial valuation on Rational field

julia> nu(1)
(0)

julia> nu(0)
infty

julia> nu = tropical_semiring_map(QQ,max) # arbitrary rings possible
Map into Max tropical semiring encoding the trivial valuation on Rational field

julia> nu(1)
(0)

julia> nu(0)
-infty

```
"""
function tropical_semiring_map(R::PuiseuxMPolyRing, t::PuiseuxMPolyRingElem, minOrMax::Union{typeof(min),typeof(max)}=min)
    @req is_univariate(R) "The tropical_semiring_map is only implemented for Puiseux polynomial rings in one variable."
    @req is_gen(t) "The uniformizer must be the generator of the Puiseux polynomial ring."

    valuedField = nothing
    uniformizerField = nothing
    valuedRing = R
    uniformizerRing = t
    residueField = base_ring(R)
    tropicalSemiring = tropical_semiring(minOrMax)
    return TropicalSemiringMap{typeof(R),typeof(t),typeof(minOrMax)}(
        valuedField,
        uniformizerField,
        valuedRing,
        uniformizerRing,
        residueField,
        tropicalSemiring)
end

# Print string
function Base.show(io::IO, nu::TropicalSemiringMap{<:PuiseuxMPolyRing,<:PuiseuxMPolyRingElem,MinOrMax}) where {MinOrMax<:Union{typeof(min),typeof(max)}}
    print(io, "Map into $(tropical_semiring(nu)) encoding the $(uniformizer(nu))-adic valuation on $(Oscar.valued_ring(nu))")
end

# Mapping an element of the valued field or ring to the tropical semiring
function (nu::TropicalSemiringMap{<:PuiseuxMPolyRing,<:PuiseuxMPolyRingElem,MinOrMax})(c::PuiseuxMPolyRingElem) where MinOrMax<:Union{typeof(min),typeof(max)}
    # return tropical zero if c is zero and valuation otherwise
    iszero(c) && return zero(tropical_semiring(nu))
    # preserve_ordering ensures that valuation is negated if convention(nu)==max
    return tropical_semiring(nu)(valuation(c); preserve_ordering=true)
end

# Mapping an element of the valued field or ring to the residue field
function initial(c::PuiseuxMPolyRingElem, nu::TropicalSemiringMap{<:PuiseuxMPolyRing,<:PuiseuxMPolyRingElem,MinOrMax}) where MinOrMax<:Union{typeof(min),typeof(max)}
    # return residue field zero if c is zero and the correct non-zero residue otherwise
    iszero(c) && return zero(residue_field(nu))
    C = collect(coefficients(c))
    return last(C)
end
