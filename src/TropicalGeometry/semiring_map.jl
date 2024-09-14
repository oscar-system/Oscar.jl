################################################################################
#
#  Tropical semiring maps
#  ======================
#  maps from a field K to a tropical semiring T with the purpose of encoding
#    - a valuation on K
#    - a choice of min- or max-convention
#
#  also collects dictionaries that assign to a polynomial ring over K
#    - a polynomial ring over the ring of integers with an extra variable
#      (for the computation of tropical groebner bases)
#    - a polynomial ring over the residue field
#      (for initial forms and initial ideals)
#
################################################################################
struct TropicalSemiringMap{typeofValuedField,typeofUniformizer,minOrMax}
    valued_field::typeofValuedField
    uniformizer_field::Union{Nothing,FieldElem}
    valued_ring::Ring
    uniformizer_ring::typeofUniformizer
    residue_field::Field
    tropical_semiring::TropicalSemiring{minOrMax}
    ###
    # Dictionaries for polynomial ring caches consisting of pairs
    #   R => S
    # where
    # - R is a polynomial ring over valued_field
    # - S is a polynomial ring over valued_ring or residue_field
    #   (with the same variables as R)
    ###
    polynomial_rings_for_groebner::Dict{MPolyRing,MPolyRing}
    polynomial_rings_for_initial::Dict{MPolyRing,MPolyRing}

    # Constructor with empty dictionaries
    function TropicalSemiringMap{typeofValuedField,typeofUniformizer,minOrMax}(valuedField::typeofValuedField,uniformizerField::Union{Nothing,FieldElem},valuedRing::Ring,uniformizerRing::typeofUniformizer,residueField::Field,tropicalSemiring::TropicalSemiring{minOrMax}) where {typeofValuedField<:Field,typeofUniformizer<:Union{Nothing,RingElem},minOrMax<:Union{typeof(min),typeof(max)}}
        return new{typeofValuedField,typeofUniformizer,minOrMax}(valuedField,uniformizerField,valuedRing,uniformizerRing,residueField,tropicalSemiring,Dict{MPolyRing,MPolyRing}(),Dict{MPolyRing,MPolyRing}())
    end
end



################################################################################
#
#  Properties
#
################################################################################
valued_field(nu::TropicalSemiringMap) = nu.valued_field
uniformizer_field(nu::TropicalSemiringMap) = nu.uniformizer_field
valued_ring(nu::TropicalSemiringMap) = nu.valued_ring
uniformizer_ring(nu::TropicalSemiringMap) = nu.uniformizer_ring
uniformizer(nu::TropicalSemiringMap) = uniformizer_ring(nu)
residue_field(nu::TropicalSemiringMap) = nu.residue_field
tropical_semiring(nu::TropicalSemiringMap) = nu.tropical_semiring

convention(::TropicalSemiringMap{typeofValuedField,typeofUniformizer,typeof(min)}) where {typeofValuedField,typeofUniformizer} = min
convention(::TropicalSemiringMap{typeofValuedField,typeofUniformizer,typeof(max)}) where {typeofValuedField,typeofUniformizer} = max

is_trivial(::TropicalSemiringMap{K,Nothing,minOrMax}) where {K,minOrMax<:Union{typeof(min),typeof(max)}} = true
is_trivial(::TropicalSemiringMap) = false

polynomial_rings_for_groebner(nu::TropicalSemiringMap) = nu.polynomial_rings_for_groebner
polynomial_rings_for_initial(nu::TropicalSemiringMap) = nu.polynomial_rings_for_initial



################################################################################
#
#  Trivial valuation on any field
#
################################################################################

# Constructor:
@doc raw"""
    tropical_semiring_map(K::Field, minOrMax::Union{typeof(min),typeof(max)}=min)

Return a map `nu` from `K` to the min (default) or max tropical semiring `T` such that `nu(0)=zero(T)` and `nu(c)=one(T)` for `c` non-zero.  In other words, `nu` extends the trivial valuation on `K`.

# Example
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
function tropical_semiring_map(K::Field, minOrMax::Union{typeof(min),typeof(max)}=min)
    valuedField = K
    uniformizerField = nothing
    valuedRing = K
    uniformizerRing = nothing
    residueField = K
    tropicalSemiring = tropical_semiring(minOrMax)
    return TropicalSemiringMap{typeof(valuedField),typeof(uniformizerRing),typeof(minOrMax)}(valuedField,uniformizerField,valuedRing,uniformizerRing,residueField,tropicalSemiring)
end

# Print string
function Base.show(io::IO, nu::TropicalSemiringMap{K,Nothing,minOrMax} where {K<:Ring, minOrMax<:Union{typeof(min),typeof(max)}})
    print(io, "Map into $(tropical_semiring(nu)) encoding the trivial valuation on $(valued_field(nu))")
end

# Mapping an element of the valued field or ring to the tropical semiring
function (nu::TropicalSemiringMap{K,Nothing,minOrMax})(c::Union{RingElem,Integer,Rational}) where {K<:Field, minOrMax<:Union{typeof(min),typeof(max)}}
    # return tropical zero if c is zero
    #   and tropical one otherwise
    return (iszero(valued_field(nu)(c)) ? zero(tropical_semiring(nu)) : one(tropical_semiring(nu)))
end

# Mapping an element of the valued field or ring to the residue field
function initial(c::Union{RingElem,Integer,Rational}, nu::TropicalSemiringMap{K,Nothing,minOrMax}) where {K<:Field, minOrMax<:Union{typeof(min),typeof(max)}}
    return residue_field(nu)(c)
end



################################################################################
#
#  p-adic valuation on QQ
#
################################################################################

# Constructor:
@doc raw"""
    tropical_semiring_map(QQ::QQField, p::QQFieldElem, minOrMax::Union{typeof(min),typeof(max)}=min)

Return a map `nu` from `QQ` to the min (default) or max tropical semiring `T` such that `nu(0)=zero(T)` and `nu(c)=+/-val(c)` for `c` non-zero, where `val` denotes the `p`-adic valuation.  Requires `p` to be a prime.

# Example
```jldoctest
julia> nu_2 = tropical_semiring_map(QQ,2)
Map into Min tropical semiring encoding the 2-adic valuation on Rational field

julia> nu_2(4)
(2)

julia> nu_2(1//4)
(-2)

julia> nu_2 = tropical_semiring_map(QQ,2,max);

julia> nu_2(4)
(-2)

julia> nu_2(1//4)
(2)

```
"""
function tropical_semiring_map(::QQField, p::Union{RingElem,Integer,Rational}, minOrMax::Union{typeof(min),typeof(max)}=min)
    valuedField = QQ
    uniformizerField = QQ(p)
    valuedRing = ZZ
    uniformizerRing = ZZ(p)
    residueField = GF(uniformizerRing)
    tropicalSemiring = tropical_semiring(minOrMax)
    return TropicalSemiringMap{typeof(valuedField),typeof(uniformizerRing),typeof(minOrMax)}(valuedField,uniformizerField,valuedRing,uniformizerRing,residueField,tropicalSemiring)
end

# Print string
function Base.show(io::IO, nu::TropicalSemiringMap{QQField,ZZRingElem,minOrMax}) where {minOrMax<:Union{typeof(min),typeof(max)}}
    print(io, "Map into $(tropical_semiring(nu)) encoding the $(uniformizer(nu))-adic valuation on $(valued_field(nu))")
end

# Mapping an element of the valued field or ring to the tropical semiring
function (nu::TropicalSemiringMap{QQField,ZZRingElem,minOrMax})(c::Union{RingElem,Integer,Rational}) where minOrMax<:Union{typeof(min),typeof(max)}
    c = valued_field(nu)(c)
    # return tropical zero if c is zero
    #   and p-adic valuation otherwise
    # preserve_ordering ensures that valuation is negated if convention(nu)==max
    iszero(c) && return zero(tropical_semiring(nu))
    return tropical_semiring(nu)(valuation(c,uniformizer(nu)); preserve_ordering=true)
end

# Mapping an element of the valued field or ring to the residue field
function initial(c::Union{RingElem,Integer,Rational}, nu::TropicalSemiringMap{QQField,ZZRingElem,minOrMax}) where minOrMax<:Union{typeof(min),typeof(max)}
    c = valued_field(nu)(c)
    # return residue field zero if c is zero
    #   and the correct non-zero residue otherwise
    iszero(c) && return zero(residue_field(nu)) # if c is zero, return 0
    c *= uniformizer_field(nu)^(-valuation(c,uniformizer(nu)))
    return residue_field(nu)(c)
end



################################################################################
#
#  t-adic valuation on K(t)
#
################################################################################

# Constructor:
@doc raw"""
    tropical_semiring_map(Kt::Generic.RationalFunctionField, t::Generic.RationalFunctionFieldElem, minOrMax::Union{typeof(min),typeof(max)}=min)

Return a map `nu` from rational function field `Kt` to the min (default) or max tropical semiring `T` such that `nu(0)=zero(T)` and `nu(c)=+/-val(c)` for `c` non-zero, where `val` denotes the `t`-adic valuation with uniformizer `t`.  Requires `t` to be non-constant and have denominator `1`.

# Example
```jldoctest
julia> Kt,t = rational_function_field(QQ,"t");

julia> nu_t = tropical_semiring_map(Kt,t)
Map into Min tropical semiring encoding the t-adic valuation on Rational function field over QQ

julia> nu_t(t^2)
(2)

julia> nu_t(1//t^2)
(-2)

julia> nu_t = tropical_semiring_map(Kt,t,max)
Map into Max tropical semiring encoding the t-adic valuation on Rational function field over QQ

julia> nu_t(t^2)
(-2)

julia> nu_t(1//t^2)
(2)

```
"""
function tropical_semiring_map(Kt::Generic.RationalFunctionField, t::Generic.RationalFunctionFieldElem, minOrMax::Union{typeof(min),typeof(max)}=min)
    valuedField = Kt
    uniformizerField = t
    uniformizerRing = numerator(t)
    valuedRing = parent(uniformizerRing)
    residueField = base_ring(valuedRing)
    tropicalSemiring = tropical_semiring(minOrMax)
    @req is_unit(denominator(uniformizerField)) "input uniformizer denominator not unit"
    @req degree(uniformizerRing)>0 "input uniformizer constant"
    return TropicalSemiringMap{typeof(valuedField),typeof(uniformizerRing),typeof(minOrMax)}(valuedField,uniformizerField,valuedRing,uniformizerRing,residueField,tropicalSemiring)
end

# Print String
function Base.show(io::IO, nu::TropicalSemiringMap{Kt,t,minOrMax}) where {Kt<:Generic.RationalFunctionField, t<:PolyRingElem, minOrMax<:Union{typeof(min),typeof(max)}}
    print(io, "Map into $(tropical_semiring(nu)) encoding the $(uniformizer(nu))-adic valuation on $(valued_field(nu))")
end

# Mapping an element of the valued field or ring to the tropical semiring
function t_adic_valuation(c::Generic.RationalFunctionFieldElem, t::PolyRingElem)
    return valuation(numerator(c),t)-valuation(denominator(c),t)
end

function (nu::TropicalSemiringMap{Kt,t,minOrMax})(c::Union{RingElem,Integer,Rational}) where {Kt<:Generic.RationalFunctionField, t<:PolyRingElem, minOrMax<:Union{typeof(min),typeof(max)}}
    c = valued_field(nu)(c)
    # return tropical zero if c is zero
    #   and p-adic valuation otherwise
    # preserve_ordering ensures that valuation is negated if convention(nu)==max
    iszero(c) && return zero(tropical_semiring(nu)) # if c zero, return tropical zero
    return tropical_semiring(nu)(t_adic_valuation(c,uniformizer(nu)); preserve_ordering=true)
end

# Mapping an element of the valued field or ring to the residue field
function initial(c::Union{RingElem,Integer,Rational}, nu::TropicalSemiringMap{Kt,t,minOrMax}) where {Kt<:Generic.RationalFunctionField, t<:PolyRingElem, minOrMax<:Union{typeof(min),typeof(max)}}
    c = valued_field(nu)(c)
    # return residue field zero if c is zero
    #   and the correct non-zero residue otherwise
    iszero(c) && return zero(residue_field(nu)) # if c is zero, return 0
    c *= uniformizer_field(nu)^(-t_adic_valuation(c,uniformizer(nu)))
    return evaluate(numerator(c),0)
end



################################################################################
#
#  constructing and caching polynomial rings for groebner and initial
#  (see groebner_basis.jl and initial.jl)
#
################################################################################

function get_polynomial_ring_for_groebner_simulation(R::MPolyRing, nu::TropicalSemiringMap)
    @req coefficient_ring(R)==valued_field(nu) "coefficient ring is not valued field"
    # return cached polynomial ring if available, create and cache it otherwise
    return get!(polynomial_rings_for_groebner(nu), R, first(polynomial_ring(valued_ring(nu),vcat([:tsim],symbols(R)); cached=false)))
end

# special function for trivial valuation to ensure reusing original ring
function get_polynomial_ring_for_groebner_simulation(R::MPolyRing, nu::TropicalSemiringMap{K,Nothing,minOrMax}) where {K<:Field, minOrMax<:Union{typeof(min),typeof(max)}}
    @req coefficient_ring(R)==valued_field(nu) "coefficient ring is not valued field"
    # return cached polynomial ring if available, create and cache it otherwise
    return get!(polynomial_rings_for_groebner(nu), R, R)
end

function get_polynomial_ring_for_groebner_desimulation(S::MPolyRing, nu::TropicalSemiringMap)
    @req coefficient_ring(S)==valued_ring(nu) "coefficient ring is not valued ring"
    # return cached polynomial ring if available, raise error otherwise

    R = findfirst(isequal(S), polynomial_rings_for_groebner(nu))
    @req !isnothing(R) error("no polynomial ring for groebner basis desimulation found")
    return R
end

# special function for trivial valuation to ensure reusing original ring
function get_polynomial_ring_for_groebner_desimulation(S::MPolyRing, nu::TropicalSemiringMap{K,Nothing,minOrMax}) where {K<:Field, minOrMax<:Union{typeof(min),typeof(max)}}
    @req coefficient_ring(S)==valued_ring(nu) "coefficient ring is not valued ring"
    # return cached polynomial ring if available, raise error otherwise
    return S
end


function get_polynomial_ring_for_initial(R::MPolyRing, nu::TropicalSemiringMap)
    @req coefficient_ring(R)==valued_field(nu) "coefficient ring is not valued field"
    # return cached polynomial ring if available, create and cache it otherwise
    return get!(polynomial_rings_for_initial(nu), R, first(polynomial_ring(residue_field(nu),symbols(R); cached=false)))
end

# special function for trivial valuation to ensure reusing original ring
function get_polynomial_ring_for_initial(R::MPolyRing, nu::TropicalSemiringMap{K,Nothing,minOrMax}) where {K<:Field, minOrMax<:Union{typeof(min),typeof(max)}}
    @req coefficient_ring(R)==valued_field(nu) "coefficient ring is not valued field"
    # return cached polynomial ring if available, create and cache it otherwise
    return get!(polynomial_rings_for_initial(nu), R, R)
end
