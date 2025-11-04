################################################################################
#
#  Tropical semiring map
#
################################################################################

@doc raw"""
    TropicalSemiringMap{DomainType, UniformizerType, MinOrMax}

A structure representing a map from a valued field or ring to a tropical
semiring, encoding a valuation and a choice of min- or max-convention.

# Parameters
- `DomainType`: type of valued field (if specified)
    or type of valued ring (if valued field unspecified)
- `UniformizerType`: type of uniformizer in valued ring
- `MinOrMax`: either `typeof(min)` or `typeof(max)` depending on convention

# Fields
- `valued_field`: the valued field. can be `nothing` if valued field unspecified
- `uniformizer_in_field`: the uniformizer as element in the valued field
    can be `nothing` if valued field unspecified or valuation trivial
- `valued_ring`: the valued ring, either specified by user or
    inferred from the valued field specified by user
- `uniformizer_in_ring`: the uniformizer as element in the valued ring specified by user
    can be `nothing` if valuation trivial
- `residue_field`: the residue field, inferred from fields above
- `tropical_semiring`: either `tropical_semiring(min)` or `tropical_semiring(max)`
    depending on convention
"""
struct TropicalSemiringMap{DomainType, UniformizerType, MinOrMax}
    valued_field::Union{Nothing,Field}
    uniformizer_in_field::Union{Nothing,FieldElem}
    valued_ring::Ring
    uniformizer_in_ring::Union{Nothing,RingElem}
    residue_field::Field
    tropical_semiring::TropicalSemiring{MinOrMax}

    function TropicalSemiringMap{DomainType, UniformizerType, MinOrMax}(
        valuedField::Union{Nothing,Field},
        uniformizerField::Union{Nothing,FieldElem},
        valuedRing::Ring,
        uniformizerRing::Union{Nothing,RingElem},
        residueField::Field,
        tropicalSemiring::TropicalSemiring{MinOrMax}
        ) where {
            DomainType <: Ring,
            UniformizerType <: Union{Nothing,RingElem},
            MinOrMax <: Union{typeof(min),typeof(max)}
        }

        # if valuedField is unspecified, check that uniformizerField is unspecified
        # if uniformizerRing is unspecified, check that uniformizerField is unspecified
        @req !isnothing(valuedField) || isnothing(uniformizerField) "valuedField / uniformizerField mismatch"
        @req !isnothing(uniformizerRing) || isnothing(uniformizerField) "uniformizerRing / uniformizerField mismatch"

        # if no valued field specified, first parameter captures the valued ring
        if isnothing(valuedField)
            return new{typeof(valuedRing),typeof(uniformizerRing),MinOrMax}(
                valuedField,
                uniformizerField,
                valuedRing,
                uniformizerRing,
                residueField,
                tropicalSemiring)
        end

        # otherwise, first parameter captures the valued field
        return new{typeof(valuedField),typeof(uniformizerRing),MinOrMax}(
            valuedField,
            uniformizerField,
            valuedRing,
            uniformizerRing,
            residueField,
            tropicalSemiring)
    end
end



################################################################################
#
#  Properties
#
################################################################################
valued_field(nu::TropicalSemiringMap) = nu.valued_field
uniformizer_in_field(nu::TropicalSemiringMap) = nu.uniformizer_in_field
valued_ring(nu::TropicalSemiringMap) = nu.valued_ring
uniformizer_in_ring(nu::TropicalSemiringMap) = nu.uniformizer_in_ring
residue_field(nu::TropicalSemiringMap) = nu.residue_field
tropical_semiring(nu::TropicalSemiringMap) = nu.tropical_semiring

domain(nu::TropicalSemiringMap) = isnothing(valued_field(nu)) ? valued_ring(nu) : valued_field(nu)
codomain(nu::TropicalSemiringMap) = tropical_semiring(nu)
uniformizer(nu::TropicalSemiringMap) = uniformizer_in_ring(nu)

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

Return a map `nu` from `K` to the min (default) or max tropical semiring `T`
such that `nu(0)=zero(T)` and `nu(c)=one(T)` for `c` non-zero.  In other words,
`nu` represents the trivial valuation on `K`.

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
    print(io, "Map into $(tropical_semiring(nu)) encoding the trivial valuation on $(domain(nu))")
end

# Mapping an element of the valued field or ring to the tropical semiring
function (nu::TropicalSemiringMap{K,Nothing,minOrMax})(c::Union{RingElem,Integer,Rational}) where {K<:Field, minOrMax<:Union{typeof(min),typeof(max)}}
    # return tropical zero if c is zero
    #   and tropical one otherwise
    return (iszero(domain(nu)(c)) ? zero(tropical_semiring(nu)) : one(tropical_semiring(nu)))
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

# Examples
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
    print(io, "Map into $(tropical_semiring(nu)) encoding the $(uniformizer(nu))-adic valuation on $(domain(nu))")
end

# Mapping an element of the valued field or ring to the tropical semiring
function (nu::TropicalSemiringMap{QQField,ZZRingElem,minOrMax})(c::Union{RingElem,Integer,Rational}) where minOrMax<:Union{typeof(min),typeof(max)}
    c = domain(nu)(c)
    # return tropical zero if c is zero
    #   and p-adic valuation otherwise
    # preserve_ordering ensures that valuation is negated if convention(nu)==max
    iszero(c) && return zero(tropical_semiring(nu))
    return tropical_semiring(nu)(valuation(c,uniformizer(nu)); preserve_ordering=true)
end

# Mapping an element of the valued field or ring to the residue field
function initial(c::Union{RingElem,Integer,Rational}, nu::TropicalSemiringMap{QQField,ZZRingElem,minOrMax}) where minOrMax<:Union{typeof(min),typeof(max)}
    c = domain(nu)(c)
    # return residue field zero if c is zero
    #   and the correct non-zero residue otherwise
    iszero(c) && return zero(residue_field(nu)) # if c is zero, return 0
    c *= uniformizer_in_field(nu)^(-valuation(c,uniformizer(nu)))
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

# Examples
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
    print(io, "Map into $(tropical_semiring(nu)) encoding the $(uniformizer(nu))-adic valuation on $(domain(nu))")
end

# Mapping an element of the valued field or ring to the tropical semiring
function t_adic_valuation(c::Generic.RationalFunctionFieldElem, t::PolyRingElem)
    return valuation(numerator(c),t)-valuation(denominator(c),t)
end

function (nu::TropicalSemiringMap{Kt,t,minOrMax})(c::Union{RingElem,Integer,Rational}) where {Kt<:Generic.RationalFunctionField, t<:PolyRingElem, minOrMax<:Union{typeof(min),typeof(max)}}
    c = domain(nu)(c)
    # return tropical zero if c is zero
    #   and p-adic valuation otherwise
    # preserve_ordering ensures that valuation is negated if convention(nu)==max
    iszero(c) && return zero(tropical_semiring(nu)) # if c zero, return tropical zero
    return tropical_semiring(nu)(t_adic_valuation(c,uniformizer(nu)); preserve_ordering=true)
end

# Mapping an element of the valued field or ring to the residue field
function initial(c::Union{RingElem,Integer,Rational}, nu::TropicalSemiringMap{Kt,t,minOrMax}) where {Kt<:Generic.RationalFunctionField, t<:PolyRingElem, minOrMax<:Union{typeof(min),typeof(max)}}
    c = domain(nu)(c)
    # return residue field zero if c is zero
    #   and the correct non-zero residue otherwise
    iszero(c) && return zero(residue_field(nu)) # if c is zero, return 0
    c *= uniformizer_in_field(nu)^(-t_adic_valuation(c,uniformizer(nu)))
    return evaluate(numerator(c),0)
end
