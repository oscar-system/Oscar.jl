################################################################################
#
#  Tropical semiring maps
#  ======================
#  maps from a field K to a tropical semiring T with the purpose of encoding
#    - a valuation on K
#    - a choice of min- or max-convention
#
################################################################################
struct TropicalSemiringMap{typeofValuedField,typeofUniformizer,minOrMax}
    valued_field::typeofValuedField
    uniformizer::typeofUniformizer
    tropical_semiring::TropicalSemiring{minOrMax}
end



################################################################################
#
#  Properties
#
################################################################################
valued_field(nu::TropicalSemiringMap) = nu.valued_field
uniformizer(nu::TropicalSemiringMap) = nu.uniformizer
tropical_semiring(nu::TropicalSemiringMap) = nu.tropical_semiring

convention(nu::TropicalSemiringMap{typeofValuedField,typeofUniformizer,typeof(min)}) where {typeofValuedField,typeofUniformizer} = min
convention(nu::TropicalSemiringMap{typeofValuedField,typeofUniformizer,typeof(max)}) where {typeofValuedField,typeofUniformizer} = max

is_trivial(nu::TropicalSemiringMap{K,Nothing,minOrMax}) where {K,minOrMax<:Union{typeof(min),typeof(max)}} = true
is_trivial(nu::TropicalSemiringMap) = false

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
    return TropicalSemiringMap{typeof(K),Nothing,typeof(minOrMax)}(K,nothing,tropical_semiring(minOrMax))
end

# display:
function Base.show(io::IO, nu::TropicalSemiringMap{K,Nothing,minOrMax} where {K<:Ring, minOrMax<:Union{typeof(min),typeof(max)}})
    print(io, "Map into $(tropical_semiring(nu)) encoding the trivial valuation on $(valued_field(nu))")
end

# evaluation:
function (nu::TropicalSemiringMap{K,Nothing,minOrMax})(c::Union{RingElem,Integer,Rational}) where {K<:Field, minOrMax<:Union{typeof(min),typeof(max)}}
    return (iszero(valued_field(nu)(c)) ? inf(tropical_semiring(nu)) : one(tropical_semiring(nu)))
end

# valued ring:
function valued_ring(nu::TropicalSemiringMap{K,Nothing,minOrMax}) where {K,minOrMax<:Union{typeof(min),typeof(max)}}
    return valued_field(nu)
end

# residue field:
function residue_field(nu::TropicalSemiringMap{K,Nothing,minOrMax}) where {K,minOrMax<:Union{typeof(min),typeof(max)}}
    return valued_field(nu)
end

# initial:
function initial(c::RingElem, nu::TropicalSemiringMap{K,Nothing,minOrMax}) where {K,minOrMax<:Union{typeof(min),typeof(max)}}
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
function tropical_semiring_map(K::QQField, p::Union{RingElem,Integer,Rational}, minOrMax::Union{typeof(min),typeof(max)}=min)
    p = ZZ(p)
    @req is_prime(ZZ(p)) "input p not prime"
    @req p < 2^63-1 "input p may not exceed 2^63"
    return TropicalSemiringMap{typeof(K),typeof(p),typeof(minOrMax)}(K,p,tropical_semiring(minOrMax))
end

# Display:
function Base.show(io::IO, nu::TropicalSemiringMap{QQField,ZZRingElem,minOrMax}) where {minOrMax<:Union{typeof(min),typeof(max)}}
    print(io, "Map into $(tropical_semiring(nu)) encoding the $(uniformizer(nu))-adic valuation on $(valued_field(nu))")
end

# Evaluation:
function (nu::TropicalSemiringMap{QQField,ZZRingElem,typeof(min)})(c::Union{RingElem,Integer,Rational})
    c = valued_field(nu)(c)
    iszero(c) && return zero(tropical_semiring(nu)) # if input is zero, return tropical zero
    return tropical_semiring(nu)(valuation(c,ZZ(uniformizer(nu))))
end
function (nu::TropicalSemiringMap{QQField,ZZRingElem,typeof(max)})(c::Union{RingElem,Integer,Rational})
    c = valued_field(nu)(c)
    iszero(c) && return zero(tropical_semiring(nu)) # if input is zero, return tropical zero
    return tropical_semiring(nu)(-valuation(c,ZZ(uniformizer(nu))))
end

# valued ring:
function valued_ring(nu::TropicalSemiringMap{QQField,ZZRingElem,minOrMax}) where {minOrMax<:Union{typeof(min),typeof(max)}}
    return ZZ
end

# residue field:
function residue_field(nu::TropicalSemiringMap{QQField,ZZRingElem,minOrMax}) where {minOrMax<:Union{typeof(min),typeof(max)}}
    return GF(Int(uniformizer(nu)))
end

# initial:
function initial(c::Union{RingElem,Integer,Rational}, nu::TropicalSemiringMap{QQField,ZZRingElem,minOrMax}) where {minOrMax<:Union{typeof(min),typeof(max)}}
    c = valued_field(nu)(c)
    iszero(c) && return zero(residue_field(nu)) # if c is zero, return 0
    c *= QQ(uniformizer(nu))^(-valuation(c,ZZ(uniformizer(nu))))
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
    @req isone(denominator(t)) "input uniformizer denominator not 1"
    t = numerator(t)
    @req degree(t)>0 "input uniformizer constant"
    return TropicalSemiringMap{typeof(Kt),typeof(t),typeof(minOrMax)}(Kt,t,tropical_semiring(minOrMax))
end

# Display:
function Base.show(io::IO, nu::TropicalSemiringMap{Kt,t,minOrMax}) where {Kt<:Generic.RationalFunctionField, t<:PolyRingElem, minOrMax<:Union{typeof(min),typeof(max)}}
    print(io, "Map into $(tropical_semiring(nu)) encoding the $(uniformizer(nu))-adic valuation on $(valued_field(nu))")
end

# Evaluation:
function t_adic_valuation(c::Generic.RationalFunctionFieldElem, t::PolyRingElem)
    c_num = numerator(c)
    c_nom = denominator(c)
    return valuation(c_num,t)-valuation(c_nom,t)
end

function (nu::TropicalSemiringMap{Kt,t,typeof(min)})(c::Union{RingElem,Integer,Rational}) where {Kt<:Generic.RationalFunctionField, t<:PolyRingElem}
    c = valued_field(nu)(c)
    iszero(c) && return zero(tropical_semiring(nu)) # if c zero, return tropical zero
    return tropical_semiring(nu)(t_adic_valuation(c,uniformizer(nu)))
end

function (nu::TropicalSemiringMap{Kt,t,typeof(max)})(c::Union{RingElem,Integer,Rational}) where {Kt<:Generic.RationalFunctionField, t<:PolyRingElem}
    c = valued_field(nu)(c)
    iszero(c) && return zero(tropical_semiring(nu)) # if c zero, return tropical zero
    return tropical_semiring(nu)(-t_adic_valuation(c,uniformizer(nu)))
end

# valued ring:
function valued_ring(nu::TropicalSemiringMap{Kt,t,minOrMax}) where {Kt<:Generic.RationalFunctionField, t<:PolyRingElem, minOrMax<:Union{typeof(min),typeof(max)}}
    return polynomial_ring(base_ring(valued_field(nu)),symbols(valued_field(nu)))[1]
end

# residue field:
function residue_field(nu::TropicalSemiringMap{Kt,t,minOrMax}) where {Kt<:Generic.RationalFunctionField, t<:PolyRingElem, minOrMax<:Union{typeof(min),typeof(max)}}
    return base_ring(valued_field(nu))
end

# initial:
function initial(c::Union{RingElem,Integer,Rational}, nu::TropicalSemiringMap{Kt,t,minOrMax}) where {Kt<:Generic.RationalFunctionField, t<:PolyRingElem, minOrMax<:Union{typeof(min),typeof(max)}}
    c = valued_field(nu)(c)
    iszero(c) && return zero(residue_field(nu)) # if c is zero, return 0
    c *= valued_field(nu)(uniformizer(nu))^(-t_adic_valuation(c,uniformizer(nu)))
    return residue_field(nu)(c)
end
