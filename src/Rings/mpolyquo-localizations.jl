import AbstractAlgebra: Ring, RingElem, Generic.FracFieldElem
import Base: issubset


########################################################################
# Localizations of polynomial algebras                                 #
########################################################################
# 
# Let R = 𝕜[x₁,…,xₘ] be a polynomial ring, I ⊂ R some ideal 
# and P = R/I its quotient. Then P is naturally an R-module 
# and localization of P as a ring coincides with localization 
# as an R-module in the sense that for every multiplicative 
# set T ⊂ R there is a commutative diagram 
#
#         R    →  P = R/I
#         ↓       ↓
#   W = R[T⁻¹] → P[T⁻¹].
#
# Observe that, moreover, for every multiplicative set 
# T' ⊂ P the preimage T of T' in R is also a multiplicative set. 
#
# We may therefore treat localizations of polynomial algebras 
# as localizations of modules over free polynomial rings and 
# exclusively use the types of multiplicative sets which are 
# available for the latter.
#
# Note the following differences compared to the standard usage 
# of the localization interface:
#
#  * The `base_ring` returns neither P, nor W, but R.
#  * The `BaseRingType` is the type of R and similar for 
#    the other ring-based type parameters.
#
# This is to make the data structure most accessible for 
# the computational backends.
#
#  * The type returned by `numerator` and `denominator` 
#    on an element of type `MPolyQuoLocRingElem` is 
#    not `RingElemType`, but the type of `P`. 
#
# This is to comply with the purely mathematical viewpoint
# where elements of localized rings are fractions of 
# residue classes rather than residue classes of fractions. 
#

@doc raw"""
    MPolyQuoLocRing{
        BaseRingType,
        BaseRingElemType,
        RingType,
        RingElemType,
        MultSetType <: AbsMultSet{RingType, RingElemType}
      } <: AbsLocalizedRing{
        RingType,
        RingElemType,
        MultSetType
      }
    
Localization ``L = (𝕜[x₁,…,xₙ]/I)[S⁻¹]`` of a quotient
``𝕜[x₁,…,xₙ]/I`` of a polynomial ring ``P = 𝕜[x₁,…,xₙ]``
of type `RingType` over a base ring ``𝕜`` of type `BaseRingType` at a
multiplicative set ``S ⊂ P`` of type `MultSetType`.
"""
@attributes mutable struct MPolyQuoLocRing{
    BaseRingType,
    BaseRingElemType,
    RingType,
    RingElemType,
    MultSetType <: AbsMultSet{RingType, RingElemType}
  } <: AbsLocalizedRing{
    RingType,
    RingElemType,
    MultSetType
  }
  R::RingType
  I::MPolyIdeal{RingElemType}
  S::MultSetType
  Q::MPolyQuoRing{RingElemType}
  W::MPolyLocRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}

  function MPolyQuoLocRing(
      R::RingType,
      I::MPolyIdeal{RingElemType},
      S::MultSetType,
      Q::MPolyQuoRing{RingElemType},
      W::MPolyLocRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}
    ) where {
      BaseRingType<:Ring, 
      BaseRingElemType<:RingElem, 
      RingType<:MPolyRing, 
      RingElemType<:MPolyRingElem, 
      MultSetType<:AbsMultSet{RingType, RingElemType}
    }
    base_ring(I) === R || error("Ideal does not belong to the ring")
    base_ring(Q) === R || error("The quotient ring does not come from the given ring")
    # The following line throws obscure error messages that might yield a bug for MPolyIdeals.
    # So it's commented out for now.
    #modulus(Q) == I || error("the modulus of the quotient ring does not coincide with the ideal")
    S === inverted_set(W) || error("the multiplicative set does not coincide with the inverted set of the localized ring")
    base_ring(W) === R || error("the localization does not come from the given ring")
    ring(S) === R || error("Multiplicative set does not belong to the ring")
    k = coefficient_ring(R)
    L = new{typeof(k), elem_type(k), typeof(R), RingElemType, MultSetType}(R, I, S, Q, W)
    return L
  end
end

### for convenience of later use
const MPAnyQuoRing = Union{MPolyQuoLocRing, MPolyQuoRing}

const MPAnyNonQuoRing = Union{MPolyRing, MPolyLocRing}

const MPolyAnyRing = Union{MPolyRing, MPolyQuoRing,
                MPolyLocRing,MPolyQuoLocRing
               }


### type getters 
coefficient_ring_type(::Type{MPolyQuoLocRing{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = BRT
coefficient_ring_type(L::MPolyQuoLocRing{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST} = coefficient_ring_type(typeof(L))

base_ring_type(::Type{MPolyQuoLocRing{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = RT

mult_set_type(::Type{MPolyQuoLocRing{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = MST
mult_set_type(L::MPolyQuoLocRing{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST} = mult_set_type(typeof(L))

localized_ring_type(::Type{MPolyQuoLocRing{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = MPolyLocRing{BRT, BRET, RT, RET, MST}
localized_ring_type(L::MPolyQuoLocRing) = localized_ring_type(typeof(L))

ideal_type(::Type{MPolyQuoLocalizedRingType}) where {MPolyQuoLocalizedRingType<:MPolyQuoLocRing} = MPolyQuoLocalizedIdeal{MPolyQuoLocalizedRingType, elem_type(MPolyQuoLocalizedRingType), ideal_type(localized_ring_type(MPolyQuoLocalizedRingType))}



### required getter functions 
base_ring(L::MPolyQuoLocRing) = L.R
inverted_set(L::MPolyQuoLocRing) = L.S

### additional getter functions
@doc raw"""
    modulus(L::MPolyQuoLocRing)

Given ``L = (𝕜[x₁,…,xₙ]/I)[S⁻¹]``, return ``IS⁻¹``.
"""
@attr ideal_type(localized_ring_type(L)) function modulus(L::MPolyQuoLocRing) 
  return localized_ring(L)(L.I)
end

### for compatibility -- also provide modulus in the trivial case
modulus(R::MPAnyNonQuoRing)=ideal(R, elem_type(R)[])


@doc raw"""
    underlying_quotient(L::MPolyQuoLocRing)

Given ``L = (𝕜[x₁,…,xₙ]/I)[S⁻¹]``, return ``𝕜[x₁,…,xₙ]/I``.
"""
underlying_quotient(L::MPolyQuoLocRing) = L.Q

## 3 more signatures for compatibility to make quotient_ring agnostic
underlying_quotient(L::MPolyQuoRing) = L

@attr MPolyQuoRing function underlying_quotient(L::MPolyRing)
   return quo(L,ideal(L,[zero(L)]))[1]
end

@attr MPolyQuoRing function underlying_quotient(L::MPolyLocRing)
   P = base_ring(L)
   return quo(P,ideal(P,[zero(P)]))[1]
end

@doc raw"""
    localized_ring(L::MPolyQuoLocRing)

Given ``L = (𝕜[x₁,…,xₙ]/I)[S⁻¹]``, return ``𝕜[x₁,…,xₙ][S⁻¹]``.
"""
localized_ring(L::MPolyQuoLocRing) = L.W

## 3 more signatures for compatibility to make localized_ring agnostic
localized_ring(L::MPolyLocRing) = L

@attr MPolyLocRing function localized_ring(L::MPolyQuoRing)
   P = base_ring(L)
   return localization(P, units_of(P))[1]
end

@attr MPolyLocRing function localized_ring(L::MPolyRing)
   return localization(L, units_of(L))[1]
end

@doc raw"""
    gens(L::MPolyQuoLocRing)

Given ``L = (𝕜[x₁,…,xₙ]/I)[S⁻¹]``, return the vector ``[x₁//1,…,xₙ//1]∈ Lⁿ``.
"""
gens(L::MPolyQuoLocRing) = L.(gens(base_ring(L)))

gen(L::MPolyQuoLocRing, i::Int) = L(gen(base_ring(L), i))

### printing

function Base.show(io::IO, ::MIME"text/plain", L::MPolyQuoLocRing)
  io = pretty(io)
  println(io, "Localization")
  print(io, Indent())
  print(io, "of ", Lowercase())
  show(io, MIME("text/plain"), underlying_quotient(L))
  println(io)
  print(io, "at ", Lowercase(), inverted_set(L))
  print(io, Dedent())
end

function Base.show(io::IO, L::MPolyQuoLocRing)
  if is_terse(io)
    print(io, "Localized quotient of multivariate polynomial ring")
  else
    io = terse(pretty(io))
    print(io, "Localization of ")
    print(io, Lowercase(), underlying_quotient(L))
    print(io, " at ", Lowercase(), inverted_set(L))
  end
end

### additional constructors
function quo(
    W::MPolyLocRing,
    I::MPolyLocalizedIdeal
  )
  R = base_ring(W)
  S = inverted_set(W)
  J = ideal(R, numerator.(gens(I)))
  L = MPolyQuoLocRing(R, J, S, quo(R, J)[1], W)
  return L, hom(W, L, hom(R, L, gens(L), check=false), check=false)
end

function quo(
    L::MPolyQuoLocRing,
    I::MPolyLocalizedIdeal
  )
  R = base_ring(L)
  S = inverted_set(L)
  base_ring(I) = localized_ring(L) || error("ideal does not belong to the correct ring")
  J = pre_saturated_ideal(I)
  W = MPolyQuoLocRing(R, J, S, quo(R, J)[1], localized_ring(L))
  return W, hom(L, W, gens(W), check=false)
end

function quo(
    L::MPolyQuoLocRing{BRT, BRET, RT, RET, MST},
    J::MPolyIdeal{RET}
  ) where {BRT, BRET, RT, RET, MST}
  R = base_ring(L)
  S = inverted_set(L)
  W = localized_ring(L) 
  J = J + modulus(underlying_quotient(L))
  P = MPolyQuoLocRing(R, J, S, quo(R, J)[1], W)
  return P, hom(L, P, gens(P), check=false)
end

@doc raw"""
    localization(RQ::MPolyQuoRing, U::AbsMPolyMultSet)

Given a quotient `RQ` of a multivariate polynomial ring `R` with projection map
`p : R -> RQ`, say, and given a multiplicatively closed subset `U` of `R`, return the 
localization of `RQ` at `p(U)`, together with the localization map.

# Examples

```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> I = ideal(R, [2*x^2-y^3, 2*x^2-y^5])
Ideal generated by
  2*x^2 - y^3
  2*x^2 - y^5

julia> U = complement_of_point_ideal(R, [0, 0])
Complement
  of maximal ideal corresponding to rational point with coordinates (0, 0)
  in multivariate polynomial ring in 2 variables over QQ

julia> RQ, _ = quo(R, I);

julia> RQL, iota = localization(RQ, U);

julia> RQL
Localization
  of quotient
    of multivariate polynomial ring in 2 variables x, y
      over rational field
    by ideal (2*x^2 - y^3, 2*x^2 - y^5)
  at complement of maximal ideal of point (0, 0)

julia> iota
Map defined by a julia-function
  from quotient of multivariate polynomial ring by ideal (2*x^2 - y^3, 2*x^2 - y^5)
  to localization of RQ at complement of maximal ideal

```
""" localization(A::MPolyQuoRing, U::AbsMPolyMultSet)

###localization is an Abstract Algebra alias for Localization

function localization(Q::MPolyQuoRing{RET}, S::MultSetType) where {RET <: RingElem, MultSetType <: AbsMultSet}
  L = MPolyQuoLocRing(base_ring(Q), modulus(Q), S, Q, localization(S)[1])
  return L, MapFromFunc(Q, L, (x->L(lift(x))))
end

function localization(
    L::MPolyQuoLocRing{BRT, BRET, RT, RET, MST}, 
    S::AbsMPolyMultSet{BRT, BRET, RT, RET}
  ) where {BRT, BRET, RT, RET, MST}
  ring(S) === base_ring(L) || error("multiplicative set does not belong to the correct ring")
  #issubset(S, inverted_set(L)) && return L, MapFromFunc(L, L, identity)
  U = inverted_set(L)*S
  W = MPolyQuoLocRing(base_ring(L), modulus(underlying_quotient(L)), U, underlying_quotient(L), localization(U)[1])
  return W, MapFromFunc(L, W, (x->W(lifted_numerator(x), lifted_denominator(x), check=false)))
end

function MPolyQuoLocRing(R::RT, I::Ideal{RET}, T::MultSetType) where {RT<:MPolyRing, RET<:MPolyRingElem, MultSetType<:AbsMultSet} 
  return MPolyQuoLocRing(R, I, T, quo(R, I)[1], localization(T)[1])
end

function MPolyQuoLocRing(R::RT) where {RT<:MPolyRing} 
  I = ideal(R, zero(R))
  Q, _ = quo(R, I)
  U = units_of(R)
  W, _ = localization(U)
  return MPolyQuoLocRing(R, I, U, Q, W)
end

function MPolyQuoLocRing(Q::RT) where {RT<:MPolyQuoRing}
  R = base_ring(Q)
  I = modulus(Q)
  U = units_of(R)
  W, _ = localization(U)
  return MPolyQuoLocRing(R, I, U, Q, W)
end

function MPolyQuoLocRing(W::MPolyLocRing)
  R = base_ring(W)
  I = ideal(R, zero(R))
  Q, _ = quo(R, I)
  U = inverted_set(W)
  return MPolyQuoLocRing(R, I, U, Q, W)
end

function Base.in(f::AbstractAlgebra.Generic.FracFieldElem{RET}, L::MPolyQuoLocRing{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST}
  R = base_ring(L)
  R === parent(numerator(f)) || error("element does not belong to the correct ring")
  denominator(f) in inverted_set(L) && return true
  return numerator(f) in ideal(L, denominator(f))
end


### generation of random elements 
function rand(W::MPolyQuoLocRing, v1::AbstractUnitRange{Int}, v2::AbstractUnitRange{Int}, v3::AbstractUnitRange{Int})
  return W(rand(localized_ring(W), v1, v2, v3))
end

########################################################################
# Elements of localizations of polynomial algebras                     #
########################################################################

@doc raw"""
    MPolyQuoLocRingElem{
      BaseRingType, 
      BaseRingElemType,
      RingType,
      RingElemType, 
      MultSetType
    } <: AbsLocalizedRingElem{
      RingType,
      RingElemType, 
      MultSetType
    } 

Elements ``a//b`` of localizations ``L = (𝕜[x₁,…,xₙ]/I)[S⁻¹]`` of type 
`MPolyQuoLocRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}`.
"""
mutable struct MPolyQuoLocRingElem{
    BaseRingType, 
    BaseRingElemType,
    RingType,
    RingElemType, 
    MultSetType
  } <: AbsLocalizedRingElem{
    RingType,
    RingElemType, 
    MultSetType
  } 

  # the parent ring
  L::MPolyQuoLocRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}
  # representatives of numerator and denominator
  numerator::RingElemType
  denominator::RingElemType
  is_reduced::Bool

  function MPolyQuoLocRingElem(
      L::MPolyQuoLocRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}, 
      a::RingElemType,
      b::RingElemType;
      check::Bool=true,
      is_reduced::Bool=false
    ) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}

    S = inverted_set(L)
    R = base_ring(L)
    parent(a) === parent(b) === R || error("elements do not belong to the correct ring")
    @check b in S || error("denominator is not admissible")
    return new{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}(L, a, b, is_reduced)
  end
end

### type getters
coefficient_ring_type(::Type{MPolyQuoLocRingElem{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = BRT
coefficient_ring_type(f::MPolyQuoLocRingElem{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST} = coefficient_ring_type(typeof(f))

base_ring_type(::Type{MPolyQuoLocRingElem{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = RT

mult_set_type(::Type{MPolyQuoLocRingElem{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = MST
mult_set_type(f::MPolyQuoLocRingElem{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST} = mult_set_type(typeof(f))

### required getter functions 
parent(a::MPolyQuoLocRingElem) = a.L
numerator(a::MPolyQuoLocRingElem) = underlying_quotient(parent(a))(a.numerator) 
denominator(a::MPolyQuoLocRingElem) = underlying_quotient(parent(a))(a.denominator) 

### additional getter functions
underlying_quotient(a::MPolyQuoLocRingElem) = underlying_quotient(parent(a))
localized_ring(a::MPolyQuoLocRingElem) = localized_ring(parent(a))
is_reduced(a::MPolyQuoLocRingElem) = a.is_reduced

@doc raw"""
    lifted_numerator(a::MPolyQuoLocRingElem)

Given ``A//B ∈ (𝕜[x₁,…,xₙ]/I)[S⁻¹]``, return a representative
``a ∈ 𝕜[x₁,…,xₙ]`` of the numerator. 
"""
lifted_numerator(a::MPolyQuoLocRingElem) = a.numerator

@doc raw"""
    lifted_denominator(a::MPolyQuoLocRingElem)

Given ``A//B ∈ (𝕜[x₁,…,xₙ]/I)[S⁻¹]``, return a representative
``b ∈  𝕜[x₁,…,xₙ]`` of the denominator.
"""
lifted_denominator(a::MPolyQuoLocRingElem) = a.denominator

@doc raw"""
    fraction(a::MPolyQuoLocRingElem)

Given ``A//B ∈ (𝕜[x₁,…,xₙ]/I)[S⁻¹]``, return a representative
``a//b ∈ Quot(𝕜[x₁,…,xₙ])`` of the fraction. 
"""
fraction(a::MPolyQuoLocRingElem) = lifted_numerator(a)//lifted_denominator(a)

### copying of elements
function Base.deepcopy_internal(f::MPolyQuoLocRingElem, dict::IdDict)
  return parent(f)(f, check=false, is_reduced=is_reduced(f))
end

### required conversions
(L::MPolyQuoLocRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType})(f::RingElemType; is_reduced::Bool=false, check::Bool=true) where {BaseRingType, BaseRingElemType, RingType, RingElemType<:RingElem, MultSetType} = MPolyQuoLocRingElem(L, f, one(f), check=false, is_reduced=is_reduced)

function (L::MPolyQuoLocRing{
                                   BaseRingType, 
                                   BaseRingElemType, 
                                   RingType, 
                                   RingElemType, 
                                   MultSetType
                                  })(
                                     a::RingElemType, 
                                     b::RingElemType;
                                     check::Bool=true,
                                     is_reduced::Bool=false
                                    ) where {
                                             BaseRingType, 
                                             BaseRingElemType, 
                                             RingType, 
                                             RingElemType, 
                                             MultSetType
                                            } 
  @check begin
    b in inverted_set(L) || return convert(L, a//b)
  end
  return MPolyQuoLocRingElem(L, a, b, check=false, is_reduced=is_reduced)
end

### additional conversions
function (L::MPolyQuoLocRing{BRT, BRET, RT, RET, MST})(f::FracFieldElem{RET}; check::Bool=true, is_reduced::Bool=false) where {BRT, BRET, RT, RET, MST}
  R = base_ring(L)
  return L(R(numerator(f)), R(denominator(f)), check=check, is_reduced=is_reduced)
end

(L::MPolyQuoLocRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType})(a::T, b::T; check::Bool=true, is_reduced::Bool=false) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType, T<:MPolyQuoRingElem{RingElemType}} = L(lift(a), lift(b), check=check, is_reduced=is_reduced)

function (L::MPolyQuoLocRing{BRT, BRET, RT, RET, MST})(f::MPolyQuoLocRingElem{BRT, BRET, RT, RET, MST}; check::Bool=true, is_reduced::Bool=false) where {BRT, BRET, RT, RET, MST}
  parent(f) === L && return f
  return L(lifted_numerator(f), lifted_denominator(f), check=check, is_reduced=is_reduced)
end

function (L::MPolyQuoLocRing{BRT, BRET, RT, RET, MST})(f::MPolyLocRingElem{BRT, BRET, RT, RET, MST}; check::Bool=true, is_reduced::Bool=false) where {BRT, BRET, RT, RET, MST}
  parent(f) === localized_ring(L) && return L(numerator(f), denominator(f), check=false, is_reduced=is_reduced)
  return L(numerator(f), denominator(f), check=check, is_reduced=is_reduced)
end

function (L::MPolyQuoLocRing{BRT, BRET, RT, RET, MST})(f::MPolyLocRingElem{BRT, BRET, RT, RET, MST}; check::Bool=true, is_reduced::Bool=false) where {BRT, BRET, RT, RET, MST<:MPolyComplementOfKPointIdeal}
  return L(numerator(f), denominator(f), check=check, is_reduced=is_reduced)
end

function (L::MPolyQuoLocRing{BRT, BRET, RT, RET, MST})(f::MPolyQuoRingElem{RET}; check::Bool=true, is_reduced::Bool=false) where {BRT, BRET, RT, RET, MST} 
  base_ring(parent(f)) === base_ring(L) || error("the given element does not belong to the correct ring") 
  if parent(f) !== underlying_quotient(L) 
    @check all(x->(iszero(L(x))), gens(modulus(parent(f)))) "coercion is not well defined"
  end
  return L(lift(f), check=check)
end

### additional functionality
@doc raw"""
    lift(f::MPolyQuoLocRingElem)

Given ``f = A//B ∈ (𝕜[x₁,…,xₙ]/I)[S⁻¹]``, return a representative
``a//b ∈  𝕜[x₁,…,xₙ][S⁻¹]`` of the fraction. 
"""
lift(f::MPolyQuoLocRingElem) = localized_ring(f)(lifted_numerator(f), lifted_denominator(f), check=false)


@doc raw"""
    is_unit(f::MPolyQuoLocRingElem) 

Return `true`, if `f` is a unit of `parent(f)`, `true` otherwise.

# Examples

```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> I = ideal(R, [2*x^2-y^3, 2*x^2-y^5])
Ideal generated by
  2*x^2 - y^3
  2*x^2 - y^5

julia> U = complement_of_point_ideal(R, [0, 0])
Complement
  of maximal ideal corresponding to rational point with coordinates (0, 0)
  in multivariate polynomial ring in 2 variables over QQ

julia> RQ, p = quo(R, I);

julia> RQL, iota = localization(RQ, U);

julia> phi = compose(p, iota);

julia> is_unit(phi(x-1))
true

```
""" is_unit(f::MPolyQuoLocRingElem)

function is_unit(f::MPolyQuoLocRingElem) 
  lifted_numerator(f) in inverted_set(parent(f)) && return true
  R=localized_ring(parent(f))
  return one(R) in modulus(parent(f)) + ideal(R, lift(f))
end

function is_unit(L::MPolyQuoLocRing, f::MPolyLocRingElem) 
  parent(f) === localized_ring(L) || error("element does not belong to the correct ring")
  numerator(f) in inverted_set(L) && return true
  one(localized_ring(L)) in modulus(L) + ideal(localized_ring(L), f)
end

function is_unit(L::MPolyQuoLocRing{BRT, BRET, RT, RET, MST}, f::RET) where {BRT, BRET, RT, RET, MST}
  parent(f) === base_ring(L) || error("element does not belong to the correct ring")
  f in inverted_set(L) && return true
  return one(localized_ring(L)) in modulus(L) + ideal(localized_ring(L), localized_ring(L)(f))
end

function is_unit(L::MPolyQuoLocRing{BRT, BRET, RT, RET, MST}, f::MPolyQuoRingElem{RET}) where {BRT, BRET, RT, RET, MST}
  parent(f) === underlying_quotient(L) || error("element does not belong to the correct ring")
  lift(f) in inverted_set(L) && return true
  one(localized_ring(L)) in modulus(L) + ideal(localized_ring(L), localized_ring(L)(f))
end

function is_zero_divisor(f::MPolyQuoLocRingElem{<:Field})
  iszero(f) && return true
  # The next block is basically useless when the coefficient ring is
  # a field, because it is merely another `is_zero`-check. However,
  # once more functionality is working, it will actually do stuff and
  # the above signature can be widened.
  if is_constant(lifted_numerator(f)) && is_constant(lifted_denominator(f))
    c = first(AbstractAlgebra.coefficients(lift(numerator(f))))
    return is_zero_divisor(c)
  end
  return !is_zero(quotient(ideal(parent(f), zero(f)), ideal(parent(f), f)))
end

# WARNING: This routine runs forever if f is not a unit in L. 
# So this needs to be checked first!
function inv(L::MPolyQuoLocRing{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}}, 
    f::MPolyQuoRingElem{RET}) where {BRT, BRET, RT, RET}
  Q = underlying_quotient(L)
  parent(f) === underlying_quotient(L) || error("element does not belong to the correct ring")
  W = localized_ring(L)
  R = base_ring(L)
  I = saturated_ideal(modulus(L))
  d = prod(denominators(inverted_set(W)); init=one(R))
  powers_of_d = [d]
  ### apply logarithmic bisection to find a power dᵏ ≡  c ⋅ f mod I
  (result, coefficient) = divides(one(Q), f)
  # check whether f is already a unit
  result && return L(coefficient)
  push!(powers_of_d, d)
  abort = false
  # find some power which works
  while !abort
    (abort, coefficient) = divides(Q(last(powers_of_d)), f)
    if !abort
      push!(powers_of_d, last(powers_of_d)^2)
    end
  end
  # find the minimal power that works
  upper = pop!(powers_of_d)
  lower = pop!(powers_of_d)
  while length(powers_of_d) > 0
    middle = lower*pop!(powers_of_d)
    (result, coefficient) = divides(Q(middle), f)
    if result 
      upper = middle
    else 
      lower = middle
    end
  end
  return L(lift(coefficient), upper)
end

function inv(f::MPolyQuoLocRingElem{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}}) where {BRT, BRET, RT, RET}
  isone(f) && return f
  
  if length(lifted_numerator(f)) > 10000
    L = parent(f)
    id, id_inv = _as_affine_algebra_with_many_variables(L)
    aa = simplify(id(L(numerator(f))))
    bb = simplify(id(L(denominator(f))))
    success, cc = _divides_hack(aa, bb)
    !success && error("element can not be converted to localization")
    return id_inv(simplify(cc))
  end

  lifted_numerator(f) in inverted_set(parent(f)) && return parent(f)(denominator(f), numerator(f), check=false)

  return convert(parent(f), lifted_denominator(f)//lifted_numerator(f))
  return parent(f)(denominator(f), numerator(f))
  # The following was the original line:
  return parent(f)(denominator(f), numerator(f))
end

### 
# Assume that [f] = [a]//[b] is an admissible element of L = (R/I)[S⁻¹] and bring it 
# to the form [f] = [c]//[dᵏ] with d∈ S.
function convert(
    L::MPolyQuoLocRing{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}}, 
    f::AbstractAlgebra.Generic.FracFieldElem{RET}
  ) where {BRT, BRET, RT, RET}
  a = numerator(f)
  b = denominator(f)
  Q = underlying_quotient(L)
  parent(a) === base_ring(L) || error("element does not belong to the correct ring")
  W = localized_ring(L)
  R = base_ring(L)
  I = saturated_ideal(modulus(L))
  isone(I) && return zero(L)
  denoms = denominators(inverted_set(W))
  if iszero(length(denoms)) || all(is_one, denoms)
    success, q = divides(Q(a), Q(b))
    success || error("element can not be converted")
    return L(q)
  end
  d = prod(denominators(inverted_set(W)); init=one(R))
  powers_of_d = [d]
  ### apply logarithmic bisection to find a power a ⋅dᵏ ≡  c ⋅ b mod I
  (result, coefficient) = _divides_hack(Q(a), Q(b))
  # check whether f is already a unit
  result && return L(coefficient)
  # If we have localized at the trivial set, then this is the end.
  isone(d) && error("element can not be converted to the localization")
  push!(powers_of_d, d)
  abort = false
  # find some power which works
  while !abort
   if length(terms(last(powers_of_d))) > 100
      id, id_inv = _as_affine_algebra_with_many_variables(L)
      cc = id(L(a)) # the result
      fac_b = factor(b)
      for (f, k) in fac_b
        ff = id(L(f))
        for i in 1:k
          success, cc = _divides_hack(cc, ff)
          @assert success "element can not be converted to localization"
        end
      end
      return id_inv(simplify(cc))*inv(unit(fac_b))
    end
    (abort, coefficient) = _divides_hack(Q(a*last(powers_of_d)), Q(b))
    if !abort
      push!(powers_of_d, last(powers_of_d)^2)
    end
  end
  # find the minimal power that works
  upper = pop!(powers_of_d)
  lower = pop!(powers_of_d)
  while length(powers_of_d) > 0
    middle = lower*pop!(powers_of_d)
    (result, coefficient) = _divides_hack(Q(a*middle), Q(b))
    if result 
      upper = middle
    else 
      lower = middle
    end
  end
  return L(lift(coefficient), upper, check=false)
end

function _is_regular_fraction(R::RingType, p::MPolyRingElem, q::MPolyRingElem) where {RingType<:Union{MPolyQuoLocRing, MPolyLocRing, MPolyRing, MPolyQuoRing}}
  return divides(R(p), R(q))[1]
end

### Extensions for coherence
#function _is_regular_fraction(R::MPolyLocRing, p::MPolyRingElem, q::MPolyRingElem)
#  return divides(R(p), R(q))[1]
#end
#
#function _is_regular_fraction(R::MPolyRing, p::MPolyRingElem, q::MPolyRingElem)
#  return divides(p, q)[1]
#end
#
#function _is_regular_fraction(R::MPolyQuoRing, p::MPolyRingElem, q::MPolyRingElem)
#  return divides(R(p), R(q))[1]
#end

### arithmetic #########################################################
function +(a::T, b::T) where {T<:MPolyQuoLocRingElem}
  parent(a) === parent(b) || error("the arguments do not have the same parent ring")
  if lifted_denominator(a) == lifted_denominator(b) 
    return (parent(a))(lifted_numerator(a) + lifted_numerator(b), lifted_denominator(a), check=false)
  end
  gcd_ab = gcd([lifted_denominator(b), lifted_denominator(a)])
  p = divexact(lifted_denominator(a), gcd_ab)
  q = divexact(lifted_denominator(b), gcd_ab)
  new_den = p*lifted_denominator(b)
  return (parent(a))(lifted_numerator(a)*q + lifted_numerator(b)*p, new_den, check=false)
end


function -(a::T, b::T) where {T<:MPolyQuoLocRingElem}
  return a + (-b)
end

function *(a::T, b::T) where {T<:MPolyQuoLocRingElem}
  parent(a) === parent(b) || error("the arguments do not have the same parent ring")
  p = gcd([lifted_numerator(a), lifted_denominator(b)])
  q = gcd([lifted_numerator(b), lifted_denominator(a)])
  aa = divexact(lifted_numerator(a), p)
  bb = divexact(lifted_numerator(b), q)
  da = divexact(lifted_denominator(a), q)
  db = divexact(lifted_denominator(b), p)
  return (parent(a))(aa*bb, da*db, check=false)
end

function *(a::RET, b::MPolyQuoLocRingElem{BRT, BRET, RT, RET, MST}) where {BRT<:Ring, BRET<:RingElem, RT<:Ring, RET <: RingElem, MST}
  return (parent(b))(a*lifted_numerator(b), lifted_denominator(b), check=false)
end

function *(a::MPolyQuoLocRingElem{BRT, BRET, RT, RET, MST}, b::RET) where {BRT<:Ring, BRET<:RingElem, RT<:Ring, RET <: RingElem, MST}
  return b*a
end

function *(a::BRET, b::MPolyQuoLocRingElem{BRT, BRET, RT, RET, MST}) where {BRT<:Ring, BRET<:RingElem, RT<:Ring, RET <: RingElem, MST}
  return (parent(b))(a*lifted_numerator(b), lifted_denominator(b), check=false)
end

function *(a::MPolyQuoLocRingElem{BRT, BRET, RT, RET, MST}, b::BRET) where {BRT<:Ring, BRET<:RingElem, RT<:Ring, RET <: RingElem, MST}
  return b*a
end

### Why are the `//`-methods not implemented?
# Since a quotient ring Q = R/I of a polynomial ring R is not necessarily 
# factorial, it is difficult to decide, whether or not a and b have a 
# common factor g that can be canceled so that b'= b/g ∈  Q belongs 
# to the multiplicative set. Moreover, this would be the case if any 
# lift of b' belonged to S + I where S ⊂ R is the original multiplicative 
# set. Such containment can not easily be checked based only on the 
# functionality provided for S: Depending on the concrete type of 
# S, this task is algorithmically difficult, if not impossible.
#
# To remedy for this, we pursue the following pattern: 
#
# * Creation of elements [a]/[b] ∈ Q[S⁻¹] is possible only from 
#   representatives a/b ∈ R[S⁻¹] with b ∈ S.
# * The ring arithmetic attempts to cancel fractions which includes 
#   reduction modulo I of both the numerator and the denominator. 
#   This leads to representatives which would not be admissible 
#   for creation of elements in Q[S⁻¹].
# * Division routines can be used for the ring R[S⁻¹] with subsequent
#   conversion. 

function Base.:(/)(a::Oscar.IntegerUnion, b::MPolyQuoLocRingElem)
  success, c = divides(parent(b), b)
  !success && error("$b does not divide $a")
  return c
end

function Base.:(/)(a::T, b::T) where {T<:MPolyQuoLocRingElem}
  success, c = divides(a, b)
  !success && error("$b does not divide $a")
  return c
end

function divexact(a::Oscar.IntegerUnion, b::MPolyQuoLocRingElem; check::Bool=true)
  return a/b
end

function divexact(a::T, b::T; check::Bool=true) where {T<:MPolyQuoLocRingElem}
  return a/b
end

function ==(a::T, b::T) where {T<:MPolyQuoLocRingElem}
  parent(a) === parent(b) || error("the arguments do not have the same parent ring")
  a === b && return true
  return lifted_numerator(a)*lifted_denominator(b) - lifted_numerator(b)*lifted_denominator(a) in modulus(parent(a))
end

function ^(a::MPolyQuoLocRingElem, i::ZZRingElem)
  return parent(a)(lifted_numerator(a)^i, lifted_denominator(a)^i, check=false)
end

function ^(a::MPolyQuoLocRingElem, i::Integer)
  return parent(a)(lifted_numerator(a)^i, lifted_denominator(a)^i, check=false)
end

function isone(a::MPolyQuoLocRingElem) 
  return lifted_numerator(a) - lifted_denominator(a) in modulus(parent(a))
end

function iszero(a::MPolyQuoLocRingElem)
  return lift(a) in modulus(parent(a))
end

function iszero(a::MPolyQuoLocRingElem{<:Any, <:Any, <:Any, <:Any, <:MPolyComplementOfPrimeIdeal})
  # In case that the original quotient ring A is an integral domain
  # the localization map is injective and a is zero iff its numerator is zero.
  I = modulus(underlying_quotient(parent(a)))
  if get_attribute(I, :is_prime, false)
    return lifted_numerator(a) in I
  end
  return lift(a) in modulus(parent(a))
end

### enhancement of the arithmetic
function reduce_fraction(f::MPolyQuoLocRingElem{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST<:MPolyPowersOfElement}
  return f # Disable reduction here, because it slows down arithmetic.
  # return parent(f)(lift(simplify(numerator(f))), lifted_denominator(f), check=false)
end

# for local orderings, reduction does not give the correct result.
function reduce_fraction(f::MPolyQuoLocRingElem{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST<:MPolyComplementOfKPointIdeal}
  is_reduced(f) && return f
  return f
end

### implementation of Oscar's general ring interface
one(W::MPolyQuoLocRing) = W(one(base_ring(W)))
zero(W::MPolyQuoLocRing)= W(zero(base_ring(W)))

elem_type(::Type{MPolyQuoLocRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyQuoLocRingElem{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}
parent_type(::Type{MPolyQuoLocRingElem{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyQuoLocRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}


@doc raw"""
    bring_to_common_denominator(f::Vector{T}) where {T<:MPolyQuoLocRingElem}

Given a vector of fractions ``[a₁//b₁,…,aₙ//bₙ]`` return a pair 
``(d, λ)`` consisting of a common denominator ``d`` and a vector 
``λ = [λ₁,…,λₙ]`` such that ``aᵢ//bᵢ = λᵢ⋅aᵢ//d``.
"""
function bring_to_common_denominator(f::Vector{T}) where {T<:MPolyQuoLocRingElem}
  length(f) == 0 && error("need at least one argument to determine the return type")
  R = base_ring(parent(f[1]))
  for a in f
    R === base_ring(parent(a)) || error("elements do not belong to the same ring")
  end
  d = one(R)
  a = Vector{elem_type(R)}()
  for den in lifted_denominator.(f)
    b = gcd(d, den)
    c = divexact(den, b)
    e = divexact(d, b)
    d = d*c
    a = [c*k for k in a]
    push!(a, e)
  end
  return d, a
end

@doc raw"""
    write_as_linear_combination(f::T, g::Vector{T}) where {T<:MPolyLocRingElem} 

Write ``f = ∑ᵢ λᵢ⋅gᵢ`` for some ``λᵢ`` and return the vector ``[λ₁,…,λₙ]``.
"""
function write_as_linear_combination(
    f::RingElemType,
    g::Vector{RingElemType}
  ) where {RingElemType<:MPolyQuoLocRingElem}
  n = length(g)
  L = parent(f)
  W = localized_ring(L)
  for a in g 
    parent(a) === L || error("elements do not belong to the same ring")
  end
  return L.(_vec(coordinates(lift(f), ideal(L, g)))[1:length(g)]) # temporary hack; to be replaced.
end

write_as_linear_combination(f::MPolyQuoLocRingElem, g::Vector) = write_as_linear_combination(f, parent(f).(g))


########################################################################
# Homomorphisms of quotients of localized polynomial algebras          #
########################################################################
# 
# Suppose we are given two localizations of polynomial algebras 
# by means of commutative diagrams 
#
#       R   →    P = R/I
#       ↓        ↓ 
# V = R[T⁻¹] →  P[T⁻¹]
#
# and 
#
#       S   →    Q = S/J
#       ↓        ↓ 
# W = S[U⁻¹] →  Q[U⁻¹].
#
# Lemma:
# For any homomorphism φ : P[T⁻¹] → Q[U⁻¹] the following holds. 
#
#             φ
#     P[T⁻¹]  →  Q[U⁻¹]
#       ↑          ↑
#     R[T⁻¹] --> S[U⁻¹]
#       ↑    ↗ ψ   ↑ ι
#       R     →  S[c⁻¹]
#             η    ↑ κ
#                  S
#
# a) The composition of maps R → Q[U⁻¹] completely determines φ by 
#    the images xᵢ ↦ [aᵢ]/[bᵢ] with aᵢ ∈ S, bᵢ ∈ U.
# b) Let ψ : R → S[U⁻¹] be the map determined by some choice of 
#    the images xᵢ↦ aᵢ/bᵢ as above. Then ψ extends to a map 
#    R[T⁻¹] → S[U⁻¹] if and only if 
#    
#       for all t ∈ T : ψ(t) ∈ U.
#
#    This is not necessarily the case as the lift of images 
#    φ(t) ∈ Q[U⁻¹] in S[U⁻¹] need only be elements of U + J.
# c) Choosing a common denominator c for all ψ(xᵢ), we obtain a 
#    ring homomorphism η : R → S[c⁻¹] such that ψ = ι ∘ η.
#
# Upshot: In order to describe φ, we may store some homomorphism 
#     
#       ψ : R → S[U⁻¹] 
#
# lifting it and keep in mind the ambiguity of choices for such ψ.
# The latter point c) will be useful for reducing to a homomorphism 
# of finitely generated algebras.

@doc raw"""
    MPolyQuoLocalizedRingHom{
         DomainType<:MPolyQuoLocRing, 
         CodomainType<:Ring, 
         RestrictedMapType<:Map
        } <: AbsLocalizedRingHom{
                                 DomainType, 
                                 CodomainType, 
                                 RestrictedMapType
                                }

Homomorphisms of localizations of affine algebras 

  ``ϕ : (𝕜[x₁,…,xₘ]/I)[S⁻¹] → (𝕜[y₁,…,yₙ]/J)[T⁻¹]``

of types `MPolyQuoLocRing{BaseRingType, BaseRingElemType, RingType, RingElemType, DomainMultSetType}` and `MPolyQuoLocRing{BaseRingType, BaseRingElemType, RingType, RingElemType, CodomainMultSetType}`.
These are completely determined by the images of the 
variables ``ϕ(xᵢ) ∈ (𝕜[y₁,…,yₙ]/J)[T⁻¹]`` so that the 
constructor takes as input the triple 
``((𝕜[x₁,…,xₘ]/I)[S⁻¹], (𝕜[y₁,…,yₙ]/J)[T⁻¹], [ϕ(x₁),…,ϕ(xₘ)])``.
"""
@attributes mutable struct MPolyQuoLocalizedRingHom{
                                     DomainType<:MPolyQuoLocRing, 
                                     CodomainType<:Ring, 
                                     RestrictedMapType<:Map
                                    } <: AbsLocalizedRingHom{
                                                             DomainType, 
                                                             CodomainType, 
                                                             RestrictedMapType
                                                            }
  domain::DomainType
  codomain::CodomainType
  res::RestrictedMapType

  function MPolyQuoLocalizedRingHom(
      L::DomainType,
      S::CodomainType,
      res::RestrictedMapType;
      check::Bool=true
    ) where {DomainType<:MPolyQuoLocRing, CodomainType<:Ring, RestrictedMapType<:Map}
    R = base_ring(L)
    R === domain(res) || error("restriction map is not compatible")
    U = inverted_set(L)
    @check begin
      for f in U
        is_unit(S(res(f))) || error("map is not well defined")
      end
      for g in gens(modulus(underlying_quotient(L)))
        iszero(S(res(g))) || error("map is not well defined")
      end
    end
    return new{DomainType, CodomainType, RestrictedMapType}(L, S, res)
  end
end

### type getters 
domain_type(::Type{MPolyQuoLocalizedRingHom{D, C, M}}) where {D, C, M} = D
domain_type(f::MPolyQuoLocalizedRingHom) = domain_type(typeof(f))
codomain_type(::Type{MPolyQuoLocalizedRingHom{D, C, M}}) where {D, C, M} = C
codomain_type(f::MPolyQuoLocalizedRingHom) = domain_type(typeof(f))
restricted_map_type(::Type{MPolyQuoLocalizedRingHom{D, C, M}}) where {D, C, M} = M
restricted_map_type(f::MPolyQuoLocalizedRingHom) = domain_type(typeof(f))

morphism_type(::Type{R}, ::Type{S}) where {R<:MPolyQuoLocRing, S<:Ring} = MPolyQuoLocalizedRingHom{R, S, morphism_type(base_ring_type(R), S)}
morphism_type(L::MPolyQuoLocRing, S::Ring) = morphism_type(typeof(L), typeof(S))

### TODO: Move to other file
morphism_type(::Type{R}, ::Type{S}) where {R<:MPolyRing, S<:Ring} = Oscar.MPolyAnyMap{R, S, Nothing, elem_type(S)}
morphism_type(R::MPolyRing, S::Ring) = morphism_type(typeof(R), typeof(S))

### required getter functions
domain(f::MPolyQuoLocalizedRingHom) = f.domain
codomain(f::MPolyQuoLocalizedRingHom) = f.codomain

@doc raw"""
    restricted_map(f::MPolyQuoLocalizedRingHom)

Given a homomorphism ``ϕ : (𝕜[x₁,…,xₘ]/I)[U⁻¹] → S``, return
the canonically associated map ``ϕ' : 𝕜[x₁,…,xₘ] → S``.
"""
restricted_map(f::MPolyQuoLocalizedRingHom) = f.res

function images(f::MPolyQuoLocalizedRingHom{<:Any, <:MPolyQuoLocRing})
  return lift.((codomain(f)).(restricted_map(f).(gens(base_ring(domain(f))))))
end

### additional constructors
function MPolyQuoLocalizedRingHom(
    L::MPolyQuoLocRing,
    S::Ring,
    a::Vector{T};
    check::Bool=true
  ) where {T<:RingElem}
  return MPolyQuoLocalizedRingHom(L, S, hom(base_ring(L), S, a, check=check), check=check)
end

hom(L::MPolyQuoLocRing, S::Ring, res::Map; check::Bool=true) = MPolyQuoLocalizedRingHom(L, S, res, check=check)

function hom(L::MPolyQuoLocRing, S::Ring, a::Vector{T}; check::Bool=true) where {T<:RingElem}
  R = base_ring(L)
  res = hom(R, S, a, check=check)
  MPolyQuoLocalizedRingHom(L, S, res, check=check)
end

### implementing the Oscar map interface
function identity_map(W::T) where {T<:MPolyQuoLocRing} 
  MPolyQuoLocalizedRingHom(W, W, identity_map(base_ring(W)))
end

function simplify(a::MPolyQuoLocRingElem)
  p = simplify(numerator(a))
  return parent(a)(lift(p), lifted_denominator(a); check=false)
end

### we need to overwrite the following method because of the 
# uncommon implementation of the numerator and denominator methods
function (f::MPolyQuoLocalizedRingHom)(a::AbsLocalizedRingElem)
  parent(a) === domain(f) || return f(domain(f)(a))
  isone(lifted_denominator(a)) && return codomain(f)(restricted_map(f)(lifted_numerator(a)))
  if total_degree(lifted_denominator(a)) > 10
    res = restricted_map(f)
    img_num = res(lifted_numerator(a))
    den = lifted_denominator(a)
    img_den = one(img_num)
    fac_den = factor(den)
    for (a, k) in fac_den
      img_den = img_den * inv(res(a))^k
    end
    img_den = img_den * inv(res(unit(fac_den)))
    return img_num * img_den
  end
  b = a #simplify(a)
  return codomain(f)(restricted_map(f)(lifted_numerator(b)))*inv(codomain(f)(restricted_map(f)(lifted_denominator(b))))
end

function compose(
    f::MPolyQuoLocalizedRingHom, 
    g::Map{<:Ring, <:Ring}
  )
  codomain(f) === domain(g) || error("maps are not compatible")

  res = restricted_map(f)
  b = res(zero(base_ring(domain(f))))
  if parent(b) === domain(g)
    return MPolyQuoLocalizedRingHom(domain(f), codomain(g), compose(res, g), check=false)
  else
    new_res = MapFromFunc(base_ring(domain(f)), codomain(g), x->g(domain(g)(res(x))))
    return MPolyQuoLocalizedRingHom(domain(f), codomain(g), new_res, check=false)
  end
#  if codomain(restricted_map(f)) === domain(g)
#    return MPolyQuoLocalizedRingHom(domain(f), codomain(g), compose(restricted_map(f), g))
#  elseif codomain(restricted_map(f)) === base_ring(domain(g)) 
#    h = hom(base_ring(domain(g)), domain(g), domain(g).(gens(base_ring(domain(g)))))
#    return MPolyQuoLocalizedRingHom(domain(f), codomain(g), compose(compose(restricted_map(f), h), g))
#  end
  ### The fallback version. Careful: This might not carry over maps on the coefficient rings!
  R = base_ring(domain(f))
  return MPolyQuoLocalizedRingHom(domain(f), codomain(g), compose(res, g), check=false)
  #return MPolyQuoLocalizedRingHom(domain(f), codomain(g), hom(R, codomain(g), [g(f(x)) for x in gens(R)], check=false), check=false)
end

(f::MPolyQuoLocalizedRingHom)(I::Ideal) = ideal(codomain(f), f.(domain(f).(gens(I))))

function ==(f::MPolyQuoLocalizedRingHom, g::MPolyQuoLocalizedRingHom) 
  f === g && return true
  domain(f) === domain(g) || return false
  codomain(f) === codomain(g) || return false
  for x in gens(base_ring(domain(f)))
    f(x) == g(x) || return false
  end
  return true
end

function Base.hash(f::MPolyQuoLocalizedRingHom, h::UInt)
  b = 0xd6d389598ad28724  % UInt
  h = hash(domain(f), h)
  h = hash(codomain(f), h)
  for x in gens(base_ring(domain(f)))
    h = hash(f(x), h)
  end
  return xor(h, b)
end

### printing
function Base.show(io::IO, ::MIME"text/plain", phi::MPolyQuoLocalizedRingHom)
  io = pretty(io)
  println(terse(io), phi)
  print(io, Indent())
  println(io, "from ", Lowercase(), domain(phi))
  println(io, "to ", Lowercase(), codomain(phi))
  println(io, Dedent(), "defined by", Indent())
  R = base_ring(domain(phi))
  psi = restricted_map(phi)
  to = is_unicode_allowed() ? " ↦ " : " -> "
  for i in 1:ngens(R)-1
    println(io, R[i], to, psi(R[i]))
  end
  n = ngens(R)
  print(io, R[n], to, psi(R[n]))
  print(io, Dedent())
end

function Base.show(io::IO, phi::MPolyQuoLocalizedRingHom)
  if is_terse(io)
    print(io, "Ring homomorphism")
  else
    R = base_ring(domain(phi))
    psi = restricted_map(phi)
    io = pretty(io)
    io = terse(io)
    print(io, "hom: ", domain(phi))
    if is_unicode_allowed()
      print(io, " → ")
    else
      print(io, " -> ")
    end
    print(io, codomain(phi))
  end
end

### helper_ring
# Sets up the ring S[c⁻¹] from the Lemma.
function helper_ring(f::MPolyQuoLocalizedRingHom{<:Any, <:MPolyQuoLocRing})
  if !has_attribute(f, :helper_ring)
    minimal_denominators = Vector{elem_type(base_ring_type(domain(f)))}()
    R = base_ring(domain(f))
    S = base_ring(codomain(f))
    p = one(S)

    for d in [denominator(y) for y in images(f)]
      g = gcd(d, p)
      d_min = divexact(d, g)
      push!(minimal_denominators, d)
      p = p*d_min
    end
    set_attribute!(f, :minimal_denominators, minimal_denominators)

    help_ring, help_kappa, theta = _add_variables(S, [:θ])
    set_attribute!(f, :helper_ring, help_ring)
    kappa = help_kappa
    set_attribute!(f, :kappa, help_kappa)
    c_inv = theta[1]
    helper_images = [kappa(numerator(y))*c_inv*kappa(divexact(p, denominator(y))) for y in images(f)]
    set_attribute!(f, :helper_images, helper_images)
    eta = hom(R, help_ring, helper_images, check=false)
    set_attribute!(f, :eta, eta)
  end
  return get_attribute(f, :helper_ring)::base_ring_type(domain(f))
end

function helper_images(
    f::MPolyQuoLocalizedRingHom{<:Any, <:MPolyQuoLocRing}
  )
  if !has_attribute(f, :helper_images) 
    helper_ring(f)
  end
  return get_attribute(f, :helper_images)::Vector{elem_type(base_ring_type(domain(f)))}
end

function minimal_denominators(
    f::MPolyQuoLocalizedRingHom{<:Any, <:MPolyQuoLocRing}
  )
  if !has_attribute(f, :minimal_denominators) 
    helper_ring(f)
  end
  return get_attribute!(f, :minimal_denominators)::Vector{elem_type(base_ring_type(domain(f)))}
end

function helper_eta(
    f::MPolyQuoLocalizedRingHom{<:Any, <:MPolyQuoLocRing}
  )
  if !has_attribute(f, :eta) 
    helper_ring(f)
  end
  return get_attribute(f, :eta)::morphism_type(base_ring_type(domain(f)), base_ring_type(domain(f)))
end

function helper_kappa(
    f::MPolyQuoLocalizedRingHom{<:Any, <:MPolyQuoLocRing}
  )
  if !has_attribute(f, :kappa) 
    helper_ring(f)
  end
  return get_attribute(f, :kappa)::morphism_type(base_ring_type(domain(f)), base_ring_type(domain(f)))
end

function common_denominator(
    f::MPolyQuoLocalizedRingHom{<:Any, <:MPolyQuoLocRing}
  )
  if !has_attribute(f, :minimal_denominators) 
    helper_ring(f)
  end
  d = get_attribute(f, :minimal_denominators)::Vector{elem_type(base_ring_type(domain(f)))}
  return (length(d) == 0 ? one(base_ring(codomain(f))) : prod(d))
end

function helper_ideal(
    f::MPolyQuoLocalizedRingHom{<:Any, <:MPolyQuoLocRing}
  )
  Sc = helper_ring(f)
  return ideal(Sc, one(Sc)-last(gens(Sc))*helper_kappa(f)(common_denominator(f)))
end

# return the localized ring as a quotient of a polynomial ring using Rabinowitsch's trick.
function as_affine_algebra(
    L::MPolyQuoLocRing{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement};
    inverse_name::VarName=:_0
  )
  R = base_ring(L)
  A, phi, t = _add_variables_first(R, [Symbol(inverse_name)])
  theta = t[1]
  f = prod(denominators(inverted_set(L)); init=one(R))
  I = ideal(A, [phi(g) for g in gens(modulus(underlying_quotient(L)))]) + ideal(A, [one(A)-theta*phi(f)])
  return A, I, f, phi, theta
end

# needed for instance to compute kernels
# adds a single extra variable to turn the localization into an affine_algebra
# return the isomorphism L -> SomeAffineAlgebra
function _as_affine_algebra(
    L::MPolyQuoLocRing{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement};
    inverse_name::VarName=:_0
  )
  R = base_ring(L)
  A, phi, t = _add_variables_first(R, [Symbol(inverse_name)])
  theta = t[1]
  f = prod(denominators(inverted_set(L)); init=one(R))
  I = ideal(A, [phi(g) for g in gens(modulus(underlying_quotient(L)))]) + ideal(A, [one(A)-theta*phi(f)])
  Q, _ = quo(A, I)
  id = hom(L, Q, gens(A)[2:end], check=false)
  id_inv = hom(Q, L, pushfirst!(gens(L), inv(L(f))), check=false)
  set_attribute!(id, :inverse, id_inv)
  set_attribute!(id_inv, :inverse, id)
  return id
end

@attr MPolyIdeal function kernel(
    f::MPolyAnyMap{<:MPolyRing, 
                   <:MPolyQuoLocRing{<:Any, <:Any, <:Any, <:Any,
                                     <:MPolyPowersOfElement}})
  R = domain(f)
  W = codomain(f)
  I = saturated_ideal(modulus(W))
  P = base_ring(W)
  img_gens = f.(gens(R))
  nums = lifted_numerator.(img_gens)
  denoms = lifted_denominator.(img_gens)

  # Build up a helper ring for the graph of f using Rabinowitschs trick.
  inverse_name=:_0
  r = length(denoms)
  kk = coefficient_ring(R)
  A, t = polynomial_ring(kk, vcat([Symbol(inverse_name,k) for k in 1:r],
                                  symbols(P), symbols(R)); cached=false)
  r = length(denoms)
  theta = t[1:r]
  n = ngens(P)
  imgs_y = t[r+1:(r+n)]
  imgs_x = t[r+n+1:end]
  # Sometimes for unnecessarily complicated sets of generators for I the computation 
  # wouldn't finish. We try to pass to a `small_generating_set` to hopefully reduce the dependency 
  # on a particular set of generators. 
  J = ideal(A, vcat([one(A) - theta[i]*evaluate(den, imgs_y) for (i, den) in enumerate(denoms)], # Rabinowitsch relations
                    [theta[i]*evaluate(num, imgs_y) - imgs_x[i] for (i, num) in enumerate(nums)], # Graph relations
                    [evaluate(g, imgs_y) for g in small_generating_set(I)])) # codomain's modulus
  # We eliminate the Rabinowitsch variables first, the codomain variables second, 
  # and finally get to the domain variables. This elimination should be quicker 
  # than one which does not know the Rabinowitsch property.
  oo = degrevlex(theta)*degrevlex(imgs_y)*degrevlex(imgs_x)
  #oo = lex(theta)*lex(imgs_y)*lex(imgs_x)
  gb = groebner_basis(J, ordering=oo)

  # TODO: Speed up and use build context.
  res_gens = elem_type(A)[f for f in gb if all(e -> is_zero(view(e, 1:(n+r))), AbstractAlgebra.exponent_vectors(f))]
  img_gens2 = vcat([zero(R) for i in 1:(n+r)], gens(R))
  result = ideal(R, elem_type(R)[evaluate(g, img_gens2) for g in res_gens])
  return result
  
  # deprecated code below
  id, _ = _as_affine_algebra_with_many_variables(codomain(f))
  g = hom(domain(f), codomain(id), id.(f.(gens(domain(f)))))
  return K
end

@attr MPolyQuoIdeal function kernel(
    f::MPolyAnyMap{<:MPolyQuoRing, 
                   <:MPolyQuoLocRing{<:Any, <:Any, <:Any, <:Any,
                                     <:MPolyPowersOfElement}})
  A = domain(f)
  R = base_ring(A)
  g = hom(R, codomain(f), f.(gens(A)); check=false)
  K = kernel(g)
  return ideal(A, elem_type(A)[h for h in A.(gens(K)) if !is_zero(h)])
  id, _ = _as_affine_algebra_with_many_variables(codomain(f))
  g = hom(domain(f), codomain(id), id.(f.(gens(domain(f)))); check=false)
  return kernel(g)
end

### The following method is also required for the internals of the generic 
# kernel routine for localized rings.
@attr MPolyIdeal function kernel(f::MPolyAnyMap{<:MPolyRing, <:MPolyQuoLocRing})
  P = domain(f)
  L = codomain(f)
  I = ideal(L, zero(L))
  R = base_ring(L)
  J = saturated_ideal(I)
  d = [lifted_denominator(g) for g in f.(gens(domain(f)))]
  S = simplify(MPolyPowersOfElement(R, d))
  W = MPolyQuoLocRing(R, modulus(underlying_quotient(L)), S)
  id, _ =  _as_affine_algebra_with_many_variables(W)
  A = codomain(id)
  h = hom(P, A, elem_type(A)[id(W(f(x), check=false)) for x in gens(P)], check=false)
  gg = Vector{elem_type(A)}(id.(W.(gens(J))))
  return preimage(h, ideal(A, gg))
end

@attr MPolyQuoIdeal function kernel(f::MPolyAnyMap{<:MPolyQuoRing, <:MPolyQuoLocRing})
  A = domain(f)
  R = base_ring(A)
  ff = hom(R, codomain(f), f.(gens(A)), check=false)
  K = kernel(ff)
  return ideal(A, [g for g in A.(gens(K)) if !iszero(g)])
end

function is_isomorphism(
    phi::MPolyQuoLocalizedRingHom{T, T}
  ) where {T<:MPolyQuoLocRing}
  if has_attribute(phi, :inverse)
    return true
  end
  K = domain(phi)
  L = codomain(phi)
  A, I, d1, inc1, theta1 = as_affine_algebra(K, inverse_name=:s)
  B, J, d2, inc2, theta2 = as_affine_algebra(L, inverse_name=:t)

  # write the denominators of the images of the variables xᵢ of K as 
  # polynomials in the variables yⱼ of L and the localization variable t:
  #
  # TODO: This can most probably be fine tuned as follows. Make sure that 
  # the ideal J is of the form J' + ⟨1-t⋅g(y)⟩ where J' is generated by a 
  # Groebner basis for the saturated ideal of L in 𝕜[y]. Then choose an 
  # elimination ordering on t on 𝕜[t,y] which coincides on the y-variables 
  # with the ordering that had been used to compute the Groebner-basis of 
  # J'. Then the generators above for J are also already a Groebner basis 
  # since J' was already saturated. 
  #
  # Now a Groebner-basis computation for the ideal J_ext should proceed 
  # quickly by multiplying q with subsequently higher powers of t and reduction 
  # modulo J until a constant polynomial is attained. 
  #
  # With this strategy, the Groebner-basis computation should not even be 
  # expensive. 
  denoms = Vector{elem_type(B)}()
  for f in images(phi)
    J_ext = ideal(B, push!(gens(J), inc2(denominator(f))))
    G, M = standard_basis_with_transformation_matrix(J_ext)
    gen(G, 1)==one(B) || error("the denominator is not a unit in the target ring")
    push!(denoms, last(collect(M)))
  end

  # write the images as polynomials in B 
  imagesB = [ inc2(numerator(images(phi)[i]))*denoms[i] for i in (1:length(denoms))]

  # expand the localization variable s in A as a polynomial in B
  denoms = Vector{elem_type(B)}()
  for h in denominators(inverted_set(K))
    phi_h = phi(h)
    p = lifted_numerator(phi_h)
    q = lifted_denominator(phi_h)
    J_ext = ideal(B, push!(gens(J), inc2(p)))
    G, M = standard_basis_with_transformation_matrix(J_ext)
    gen(G, 1)==one(B) || error("the denominator is not a unit in the target ring")
    push!(denoms, inc2(q)*last(collect(M)))
  end
  pushfirst!(imagesB, prod(denoms; init=one(B)))

  # perform a sanity check
  phiAB = hom(A, B, imagesB, check=false)
  issubset(ideal(B, [phiAB(g) for g in gens(I)]), J) || error("the homomorphism is not well defined")

  # assemble a common ring in which the equations for the graph of phi can 
  # be realized.
  C, j1, B_vars = _add_variables_first(A, symbols(B))
  j2 = hom(B, C, B_vars, check=false)
  G = ideal(C, [j1(gen(A, i)) - j2(imagesB[i]) for i in 1:ngens(A)]) + ideal(C, j2.(gens(J))) + ideal(C, j1.(gens(I)))
  singC, _ = Singular.polynomial_ring(Oscar.singular_coeff_ring(base_ring(C)), 
            symbols(C),
            ordering=Singular.ordering_dp(1)
                    *Singular.ordering_dp(nvars(B)-1)
                    *Singular.ordering_dp(1)
                    *Singular.ordering_dp(nvars(A)-1))
  # TODO: adjust this to the orderings used for the previous groebner basis 
  # computations in A and B once such things are respected. 
  singG = Singular.Ideal(singC, singC.(gens(G)))
  stdG = Singular.std(singG) 
  
  # Compute the inverse images of the variables of L
  m = nvars(A)-1
  n = nvars(B)-1
  R = base_ring(K)
  V = localized_ring(K)
  pre_images = Vector{elem_type(V)}()
  #pre_imagesA = Vector{elem_type(A)}()
  # the first variable needs special treatment
  #singp = Singular.reduce(gen(singC, 1), stdG)
  #singp < gen(singC, n+1) || return false
  #p = C(singp)
  #push!(pre_imagesA, evaluate(p, vcat([zero(A) for i in 0:n], gens(A))))
  for i in 1:n
    singp = Singular.reduce(gen(singC, i+1), stdG)
    singp < gen(singC, n+1) || return false
    p = C(singp)
    # Write p as an element in the very original ring
    #push!(pre_imagesA, evaluate(p, vcat([zero(A) for i in 0:n], gens(A))))
    push!(pre_images, evaluate(p, vcat([zero(V) for i in 0:n], [V(one(R), d1)], V.(gens(R)))))
  end

  invJ = ideal(A, [(p < gen(singC, n+1) ? evaluate(C(p), vcat([zero(A) for i in 0:n], gens(A))) : zero(A)) for p in gens(stdG)])
  # TODO: invJ is already a Groebner basis, but only for the ordering used 
  # in the above elimination.
  # Make sure, this ordering is used again for the sanity check below!
  invJ == I || return false

  set_attribute!(phi, :inverse, MPolyQuoLocalizedRingHom(L, K, pre_images))
  psi = get_attribute(phi, :inverse)
  set_attribute!(psi, :inverse, phi)
  return true
end

function inverse(
    f::MPolyQuoLocalizedRingHom{
                                <:MPolyQuoLocRing{<:Any, <:Any, <:Any, <:Any, 
                                                        <:MPolyPowersOfElement
                                                       },
                                <:MPolyQuoLocRing{<:Any, <:Any, <:Any, <:Any, 
                                                        <:MPolyPowersOfElement
                                                       }
                               }
  )
  is_isomorphism(f) || error("the given morphism is not an isomorphism")
  return get_attribute(f, :inverse)::morphism_type(codomain(f), domain(f))
end



### adds the variables with names specified in v to the polynomial 
# ring R and returns a triple consisting of the new ring, the embedding 
# of the original one, and a list of the new variables. 
function _add_variables(R::RingType, v::Vector{<:VarName}) where {RingType<:MPolyRing}
  ext_R, _ = polynomial_ring(coefficient_ring(R), vcat(symbols(R), Symbol.(v)); cached = false)
  n = ngens(R)
  phi = hom(R, ext_R, gens(ext_R)[1:n], check=false)
  return ext_R, phi, gens(ext_R)[(n+1):ngens(ext_R)]
end

function _add_variables_first(R::RingType, v::Vector{<:VarName}) where {RingType<:MPolyRing}
  ext_R, _ = polynomial_ring(coefficient_ring(R), vcat(Symbol.(v), symbols(R)); cached = false)
  n = ngens(R)
  phi = hom(R, ext_R, gens(ext_R)[1+length(v):n+length(v)], check=false)
  return ext_R, phi, gens(ext_R)[(1:length(v))]
end

@doc raw"""
    simplify(L::MPolyQuoLocRing{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement})

Use `elimpart` from the `Singular` library `Presolve.lib` to simplify the presentation 
of `L` by eliminating superfluous variables; return a triple ``(L', f, g)`` where 
``L' ≅ L`` and ``f : L ↔ L' : g`` are the identifying isomorphisms.
"""
function simplify(L::MPolyQuoLocRing{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement})
  W = localized_ring(L)
  I = modulus(L)
  J = modulus(underlying_quotient(L))
  R = base_ring(L)
  SR = singular_poly_ring(R)
  SJ = singular_generators(J)

  # collect the output from elimpart in Singular
  l = Singular.LibPresolve.elimpart(SJ)

  # set up the ring with the fewer variables 
  kept_var_symb = [symbols(R)[i] for i in 1:ngens(R) if !iszero(l[4][i])]
  Rnew, new_vars = polynomial_ring(coefficient_ring(R), kept_var_symb; cached = false)

  # and the maps to go back and forth
  subst_map_R = hom(R, R, R.(gens(l[5])), check=false)
  imgs = Vector{elem_type(Rnew)}()
  j = 1
  for i in 1:ngens(R)
    if !iszero(l[4][i])
      push!(imgs, gen(Rnew, j))
      j = j+1
    else
      push!(imgs, zero(Rnew))
    end
  end
  proj_map = hom(R, Rnew, imgs, check=false)

  # the full substitution map 
  f = compose(subst_map_R, proj_map)

  # the transformed ideal
  Jnew = ideal(Rnew, f.(gens(J)))

  # the translated inverted set
  U = inverted_set(L)
  Unew = MPolyPowersOfElement(Rnew, f.(denominators(U)))

  # the new localized ring
  Lnew = MPolyQuoLocRing(Rnew, Jnew, Unew)

  # the localized map and its inverse
  floc = hom(L, Lnew, Lnew.(f.(gens(R))), check=false)
  flocinv = hom(Lnew, L, [L(R(a)) for a in gens(l[4]) if !iszero(a)], check=false)

  return Lnew, floc, flocinv
end

function simplify(L::MPolyLocRing{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement})
  Lnew = MPolyLocRing(base_ring(L), inverted_set(L))
  R = base_ring(L)
  finv = hom(R, L, gens(L), check=false)
  f = hom(R, Lnew, gens(Lnew), check=false)
  return Lnew, hom(L, Lnew, f), hom(Lnew, L, finv, check=false)
end

function simplify(L::MPolyQuoRing)
  J = modulus(L)
  R = base_ring(L)
  is_zero(ngens(R)) && return L, identity_map(L), identity_map(L)
  SR = singular_poly_ring(R)
  SJ = singular_generators(J)

  # collect the output from elimpart in Singular
  l = Singular.LibPresolve.elimpart(SJ)

  # set up the ring with the fewer variables 
  kept_var_symb = [symbols(R)[i] for i in 1:ngens(R) if !iszero(l[4][i])]
  Rnew, new_vars = polynomial_ring(coefficient_ring(R), kept_var_symb; cached=false)

  # and the maps to go back and forth
  subst_map_R = hom(R, R, R.(gens(l[5])), check=false)
  imgs = Vector{elem_type(Rnew)}()
  j = 1
  for i in 1:ngens(R)
    if !iszero(l[4][i])
      push!(imgs, gen(Rnew, j))
      j = j+1
    else
      push!(imgs, zero(Rnew))
    end
  end
  proj_map = hom(R, Rnew, imgs, check=false)

  # the full substitution map 
  f = compose(subst_map_R, proj_map)

  # the transformed ideal
  Jnew = ideal(Rnew, f.(gens(J)))
  Lnew, _ = quo(Rnew, Jnew)

  # the inverse of the identification map
  fres = hom(L, Lnew, Lnew.(f.(gens(R))), check=false)
  fresinv = hom(Lnew, L, [L(R(a)) for a in gens(l[4]) if !iszero(a)], check=false)

  return Lnew, fres, fresinv
end

function simplify(R::MPolyRing)
  Rnew, new_vars = polynomial_ring(coefficient_ring(R), symbols(R); cached=false)
  f = hom(R, Rnew, gens(Rnew), check=false)
  finv = hom(Rnew, R, gens(R), check=false)
  return Rnew, f, finv
end

@doc raw"""
    MPolyQuoLocalizedIdeal{
        LocRingType<:MPolyQuoLocRing, 
        LocRingElemType<:MPolyQuoLocRingElem
      } <: AbsLocalizedIdeal{LocRingElemType}

Ideals in localizations of affine algebras.
"""
@attributes mutable struct MPolyQuoLocalizedIdeal{
     LocRingType<:MPolyQuoLocRing, 
     LocRingElemType<:MPolyQuoLocRingElem, 
     MPolyLocalizedIdealType<:MPolyLocalizedIdeal
    } <: AbsLocalizedIdeal{LocRingElemType}
  # the initial set of generators, not to be changed ever!
  gens::Vector{LocRingElemType}
  # the ambient ring for this ideal
  W::LocRingType

  # fields for caching 
  map_from_base_ring::Map

  J::MPolyLocalizedIdealType
 
  function MPolyQuoLocalizedIdeal(
      W::MPolyQuoLocRing, 
      g::Vector{LocRingElemType};
      map_from_base_ring::Map = MapFromFunc(
          base_ring(W), 
          W,
          W,
          y->(isone(lifted_denominator(y)) ? lifted_numerator(y) : divexact(lifted_numerator(y), lifted_denominator(y))),
        )
    ) where {LocRingElemType<:MPolyQuoLocRingElem}
    for f in g
      parent(f) === W || error("generator is not an element of the given ring")
    end

    L = localized_ring(W)
    J = ideal(L, vcat(lift.(g), gens(modulus(W))))
    I = new{typeof(W), LocRingElemType, typeof(J)}()
    I.gens = g
    I.W = W
    I.map_from_base_ring = map_from_base_ring
    I.J = J
    return I
  end
end
 
### required getter functions
gens(I::MPolyQuoLocalizedIdeal) = copy(I.gens)
gen(I::MPolyQuoLocalizedIdeal, i::Int) = I.gens[i]
base_ring(I::MPolyQuoLocalizedIdeal) = I.W
base_ring_type(::Type{MPolyQuoLocalizedIdeal{LRT, LET, MPT}}) where {LRT, LET, MPT} = LRT

### additional getter functions 
map_from_base_ring(I::MPolyQuoLocalizedIdeal) = I.map_from_base_ring
pre_image_ideal(I::MPolyQuoLocalizedIdeal) = I.J
number_of_generators(I::MPolyQuoLocalizedIdeal) = length(I.gens)

### a shorthand notation for any MPolyIdeal 
const MPolyAnyIdeal = Union{MPolyIdeal, MPolyQuoIdeal,
                 MPolyLocalizedIdeal, MPolyQuoLocalizedIdeal
                }


### Additional constructors
@doc raw"""
    intersect(I::MPolyQuoLocalizedIdeal, J::MPolyQuoLocalizedIdeal)
    intersect(I::MPolyQuoLocalizedIdeal, J::MPolyQuoLocalizedIdeal...)
    intersect(VI::Vector{<:MPolyQuoLocalizedIdeal{T}}) where T
Return the intersection of two or more ideals.

# Examples
```jldoctest
julia> R,(x,y,z,w) = QQ[:x, :y, :z, :w];

julia> Q = ideal(R,[x*y-z*w]);

julia> RQ,phiQ = quo(R,Q);

julia> T = complement_of_point_ideal(R,[0,0,0,0]);

julia> RQL, phiQL = localization(RQ,T);

julia> I = ideal(RQL,RQL.([x,z]))
Ideal generated by
  x
  z

julia> J = ideal(RQL,RQL.([y]))
Ideal generated by
  y

julia> intersect(I,J)
Ideal generated by
  z*w
  y*z
  x*y

julia> K = intersect([I,J])
Ideal generated by
  z*w
  y*z
  x*y

julia> (I,J,K)
(Ideal (x, z), Ideal (y), Ideal (z*w, y*z, x*y))

```
"""
function intersect(I::MPolyQuoLocalizedIdeal, J::MPolyQuoLocalizedIdeal)
  L = base_ring(I)
  L === base_ring(J) || error("ideals must be defined in the same ring")
  preI = pre_image_ideal(I)
  preJ = pre_image_ideal(J)
  R = base_ring(L)
  K = intersect(preI, preJ)
  return L(K)
end

function intersect(I::MPolyQuoLocalizedIdeal, J::MPolyQuoLocalizedIdeal...)
  L = base_ring(I)
  erg = pre_image_ideal(I)
  for K in J
    base_ring(K) === L || error("base rings must match")
    erg = intersect(erg,pre_image_ideal(K))
  end
  return L(erg)
end

function intersect(VI::Vector{<:MPolyQuoLocalizedIdeal{T}}) where T
  @assert length(VI) != 0
  L = base_ring(VI[1])
  all(J -> base_ring(J) === L,VI) || error("base rings must match")
  VIpre = [pre_image_ideal(J) for J in VI]
  erg = intersect(VIpre)
  return L(erg)
end

### Basic functionality
function ideal_membership(a::RingElem, I::MPolyQuoLocalizedIdeal)
  L = base_ring(I)
  parent(a) === L || return L(a) in I
  return lift(a) in pre_image_ideal(I)
end

function coordinates(a::RingElem, I::MPolyQuoLocalizedIdeal)
  L = base_ring(I)
  parent(a) === L || return coordinates(L(a), I)
  a in I || error("the given element is not in the ideal")
  x = coordinates(lift(a), pre_image_ideal(I), check=false)
  return map_entries(L, x[1:1, 1:ngens(I)])
end

function saturated_ideal(I::MPolyQuoLocalizedIdeal)
  return saturated_ideal(pre_image_ideal(I))
end

function saturated_ideal(I::MPolyQuoLocalizedIdeal{LRT}) where {LRT<:MPolyQuoLocRing{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement}}
  return saturated_ideal(pre_image_ideal(I),strategy=:iterative_saturation,with_generator_transition=false)
end 

function vector_space_dimension(R::MPolyQuoLocRing{<:Field, <:Any,<:Any, <:Any,
                                 <:MPolyComplementOfKPointIdeal})
  I = shifted_ideal(modulus(R))
  o = negdegrevlex(gens(base_ring(R)))
  LI=leading_ideal(standard_basis(I, ordering = o))
  return vector_space_dimension(quo(base_ring(R),ideal(base_ring(R),gens(LI)))[1])
end

@doc raw"""
  _monomial_basis(L::MPolyLocRing{<:Field, <:Any, <:Any, <:Any, <:MPolyComplementOfKPointIdeal}, I::MPolyLocalizedIdeal)

If, say, `A = L/I`, where `L` is a localization of multivariate polynomial ring over a field
`K` at a point `p`, and `I` is an ideal of `L`, return a vector of monomials of `L`
such that the residue classes of these monomials form a basis of `A` as a `K`-vector
space.
!!! note 
    This is an internal method for computing a monomial basis without creating the quotient. 
"""
function _monomial_basis(L::MPolyLocRing{<:Field, <:Any, <:Any, <:Any, <:MPolyComplementOfKPointIdeal}, I::MPolyLocalizedIdeal)
  base_ring(I) == L || error("ideal does not belong to the correct ring")
  G = numerator.(gens(I))
  shift, _ = base_ring_shifts(L)
  G_0 = shift.(G) 
  R = base_ring(L)
  LI = leading_ideal(ideal(R, G_0), ordering = negdeglex(R))
  return L.(monomial_basis(quo(R, LI)[1]))
end

@doc raw"""
    monomial_basis(RQL::MPolyQuoLocRing{<:Field, <:Any, <:Any, <:Any, <:MPolyComplementOfKPointIdeal})

Given a localized ring $RQL$ of type $RQL = (R/I)_P$, where
- ``R`` is a multivariate polynomial ring over a field $K$, say $R = K[x_1,\dots, x_n]$,
- ``I`` is an ideal of $R$ such that $R/I$ is finite-dimensional as a $K$-vector space, and
- ``P`` is of type $P = \langle x_1-a_1,...,x_n-a_n\rangle$, with $a_1,\dots, a_n\in K$,
return a vector of monomials of $R$ such that the images of these monomials under the
composition of the localization map with the projection map form a $K$-basis of $RQL$.

!!! note
    The conditions on $RQL$ are automatically checked.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> f = (y-1)^3-(x)^2;

julia> g = (y-1)^5-(x)^2;

julia> I = ideal(R, [f, g]);  # tangential cusps at  [0, 1] and transversal
                              # intersections at [1, 2], [-1, 2], [1, i], [1, -i]

julia> RQ, _ = quo(R, I);

julia> monomial_basis(RQ)
10-element Vector{QQMPolyRingElem}:
 x*y^2
 y^2
 x^3*y
 x^2*y
 x*y
 y
 x^3
 x^2
 x
 1

julia> P = [0, 1]
2-element Vector{Int64}:
 0
 1

julia> U = complement_of_point_ideal(R, P)
Complement
  of maximal ideal corresponding to rational point with coordinates (0, 1)
  in multivariate polynomial ring in 2 variables over QQ

julia> RQL, _ = localization(RQ, U);

julia> monomial_basis(RQL)
6-element Vector{MPolyLocRingElem{QQField, QQFieldElem, QQMPolyRing, QQMPolyRingElem, MPolyComplementOfKPointIdeal{QQField, QQFieldElem, QQMPolyRing, QQMPolyRingElem}}}:
 x*y^2
 y^2
 x*y
 y
 x
 1

julia> C = plane_curve(f);

julia> D = plane_curve(g);

julia> intersection_multiplicity(C, D, C(P))
6

```

```jldoctest
julia> R, (x,y) = QQ["x","y"];

julia> f = (x^2-y^3)*(y-1);   # 3 singularities, a cusp and two nodes

julia> Tf = tjurina_algebra(f)
Quotient
  of multivariate polynomial ring in 2 variables x, y
    over rational field
  by ideal (x^2*y - x^2 - y^4 + y^3, 2*x*y - 2*x, x^2 - 4*y^3 + 3*y^2)

julia> monomial_basis(Tf)
4-element Vector{QQMPolyRingElem}:
 y^2
 y
 x
 1

julia> TfL0,_ = localization(Tf, complement_of_point_ideal(R, [0,0]));

julia> monomial_basis(TfL0)
2-element Vector{MPolyLocRingElem{QQField, QQFieldElem, QQMPolyRing, QQMPolyRingElem, MPolyComplementOfKPointIdeal{QQField, QQFieldElem, QQMPolyRing, QQMPolyRingElem}}}:
 y
 1

julia> TfL1,_ = localization(Tf, complement_of_point_ideal(R, [1,1]));

julia> monomial_basis(TfL1)
1-element Vector{MPolyLocRingElem{QQField, QQFieldElem, QQMPolyRing, QQMPolyRingElem, MPolyComplementOfKPointIdeal{QQField, QQFieldElem, QQMPolyRing, QQMPolyRingElem}}}:
 1

julia> TfL2,_ = localization(Tf, complement_of_point_ideal(R, [-1,1]));

julia> monomial_basis(TfL2)
1-element Vector{MPolyLocRingElem{QQField, QQFieldElem, QQMPolyRing, QQMPolyRingElem, MPolyComplementOfKPointIdeal{QQField, QQFieldElem, QQMPolyRing, QQMPolyRingElem}}}:
 1

julia> vector_space_dimension(Tf) == vector_space_dimension(TfL0) + vector_space_dimension(TfL1) + vector_space_dimension(TfL2)
true
```
"""
function monomial_basis(A::MPolyQuoLocRing{<:Field, <:Any, <:Any, <:Any, <:MPolyComplementOfKPointIdeal})
  return _monomial_basis(localized_ring(A), modulus(A))
end

function is_finite_dimensional_vector_space(R::MPolyQuoLocRing)
  throw(NotImplementedError(:is_finite_dimensional_vector_space, R))
end

function is_finite_dimensional_vector_space(R::MPAnyNonQuoRing)
  return false
end

### Conversion of ideals in the original ring to localized ideals
function (W::MPolyQuoLocRing{BRT, BRET, RT, RET, MST})(I::MPolyIdeal{RET}) where {BRT, BRET, RT, RET, MST}
  return MPolyQuoLocalizedIdeal(W, W.(gens(I)))
end

### required constructors 
function ideal(
    W::MPolyQuoLocRing, f
  )
  return MPolyQuoLocalizedIdeal(W, [W(f)])
end

function ideal(
    W::MPolyQuoLocRing, gens::Vector
  )
  return MPolyQuoLocalizedIdeal(W, Vector{elem_type(W)}(W.(gens)))
end

function ideal(
    W::MPolyQuoLocRing,
    I::MPolyLocalizedIdeal
  )
  return MPolyQuoLocalizedIdeal(W, W.(gens(I)))
end

function ideal(
    W::MPolyQuoLocRing,
    I::MPolyIdeal
  )
  return MPolyQuoLocalizedIdeal(W, W.(gens(I)))
end

### Printing for all mpoly ideal types
"""
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> S = complement_of_point_ideal(R,[1,2]);

julia> SinvR = localization(R,S)[1]
Localization
  of multivariate polynomial ring in 2 variables x, y
    over rational field
  at complement of maximal ideal of point (1, 2)

julia> I = ideal(SinvR, gens(SinvR))
Ideal generated by
  x
  y

julia> (I,)
(Ideal (x, y),)

```
"""
function Base.show(io::IO, ::MIME"text/plain", I::MPolyAnyIdeal)
  io = pretty(io)
  print(io, "Ideal ")
  if ngens(I) == 0
    print(io, "with ", ItemQuantity(ngens(I), "generator"))
  else
    print(io, "generated by")
    for f in gens(I)
      print(io, "\n", Indent(), f, Dedent())
    end
  end
end

function _get_generators_string_one_line(I::MPolyAnyIdeal, character_limit::Int = 100)
  # Try a full list of generators if it fits $character_limit characters, otherwise
  # print `default`
  default = "with $(ItemQuantity(ngens(I), "generator"))"

  if ngens(I)*3 > character_limit
    # We need at least 3 characters (generator, comma, space) per generator, so
    # we don't need to build the whole string
    return default
  end

  # Generate the full string
  gen_string = "("*join(gens(I), ", ")*")"
  if length(gen_string) <= character_limit
    return gen_string
  end

  return default
end

function Base.show(io::IO, I::MPolyAnyIdeal)
  if is_terse(io)
    print(io, "Ideal")
  else
    print(io, "Ideal ", _get_generators_string_one_line(I))
  end
end

### Further constructors for quotient rings
function quo(
    L::MPolyQuoLocRing,
    I::MPolyQuoLocalizedIdeal
  )
  base_ring(I) === L || error("ideal does not belong to the correct ring")
  R = base_ring(L)
  J = modulus(underlying_quotient(L))
  J = ideal(R, vcat([g for g in gens(J) if !iszero(g)], 
                    [g for g in lifted_numerator.(gens(pre_image_ideal(I))) if !(g in J)]))
  W, _ = quo(localized_ring(L), localized_ring(L)(J))
  return W, hom(L, W, gens(W), check=false)
end

function quo(A::MPolyQuoRing, I::MPolyQuoIdeal)
  base_ring(I) === A || error("ideal does not belong to the correct ring")
  R = base_ring(A)
  Q, _ = quo(R, modulus(A) + ideal(R, lift.(gens(I))))
  return Q, hom(A, Q, Q.(gens(R)), check=false)
end

function divides(a::MPolyQuoLocRingElem, b::MPolyQuoLocRingElem)
  W = parent(a)
  W === parent(b) || error("elements do not belong to the same ring")
  F = FreeMod(W, 1)
  A = matrix(W, 1, 1, [b])
  M, _ = sub(F, A)
  represents_element(a*F[1], M) || return (false, zero(W))
  x = coordinates(a*F[1], M)
  return true, W(x[1])
end

function quotient(I::IdealType, J::IdealType) where {IdealType<:MPolyQuoLocalizedIdeal}
  return base_ring(I)(quotient(pre_image_ideal(I), pre_image_ideal(J)))
end

##############################################################################
# Further functionality for elements of localized quotients of rings
##############################################################################

function derivative(f::MPolyQuoLocRingElem, i::Int)
  num = derivative(lifted_numerator(f), i)*lifted_denominator(f) - derivative(lifted_denominator(f), i)*lifted_numerator(f)
  den = lifted_denominator(f)^2
  g = gcd(num, den)
  return parent(f)(divexact(num, g), divexact(den, g), check=false)
end

function jacobian_matrix(f::MPolyQuoLocRingElem)
  L = parent(f)
  n = nvars(base_ring(L))
  return matrix(L, n, 1, [derivative(f, i) for i=1:n])
end

function jacobian_matrix(g::Vector{<:MPolyQuoLocRingElem})
  L = parent(g[1])
  n = nvars(base_ring(L))
  @assert all(x->parent(x) === L, g)
  return matrix(L, n, length(g), [derivative(x, i) for i=1:n for x = g])
end

@attr Bool function is_prime(I::MPolyQuoLocalizedIdeal)
  return is_prime(saturated_ideal(I))
end

@attr Bool function _is_integral_domain(W::MPolyQuoLocRing)
  return is_prime(modulus(W))
end

@doc raw"""
    is_integral_domain(R::Ring)

Return whether or not `R` is an integral domain.
"""
function _is_integral_domain(R::Ring)
  is_domain_type(typeof(R)) && return true
  error("method not implemented for rings of type $(typeof(R))")
end

### Some auxiliary functions

@attr T function radical(I::T) where {T<:MPolyQuoLocalizedIdeal}
  get_attribute(I, :is_prime, false) && return I
  get_attribute(I, :is_radical, false) && return I
  R = base_ring(I)
  R_simp, iso, iso_inv = simplify(R) # This usually does not cost much
  I_simp = ideal(R_simp, restricted_map(iso).(lifted_numerator.(gens(I))))
  J = pre_image_ideal(I_simp)
  pre_result = ideal(R_simp, [R_simp(g) for g in gens(radical(J)) if !iszero(g)])
  return ideal(R, restricted_map(iso_inv).(lifted_numerator.(gens(pre_result))))
end

@attr Union{Int, NegInf} function dim(I::MPolyQuoLocalizedIdeal)
  return dim(pre_image_ideal(I))
end

#############################################################################
## compatibility functions to allow user to not care about which type of
## MPolyAnyIdeal they are handling
############################################################################# 

@doc raw"""
    primary_decomposition(I::Union{<:MPolyQuoIdeal, <:MPolyQuoLocalizedIdeal, <:MPolyLocalizedIdeal})

Return the primary decomposition of ``I``.
"""
function primary_decomposition(
    I::Union{<:MPolyQuoIdeal, <:MPolyQuoLocalizedIdeal, <:MPolyLocalizedIdeal};
    algorithm::Symbol=:GTZ, cache::Bool=true
  )
  if has_attribute(I, :primary_decomposition)
    return get_attribute(I, :primary_decomposition)::Vector{Tuple{typeof(I), typeof(I)}}
  end
  Q = base_ring(I)
  R = base_ring(Q)
  decomp = primary_decomposition(saturated_ideal(I); algorithm, cache)
  result = [(ideal(Q, Q.(gens(a))), ideal(Q, Q.(gens(b)))) for (a, b) in decomp]
  erase = Int[]
  for i in 1:length(result)
    # if component is trivial, erase it
    is_one(result[i][2]) || continue
    push!(erase,i)
  end
  deleteat!(result, erase)
  
  for (Q,P) in result
    set_attribute!(P, :is_prime=>true)
    set_attribute!(Q, :is_primary=>true)
  end

  cache && set_attribute!(I, :primary_decomposition=>result)
  return result
end

@doc raw"""
    minimal_primes(I::Union{<:MPolyQuoIdeal, <:MPolyQuoLocalizedIdeal, <:MPolyLocalizedIdeal})

Return the minimal associated primes of I.
"""
function minimal_primes(
    I::Union{<:MPolyQuoIdeal, <:MPolyQuoLocalizedIdeal, <:MPolyLocalizedIdeal};
    algorithm::Symbol=:GTZ
  )
  Q = base_ring(I)
  R = base_ring(Q)
  decomp = minimal_primes(saturated_ideal(I); algorithm)
  result = [ideal(Q, Q.(gens(b))) for b in decomp]
  erase = Int[]
  for i in 1:length(result)
    # if a component is trivial, erase it
    is_one(result[i]) || continue
    push!(erase,i)
  end
  deleteat!(result, erase)

  for Ptemp in result
    set_attribute!(Ptemp, :is_prime=>true)
  end

  return result
end

@doc raw"""
    saturation(I::T, J::T) where T <: Union{ MPolyQuoIdeal, MPolyLocalizedIdeal, MPolyQuoLocalizedIdeal}

Return ``I:J^\infty``.
"""
function saturation(I::IdealType, J::IdealType) where {IdealType<:Union{MPolyQuoIdeal, MPolyLocalizedIdeal, MPolyQuoLocalizedIdeal}}
  A = base_ring(I)
  A === base_ring(J) || error("ideals must lie in the same ring")
  R = base_ring(A)

  I_sat = saturated_ideal(I)
  J_sat = saturated_ideal(J)
  K = saturation(I_sat, J_sat)
  return ideal(A, [g for g in A.(gens(K)) if !iszero(g)])
end

@doc raw"""
    saturation_with_index(I::T, J::T) where T <: Union{ MPolyQuoIdeal, MPolyLocalizedIdeal, MPolyQuoLocalizedIdeal}

Return ``I:J^{\infty}`` together with the smallest integer ``m`` such that ``I:J^m = I:J^{\infty}``.
"""
function saturation_with_index(I::T,J::T) where T <: Union{ MPolyQuoIdeal, MPolyLocalizedIdeal, MPolyQuoLocalizedIdeal}
  R = base_ring(I)
  R === base_ring(J) || error("Ideals do not live in the same ring.")

  I_sat = saturated_ideal(I)
  J_sat = saturated_ideal(J)

  I_result,k = saturation_with_index(I_sat,J_sat)
  return (ideal(R,gens(I_result)),k)
end

@doc raw"""
    iterated_quotients(I::T, J::T) where T <: MPolyAnyIdeal
    iterated_quotients(I::T, J::T, b::Int) where T <: MPolyAnyIdeal

Return ``I:J^m`` and maximal ``m`` such that ``J(I:J^m)== (I:J^(m-1))``, if no ``b`` has been specified
Return ``I:J^b`` and ``b``, for the given natural number ``b``.

Internal function for weak and controlled transform.
"""
function iterated_quotients(I::T, J::T, b::Int=0) where T <: MPolyAnyIdeal
  R = base_ring(I)
  R === base_ring(J) || error("Ideals do not live in the same ring.")
  b > -1 || error("negative multiplicity not allowed")

  Itemp = I
  k = 0

  ## iterate ideal quotients b times -- 
  ## or by default (i.e. for b=0) a maximal number of times
  while (b == 0 || k < b)
    Itemp2 = quotient(Itemp, J)
    if !issubset(Itemp, Itemp2 * J)
       b == 0 || error("cannot extract J from I with multiplicity b")
       break
    end
    Itemp = Itemp2
    k = k+1
  end

  return Itemp,k
end

@attr Bool function is_one(I::Union{MPolyQuoIdeal, MPolyLocalizedIdeal, MPolyQuoLocalizedIdeal})
  return is_one(saturated_ideal(I))
end

### Hack for a detour to speed up mapping of elements 
# This is terribly slow in all kinds of quotient rings 
# because of massive checks for `iszero` due to memory 
# management.
function (f::Oscar.MPolyAnyMap{<:MPolyRing, <:MPolyQuoLocRing, <:Nothing})(a::MPolyRingElem)
  g = get_attribute!(f, :lifted_map) do
    S = domain(f)
    W = codomain(f)
    L = localized_ring(W)
    hom(S, L, lift.(f.img_gens), check=false)
  end::Map{typeof(domain(f)), typeof(localized_ring(codomain(f)))}
  b = g(a)::MPolyLocRingElem
  return codomain(f)(numerator(b), denominator(b), check=false)
end

function (f::Oscar.MPolyAnyMap{<:MPolyRing, <:MPolyQuoLocRing, <:MPolyQuoLocalizedRingHom})(a::MPolyRingElem)
  g = get_attribute!(f, :lifted_map) do
    S = domain(f)
    W = codomain(f)
    L = localized_ring(W)
    g = hom(S, L, x -> lift(f.coeff_map(x)), lift.(f.img_gens), check=false)
    set_attribute!(f, :lifted_map, g)
    g
  end::Map{typeof(domain(f)), typeof(localized_ring(codomain(f)))}

  b = g(a)::MPolyLocRingElem
  return codomain(f)(numerator(b), denominator(b), check=false)
end

function vector_space(kk::Field, W::MPolyQuoLocRing;
    ordering::MonomialOrdering=degrevlex(gens(base_ring(W)))
  )
  R = base_ring(W)::MPolyRing
  kk === coefficient_ring(R)::Field || error("change of base field not implemented")
  I = modulus(W)::MPolyLocalizedIdeal
  I_sat = saturated_ideal(modulus(W))::MPolyIdeal
  @assert iszero(dim(I_sat)) "algebra must be zero dimensional"
  o = ordering
  # We set up an algebra isomorphic to the given one. Since we do not know anything about the localization, this is the best we can do.
  A, pr = quo(R, I_sat) 
  f = hom(W, A, gens(A))
  g = hom(A, W, gens(W))
  set_attribute!(f, :inverse, g)
  set_attribute!(g, :inverse, f)
  V, id = vector_space(kk, A)
  return V, MapFromFunc(V, W, v->g(id(v)), a->preimage(id, f(a)))
end

function vector_space(kk::Field, W::MPolyQuoLocRing{<:Field, <:FieldElem, 
                                                    <:MPolyRing, <:MPolyRingElem,
                                                    <:MPolyComplementOfKPointIdeal
                                                   };
    ordering::MonomialOrdering=negdegrevlex(gens(base_ring(W)))
  )
  R = base_ring(W)::MPolyRing
  kk === coefficient_ring(R)::Field || error("change of base field not implemented")
  I = modulus(W)::MPolyLocalizedIdeal
  I_shift = shifted_ideal(I)::MPolyIdeal

  # Collect all monomials which are not in the leading ideal as representatives 
  # of a basis over kk.
  lead_I = leading_ideal(I_shift, ordering=ordering)
  @assert iszero(dim(lead_I)) "quotient must be zero dimensional"
  V_gens = elem_type(R)[]
  done = false
  d = 0
  while !done
    inc = [m for m in monomials_of_degree(R, d) if !(m in lead_I)]
    if iszero(length(inc))
      done = true
      break
    end
    V_gens = vcat(V_gens, [m for m in monomials_of_degree(R, d) if !(m in lead_I)])
    d = d + 1
  end

  n = length(V_gens)
  V = free_module(kk, n)
  L = localized_ring(W)::MPolyLocRing
  shift, backshift = base_ring_shifts(L)

  function im(a::Generic.FreeModuleElem)
    @assert parent(a) === V
    b = R(0)
    for k=1:length(V_gens)
      c = a[k]
      if !iszero(c)
        b += c*V_gens[k]
      end
    end
    return W(backshift(b))
  end

  # The inverse function. We use the fact that for a chosen monomial ordering 
  # the monomials which are not in the leading ideal, form a basis for the 
  # quotient; see Greuel/Pfister "A singular introduction to Commutative Algebra".
  function prim(f::MPolyQuoLocRingElem)
    @assert parent(f) === W
    isone(denominator(f)) || error("handling of non-trivial denominators not implemented")
    b = lifted_numerator(f)
    b = shift(b)
    # TODO: The normal form command seems to do something unexpected here.
    # Please investigate!
    #b = normal_form(b, I_shift, ordering=ordering)
    result = zero(V)
    # The following is an ugly hack, because normal_form is currently broken 
    # for local orderings.
    while !iszero(b)
      m = leading_monomial(b, ordering=ordering)
      c = leading_coefficient(b, ordering=ordering)
      t = normal_form(c*m, I_shift, ordering=ordering)
      if t == c*m 
        j = findfirst(==(m), V_gens)
        result = result + c * V[j]
        b = b - c * m
      else
        b = b - c*m + t
      end
    end
    return result
  end
  return V, MapFromFunc(V, W, im, prim)
end


# disambiguate some conversions
# ... but this needs to be after some type declarations... so here it is
function (W::MPolyDecRing)(f::MPolyQuoRingElem)
  return W(forget_decoration(W)(f))
end

function (W::MPolyDecRing)(f::MPolyQuoLocRingElem)
  return W(forget_decoration(W)(f))
end

function (W::MPolyDecRing)(f::MPolyLocRingElem)
  return W(forget_decoration(W)(f))
end

@attr Tuple{<:Map, <:Map} function _as_affine_algebra_with_many_variables(
    L::MPolyLocRing{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement}
  )
  inverse_name=:_0
  R = base_ring(L)
  f = denominators(inverted_set(L))
  f = sort(f; by=total_degree, rev=true)
  r = length(f)
  A, phi, t = _add_variables_first(R, [Symbol(inverse_name,k) for k in 1:r])
  theta = t[1:r]
  I = ideal(A, [one(A)-theta[k]*phi(f[k]) for k in 1:r])
  ordering = degrevlex(gens(A)[r+1:end])
  if r > 0 
    ordering = deglex(theta)*ordering
  end
  Q = MPolyQuoRing(A, I, ordering)
  function my_fun(g)
    a = Q(phi(lifted_numerator(g)))
    isone(lifted_denominator(g)) && return a
    b = Q(phi(lifted_denominator(g)))
    #success, c = divides(a, b)
    success, c = _divides_hack(a, b)
    success || error("element can not be mapped")
    return c
  end
  id = MapFromFunc(L, Q, my_fun)
  #id = hom(L, Q, gens(A)[r+1:end], check=false)
  id_inv = hom(Q, L, vcat([L(one(R), b, check=false) for b in f], gens(L)), check=false)
  set_attribute!(id, :inverse, id_inv)
  set_attribute!(id_inv, :inverse, id)
  return id, id_inv
end

@attr Tuple{<:Map, <:Map} function _as_affine_algebra_with_many_variables(
    L::MPolyQuoLocRing{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement}
  )
  inverse_name=:_0
  R = base_ring(L)
  f = denominators(inverted_set(L))
  f = sort(f; by=total_degree, rev=true)
  r = length(f)
  A, phi, t = _add_variables_first(R, [Symbol(inverse_name,k) for k in 1:r])
  theta = t[1:r]
  I = ideal(A, [phi(g) for g in gens(modulus(underlying_quotient(L)))]) + ideal(A, [one(A)-theta[k]*phi(f[k]) for k in 1:r])
  ordering = degrevlex(gens(A)[r+1:end])
  if r > 0 
    ordering = deglex(theta)*ordering
  end
  Q = MPolyQuoRing(A, I, ordering)
  function my_fun(g)
    a = Q(phi(lifted_numerator(g)))
    isone(lifted_denominator(g)) && return a
    b = Q(phi(lifted_denominator(g)))
    #success, c = divides(a, b)
    success, c = _divides_hack(a, b)
    success || error("element can not be mapped")
    return c
  end
  id = MapFromFunc(L, Q, my_fun)
  #id = hom(L, Q, gens(A)[r+1:end], check=false)
  id_inv = hom(Q, L, vcat([L(one(R), b, check=false) for b in f], gens(L)), check=false)
  set_attribute!(id, :inverse, id_inv)
  set_attribute!(id_inv, :inverse, id)
  return id, id_inv
end

function inverse(phi::MPolyQuoLocalizedRingHom)
  has_attribute(phi, :inverse) && return get_attribute(phi, :inverse)
  error("computation of inverse not implemented")
end

@doc raw"""
    minimal_generating_set(I::MPolyQuoLocalizedIdeal{<:MPolyQuoLocRing{<:Field, <:Any, <:Any, <:Any, <:MPolyComplementOfKPointIdeal}, <:Any,<:Any})

Given an ideal $I$ in a localized ring $RQL$ of type $RQL = (R/J)_P$, where
- ``R`` is a multivariate polynomial ring over a field $K$, say $R = K[x_1,\dots, x_n]$, and
- ``P`` is of type $P = \langle x_1-a_1,...,x_n-a_n\rangle$, with $a_1,\dots, a_n\in K$,
return a vector containing a minimal set of generators of `I`.

If `I` is the zero ideal, then an empty list is returned.
"""
@attr Vector{<:MPolyQuoLocRingElem} function minimal_generating_set(
    I::MPolyQuoLocalizedIdeal{<:MPolyQuoLocRing{<:Field, <:Any,
                                          <:Any, <:Any,
                                          <:MPolyComplementOfKPointIdeal},
                              <:Any,<:Any}
  )
  Q = base_ring(I)
  !is_zero(I) || return typeof(zero(Q))[]
  L = localized_ring(Q)

  ## list of generators in the localized ring, append modulus
  Jlist = lift.(gens(I))
  nJlist = length(Jlist)
  append!(Jlist,gens(modulus(Q)))

  ## move to origin
  shift, back_shift = base_ring_shifts(L)
  I_shift = shifted_ideal(ideal(L,Jlist))
  R = base_ring(I_shift)
  oL = negdegrevlex(R)                    # the default local ordering


  ## determine the relations
  gensSord_shift = singular_generators(I_shift, oL)
  syz_mod = Singular.syz(gensSord_shift)

  ## prepare Nakayama-check for minimal generating system
  F = free_module(R, length(Jlist))
  Imax = ideal(R,gens(R))
  M = sub(F,[F(syz_mod[i]) for i=1:Singular.ngens(syz_mod)])[1] + (Imax*F)[1]
  oF =  negdegrevlex(R)*invlex(F)
  res_vec = typeof(gen(I,1))[]

  ## select by Nakayama
  for i in 1:nJlist
    if !(gen(F,i) in leading_module(M,oF))
      append!(res_vec,[Q.(gen(I,i))])
    end
  end

  return res_vec
end

@doc raw"""
    small_generating_set(I::MPolyQuoLocalizedIdeal)

Given an ideal `I` in a quotient of a localization of a multivariate
polynomial ring over a field, return an array containing a set of
generators of `I`, which is usually smaller than the original one.

If `I` is the zero ideal an empty list is returned.

If the localization is at a point, a minimal set of generators is returned.
"""
@attr Vector{elem_type(base_ring(I))} function small_generating_set(
      I::MPolyQuoLocalizedIdeal{<:MPolyQuoLocRing{<:Field, <:FieldElem,
                                          <:MPolyRing, <:MPolyRingElem,
                                          <:MPolyComplementOfKPointIdeal},
                              <:Any,<:Any};
      algorithm::Symbol=:simple
  )
  Q = base_ring(I)
  L = localized_ring(Q)
  J = pre_image_ideal(I)
  return unique!(filter(!iszero, Q.(small_generating_set(J; algorithm))))
end

@attr Vector{elem_type(base_ring(I))} function small_generating_set(
    I::MPolyQuoLocalizedIdeal{<:MPolyQuoLocRing{<:Field, <:FieldElem,
                                          <:MPolyRing, <:MPolyRingElem,
                                          <:MPolyPowersOfElement}
                          };
      algorithm::Symbol=:simple
  )
  Q = base_ring(I)
  L = localized_ring(Q)

  J = pre_image_ideal(I)
  return unique!(filter(!iszero, Q.(small_generating_set(J; algorithm))))
end

@attr Int function dim(R::MPolyLocRing)
  error("Not implemented")
end

@attr Union{Int, NegInf} function dim(R::MPolyQuoLocRing{<:Any, <:Any, <:MPolyRing, <:MPolyRingElem, <:Union{MPolyComplementOfPrimeIdeal, MPolyComplementOfKPointIdeal}})
  P = prime_ideal(inverted_set(R))
  I = saturated_ideal(modulus(R))
  dim(I) === -inf && return -inf
  return dim(I) - dim(P)
end

@attr Union{Int, NegInf} function dim(R::MPolyQuoLocRing{<:Any, <:Any, <:MPolyRing, <:MPolyRingElem, <:MPolyPowersOfElement})
  return dim(saturated_ideal(modulus(R)))
end

@attr Union{Int, NegInf} function dim(R::MPolyLocRing{<:Any,<:Any,<:MPolyRing,<:MPolyRingElem, <:MPolyPowersOfElement})
  # zariski open subset of A^n
  return dim(base_ring(R))
end

@attr Union{Int, NegInf} function dim(R::MPolyLocRing{<:Any,<:Any,<:MPolyRing,<:MPolyRingElem, <:MPolyComplementOfPrimeIdeal})
  P = prime_ideal(inverted_set(R))
  return codim(P)
end


@attr Union{Int, NegInf} function dim(R::MPolyLocRing{<:Field,<:Any,<:MPolyRing,<:MPolyRingElem, <:MPolyComplementOfKPointIdeal})
  # localization of a polynomial ring over a field at a maximal ideal does not change the dimension
  # because all maximal ideals have the same dimension in this case. 
  return dim(base_ring(R))
end

########################################################################
# Localizations of graded rings                                        #
########################################################################
function is_graded(L::MPolyQuoLocRing{<:Ring, <:RingElem, <:MPolyDecRing})
  return true
end

function grading_group(L::MPolyQuoLocRing{<:Ring, <:RingElem, <:MPolyDecRing})
  return grading_group(base_ring(L))
end

function degree(a::MPolyQuoLocRingElem{<:Ring, <:RingElem, <:MPolyDecRing}; check::Bool=true)
  return degree(numerator(a); check) - degree(denominator(a); check)
end

function degree(::Type{Int}, a::MPolyQuoLocRingElem{<:Ring, <:RingElem, <:MPolyDecRing}; check::Bool=true)
  return degree(Int, numerator(a); check) - degree(Int, denominator(a); check)
end

function degree(::Type{Vector{Int}}, a::MPolyQuoLocRingElem{<:Ring, <:RingElem, <:MPolyDecRing}; check::Bool=true)
  return degree(Vector{Int}, numerator(a); check) - degree(Vector{Int}, denominator(a); check)
end

function _degree_fast(a::MPolyQuoLocRingElem{<:Ring, <:RingElem, <:MPolyDecRing})
  return _degree_fast(numerator(a)) - _degree_fast(denominator(a))
end

function is_homogeneous(a::MPolyQuoLocRingElem{<:Ring, <:RingElem, <:MPolyDecRing})
  return is_homogeneous(numerator(a)) && is_homogeneous(denominator(a))
end

########################################################################
# Inverses of homomorphisms                                            #
########################################################################

function inverse(
    f::Map{<:MPolyAnyRing, <:MPolyAnyRing}
  )
  has_attribute(f, :inverse) && return get_attribute(f, :inverse)::Map
  R = domain(f)
  S = codomain(f)
  RR, iso_R, iso_R_inv = _as_localized_quotient(R)
  SS, iso_S, iso_S_inv = _as_localized_quotient(S)

  ff = hom(RR, SS, iso_S.(f.(gens(R))))

  inv = inverse(ff)
  result = compose(iso_S, compose(inv, iso_R_inv))
  set_attribute!(f, :inverse=>result)
  return result
end

# It seems tedious to implement special methods for all constellations of 
# rings that can appear. But all rings considered can be realized 
# as a localized quotient and the functionality exists for those. 
# Since this does not produce significant overhead, we reroute everything 
# to those. 
function is_isomorphism(
    f::Map{<:MPolyAnyRing, <:MPolyAnyRing}
  )
  has_attribute(f, :inverse) && return true
  R = domain(f)
  S = codomain(f)
  RR, iso_R, iso_R_inv = _as_localized_quotient(R)
  SS, iso_S, iso_S_inv = _as_localized_quotient(S)

  ff = hom(RR, SS, iso_S.(f.(gens(R))))

  result = is_isomorphism(ff)
  if result
    # Supposedly the inverse has been computed and cached by the call to 
    # is_isomorphism
    inv = compose(iso_S, compose(inverse(ff), iso_R_inv))
    set_attribute!(f, :inverse=>inv)
  end
  return result
end

function _as_localized_quotient(R::MPolyRing)
  result = MPolyQuoLocRing(R, ideal(R, elem_type(R)[]), powers_of_element(one(R)))
  iso = hom(R, result, gens(result), check=false)
  iso_inv = hom(result, R, gens(R), check=false)
  return result, iso, iso_inv
end

function _as_localized_quotient(A::MPolyQuoRing)
  R = base_ring(A)
  result = MPolyQuoLocRing(R, modulus(A), powers_of_element(one(R)))
  iso = hom(A, result, gens(result), check=false)
  iso_inv = hom(result, A, gens(R), check=false)
  return result, iso, iso_inv
end

function _as_localized_quotient(L::MPolyLocRing)
  R = base_ring(L)
  result = MPolyQuoLocRing(R, ideal(R, elem_type(R)[]), inverted_set(L))
  iso = hom(L, result, gens(result), check=false)
  iso_inv = hom(result, L, gens(R), check=false)
  return result, iso, iso_inv
end


function _as_localized_quotient(W::MPolyQuoLocRing)
  return W, identity_map(W), identity_map(W)
end

# Problems arise with comparison for non-trivial coefficient maps 
# in all constellations. This calls for a streamlined set of commands 
# to check whether coefficient maps exist and if so, what they are, 
# so that they can be compared. 
function Base.:(==)(
    f::Map{<:MPolyAnyRing, <:MPolyAnyRing},
    g::Map{<:MPolyAnyRing, <:MPolyAnyRing}
   )
  domain(f) === domain(g) || return false
  codomain(f) === codomain(g) || return false

  if _has_coefficient_map(f)
    if _has_coefficient_map(g) 
      return is_equal_as_morphism(coefficient_map(f), coefficient_map(g))
    else
      # In this case, type casting needs to have the same result 
      # as the coefficient map on the other end
      return is_equal_as_morphism(coefficient_ring(codomain(g)), coefficient_map(f))
    end
  elseif _has_coefficient_map(g)
    if _has_coefficient_map(f) 
      return is_equal_as_morphism(coefficient_map(f), coefficient_map(g))
    else
      return is_equal_as_morphism(coefficient_ring(codomain(f)), coefficient_map(g))
    end
  end

  return all(f(x) == g(x) for x in gens(domain(f)))
end

# Accompanying implementation of hash
function Base.hash(f::Map{<:MPolyAnyRing, <:MPolyAnyRing}, h::UInt)
  h = hash(domain(f), h)
  h = hash(codomain(f), h)
  # TODO: add in coefficient_map if available
  h = hash(f.(gens(domain(f))), h)
  return h
end

# This rerouting is necessary as internally the function 
# `is_equal_as_morphism` is called in several places 
# and we wish to have the above implementation for such 
# signatures. 
function is_equal_as_morphism(
    f::Map{<:MPolyAnyRing, <:MPolyAnyRing},
    g::Map{<:MPolyAnyRing, <:MPolyAnyRing}
   )
  return f == g
end

coefficient_map(f::MPolyLocalizedRingHom) = coefficient_map(restricted_map(f))
coefficient_map(f::MPolyQuoLocalizedRingHom) = coefficient_map(restricted_map(f))

_has_coefficient_map(::Type{T}) where {U, T<:MPolyLocalizedRingHom{<:Any, <:Any, U}} = _has_coefficient_map(U)
_has_coefficient_map(f::MPolyLocalizedRingHom) = _has_coefficient_map(typeof(f))

_has_coefficient_map(::Type{T}) where {U, T<:MPolyQuoLocalizedRingHom{<:Any, <:Any, U}} = _has_coefficient_map(U)
_has_coefficient_map(f::MPolyQuoLocalizedRingHom) = _has_coefficient_map(typeof(f))

_has_coefficient_map(::Type{T}) where {U, T<:MPolyAnyMap{<:Any, <:Any, U}} = (U !== Nothing)
_has_coefficient_map(f::MPolyAnyMap) = _has_coefficient_map(typeof(f))

# for the default, we have to assume that there is no coefficient map.
_has_coefficient_map(::Type{T}) where {T<:Map} = false
_has_coefficient_map(f::Map) = _has_coefficient_map(typeof(f))


# Composition of maps is complicated, if there are coefficient maps, too.
# By convention, the coefficient maps are allowed to take values in the 
# coefficient ring of the codomain, or the codomain itself. Depending 
# on the case, we have to select the correct way to compose maps.
function compose(f::MPolyAnyMap, 
    g::Union{<:MPolyQuoLocalizedRingHom, <:MPolyLocalizedRingHom}
  )
  @assert codomain(f) === domain(g)
  if _has_coefficient_map(f)
    b = coefficient_map(f)(one(coefficient_ring(domain(f))))
    if _has_coefficient_map(g)
      if parent(b) === coefficient_ring(domain(g))
        c = coefficient_map(g)(b)
        return hom(domain(f), codomain(g), 
                   _compose(coefficient_map(f), coefficient_map(g), coefficient_ring(domain(f)), parent(c)),
                   g.(f.(gens(domain(f)))); check=false)
      elseif parent(b) === domain(g)
        return hom(domain(f), codomain(g), 
                   _compose(coefficient_map(f), g, coefficient_ring(domain(f)), codomain(g)), 
                   g.(f.(gens(domain(f)))); check=false)
      else 
        # assume that things can be coerced 
        c = g(domain(g)(b))
        return hom(domain(f), codomain(g), 
                   MapFromFunc(coefficient_ring(domain(f)), parent(c), x-> g(domain(g)(coefficient_map(f)(x)))),
                   g.(f.(gens(domain(f)))); check=false)
      end
    else
      if parent(b) === coefficient_ring(domain(g))
        return hom(domain(f), codomain(g), coefficient_map(f), 
                   g.(f.(gens(domain(f)))); check=false)
      elseif parent(b) === domain(g)
        return hom(domain(f), codomain(g), 
                   _compose(coefficient_map(f), g, coefficient_ring(domain(f)), codomain(g)),
                   g.(f.(gens(domain(f)))); check=false)
      else 
        # assume that things can be coerced 
        return hom(domain(f), codomain(g), 
                   MapFromFunc(coefficient_ring(domain(f)), codomain(g), x-> g(domain(g)(coefficient_map(f)(x)))),
                   g.(f.(gens(domain(f)))); check=false)
      end
    end
  elseif _has_coefficient_map(g)
    return hom(domain(f), codomain(g), coefficient_map(g),
               g.(f.(gens(domain(f)))); check=false)
  else
    return hom(domain(f), codomain(g), 
               g.(f.(gens(domain(f)))); check=false)
  end
end

_compose(f::Map, g::Map, dom::Any, cod::Any) = compose(f, g)
_compose(f::Any, g::Map, dom::Any, cod::Any) = MapFromFunc(dom, cod, x->g(f(x)))
_compose(f::Map, g::Any, dom::Any, cod::Any) = MapFromFunc(dom, cod, x->g(f(x)))
_compose(f::Any, g::Any, dom::Any, cod::Any) = MapFromFunc(dom, cod, x->g(f(x)))

morphism_type(::Type{DT}, ::Type{CT}) where {DT<:MPolyLocRing, CT<:Ring} = MPolyLocalizedRingHom{DT, CT, morphism_type(base_ring_type(DT), CT)}

base_ring_type(::Type{T}) where {BRT, T<:MPolyLocRing{<:Any, <:Any, BRT}} = BRT

base_ring_type(::Type{T}) where {BRT, T<:MPolyQuoLocRing{<:Any, <:Any, BRT}} = BRT

function dim(R::MPolyQuoLocRing)
  return dim(modulus(R))
end

### extra methods for speedup of mappings
# See `src/Rings/MPolyMap/MPolyRing.jl` for the original implementation and 
# the rationale. These methods make speedup available for quotient rings 
# and localizations of polynomial rings.
function _allunique(lst::Vector{T}) where {T<:MPolyLocRingElem}
  return _allunique(numerator.(lst))
end

function _allunique(lst::Vector{T}) where {T<:MPolyQuoLocRingElem}
  return _allunique(lifted_numerator.(lst))
end

_is_gen(x::MPolyLocRingElem) = is_one(denominator(x)) && _is_gen(numerator(x))
_is_gen(x::MPolyQuoLocRingElem) = is_one(lifted_denominator(x)) && _is_gen(lifted_numerator(x))

function _evaluate_plain(
    F::MPolyAnyMap{<:MPolyRing, CT}, u
  ) where {CT <: Union{<:MPolyLocRing, <:MPolyQuoLocRing}}
  W = codomain(F)::MPolyLocRing
  S = base_ring(W)
  fl, var_ind = _maps_variables_to_variables(F)
  fl && return W(_build_poly(u, var_ind, S))
  return evaluate(u, F.img_gens)
end

function _evaluate_general(
    F::MPolyAnyMap{<:MPolyRing, CT}, u
  ) where {CT <: Union{<:MPolyQuoRing, <:MPolyLocRing, <:MPolyQuoLocRing}}
  S = temp_ring(F)
  fl, var_ind = _maps_variables_to_variables(F)
  if S !== nothing
    if !fl || coefficient_ring(S) !== codomain(F)
      return evaluate(map_coefficients(coefficient_map(F), u,
                                       parent = S), F.img_gens)
    else
      tmp_poly = map_coefficients(coefficient_map(F), u, parent = S)
      return _evaluate_with_build_ctx(tmp_poly, var_ind, codomain(F))
    end
  else
    if !fl
      return evaluate(map_coefficients(coefficient_map(F), u), F.img_gens)
    else
      # For the case where we can recycle the method above, do so.
      tmp_poly = map_coefficients(coefficient_map(F), u)
      coefficient_ring(parent(tmp_poly)) === codomain(F) && return _evaluate_with_build_ctx(
                                                                                            tmp_poly,
                                                                                            var_ind,
                                                                                            codomain(F)
                                                                                           )
      # Otherwise default to the standard evaluation for the time being.
      return evaluate(tmp_poly, F.img_gens)
    end
  end
end

function _evaluate_with_build_ctx(
    p::MPolyRingElem, ind::Vector{Int}, 
    Q::Union{<:MPolyQuoRing, <:MPolyLocRing, <:MPolyQuoLocRing}
  )
  @assert Q === coefficient_ring(parent(p))
  cod_ring = base_ring(Q)
  r = ngens(cod_ring)
  kk = coefficient_ring(cod_ring)
  ctx = MPolyBuildCtx(cod_ring)
  for (q, e) in zip(AbstractAlgebra.coefficients(p), AbstractAlgebra.exponent_vectors(p))
    ee = [0 for _ in 1:r]
    for (i, k) in enumerate(e)
      ee[ind[i]] = k
    end
    for (c, d) in zip(_coefficients(q), _exponents(q))
      push_term!(ctx, kk(c), ee+d)
    end
  end
  return Q(finish(ctx))
end

# The following methods are only safe, because they are called 
# exclusively in a setup where variables map to variables 
# and no denominators are introduced. 
_coefficients(x::MPolyRingElem) = AbstractAlgebra.coefficients(x)
_coefficients(x::MPolyQuoRingElem) = AbstractAlgebra.coefficients(lift(x))
_coefficients(x::MPolyLocRingElem) = AbstractAlgebra.coefficients(numerator(x))
_coefficients(x::MPolyQuoLocRingElem) = AbstractAlgebra.coefficients(lifted_numerator(x))

_exponents(x::MPolyRingElem) = AbstractAlgebra.exponent_vectors(x)
_exponents(x::MPolyQuoRingElem) = AbstractAlgebra.exponent_vectors(lift(x))
_exponents(x::MPolyLocRingElem) = AbstractAlgebra.exponent_vectors(numerator(x))
_exponents(x::MPolyQuoLocRingElem) = AbstractAlgebra.exponent_vectors(lifted_numerator(x))

# overwriting the comparison method to avoid computing saturations and groebner bases.
_cmp_reps(a::MPolyLocRingElem) = y->(fraction(y) == fraction(a))
_cmp_reps(a::MPolyQuoLocRingElem) = y->(fraction(y) == fraction(a))
_cmp_reps(a::MPolyQuoRingElem) = y->(y.f == a.f)

