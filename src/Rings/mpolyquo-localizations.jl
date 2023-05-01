import AbstractAlgebra: Ring, RingElem, Generic.Frac
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
    base_ring(I) == R || error("Ideal does not belong to the ring")
    base_ring(Q) == R || error("The quotient ring does not come from the given ring")
    # The following line throws obscure error messages that might yield a bug for MPolyIdeals.
    # So it's commented out for now.
    #modulus(Q) == I || error("the modulus of the quotient ring does not coincide with the ideal")
    S == inverted_set(W) || error("the multiplicative set does not coincide with the inverted set of the localized ring")
    base_ring(W) == R || error("the localization does not come from the given ring")
    ambient_ring(S) == R || error("Multiplicative set does not belong to the ring")
    k = coefficient_ring(R)
    L = new{typeof(k), elem_type(k), typeof(R), RingElemType, MultSetType}(R, I, S, Q, W)
    return L
  end
end

### for convenience of later use
MPAnyQuoRing = Union{MPolyQuoLocRing, 
                MPolyQuoRing
               }

MPAnyNonQuoRing = Union{MPolyRing, MPolyLocRing
                  }

MPolyAnyRing = Union{MPolyRing, MPolyQuoRing,
                MPolyLocRing,MPolyQuoLocRing
               }


### type getters 
coefficient_ring_type(::Type{MPolyQuoLocRing{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = BRT
coefficient_ring_type(L::MPolyQuoLocRing{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST} = coefficient_ring_type(typeof(L))
coefficient_ring_elem_type(::Type{MPolyQuoLocRing{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = BRET
coefficient_ring_elem_type(L::MPolyQuoLocRing{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST} = coefficient_ring_elem_type(typeof(L))

base_ring_type(::Type{MPolyQuoLocRing{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = RT
base_ring_type(L::MPolyQuoLocRing{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST} = base_ring_type(typeof(L))
base_ring_elem_type(::Type{MPolyQuoLocRing{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = RET
base_ring_elem_type(L::MPolyQuoLocRing{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST} = base_ring_elem_type(typeof(L))

mult_set_type(::Type{MPolyQuoLocRing{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = MST
mult_set_type(L::MPolyQuoLocRing{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST} = mult_set_type(typeof(L))

localized_ring_type(::Type{MPolyQuoLocRing{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = MPolyLocRing{BRT, BRET, RT, RET, MST}
localized_ring_type(L::MPolyQuoLocRing) = localized_ring_type(typeof(L))

ideal_type(::Type{MPolyQuoLocalizedRingType}) where {MPolyQuoLocalizedRingType<:MPolyQuoLocRing} = MPolyQuoLocalizedIdeal{MPolyQuoLocalizedRingType, elem_type(MPolyQuoLocalizedRingType), ideal_type(localized_ring_type(MPolyQuoLocalizedRingType))}
ideal_type(W::MPolyQuoLocRing) = ideal_type(typeof(W))





### required getter functions 
base_ring(L::MPolyQuoLocRing) = L.R
inverted_set(L::MPolyQuoLocRing) = L.S

### additional getter functions
@doc raw"""
    modulus(L::MPolyQuoLocRing)

Given ``L = (𝕜[x₁,…,xₙ]/I)[S⁻¹]``, return ``IS⁻¹``.
"""
function modulus(L::MPolyQuoLocRing) 
  if !has_attribute(L, :modulus)
    set_attribute!(L, :modulus, localized_ring(L)(L.I))
  end
  return get_attribute(L, :modulus)::ideal_type(localized_ring_type(L))
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
function Base.show(io::IO, L::MPolyQuoLocRing)
  print(io, "Localization of $(underlying_quotient(L)) at the multiplicative set $(inverted_set(L))")
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
  return L, hom(W, L, hom(R, L, gens(L)))
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
  return W, hom(L, W, gens(W))
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
  return P, hom(L, P, gens(P))
end

@doc raw"""
    localization(RQ::MPolyQuoRing, U::AbsMPolyMultSet)

Given a quotient `RQ` of a multivariate polynomial ring `R` with projection map
`p : R -> RQ`, say, and given a multiplicatively closed subset `U` of `R`, return the 
localization of `RQ` at `p(U)`, together with the localization map.

# Examples

```jldoctest
julia> T, t = polynomial_ring(QQ, "t");

julia> K, a =  number_field(2*t^2-1, "a");

julia> R, (x, y) = polynomial_ring(K, ["x", "y"]);

julia> I = ideal(R, [2*x^2-y^3, 2*x^2-y^5])
ideal(2*x^2 - y^3, 2*x^2 - y^5)

julia> P = ideal(R, [y-1, x-a])
ideal(y - 1, x - a)

julia> U = complement_of_prime_ideal(P)
complement of ideal(y - 1, x - a)

julia> RQ, _ = quo(R, I);

julia> RQL, iota = localization(RQ, U);

julia> RQL
Localization of Quotient of Multivariate Polynomial Ring in x, y over Number field over Rational Field with defining polynomial 2*t^2 - 1 by ideal(2*x^2 - y^3, 2*x^2 - y^5) at the multiplicative set complement of ideal(y - 1, x - a)

julia> iota
Map from
Quotient of Multivariate Polynomial Ring in x, y over Number field over Rational Field with defining polynomial 2*t^2 - 1 by ideal(2*x^2 - y^3, 2*x^2 - y^5) to Localization of Quotient of Multivariate Polynomial Ring in x, y over Number field over Rational Field with defining polynomial 2*t^2 - 1 by ideal(2*x^2 - y^3, 2*x^2 - y^5) at the multiplicative set complement of ideal(y - 1, x - a) defined by a julia-function
```
""" localization(A::MPolyQuoRing, U::AbsMPolyMultSet)

###localization is an Abstract Algebra alias for Localization

function Localization(Q::MPolyQuoRing{RET}, S::MultSetType) where {RET <: RingElem, MultSetType <: AbsMultSet}
  L = MPolyQuoLocRing(base_ring(Q), modulus(Q), S, Q, Localization(S)[1])
  return L, MapFromFunc((x->L(lift(x))), Q, L)
end

function Localization(
    L::MPolyQuoLocRing{BRT, BRET, RT, RET, MST}, 
    S::AbsMPolyMultSet{BRT, BRET, RT, RET}
  ) where {BRT, BRET, RT, RET, MST}
  ambient_ring(S) == base_ring(L) || error("multiplicative set does not belong to the correct ring")
  issubset(S, inverted_set(L)) && return L, MapFromFunc(x->x, L, L)
  U = inverted_set(L)*S
  W = MPolyQuoLocRing(base_ring(L), modulus(underlying_quotient(L)), U, underlying_quotient(L), Localization(U)[1])
  return W, MapFromFunc((x->W(lifted_numerator(x), lifted_denominator(x), check=false)), L, W)
end

function MPolyQuoLocRing(R::RT, I::Ideal{RET}, T::MultSetType) where {RT<:MPolyRing, RET<:MPolyRingElem, MultSetType<:AbsMultSet} 
  return MPolyQuoLocRing(R, I, T, quo(R, I)[1], Localization(T)[1])
end

function MPolyQuoLocRing(R::RT) where {RT<:MPolyRing} 
  I = ideal(R, zero(R))
  Q, _ = quo(R, I)
  U = units_of(R)
  W, _ = Localization(U)
  return MPolyQuoLocRing(R, I, U, Q, W)
end

function MPolyQuoLocRing(Q::RT) where {RT<:MPolyQuoRing}
  R = base_ring(Q)
  I = modulus(Q)
  U = units_of(R)
  W, _ = Localization(U)
  return MPolyQuoLocRing(R, I, U, Q, W)
end

function MPolyQuoLocRing(W::MPolyLocRing)
  R = base_ring(W)
  I = ideal(R, zero(R))
  Q, _ = quo(R, I)
  U = inverted_set(W)
  return MPolyQuoLocRing(R, I, U, Q, W)
end

function Base.in(f::AbstractAlgebra.Generic.Frac{RET}, L::MPolyQuoLocRing{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST}
  R = base_ring(L)
  R == parent(numerator(f)) || error("element does not belong to the correct ring")
  denominator(f) in inverted_set(L) && return true
  return numerator(f) in ideal(L, denominator(f))
end


### generation of random elements 
function rand(W::MPolyQuoLocRing, v1::UnitRange{Int}, v2::UnitRange{Int}, v3::UnitRange{Int})
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
    parent(a) == parent(b) == R || error("elements do not belong to the correct ring")
    check && (b in S || error("denominator is not admissible"))
    return new{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}(L, a, b, is_reduced)
  end
end

### type getters
coefficient_ring_type(::Type{MPolyQuoLocRingElem{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = BRT
coefficient_ring_type(f::MPolyQuoLocRingElem{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST} = base_ring_type(typeof(f))
coefficient_ring_elem_type(::Type{MPolyQuoLocRingElem{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = BRET
coefficient_ring_elem_type(f::MPolyQuoLocRingElem{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST} = base_ring_type(typeof(f))

base_ring_type(::Type{MPolyQuoLocRingElem{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = RT
base_ring_type(f::MPolyQuoLocRingElem{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST} = base_ring_type(typeof(f))
base_ring_elem_type(::Type{MPolyQuoLocRingElem{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = RET
base_ring_elem_type(f::MPolyQuoLocRingElem{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST} = base_ring_type(typeof(f))

mult_set_type(::Type{MPolyQuoLocRingElem{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = MST
mult_set_type(f::MPolyQuoLocRingElem{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST} = base_ring_type(typeof(f))

### required getter functions 
parent(a::MPolyQuoLocRingElem) = a.L
numerator(a::MPolyQuoLocRingElem) = underlying_quotient(parent(a))(a.numerator) 
denominator(a::MPolyQuoLocRingElem) = underlying_quotient(parent(a))(a.denominator) 

### additional getter functions
underlying_quotient(a::MPolyQuoLocRingElem) = underlying_quotient(parent(a))
localized_ring(a::MPolyQuoLocRingElem) = localized_ring(parent(a))
base_ring(a::MPolyQuoLocRingElem) = base_ring(parent(a))
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
(L::MPolyQuoLocRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType})(f::RingElemType, is_reduced::Bool=false) where {BaseRingType, BaseRingElemType, RingType, RingElemType<:RingElem, MultSetType} = MPolyQuoLocRingElem(L, f, one(f), check=false, is_reduced=is_reduced)

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
  check || return MPolyQuoLocRingElem(L, a, b, check=false, is_reduced=is_reduced)
  b in inverted_set(L) || return convert(L, a//b)
  return MPolyQuoLocRingElem(L, a, b, check=false, is_reduced=is_reduced)
end

### additional conversions
function (L::MPolyQuoLocRing{BRT, BRET, RT, RET, MST})(f::Frac{RET}; check::Bool=true, is_reduced::Bool=false) where {BRT, BRET, RT, RET, MST}
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
  base_ring(parent(f)) == base_ring(L) || error("the given element does not belong to the correct ring") 
  check && (parent(f) == underlying_quotient(L) || all(x->(iszero(L(x))), gens(modulus(parent(f)))) || error("coercion is not well defined"))
  return L(lift(f))
end

### additional functionality
@doc raw"""
    lift(f::MPolyQuoLocRingElem)

Given ``f = A//B ∈ (𝕜[x₁,…,xₙ]/I)[S⁻¹]``, return a representative
``a//b ∈  𝕜[x₁,…,xₙ][S⁻¹]`` of the fraction. 
"""
lift(f::MPolyQuoLocRingElem) = localized_ring(f)(lifted_numerator(f), lifted_denominator(f))


@doc raw"""
    is_unit(f::MPolyQuoLocRingElem) 

Return `true`, if `f` is a unit of `parent(f)`, `true` otherwise.

# Examples

```jldoctest
julia> T, t = polynomial_ring(QQ, "t");

julia> K, a =  number_field(2*t^2-1, "a");

julia> R, (x, y) = polynomial_ring(K, ["x", "y"]);

julia> I = ideal(R, [2*x^2-y^3, 2*x^2-y^5])
ideal(2*x^2 - y^3, 2*x^2 - y^5)

julia> P = ideal(R, [y-1, x-a])
ideal(y - 1, x - a)

julia> U = complement_of_prime_ideal(P)
complement of ideal(y - 1, x - a)

julia> RQ, p = quo(R, I);

julia> RQL, iota = Localization(RQ, U);

julia> is_unit(iota(p(x)))
true
```
""" is_unit(f::MPolyQuoLocRingElem)

function is_unit(f::MPolyQuoLocRingElem) 
  lifted_numerator(f) in inverted_set(parent(f)) && return true
  R=localized_ring(parent(f))
  return one(R) in modulus(parent(f)) + ideal(R, lift(f))
end

function is_unit(L::MPolyQuoLocRing, f::MPolyLocRingElem) 
  parent(f) == localized_ring(L) || error("element does not belong to the correct ring")
  numerator(f) in inverted_set(L) && return true
  one(localized_ring(L)) in modulus(L) + ideal(localized_ring(L), f)
end

function is_unit(L::MPolyQuoLocRing{BRT, BRET, RT, RET, MST}, f::RET) where {BRT, BRET, RT, RET, MST}
  parent(f) == base_ring(L) || error("element does not belong to the correct ring")
  f in inverted_set(L) && return true
  return one(localized_ring(L)) in modulus(L) + ideal(localized_ring(L), localized_ring(L)(f))
end

function is_unit(L::MPolyQuoLocRing{BRT, BRET, RT, RET, MST}, f::MPolyQuoRingElem{RET}) where {BRT, BRET, RT, RET, MST}
  parent(f) == underlying_quotient(L) || error("element does not belong to the correct ring")
  lift(f) in inverted_set(L) && return true
  one(localized_ring(L)) in modulus(L) + ideal(localized_ring(L), localized_ring(L)(f))
end

# WARNING: This routine runs forever if f is not a unit in L. 
# So this needs to be checked first!
function inv(L::MPolyQuoLocRing{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}}, 
    f::MPolyQuoRingElem{RET}) where {BRT, BRET, RT, RET}
  Q = underlying_quotient(L)
  parent(f) == underlying_quotient(L) || error("element does not belong to the correct ring")
  W = localized_ring(L)
  R = base_ring(L)
  I = saturated_ideal(modulus(L))
  d = prod(denominators(inverted_set(W)))
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
  return parent(f)(denominator(f), numerator(f))
end

### 
# Assume that [f] = [a]//[b] is an admissible element of L = (R/I)[S⁻¹] and bring it 
# to the form [f] = [c]//[dᵏ] with d∈ S.
function convert(
    L::MPolyQuoLocRing{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}}, 
    f::AbstractAlgebra.Generic.Frac{RET}
  ) where {BRT, BRET, RT, RET}
  a = numerator(f)
  b = denominator(f)
  Q = underlying_quotient(L)
  parent(a) == base_ring(L) || error("element does not belong to the correct ring")
  W = localized_ring(L)
  R = base_ring(L)
  I = saturated_ideal(modulus(L))
  one(R) in I && return zero(L)
  d = prod(denominators(inverted_set(W)))
  powers_of_d = [d]
  ### apply logarithmic bisection to find a power a ⋅dᵏ ≡  c ⋅ b mod I
  (result, coefficient) = divides(Q(a), Q(b))
  # check whether f is already a unit
  result && return L(coefficient)
  push!(powers_of_d, d)
  abort = false
  # find some power which works
  while !abort
    (abort, coefficient) = divides(Q(a*last(powers_of_d)), Q(b))
    if !abort
      push!(powers_of_d, last(powers_of_d)^2)
    end
  end
  # find the minimal power that works
  upper = pop!(powers_of_d)
  lower = pop!(powers_of_d)
  while length(powers_of_d) > 0
    middle = lower*pop!(powers_of_d)
    (result, coefficient) = divides(Q(a*middle), Q(b))
    if result 
      upper = middle
    else 
      lower = middle
    end
  end
  return L(lift(coefficient), upper)
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
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  if lifted_denominator(a) == lifted_denominator(b) 
    return (parent(a))(lifted_numerator(a) + lifted_numerator(b), lifted_denominator(a), check=false)
  end
  return (parent(a))(lifted_numerator(a)*lifted_denominator(b) + lifted_numerator(b)*lifted_denominator(a), lifted_denominator(a)*lifted_denominator(b), check=false)
end

# TODO: improve this method.
function addeq!(a::T, b::T) where {T<:MPolyQuoLocRingElem}
  a = a+b
  return a
end

function -(a::T, b::T) where {T<:MPolyQuoLocRingElem}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  if lifted_denominator(a) == lifted_denominator(b) 
    return (parent(a))(lifted_numerator(a) - lifted_numerator(b), lifted_denominator(a), check=false)
  end
  return (parent(a))(lifted_numerator(a)*lifted_denominator(b) - lifted_numerator(b)*lifted_denominator(a), lifted_denominator(a)*lifted_denominator(b), check=false)
end

function *(a::T, b::T) where {T<:MPolyQuoLocRingElem}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  return (parent(a))(lifted_numerator(a)*lifted_numerator(b), lifted_denominator(a)*lifted_denominator(b), check=false)
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
# common factor g that can be cancelled so that b'= b/g ∈  Q belongs 
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

function divexact(a::Oscar.IntegerUnion, b::MPolyQuoLocRingElem)
  return a/b
end

function divexact(a::T, b::T) where {T<:MPolyQuoLocRingElem}
  return a/b
end

function ==(a::T, b::T) where {T<:MPolyQuoLocRingElem}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
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

### enhancement of the arithmetic
function reduce_fraction(f::MPolyQuoLocRingElem{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST<:MPolyPowersOfElement}
  return f # Disable reduction here, because it slows down arithmetic.
  return parent(f)(lift(simplify(numerator(f))), lifted_denominator(f), check=false)
end

# for local orderings, reduction does not give the correct result.
function reduce_fraction(f::MPolyQuoLocRingElem{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST<:MPolyComplementOfKPointIdeal}
  is_reduced(f) && return f
  return f
end

### implementation of Oscar's general ring interface
one(W::MPolyQuoLocRing) = W(one(base_ring(W)))
zero(W::MPolyQuoLocRing)= W(zero(base_ring(W)))

elem_type(W::MPolyQuoLocRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyQuoLocRingElem{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}
elem_type(T::Type{MPolyQuoLocRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyQuoLocRingElem{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}

parent_type(W::MPolyQuoLocRingElem{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyQuoLocRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}
parent_type(T::Type{MPolyQuoLocRingElem{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyQuoLocRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}



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
    R == base_ring(parent(a)) || error("elements do not belong to the same ring")
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
    parent(a) == L || error("elements do not belong to the same ring")
  end
  return L.(vec(coordinates(lift(f), ideal(L, g)))[1:length(g)]) # temporary hack; to be replaced.
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
      BaseRingType, 
      BaseRingElemType, 
      RingType, 
      RingElemType, 
      DomainMultSetType, 
      CodomainMultSetType
    } <: AbsLocalizedRingHom{
      RingType, RingElemType, DomainMultSetType, CodomainMultSetType
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
    if check
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
  return MPolyQuoLocalizedRingHom(L, S, hom(base_ring(L), S, a), check=check)
end

hom(L::MPolyQuoLocRing, S::Ring, res::Map; check::Bool=true) = MPolyQuoLocalizedRingHom(L, S, a, check=check)

function hom(L::MPolyQuoLocRing, S::Ring, a::Vector{T}; check::Bool=true) where {T<:RingElem}
  R = base_ring(L)
  res = hom(R, S, a)
  MPolyQuoLocalizedRingHom(L, S, res, check=check)
end

### implementing the Oscar map interface
function identity_map(W::T) where {T<:MPolyQuoLocRing} 
  MPolyQuoLocalizedRingHom(W, W, identity_map(base_ring(W)))
end

### we need to overwrite the following method because of the 
# uncommon implementation of the numerator and denominator methods
function (f::MPolyQuoLocalizedRingHom)(a::AbsLocalizedRingElem)
  parent(a) === domain(f) || return f(domain(f)(a))
  return codomain(f)(restricted_map(f)(lifted_numerator(a)))*inv(codomain(f)(restricted_map(f)(lifted_denominator(a))))
end

function compose(
    f::MPolyQuoLocalizedRingHom, 
    g::Hecke.Map{<:Ring, <:Ring}
  )
  codomain(f) === domain(g) || error("maps are not compatible")

  ### the following is commented out because as of now it's not type-stable!
#  if codomain(restricted_map(f)) === domain(g)
#    return MPolyQuoLocalizedRingHom(domain(f), codomain(g), compose(restricted_map(f), g))
#  elseif codomain(restricted_map(f)) === base_ring(domain(g)) 
#    h = hom(base_ring(domain(g)), domain(g), domain(g).(gens(base_ring(domain(g)))))
#    return MPolyQuoLocalizedRingHom(domain(f), codomain(g), compose(compose(restricted_map(f), h), g))
#  end
  ### The fallback version. Careful: This might not carry over maps on the coefficient rings!
  R = base_ring(domain(f))
  return MPolyQuoLocalizedRingHom(domain(f), codomain(g), hom(R, codomain(g), [g(f(x)) for x in gens(R)], check=false), check=false)
end

(f::MPolyQuoLocalizedRingHom)(I::Ideal) = ideal(codomain(f), f.(domain(f).(gens(I))))

function ==(f::MPolyQuoLocalizedRingHom, g::MPolyQuoLocalizedRingHom) 
  domain(f) === domain(g) || return false
  codomain(f) === codomain(g) || return false
  for x in gens(base_ring(domain(f)))
    f(x) == g(x) || return false
  end
  return true
end

### printing
function Base.show(io::IO, phi::MPolyQuoLocalizedRingHom)
  R = base_ring(domain(phi))
  psi = restricted_map(phi)
  println(io, "$(domain(phi)) → $(codomain(phi));")
  for i in 1:ngens(R)-1
    println(io, " $(R[i]) ↦ $(psi(R[i])),")
  end
  n = ngens(R)
  println(io, " $(R[n]) ↦ $(psi(R[n]))")
  return
end

### helper_ring
# Sets up the ring S[c⁻¹] from the Lemma.
function helper_ring(f::MPolyQuoLocalizedRingHom{<:Any, <:MPolyQuoLocRing})
  if !has_attribute(f, :helper_ring)
    minimal_denominators = Vector{base_ring_elem_type(domain(f))}()
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

    help_ring, help_kappa, theta = _add_variables(S, ["θ"])
    set_attribute!(f, :helper_ring, help_ring)
    kappa = help_kappa
    set_attribute!(f, :kappa, help_kappa)
    c_inv = theta[1]
    helper_images = [kappa(numerator(y))*c_inv*kappa(divexact(p, denominator(y))) for y in images(f)]
    set_attribute!(f, :helper_images, helper_images)
    eta = hom(R, help_ring, helper_images)
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
  return get_attribute(f, :helper_images)::Vector{base_ring_elem_type(domain(f))}
end

function minimal_denominators(
    f::MPolyQuoLocalizedRingHom{<:Any, <:MPolyQuoLocRing}
  )
  if !has_attribute(f, :minimal_denominators) 
    helper_ring(f)
  end
  return get_attribute!(f, :minimal_denominators)::Vector{base_ring_elem_type(domain(f))}
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
  d = get_attribute(f, :minimal_denominators)::Vector{base_ring_elem_type(domain(f))}
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
    inverse_name::String="θ"
  )
  R = base_ring(L)
  A, phi, t = _add_variables_first(R, [inverse_name])
  theta = t[1]
  f = prod(denominators(inverted_set(L)))
  I = ideal(A, [phi(g) for g in gens(modulus(underlying_quotient(L)))]) + ideal(A, [one(A)-theta*phi(f)])
  return A, I, f, phi, theta
end

# needed for instance to compute kernels
# adds a single extra variable to turn the localization into an affine_algebra
# return the isomorphism L -> SomeAffineAlgebra
function _as_affine_algebra(
    L::MPolyQuoLocRing{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement};
    inverse_name::String="θ"
  )
  R = base_ring(L)
  A, phi, t = _add_variables_first(R, [inverse_name])
  theta = t[1]
  f = prod(denominators(inverted_set(L)))
  I = ideal(A, [phi(g) for g in gens(modulus(underlying_quotient(L)))]) + ideal(A, [one(A)-theta*phi(f)])
  Q, _ = quo(A, I)
  id = hom(L, Q, gens(A)[2:end], check=false)
  id_inv = hom(Q, L, pushfirst!(gens(L), inv(L(f))), check=false)
  set_attribute!(id, :inverse, id_inv)
  set_attribute!(id_inv, :inverse, id)
  return id
end

### The following method is also required for the internals of the generic 
# kernel routine for localized rings.
function kernel(f::MPolyAnyMap{<:MPolyRing, <:MPolyQuoLocRing})
  P = domain(f)
  L = codomain(f)
  I = ideal(L, zero(L))
  R = base_ring(L)
  J = saturated_ideal(I)
  d = [lifted_denominator(g) for g in f.(gens(domain(f)))]
  W = MPolyQuoLocRing(R, modulus(underlying_quotient(L)), MPolyPowersOfElement(R, d))
  id =  _as_affine_algebra(W)
  A = codomain(id)
  h = hom(P, A, id.(f.(gens(P))))
  return preimage(h, ideal(A, id.(W.(gens(J)))))
end

function is_isomorphism(
    phi::MPolyQuoLocalizedRingHom{T, T}
  ) where {T<:MPolyQuoLocRing}
  if has_attribute(phi, :inverse)
    return true
  end
  K = domain(phi)
  L = codomain(phi)
  A, I, d1, inc1, theta1 = as_affine_algebra(K, inverse_name="s")
  B, J, d2, inc2, theta2 = as_affine_algebra(L, inverse_name="t")

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
  pushfirst!(imagesB, prod(denoms))

  # perform a sanity check
  phiAB = hom(A, B, imagesB)
  issubset(ideal(B, [phiAB(g) for g in gens(I)]), J) || error("the homomorphism is not well defined")

  # assemble a common ring in which the equations for the graph of phi can 
  # be realized.
  C, j1, B_vars = _add_variables_first(A, String.(symbols(B)))
  j2 = hom(B, C, B_vars)
  G = ideal(C, [j1(gen(A, i)) - j2(imagesB[i]) for i in 1:ngens(A)]) + ideal(C, j2.(gens(J))) + ideal(C, j1.(gens(I)))
  singC, _ = Singular.polynomial_ring(Oscar.singular_coeff_ring(base_ring(C)), 
				  String.(symbols(C)),  
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
function _add_variables(R::RingType, v::Vector{String}) where {RingType<:MPolyRing}
  ext_R, _ = polynomial_ring(coefficient_ring(R), vcat(symbols(R), Symbol.(v)))
  n = ngens(R)
  phi = hom(R, ext_R, gens(ext_R)[1:n])
  return ext_R, phi, gens(ext_R)[(n+1):ngens(ext_R)]
end

function _add_variables_first(R::RingType, v::Vector{String}) where {RingType<:MPolyRing}
  ext_R, _ = polynomial_ring(coefficient_ring(R), vcat(Symbol.(v), symbols(R)))
  n = ngens(R)
  phi = hom(R, ext_R, gens(ext_R)[1+length(v):n+length(v)])
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
  singular_assure(J)
  R = base_ring(L)
  SR = singular_poly_ring(R)
  SJ = J.gens.S

  # collect the output from elimpart in Singular
  l = Singular.LibPresolve.elimpart(SJ)

  # set up the ring with the fewer variables 
  kept_var_symb = [symbols(R)[i] for i in 1:ngens(R) if !iszero(l[4][i])]
  Rnew, new_vars = polynomial_ring(coefficient_ring(R), kept_var_symb)

  # and the maps to go back and forth
  subst_map_R = hom(R, R, R.(gens(l[5])))
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
  proj_map = hom(R, Rnew, imgs)

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
  finv = hom(R, L, gens(L))
  f = hom(R, Lnew, gens(Lnew))
  return Lnew, hom(L, Lnew, f), hom(Lnew, L, finv)
end

function simplify(L::MPolyQuoRing)
  J = modulus(L)
  singular_assure(J)
  R = base_ring(L)
  SR = singular_poly_ring(R)
  SJ = J.gens.S

  # collect the output from elimpart in Singular
  l = Singular.LibPresolve.elimpart(SJ)

  # set up the ring with the fewer variables 
  kept_var_symb = [symbols(R)[i] for i in 1:ngens(R) if !iszero(l[4][i])]
  Rnew, new_vars = polynomial_ring(coefficient_ring(R), kept_var_symb, cached=false)

  # and the maps to go back and forth
  subst_map_R = hom(R, R, R.(gens(l[5])))
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
  proj_map = hom(R, Rnew, imgs)

  # the full substitution map 
  f = compose(subst_map_R, proj_map)

  # the transformed ideal
  Jnew = ideal(Rnew, f.(gens(J)))
  Lnew, _ = quo(Rnew, Jnew)

  # the inverse of the identification map
  fres = hom(L, Lnew, Lnew.(f.(gens(R))), check=true)
  fresinv = hom(Lnew, L, [L(R(a)) for a in gens(l[4]) if !iszero(a)], check=false)

  return Lnew, fres, fresinv
end

function simplify(R::MPolyRing)
  Rnew, new_vars = polynomial_ring(coefficient_ring(R), symbols(R), cached=false)
  f = hom(R, Rnew, gens(Rnew))
  finv = hom(Rnew, R, gens(R))
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
  map_from_base_ring::Hecke.Map

  J::MPolyLocalizedIdealType
 
  function MPolyQuoLocalizedIdeal(
      W::MPolyQuoLocRing, 
      g::Vector{LocRingElemType};
      map_from_base_ring::Hecke.Map = MapFromFunc(
          x->W(x),
          y->(isone(lifted_denominator(y)) ? lifted_numerator(y) : divexact(lifted_numerator(y), lifted_denominator(y))),
          base_ring(W), 
          W
        )
    ) where {LocRingElemType<:MPolyQuoLocRingElem}
    for f in g
      parent(f) == W || error("generator is not an element of the given ring")
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
getindex(I::MPolyQuoLocalizedIdeal, i::Int) = I.gens[i]
base_ring(I::MPolyQuoLocalizedIdeal) = I.W

### additional getter functions 
map_from_base_ring(I::MPolyQuoLocalizedIdeal) = I.map_from_base_ring
pre_image_ideal(I::MPolyQuoLocalizedIdeal) = I.J
ngens(I::MPolyQuoLocalizedIdeal) = length(I.gens)

### a shorthand notation for any MPolyIdeal 
MPolyAnyIdeal = Union{MPolyIdeal, MPolyQuoIdeal,
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
julia> R,(x,y,z,w) = QQ["x","y","z","w"];

julia> Q = ideal(R,[x*y-z*w]);

julia> RQ,phiQ = quo(R,Q);

julia> T = MPolyComplementOfKPointIdeal(R,[0,0,0,0]);

julia> RQL, phiQL = Localization(RQ,T);

julia> I = ideal(RQL,RQL.([x,z]))
ideal in Localization of Quotient of Multivariate Polynomial Ring in x, y, z, w over Rational Field by ideal(x*y - z*w) at the multiplicative set complement of maximal ideal corresponding to point with coordinates QQFieldElem[0, 0, 0, 0] generated by [x, z]

julia> J = ideal(RQL,RQL.([y]))
ideal in Localization of Quotient of Multivariate Polynomial Ring in x, y, z, w over Rational Field by ideal(x*y - z*w) at the multiplicative set complement of maximal ideal corresponding to point with coordinates QQFieldElem[0, 0, 0, 0] generated by [y]

julia> intersect(I,J)
ideal in Localization of Quotient of Multivariate Polynomial Ring in x, y, z, w over Rational Field by ideal(x*y - z*w) at the multiplicative set complement of maximal ideal corresponding to point with coordinates QQFieldElem[0, 0, 0, 0] generated by [z*w, y*z, x*y]

julia> intersect([I,J])
ideal in Localization of Quotient of Multivariate Polynomial Ring in x, y, z, w over Rational Field by ideal(x*y - z*w) at the multiplicative set complement of maximal ideal corresponding to point with coordinates QQFieldElem[0, 0, 0, 0] generated by [z*w, y*z, x*y]
```
"""
function intersect(I::MPolyQuoLocalizedIdeal, J::MPolyQuoLocalizedIdeal)
  L = base_ring(I)
  L == base_ring(J) || error("ideals must be defined in the same ring")
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
    base_ring(K) == L || error("base rings must match")
    erg = intersect(erg,pre_image_ideal(K))
  end
  return L(erg)
end

function intersect(VI::Vector{<:MPolyQuoLocalizedIdeal{T}}) where T
  @assert length(VI) != 0
  L = base_ring(VI[1])
  all(J -> base_ring(J) == L,VI) || error("base rings must match")
  VIpre = [pre_image_ideal(J) for J in VI]
  erg = intersect(VIpre)
  return L(erg)
end

### Basic functionality
function ideal_membership(a::RingElem, I::MPolyQuoLocalizedIdeal)
  L = base_ring(I)
  parent(a) == L || return L(a) in I
  return lift(a) in pre_image_ideal(I)
end

function coordinates(a::RingElem, I::MPolyQuoLocalizedIdeal)
  L = base_ring(I)
  parent(a) == L || return coordinates(L(a), I)
  a in I || error("the given element is not in the ideal")
  x = coordinates(lift(a), pre_image_ideal(I))
  return map_entries(L, x[1, 1:ngens(I)])
end

function saturated_ideal(I::MPolyQuoLocalizedIdeal)
  return saturated_ideal(pre_image_ideal(I))
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

### printing
function Base.show(io::IO, I::MPolyQuoLocalizedIdeal)
  if ngens(I) == 0
    print(io, "zero ideal in $(base_ring(I))")
    return
  end
  str = "ideal in $(base_ring(I)) generated by [$(first(gens(I)))"
  for i in 2:ngens(I)
    str = str * ", $(gen(I, i))"
  end
  str = str * "]"
  print(io, str)
  return
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
  base_ring(I) == A || error("ideal does not belong to the correct ring")
  R = base_ring(A)
  Q, _ = quo(R, modulus(A) + ideal(R, lift.(gens(I))))
  return Q, hom(A, Q, Q.(gens(R)), check=false)
end

function divides(a::MPolyQuoLocRingElem, b::MPolyQuoLocRingElem)
  W = parent(a)
  W == parent(b) || error("elements do not belong to the same ring")
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

function jacobi_matrix(f::MPolyQuoLocRingElem)
  L = parent(f)
  n = nvars(base_ring(L))
  return matrix(L, n, 1, [derivative(f, i) for i=1:n])
end

function jacobi_matrix(g::Vector{<:MPolyQuoLocRingElem})
  L = parent(g[1])
  n = nvars(base_ring(L))
  @assert all(x->parent(x) == L, g)
  return matrix(L, n, length(g), [derivative(x, i) for i=1:n for x = g])
end

@attr function is_prime(I::MPolyQuoLocalizedIdeal)
  return is_prime(saturated_ideal(I))
end

@attr function _is_integral_domain(W::MPolyQuoLocRing)
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

@attr MPolyQuoLocalizedIdeal function radical(I::MPolyQuoLocalizedIdeal)
  W = base_ring(I)
  J = pre_image_ideal(I)
  return ideal(W, [g for g in W.(gens(radical(J))) if !iszero(g)])
end

@attr function dim(I::MPolyQuoLocalizedIdeal)
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
function primary_decomposition(I::Union{<:MPolyQuoIdeal, <:MPolyQuoLocalizedIdeal, <:MPolyLocalizedIdeal})
  Q = base_ring(I)
  R = base_ring(Q)
  decomp = primary_decomposition(saturated_ideal(I))
  result = [(ideal(Q, Q.(gens(a))), ideal(Q, Q.(gens(b)))) for (a, b) in decomp]
  erase = Int[]
  for i in 1:length(result)
    # if component is trivial, erase it
    is_one(result[i][2]) || continue
    push!(erase,i)
  end
  deleteat!(result, erase)
  return result
end

@doc raw"""
    minimal_primes(I::Union{<:MPolyQuoIdeal, <:MPolyQuoLocalizedIdeal, <:MPolyLocalizedIdeal})

Return the minimal associated primes of I.
"""
function minimal_primes(I::Union{<:MPolyQuoIdeal, <:MPolyQuoLocalizedIdeal, <:MPolyLocalizedIdeal})
  Q = base_ring(I)
  R = base_ring(Q)
  decomp = minimal_primes(saturated_ideal(I))
  result = [ideal(Q, Q.(gens(b))) for b in decomp]
  erase = Int[]
  for i in 1:length(result)
    # if a component is trivial, erase it
    is_one(result[i]) || continue
    push!(erase,i)
  end
  deleteat!(result, erase)
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
  R == base_ring(J) || error("Ideals do not live in the same ring.")

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
  R == base_ring(J) || error("Ideals do not live in the same ring.")
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
  if !has_attribute(f, :lifted_map)
    S = domain(f)
    W = codomain(f)
    L = localized_ring(W)
    g = hom(S, L, lift.(f.img_gens))
    set_attribute!(f, :lifted_map, g)
  end
  g = get_attribute(f, :lifted_map)
  return codomain(f)(g(a), check=false)
end

function (f::Oscar.MPolyAnyMap{<:MPolyRing, <:MPolyQuoLocRing, <:MPolyQuoLocalizedRingHom})(a::MPolyRingElem)
  if !has_attribute(f, :lifted_map)
    S = domain(f)
    W = codomain(f)
    L = localized_ring(W)
    g = hom(S, L, x -> lift(f.coeff_map(x)), lift.(f.img_gens))
    set_attribute!(f, :lifted_map, g)
  end
  g = get_attribute(f, :lifted_map)
  return codomain(f)(g(a), check=false)
end
