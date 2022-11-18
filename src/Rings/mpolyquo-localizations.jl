import AbstractAlgebra: Ring, RingElem, Generic.Frac
import Base: issubset

export MPolyQuoLocalizedRing
export parent, inverted_set, base_ring, underlying_quotient, localized_ring, modulus, gens
export Localization

export MPolyQuoLocalizedRingElem
export numerator, denominator, parent, lift, is_unit, inv, convert, lifted_numerator, lifted_denominator, fraction

export MPolyQuoLocalizedRingHom
export domain, codomain, images, morphism_type, domain_type, codomain_type, restricted_map_type, ideal_type
export helper_ring, helper_images, minimal_denominators, helper_eta, helper_kappa, common_denominator, helper_ideal

export MPolyQuoLocalizedIdeal

export is_isomorphism, inverse

export simplify

export mult_set_type

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
#    on an element of type `MPolyQuoLocalizedRingElem` is 
#    not `RingElemType`, but the type of `P`. 
#
# This is to comply with the purely mathematical viewpoint
# where elements of localized rings are fractions of 
# residue classes rather than residue classes of fractions. 
#

@Markdown.doc """
    MPolyQuoLocalizedRing{
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
@attributes mutable struct MPolyQuoLocalizedRing{
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
  Q::MPolyQuo{RingElemType}
  W::MPolyLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}

  function MPolyQuoLocalizedRing(
      R::RingType,
      I::MPolyIdeal{RingElemType},
      S::MultSetType,
      Q::MPolyQuo{RingElemType},
      W::MPolyLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}
    ) where {
      BaseRingType<:Ring, 
      BaseRingElemType<:RingElem, 
      RingType<:MPolyRing, 
      RingElemType<:MPolyElem, 
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
MPAnyQuoRing = Union{MPolyQuoLocalizedRing, 
                MPolyQuo
               }

MPAnyNonQuoRing = Union{MPolyRing, MPolyLocalizedRing
                  }

MPolyAnyRing = Union{MPolyRing, MPolyQuo,
                MPolyLocalizedRing,MPolyQuoLocalizedRing
               }

### type getters 
coefficient_ring_type(::Type{MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = BRT
coefficient_ring_type(L::MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST} = coefficient_ring_type(typeof(L))
coefficient_ring_elem_type(::Type{MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = BRET
coefficient_ring_elem_type(L::MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST} = coefficient_ring_elem_type(typeof(L))

base_ring_type(::Type{MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = RT
base_ring_type(L::MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST} = base_ring_type(typeof(L))
base_ring_elem_type(::Type{MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = RET
base_ring_elem_type(L::MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST} = base_ring_elem_type(typeof(L))

mult_set_type(::Type{MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = MST
mult_set_type(L::MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST} = mult_set_type(typeof(L))

localized_ring_type(::Type{MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = MPolyLocalizedRing{BRT, BRET, RT, RET, MST}
localized_ring_type(L::MPolyQuoLocalizedRing) = localized_ring_type(typeof(L))

ideal_type(::Type{MPolyQuoLocalizedRingType}) where {MPolyQuoLocalizedRingType<:MPolyQuoLocalizedRing} = MPolyQuoLocalizedIdeal{MPolyQuoLocalizedRingType, elem_type(MPolyQuoLocalizedRingType), ideal_type(localized_ring_type(MPolyQuoLocalizedRingType))}
ideal_type(W::MPolyQuoLocalizedRing) = ideal_type(typeof(W))





### required getter functions 
base_ring(L::MPolyQuoLocalizedRing) = L.R
inverted_set(L::MPolyQuoLocalizedRing) = L.S

### additional getter functions
@Markdown.doc """
    modulus(L::MPolyQuoLocalizedRing)

For ``L = (𝕜[x₁,…,xₙ]/I)[S⁻¹]`` this returns ``IS⁻¹``.
"""
function modulus(L::MPolyQuoLocalizedRing) 
  if !has_attribute(L, :modulus)
    set_attribute!(L, :modulus, localized_ring(L)(L.I))
  end
  return get_attribute(L, :modulus)::ideal_type(localized_ring_type(L))
end

### for compatibility -- also provide modulus in the trivial case
modulus(R::MPAnyNonQuoRing)=ideal(R,[zero(R)])


@Markdown.doc """
    underlying_quotient(L::MPolyQuoLocalizedRing)

For ``L = (𝕜[x₁,…,xₙ]/I)[S⁻¹]`` this returns ``𝕜[x₁,…,xₙ]/I``.
"""
underlying_quotient(L::MPolyQuoLocalizedRing) = L.Q

## 3 more signatures for compatibility to make quotient_ring agnostic
quotient_ring(L::MPolyQuo) = L

@attr MPolyQuo function quotient_ring(L::MPolyRing)
   return quo(L,ideal(L,[zero(L)]))[1]
end

@attr MPolyQuo function quotient_ring(L::MPolyLocalizedRing)
   P = base_ring(L)
   return quo(P,ideal(P,[zero(P)]))[1]
end

@Markdown.doc """
    localized_ring(L::MPolyQuoLocalizedRing)

For ``L = (𝕜[x₁,…,xₙ]/I)[S⁻¹]`` this returns ``𝕜[x₁,…,xₙ][S⁻¹]``.
"""
localized_ring(L::MPolyQuoLocalizedRing) = L.W

## 3 more signatures for compatibility to make localized_ring agnostic
localized_ring(L::MPolyLocalizedRing) = L

@attr MPolyLocalizedRing function localized_ring(L::MPolyQuo)
   P = base_ring(L)
   return localization(P, units_of(P))[1]
end

@attr MPolyLocalizedRing function localized_ring(L::MPolyRing)
   return localization(L, units_of(L))[1]
end

@Markdown.doc """
    gens(L::MPolyQuoLocalizedRing)

For ``L = (𝕜[x₁,…,xₙ]/I)[S⁻¹]`` this returns the vector ``[x₁//1,…,xₙ//1]∈ Lⁿ``.
"""
gens(L::MPolyQuoLocalizedRing) = L.(gens(base_ring(L)))

### printing
function Base.show(io::IO, L::MPolyQuoLocalizedRing)
  print(io, "Localization of $(underlying_quotient(L)) at the multiplicative set $(inverted_set(L))")
end

### additional constructors
function quo(
    W::MPolyLocalizedRing,
    I::MPolyLocalizedIdeal
  )
  R = base_ring(W)
  S = inverted_set(W)
  J = ideal(R, numerator.(gens(I)))
  L = MPolyQuoLocalizedRing(R, J, S, quo(R, J)[1], W)
  return L, hom(R, L, gens(L))
end

function quo(
    L::MPolyQuoLocalizedRing,
    I::MPolyLocalizedIdeal
  )
  R = base_ring(L)
  S = inverted_set(L)
  base_ring(I) = localized_ring(L) || error("ideal does not belong to the correct ring")
  J = pre_saturated_ideal(I)
  W = MPolyQuoLocalizedRing(R, J, S, quo(R, J)[1], localized_ring(L))
  return W, hom(L, W, gens(W))
end

function quo(
    L::MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MST},
    J::MPolyIdeal{RET}
  ) where {BRT, BRET, RT, RET, MST}
  R = base_ring(L)
  S = inverted_set(L)
  W = localized_ring(L) 
  J = J + modulus(underlying_quotient(L))
  P = MPolyQuoLocalizedRing(R, J, S, quo(R, J)[1], W)
  return P, hom(L, P, gens(P))
end

function Localization(Q::MPolyQuo{RET}, S::MultSetType) where {RET <: RingElem, MultSetType <: AbsMultSet}
  L = MPolyQuoLocalizedRing(base_ring(Q), modulus(Q), S, Q, Localization(S)[1])
  return L, MapFromFunc((x->L(lift(x))), Q, L)
end

function Localization(
    L::MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MST}, 
    S::AbsMPolyMultSet{BRT, BRET, RT, RET}
  ) where {BRT, BRET, RT, RET, MST}
  ambient_ring(S) == base_ring(L) || error("multiplicative set does not belong to the correct ring")
  issubset(S, inverted_set(L)) && return L, MapFromFunc(x->x, L, L)
  U = inverted_set(L)*S
  W = MPolyQuoLocalizedRing(base_ring(L), modulus(underlying_quotient(L)), U, underlying_quotient(L), Localization(U)[1])
  return W, MapFromFunc((x->W(lifted_numerator(x), lifted_denominator(x), check=false)), L, W)
end

function MPolyQuoLocalizedRing(R::RT, I::Ideal{RET}, T::MultSetType) where {RT<:MPolyRing, RET<:MPolyElem, MultSetType<:AbsMultSet} 
  return MPolyQuoLocalizedRing(R, I, T, quo(R, I)[1], Localization(T)[1])
end

function MPolyQuoLocalizedRing(R::RT) where {RT<:MPolyRing} 
  I = ideal(R, zero(R))
  Q, _ = quo(R, I)
  U = units_of(R)
  W, _ = Localization(U)
  return MPolyQuoLocalizedRing(R, I, U, Q, W)
end

function MPolyQuoLocalizedRing(Q::RT) where {RT<:MPolyQuo}
  R = base_ring(Q)
  I = modulus(Q)
  U = units_of(R)
  W, _ = Localization(U)
  return MPolyQuoLocalizedRing(R, I, U, Q, W)
end

function MPolyQuoLocalizedRing(W::MPolyLocalizedRing)
  R = base_ring(W)
  I = ideal(R, zero(R))
  Q, _ = quo(R, I)
  U = inverted_set(W)
  return MPolyQuoLocalizedRing(R, I, U, Q, W)
end

function Base.in(f::AbstractAlgebra.Generic.Frac{RET}, L::MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST}
  R = base_ring(L)
  R == parent(numerator(f)) || error("element does not belong to the correct ring")
  denominator(f) in inverted_set(L) && return true
  return numerator(f) in ideal(L, denominator(f))
end


### generation of random elements 
function rand(W::MPolyQuoLocalizedRing, v1::UnitRange{Int}, v2::UnitRange{Int}, v3::UnitRange{Int})
  return W(rand(localized_ring(W), v1, v2, v3))
end

########################################################################
# Elements of localizations of polynomial algebras                     #
########################################################################

@Markdown.doc """
    MPolyQuoLocalizedRingElem{
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
`MPolyQuoLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}`.
"""
mutable struct MPolyQuoLocalizedRingElem{
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
  L::MPolyQuoLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}
  # representatives of numerator and denominator
  numerator::RingElemType
  denominator::RingElemType
  is_reduced::Bool

  function MPolyQuoLocalizedRingElem(
      L::MPolyQuoLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}, 
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
coefficient_ring_type(::Type{MPolyQuoLocalizedRingElem{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = BRT
coefficient_ring_type(f::MPolyQuoLocalizedRingElem{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST} = base_ring_type(typeof(f))
coefficient_ring_elem_type(::Type{MPolyQuoLocalizedRingElem{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = BRET
coefficient_ring_elem_type(f::MPolyQuoLocalizedRingElem{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST} = base_ring_type(typeof(f))

base_ring_type(::Type{MPolyQuoLocalizedRingElem{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = RT
base_ring_type(f::MPolyQuoLocalizedRingElem{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST} = base_ring_type(typeof(f))
base_ring_elem_type(::Type{MPolyQuoLocalizedRingElem{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = RET
base_ring_elem_type(f::MPolyQuoLocalizedRingElem{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST} = base_ring_type(typeof(f))

mult_set_type(::Type{MPolyQuoLocalizedRingElem{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = MST
mult_set_type(f::MPolyQuoLocalizedRingElem{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST} = base_ring_type(typeof(f))

### required getter functions 
parent(a::MPolyQuoLocalizedRingElem) = a.L
numerator(a::MPolyQuoLocalizedRingElem) = underlying_quotient(parent(a))(a.numerator) 
denominator(a::MPolyQuoLocalizedRingElem) = underlying_quotient(parent(a))(a.denominator) 

### additional getter functions
underlying_quotient(a::MPolyQuoLocalizedRingElem) = underlying_quotient(parent(a))
localized_ring(a::MPolyQuoLocalizedRingElem) = localized_ring(parent(a))
base_ring(a::MPolyQuoLocalizedRingElem) = base_ring(parent(a))
is_reduced(a::MPolyQuoLocalizedRingElem) = a.is_reduced

@Markdown.doc """
    lifted_numerator(a::MPolyQuoLocalizedRingElem)

For ``A//B ∈ (𝕜[x₁,…,xₙ]/I)[S⁻¹]`` this returns a representative 
``a ∈ 𝕜[x₁,…,xₙ]`` of the numerator. 
"""
lifted_numerator(a::MPolyQuoLocalizedRingElem) = a.numerator

@Markdown.doc """
    lifted_denominator(a::MPolyQuoLocalizedRingElem)

For ``A//B ∈ (𝕜[x₁,…,xₙ]/I)[S⁻¹]`` this returns a representative 
``b ∈  𝕜[x₁,…,xₙ]`` of the denominator.
"""
lifted_denominator(a::MPolyQuoLocalizedRingElem) = a.denominator

@Markdown.doc """
    fraction(a::MPolyQuoLocalizedRingElem)

For ``A//B ∈ (𝕜[x₁,…,xₙ]/I)[S⁻¹]`` this returns a representative 
``a//b ∈ Quot(𝕜[x₁,…,xₙ])`` of the fraction. 
"""
fraction(a::MPolyQuoLocalizedRingElem) = lifted_numerator(a)//lifted_denominator(a)

### copying of elements
function Base.deepcopy_internal(f::MPolyQuoLocalizedRingElem, dict::IdDict)
  return parent(f)(f, check=false, is_reduced=is_reduced(f))
end

### required conversions
(L::MPolyQuoLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType})(f::RingElemType, is_reduced::Bool=false) where {BaseRingType, BaseRingElemType, RingType, RingElemType<:RingElem, MultSetType} = MPolyQuoLocalizedRingElem(L, f, one(f), check=false, is_reduced=is_reduced)

function (L::MPolyQuoLocalizedRing{
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
  check || return MPolyQuoLocalizedRingElem(L, a, b, check=false, is_reduced=is_reduced)
  b in inverted_set(L) || return convert(L, a//b)
  return MPolyQuoLocalizedRingElem(L, a, b, check=false, is_reduced=is_reduced)
end

### additional conversions
function (L::MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MST})(f::Frac{RET}; check::Bool=true, is_reduced::Bool=false) where {BRT, BRET, RT, RET, MST}
  R = base_ring(L)
  return L(R(numerator(f)), R(denominator(f)), check=check, is_reduced=is_reduced)
end

(L::MPolyQuoLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType})(a::T, b::T; check::Bool=true, is_reduced::Bool=false) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType, T<:MPolyQuoElem{RingElemType}} = L(lift(a), lift(b), check=check, is_reduced=is_reduced)

function (L::MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MST})(f::MPolyQuoLocalizedRingElem{BRT, BRET, RT, RET, MST}; check::Bool=true, is_reduced::Bool=false) where {BRT, BRET, RT, RET, MST}
  parent(f) === L && return f
  return L(lifted_numerator(f), lifted_denominator(f), check=check, is_reduced=is_reduced)
end

function (L::MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MST})(f::MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}; check::Bool=true, is_reduced::Bool=false) where {BRT, BRET, RT, RET, MST}
  parent(f) === localized_ring(L) && return L(numerator(f), denominator(f), check=false, is_reduced=is_reduced)
  return L(numerator(f), denominator(f), check=check, is_reduced=is_reduced)
end

function (L::MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MST})(f::MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}; check::Bool=true, is_reduced::Bool=false) where {BRT, BRET, RT, RET, MST<:MPolyComplementOfKPointIdeal}
  return L(numerator(f), denominator(f), check=check, is_reduced=is_reduced)
end

function (L::MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MST})(f::MPolyQuoElem{RET}; check::Bool=true, is_reduced::Bool=false) where {BRT, BRET, RT, RET, MST} 
  base_ring(parent(f)) == base_ring(L) || error("the given element does not belong to the correct ring") 
  check && (parent(f) == underlying_quotient(L) || all(x->(iszero(L(x))), gens(modulus(parent(f)))) || error("coercion is not well defined"))
  return L(lift(f))
end

### additional functionality
@Markdown.doc """
    lift(f::MPolyQuoLocalizedRingElem)

For ``f = A//B ∈ (𝕜[x₁,…,xₙ]/I)[S⁻¹]`` this returns a representative 
``a//b ∈  𝕜[x₁,…,xₙ][S⁻¹]`` of the fraction. 
"""
lift(f::MPolyQuoLocalizedRingElem) = localized_ring(f)(lifted_numerator(f), lifted_denominator(f))

function is_unit(f::MPolyQuoLocalizedRingElem) 
  lifted_numerator(f) in inverted_set(parent(f)) && return true
  R=localized_ring(parent(f))
  return one(R) in modulus(parent(f)) + ideal(R, lift(f))
end

function is_unit(L::MPolyQuoLocalizedRing, f::MPolyLocalizedRingElem) 
  parent(f) == localized_ring(L) || error("element does not belong to the correct ring")
  numerator(f) in inverted_set(L) && return true
  one(localized_ring(L)) in modulus(L) + ideal(localized_ring(L), f)
end

function is_unit(L::MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MST}, f::RET) where {BRT, BRET, RT, RET, MST}
  parent(f) == base_ring(L) || error("element does not belong to the correct ring")
  f in inverted_set(L) && return true
  return one(localized_ring(L)) in modulus(L) + ideal(localized_ring(L), localized_ring(L)(f))
end

function is_unit(L::MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MST}, f::MPolyQuoElem{RET}) where {BRT, BRET, RT, RET, MST}
  parent(f) == underlying_quotient(L) || error("element does not belong to the correct ring")
  lift(f) in inverted_set(L) && return true
  one(localized_ring(L)) in modulus(L) + ideal(localized_ring(L), localized_ring(L)(f))
end

# WARNING: This routine runs forever if f is not a unit in L. 
# So this needs to be checked first!
function inv(L::MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}}, 
    f::MPolyQuoElem{RET}) where {BRT, BRET, RT, RET}
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

function inv(f::MPolyQuoLocalizedRingElem{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}}) where {BRT, BRET, RT, RET}
  return parent(f)(denominator(f), numerator(f))
end

### 
# Assume that [f] = [a]//[b] is an admissible element of L = (R/I)[S⁻¹] and bring it 
# to the form [f] = [c]//[dᵏ] with d∈ S.
function convert(
    L::MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}}, 
    f::AbstractAlgebra.Generic.Frac{RET}
  ) where {BRT, BRET, RT, RET}
  a = numerator(f)
  b = denominator(f)
  Q = underlying_quotient(L)
  parent(a) == base_ring(L) || error("element does not belong to the correct ring")
  W = localized_ring(L)
  R = base_ring(L)
  I = saturated_ideal(modulus(L))
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


### arithmetic #########################################################
function +(a::T, b::T) where {T<:MPolyQuoLocalizedRingElem}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  if lifted_denominator(a) == lifted_denominator(b) 
    return (parent(a))(lifted_numerator(a) + lifted_numerator(b), lifted_denominator(a), check=false)
  end
  return (parent(a))(lifted_numerator(a)*lifted_denominator(b) + lifted_numerator(b)*lifted_denominator(a), lifted_denominator(a)*lifted_denominator(b), check=false)
end

# TODO: improve this method.
function addeq!(a::T, b::T) where {T<:MPolyQuoLocalizedRingElem}
  a = a+b
  return a
end

function -(a::T, b::T) where {T<:MPolyQuoLocalizedRingElem}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  if lifted_denominator(a) == lifted_denominator(b) 
    return (parent(a))(lifted_numerator(a) - lifted_numerator(b), lifted_denominator(a), check=false)
  end
  return (parent(a))(lifted_numerator(a)*lifted_denominator(b) - lifted_numerator(b)*lifted_denominator(a), lifted_denominator(a)*lifted_denominator(b), check=false)
end

function *(a::T, b::T) where {T<:MPolyQuoLocalizedRingElem}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  return (parent(a))(lifted_numerator(a)*lifted_numerator(b), lifted_denominator(a)*lifted_denominator(b), check=false)
end

function *(a::RET, b::MPolyQuoLocalizedRingElem{BRT, BRET, RT, RET, MST}) where {BRT<:Ring, BRET<:RingElem, RT<:Ring, RET <: RingElem, MST}
  return (parent(b))(a*lifted_numerator(b), lifted_denominator(b), check=false)
end

function *(a::MPolyQuoLocalizedRingElem{BRT, BRET, RT, RET, MST}, b::RET) where {BRT<:Ring, BRET<:RingElem, RT<:Ring, RET <: RingElem, MST}
  return b*a
end

function *(a::BRET, b::MPolyQuoLocalizedRingElem{BRT, BRET, RT, RET, MST}) where {BRT<:Ring, BRET<:RingElem, RT<:Ring, RET <: RingElem, MST}
  return (parent(b))(a*lifted_numerator(b), lifted_denominator(b), check=false)
end

function *(a::MPolyQuoLocalizedRingElem{BRT, BRET, RT, RET, MST}, b::BRET) where {BRT<:Ring, BRET<:RingElem, RT<:Ring, RET <: RingElem, MST}
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

function Base.:(//)(a::Oscar.IntegerUnion, b::MPolyQuoLocalizedRingElem)
  error("function `//` not implemented for elements of type $(typeof(b))")
end

function Base.:(//)(a::T, b::T) where {T<:MPolyQuoLocalizedRingElem}
  error("function `//` not implemented for elements of type $(typeof(b))")
end

function ==(a::T, b::T) where {T<:MPolyQuoLocalizedRingElem}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  return lifted_numerator(a)*lifted_denominator(b) - lifted_numerator(b)*lifted_denominator(a) in modulus(parent(a))
end

function ^(a::MPolyQuoLocalizedRingElem, i::fmpz)
  return parent(a)(lifted_numerator(a)^i, lifted_denominator(a)^i, check=false)
end

function ^(a::MPolyQuoLocalizedRingElem, i::Integer)
  return parent(a)(lifted_numerator(a)^i, lifted_denominator(a)^i, check=false)
end

function isone(a::MPolyQuoLocalizedRingElem) 
  return lifted_numerator(a) - lifted_denominator(a) in modulus(parent(a))
end

function iszero(a::MPolyQuoLocalizedRingElem)
  return lift(a) in modulus(parent(a))
end

### enhancement of the arithmetic
function reduce_fraction(f::MPolyQuoLocalizedRingElem{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST<:MPolyPowersOfElement}
  return f # Disable reduction here, because it slows down arithmetic.
  return parent(f)(lift(simplify(numerator(f))), lifted_denominator(f), check=false)
end

# for local orderings, reduction does not give the correct result.
function reduce_fraction(f::MPolyQuoLocalizedRingElem{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST<:MPolyComplementOfKPointIdeal}
  is_reduced(f) && return f
  return f
end

### implementation of Oscar's general ring interface
one(W::MPolyQuoLocalizedRing) = W(one(base_ring(W)))
zero(W::MPolyQuoLocalizedRing)= W(zero(base_ring(W)))

elem_type(W::MPolyQuoLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyQuoLocalizedRingElem{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}
elem_type(T::Type{MPolyQuoLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyQuoLocalizedRingElem{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}

parent_type(W::MPolyQuoLocalizedRingElem{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyQuoLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}
parent_type(T::Type{MPolyQuoLocalizedRingElem{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyQuoLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}



@Markdown.doc """
    bring_to_common_denominator(f::Vector{T}) where {T<:MPolyQuoLocalizedRingElem}

Given a vector of fractions ``[a₁//b₁,…,aₙ//bₙ]`` return a pair 
``(d, λ)`` consisting of a common denominator ``d`` and a vector 
``λ = [λ₁,…,λₙ]`` such that ``aᵢ//bᵢ = λᵢ⋅aᵢ//d``.
"""
function bring_to_common_denominator(f::Vector{T}) where {T<:MPolyQuoLocalizedRingElem}
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

@Markdown.doc """
    write_as_linear_combination(f::T, g::Vector{T}) where {T<:MPolyLocalizedRingElem} 

Write ``f = ∑ᵢ λᵢ⋅gᵢ`` for some ``λᵢ`` and return the vector ``[λ₁,…,λₙ]``.
"""
function write_as_linear_combination(
    f::RingElemType,
    g::Vector{RingElemType}
  ) where {RingElemType<:MPolyQuoLocalizedRingElem}
  n = length(g)
  L = parent(f)
  W = localized_ring(L)
  for a in g 
    parent(a) == L || error("elements do not belong to the same ring")
  end
  return L.(vec(coordinates(lift(f), ideal(L, g)))[1:length(g)]) # temporary hack; to be replaced.
end

write_as_linear_combination(f::MPolyQuoLocalizedRingElem, g::Vector) = write_as_linear_combination(f, parent(f).(g))


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

@Markdown.doc """
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

of types `MPolyQuoLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, DomainMultSetType}` and `MPolyQuoLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, CodomainMultSetType}`.
These are completely determined by the images of the 
variables ``ϕ(xᵢ) ∈ (𝕜[y₁,…,yₙ]/J)[T⁻¹]`` so that the 
constructor takes as input the triple 
``((𝕜[x₁,…,xₘ]/I)[S⁻¹], (𝕜[y₁,…,yₙ]/J)[T⁻¹], [ϕ(x₁),…,ϕ(xₘ)])``.
"""
@attributes mutable struct MPolyQuoLocalizedRingHom{
                                     DomainType<:MPolyQuoLocalizedRing, 
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
    ) where {DomainType<:MPolyQuoLocalizedRing, CodomainType<:Ring, RestrictedMapType<:Map}
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

morphism_type(::Type{R}, ::Type{S}) where {R<:MPolyQuoLocalizedRing, S<:Ring} = MPolyQuoLocalizedRingHom{R, S, morphism_type(base_ring_type(R), S)}
morphism_type(L::MPolyQuoLocalizedRing, S::Ring) = morphism_type(typeof(L), typeof(S))

### TODO: Move to other file
morphism_type(::Type{R}, ::Type{S}) where {R<:MPolyRing, S<:Ring} = Oscar.MPolyAnyMap{R, S, Nothing, elem_type(S)}
morphism_type(R::MPolyRing, S::Ring) = morphism_type(typeof(R), typeof(S))

### required getter functions
domain(f::MPolyQuoLocalizedRingHom) = f.domain
codomain(f::MPolyQuoLocalizedRingHom) = f.codomain

@Markdown.doc """
    restricted_map(f::MPolyQuoLocalizedRingHom)

For a homomorphism ``ϕ : (𝕜[x₁,…,xₘ]/I)[U⁻¹] → S``this returns 
the canonically associated map ``ϕ' : 𝕜[x₁,…,xₘ] → S``.
"""
restricted_map(f::MPolyQuoLocalizedRingHom) = f.res

function images(f::MPolyQuoLocalizedRingHom{<:Any, <:MPolyQuoLocalizedRing})
  return lift.((codomain(f)).(restricted_map(f).(gens(base_ring(domain(f))))))
end

### additional constructors
function MPolyQuoLocalizedRingHom(
    L::MPolyQuoLocalizedRing,
    S::Ring,
    a::Vector{T};
    check::Bool=true
  ) where {T<:RingElem}
  return MPolyQuoLocalizedRingHom(L, S, hom(base_ring(L), S, a), check=check)
end

hom(L::MPolyQuoLocalizedRing, S::Ring, a::Vector{T}; check::Bool=true) where {T<:RingElem} = MPolyQuoLocalizedRingHom(L, S, a, check=check)

### implementing the Oscar map interface
function identity_map(W::T) where {T<:MPolyQuoLocalizedRing} 
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
    g::MPolyQuoLocalizedRingHom
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
  return MPolyQuoLocalizedRingHom(domain(f), codomain(g), hom(R, codomain(g), [g(f(x)) for x in gens(R)]))
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

### helper_ring
# Sets up the ring S[c⁻¹] from the Lemma.
function helper_ring(f::MPolyQuoLocalizedRingHom{<:Any, <:MPolyQuoLocalizedRing})
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
    f::MPolyQuoLocalizedRingHom{<:Any, <:MPolyQuoLocalizedRing}
  )
  if !has_attribute(f, :helper_images) 
    helper_ring(f)
  end
  return get_attribute(f, :helper_images)::Vector{base_ring_elem_type(domain(f))}
end

function minimal_denominators(
    f::MPolyQuoLocalizedRingHom{<:Any, <:MPolyQuoLocalizedRing}
  )
  if !has_attribute(f, :minimal_denominators) 
    helper_ring(f)
  end
  return get_attribute!(f, :minimal_denominators)::Vector{base_ring_elem_type(domain(f))}
end

function helper_eta(
    f::MPolyQuoLocalizedRingHom{<:Any, <:MPolyQuoLocalizedRing}
  )
  if !has_attribute(f, :eta) 
    helper_ring(f)
  end
  return get_attribute(f, :eta)::morphism_type(base_ring_type(domain(f)), base_ring_type(domain(f)))
end

function helper_kappa(
    f::MPolyQuoLocalizedRingHom{<:Any, <:MPolyQuoLocalizedRing}
  )
  if !has_attribute(f, :kappa) 
    helper_ring(f)
  end
  return get_attribute(f, :kappa)::morphism_type(base_ring_type(domain(f)), base_ring_type(domain(f)))
end

function common_denominator(
    f::MPolyQuoLocalizedRingHom{<:Any, <:MPolyQuoLocalizedRing}
  )
  if !has_attribute(f, :minimal_denominators) 
    helper_ring(f)
  end
  d = get_attribute(f, :minimal_denominators)::Vector{base_ring_elem_type(domain(f))}
  return (length(d) == 0 ? one(base_ring(codomain(f))) : prod(d))
end

function helper_ideal(
    f::MPolyQuoLocalizedRingHom{<:Any, <:MPolyQuoLocalizedRing}
  )
  Sc = helper_ring(f)
  return ideal(Sc, one(Sc)-last(gens(Sc))*helper_kappa(f)(common_denominator(f)))
end

# return the localized ring as a quotient of a polynomial ring using Rabinowitsch's trick.
function as_affine_algebra(
    L::MPolyQuoLocalizedRing{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement};
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
    L::MPolyQuoLocalizedRing{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement};
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
function kernel(f::MPolyAnyMap{<:MPolyRing, <:MPolyQuoLocalizedRing})
  P = domain(f)
  L = codomain(f)
  I = ideal(L, zero(L))
  R = base_ring(L)
  J = saturated_ideal(I)
  d = [lifted_denominator(g) for g in f.(gens(domain(f)))]
  W = MPolyQuoLocalizedRing(R, modulus(underlying_quotient(L)), MPolyPowersOfElement(R, d))
  id =  _as_affine_algebra(W)
  A = codomain(id)
  h = hom(P, A, id.(f.(gens(P))))
  return preimage(h, ideal(A, id.(W.(gens(J)))))
end

function is_isomorphism(
    phi::MPolyQuoLocalizedRingHom{T, T}
  ) where {T<:MPolyQuoLocalizedRing}
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
    G, M = groebner_basis_with_transformation_matrix(J_ext)
    G[1]==one(B) || error("the denominator is not a unit in the target ring")
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
    G, M = groebner_basis_with_transformation_matrix(J_ext)
    G[1]==one(B) || error("the denominator is not a unit in the target ring")
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
  G = ideal(C, [j1(gens(A)[i]) - j2(imagesB[i]) for i in (1:length(gens(A)))]) + ideal(C, j2.(gens(J))) + ideal(C, j1.(gens(I)))
  singC, _ = Singular.PolynomialRing(Oscar.singular_coeff_ring(base_ring(C)), 
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
  #singp = Singular.reduce(gens(singC)[1], stdG)
  #singp < gens(singC)[n+1] || return false
  #p = C(singp)
  #push!(pre_imagesA, evaluate(p, vcat([zero(A) for i in 0:n], gens(A))))
  for i in 1:n
    singp = Singular.reduce(gens(singC)[i+1], stdG)
    singp < gens(singC)[n+1] || return false
    p = C(singp)
    # Write p as an element in the very original ring
    #push!(pre_imagesA, evaluate(p, vcat([zero(A) for i in 0:n], gens(A))))
    push!(pre_images, evaluate(p, vcat([zero(V) for i in 0:n], [V(one(R), d1)], V.(gens(R)))))
  end

  invJ = ideal(A, [(p < gens(singC)[n+1] ? evaluate(C(p), vcat([zero(A) for i in 0:n], gens(A))) : zero(A)) for p in gens(stdG)])
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
                                <:MPolyQuoLocalizedRing{<:Any, <:Any, <:Any, <:Any, 
                                                        <:MPolyPowersOfElement
                                                       },
                                <:MPolyQuoLocalizedRing{<:Any, <:Any, <:Any, <:Any, 
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
  ext_R, _ = PolynomialRing(coefficient_ring(R), vcat(symbols(R), Symbol.(v)))
  n = length(gens(R))
  phi = hom(R, ext_R, gens(ext_R)[1:n])
  return ext_R, phi, gens(ext_R)[(length(gens(R))+1):length(gens(ext_R))]
end

function _add_variables_first(R::RingType, v::Vector{String}) where {RingType<:MPolyRing}
  ext_R, _ = PolynomialRing(coefficient_ring(R), vcat(Symbol.(v), symbols(R)))
  n = length(gens(R))
  phi = hom(R, ext_R, gens(ext_R)[1+length(v):n+length(v)])
  return ext_R, phi, gens(ext_R)[(1:length(v))]
end

@Markdown.doc """
    simplify(L::MPolyQuoLocalizedRing{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement})

Use `elimpart` from the `Singular` library `Presolve.lib` to simplify the presentation 
of `L` by eliminating superfluous variables; return a triple ``(L', f, g)`` where 
``L' ≅ L`` and ``f : L ↔ L' : g`` are the identifying isomorphisms.
"""
function simplify(L::MPolyQuoLocalizedRing{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement})
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
  Rnew, new_vars = PolynomialRing(coefficient_ring(R), kept_var_symb)

  # and the maps to go back and forth
  subst_map_R = hom(R, R, R.(gens(l[5])))
  imgs = Vector{elem_type(Rnew)}()
  j = 1
  for i in 1:ngens(R)
    if !iszero(l[4][i])
      push!(imgs, gens(Rnew)[j])
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
  Lnew = MPolyQuoLocalizedRing(Rnew, Jnew, Unew)

  # the localized map and its inverse
  floc = hom(L, Lnew, Lnew.(f.(gens(R))), check=false)
  flocinv = hom(Lnew, L, [L(R(a)) for a in gens(l[4]) if !iszero(a)], check=false)

  return Lnew, floc, flocinv
end

function simplify(L::MPolyLocalizedRing{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement})
  Lnew = MPolyLocalizedRing(base_ring(L), inverted_set(L))
  R = base_ring(L)
  finv = hom(R, L, gens(L))
  f = hom(R, Lnew, gens(Lnew))
  return Lnew, hom(L, Lnew, f), hom(Lnew, L, finv)
end

function simplify(L::MPolyQuo)
  J = modulus(L)
  singular_assure(J)
  R = base_ring(L)
  SR = singular_poly_ring(R)
  SJ = J.gens.S

  # collect the output from elimpart in Singular
  l = Singular.LibPresolve.elimpart(SJ)

  # set up the ring with the fewer variables 
  kept_var_symb = [symbols(R)[i] for i in 1:ngens(R) if !iszero(l[4][i])]
  Rnew, new_vars = PolynomialRing(coefficient_ring(R), kept_var_symb, cached=false)

  # and the maps to go back and forth
  subst_map_R = hom(R, R, R.(gens(l[5])))
  imgs = Vector{elem_type(Rnew)}()
  j = 1
  for i in 1:ngens(R)
    if !iszero(l[4][i])
      push!(imgs, gens(Rnew)[j])
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
  Rnew, new_vars = PolynomialRing(coefficient_ring(R), symbols(R), cached=false)
  f = hom(R, Rnew, gens(Rnew))
  finv = hom(Rnew, R, gens(R))
  return Rnew, f, finv
end

@Markdown.doc """
    MPolyQuoLocalizedIdeal{
        LocRingType<:MPolyQuoLocalizedRing, 
        LocRingElemType<:MPolyQuoLocalizedRingElem
      } <: AbsLocalizedIdeal{LocRingElemType}

Ideals in localizations of affine algebras.
"""
@attributes mutable struct MPolyQuoLocalizedIdeal{
     LocRingType<:MPolyQuoLocalizedRing, 
     LocRingElemType<:MPolyQuoLocalizedRingElem, 
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
      W::MPolyQuoLocalizedRing, 
      g::Vector{LocRingElemType};
      map_from_base_ring::Hecke.Map = MapFromFunc(
          x->W(x),
          y->(isone(lifted_denominator(y)) ? lifted_numerator(y) : divexact(lifted_numerator(y), lifted_denominator(y))),
          base_ring(W), 
          W
        )
    ) where {LocRingElemType<:MPolyQuoLocalizedRingElem}
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
base_ring(I::MPolyQuoLocalizedIdeal) = I.W

### additional getter functions 
map_from_base_ring(I::MPolyQuoLocalizedIdeal) = I.map_from_base_ring
pre_image_ideal(I::MPolyQuoLocalizedIdeal) = I.J
ngens(I::MPolyQuoLocalizedIdeal) = length(I.gens)
getindex(I::MPolyQuoLocalizedIdeal, k::Int) = copy(I.gens[k])

### Additional constructors
function intersect(I::MPolyQuoLocalizedIdeal, J::MPolyQuoLocalizedIdeal)
  L = base_ring(I)
  L == base_ring(J) || error("ideals must be defined in the same ring")
  preI = Oscar.pre_image_ideal(I)
  preJ = Oscar.pre_image_ideal(J) 
  R = base_ring(L)
  K = intersect(preI, preJ)
  return L(K)
end

### Basic functionality
function Base.in(a::RingElem, I::MPolyQuoLocalizedIdeal)
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
function (W::MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MST})(I::MPolyIdeal{RET}) where {BRT, BRET, RT, RET, MST}
  return MPolyQuoLocalizedIdeal(W, W.(gens(I)))
end

### required constructors 
function ideal(
    W::MPolyQuoLocalizedRing, f
  )
  return MPolyQuoLocalizedIdeal(W, [W(f)])
end

function ideal(
    W::MPolyQuoLocalizedRing, gens::Vector
  )
  return MPolyQuoLocalizedIdeal(W, W.(gens))
end

function ideal(
    W::MPolyQuoLocalizedRing,
    I::MPolyLocalizedIdeal
  )
  return MPolyQuoLocalizedIdeal(W, W.(gens(I)))
end

function ideal(
    W::MPolyQuoLocalizedRing,
    I::MPolyIdeal
  )
  return MPolyQuoLocalizedIdeal(W, W.(gens(I)))
end

### Further constructors for quotient rings
function quo(
    L::MPolyQuoLocalizedRing,
    I::MPolyQuoLocalizedIdeal
  )
  base_ring(I) == L || error("ideal does not belong to the correct ring")
  W, _ = quo(localized_ring(L), modulus(L) + pre_image_ideal(I))
  return W, hom(L, W, gens(W), check=false)
end

function quo(A::MPolyQuo, I::MPolyQuoIdeal)
  base_ring(I) == A || error("ideal does not belong to the correct ring")
  R = base_ring(A)
  Q, _ = quo(R, modulus(A) + ideal(R, lift.(gens(I))))
  return Q, hom(A, Q, Q.(gens(R)))
end

function divides(a::MPolyQuoLocalizedRingElem, b::MPolyQuoLocalizedRingElem)
  W = parent(a)
  W == parent(b) || error("elements do not belong to the same ring")
  F = FreeMod(W, 1)
  A = MatrixSpace(W, 1, 1)([b])
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

function derivative(f::MPolyQuoLocalizedRingElem, i::Int)
  num = derivative(lifted_numerator(f), i)*lifted_denominator(f) - derivative(lifted_denominator(f), i)*lifted_numerator(f)
  den = lifted_denominator(f)^2
  g = gcd(num, den)
  return parent(f)(divexact(num, g), divexact(den, g), check=false)
end

function jacobi_matrix(f::MPolyQuoLocalizedRingElem)
  L = parent(f)
  n = nvars(base_ring(L))
  return matrix(L, n, 1, [derivative(f, i) for i=1:n])
end

function jacobi_matrix(g::Vector{<:MPolyQuoLocalizedRingElem})
  L = parent(g[1])
  n = nvars(base_ring(L))
  @assert all(x->parent(x) == L, g)
  return matrix(L, n, length(g), [derivative(x, i) for i=1:n for x = g])
end
