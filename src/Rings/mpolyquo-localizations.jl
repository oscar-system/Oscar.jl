import AbstractAlgebra: Ring, RingElem, Generic.Frac
import Base: issubset

export MPolyQuoLocalizedRing
export parent, inverted_set, base_ring, quotient_ring, localized_ring, modulus, localized_modulus, gens
export Localization

export MPolyQuoLocalizedRingElem
export numerator, denominator, parent, lift, isunit, inv, convert, lifted_numerator, lifted_denominator, fraction

export MPolyQuoLocalizedRingHom
export domain, codomain, images, morphism_type, domain_type, codomain_type, restricted_map_type, ideal_type
export helper_ring, helper_images, minimal_denominators, helper_eta, helper_kappa, common_denominator, helper_ideal

export is_isomorphism, inverse

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

ideal_type(::Type{MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = MPolyLocalizedIdeal{BRT, BRET, RT, RET, MST}
ideal_type(W::MPolyQuoLocalizedRing) = ideal_type(typeof(W))





### required getter functions 
base_ring(L::MPolyQuoLocalizedRing) = L.R
inverted_set(L::MPolyQuoLocalizedRing) = L.S

### additional getter functions
@Markdown.doc """
    modulus(L::MPolyQuoLocalizedRing)

For ``L = (𝕜[x₁,…,xₙ]/I)[S⁻¹]`` this returns ``I``.
"""
modulus(L::MPolyQuoLocalizedRing) = L.I

@Markdown.doc """
    localized_modulus(L::MPolyQuoLocalizedRing)

For ``L = (𝕜[x₁,…,xₙ]/I)[S⁻¹]`` this returns ``IS⁻¹``.
"""
function localized_modulus(L::MPolyQuoLocalizedRing) 
  if !has_attribute(L, :localized_modulus)
    set_attribute!(L, :localized_modulus, localized_ring(L)(modulus(L)))
  end
  return get_attribute(L, :localized_modulus)::ideal_type(L)
end

@Markdown.doc """
    quotient_ring(L::MPolyQuoLocalizedRing)

For ``L = (𝕜[x₁,…,xₙ]/I)[S⁻¹]`` this returns ``𝕜[x₁,…,xₙ]/I``.
"""
quotient_ring(L::MPolyQuoLocalizedRing) = L.Q

@Markdown.doc """
    localized_ring(L::MPolyQuoLocalizedRing)

For ``L = (𝕜[x₁,…,xₙ]/I)[S⁻¹]`` this returns ``𝕜[x₁,…,xₙ][S⁻¹]``.
"""
localized_ring(L::MPolyQuoLocalizedRing) = L.W

@Markdown.doc """
    gens(L::MPolyQuoLocalizedRing)

For ``L = (𝕜[x₁,…,xₙ]/I)[S⁻¹]`` this returns the vector ``[x₁//1,…,xₙ//1]∈ Lⁿ``.
"""
gens(L::MPolyQuoLocalizedRing) = L.(gens(base_ring(L)))

### printing
function Base.show(io::IO, L::MPolyQuoLocalizedRing)
  print(io, "Localization of $(quotient_ring(L)) at the multiplicative set $(inverted_set(L))")
end

### additional constructors
function quo(
    W::MPolyLocalizedRing{BRT, BRET, RT, RET, MST},
    I::MPolyLocalizedIdeal{BRT, BRET, RT, RET, MST}
  ) where {BRT, BRET, RT, RET, MST}
  R = base_ring(W)
  S = inverted_set(W)
  #lbpa = groebner_basis(I) # In particular, this saturates the ideal
  #J = ideal(R, numerator.(oscar_gens(lbpa))) # the preimage of I in R
  J = ideal(R, numerator.(gens(I)))
  return MPolyQuoLocalizedRing(R, J, S, quo(R, J)[1], W)
end

function quo(
    L::MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MST},
    I::MPolyLocalizedIdeal{BRT, BRET, RT, RET, MST}
  ) where {BRT, BRET, RT, RET, MST}
  R = base_ring(L)
  S = inverted_set(L)
  lbpa = groebner_basis(I) # In particular, this saturates the ideal
  J = ideal(R, numerator.(oscar_gens(lbpa))) # the preimage of I in R
  return MPolyQuoLocalizedRing(R, J, S, quo(R, J)[1], localized_ring(L))
end

function quo(
    L::MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MST},
    J::MPolyIdeal{RET}
  ) where {BRT, BRET, RT, RET, MST}
  R = base_ring(L)
  S = inverted_set(L)
  W = localized_ring(L) 
  J = J + modulus(L)
  return MPolyQuoLocalizedRing(R, J, S, quo(R, J)[1], W)
end

function Localization(Q::MPolyQuo{RET}, S::MultSetType) where {RET <: RingElem, MultSetType <: AbsMultSet}
  return MPolyQuoLocalizedRing(base_ring(Q), modulus(Q), S, Q, Localization(S))
end

function Localization(
    L::MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MST}, 
    S::AbsMPolyMultSet{BRT, BRET, RT, RET}
  ) where {BRT, BRET, RT, RET, MST}
  ambient_ring(S) == base_ring(L) || error("multiplicative set does not belong to the correct ring")
  issubset(S, inverted_set(L)) && return L
  U = inverted_set(L)*S
  return MPolyQuoLocalizedRing(base_ring(L), modulus(L), U, quotient_ring(L), Localization(U))
end

function MPolyQuoLocalizedRing(R::RT, I::Ideal{RET}, T::MultSetType) where {RT<:MPolyRing, RET<:MPolyElem, MultSetType<:AbsMultSet} 
  return MPolyQuoLocalizedRing(R, I, T, quo(R, I)[1], Localization(T))
end

function MPolyQuoLocalizedRing(R::RT) where {RT<:MPolyRing} 
  I = ideal(R, zero(R))
  Q, _ = quo(R, I)
  U = units_of(R)
  W = Localization(U)
  return MPolyQuoLocalizedRing(R, I, U, Q, W)
end

function MPolyQuoLocalizedRing(Q::RT) where {RT<:MPolyQuo}
  R = base_ring(Q)
  I = modulus(Q)
  U = units_of(R)
  W = Localization(U)
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

  function MPolyQuoLocalizedRingElem(
      L::MPolyQuoLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}, 
      a::RingElemType,
      b::RingElemType;
      check::Bool=true
    ) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}

    S = inverted_set(L)
    R = base_ring(L)
    parent(a) == parent(b) == R || error("elements do not belong to the correct ring")
    check && (b in S || error("denominator is not admissible"))
    return new{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}(L, a, b)
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
numerator(a::MPolyQuoLocalizedRingElem) = quotient_ring(parent(a))(a.numerator) 
denominator(a::MPolyQuoLocalizedRingElem) = quotient_ring(parent(a))(a.denominator) 

### additional getter functions
quotient_ring(a::MPolyQuoLocalizedRingElem) = quotient_ring(parent(a))
localized_ring(a::MPolyQuoLocalizedRingElem) = localized_ring(parent(a))
base_ring(a::MPolyQuoLocalizedRingElem) = base_ring(parent(a))

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
  return parent(f)(f, check=false)
end

### required conversions
(L::MPolyQuoLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType})(f::RingElemType) where {BaseRingType, BaseRingElemType, RingType, RingElemType<:RingElem, MultSetType} = MPolyQuoLocalizedRingElem(L, f, one(f), check=false)

function (L::MPolyQuoLocalizedRing{
                                   BaseRingType, 
                                   BaseRingElemType, 
                                   RingType, 
                                   RingElemType, 
                                   MultSetType
                                  })(
                                     a::RingElemType, 
                                     b::RingElemType;
                                     check::Bool=true
                                    ) where {
                                             BaseRingType, 
                                             BaseRingElemType, 
                                             RingType, 
                                             RingElemType, 
                                             MultSetType
                                            } 
  check || return MPolyQuoLocalizedRingElem(L, a, b, check=false)
  b in inverted_set(L) || return convert(L, a//b)
  return MPolyQuoLocalizedRingElem(L, a, b, check=false)
end

### additional conversions
function (L::MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MST})(f::Frac{RET}; check::Bool=true) where {BRT, BRET, RT, RET, MST}
  R = base_ring(L)
  return L(R(numerator(f)), R(denominator(f)), check=check)
end

(L::MPolyQuoLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType})(a::T, b::T; check::Bool=true) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType, T<:MPolyQuoElem{RingElemType}} = L(lift(a), lift(b), check=check)

function (L::MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MST})(f::MPolyQuoLocalizedRingElem{BRT, BRET, RT, RET, MST}; check::Bool=true) where {BRT, BRET, RT, RET, MST}
  parent(f) === L && return f
  return L(lifted_numerator(f), lifted_denominator(f), check=check)
end

function (L::MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MST})(f::MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}; check::Bool=true) where {BRT, BRET, RT, RET, MST}
  parent(f) === localized_ring(L) && return L(numerator(f), denominator(f), check=false)
  return L(numerator(f), denominator(f), check=check)
end

function (L::MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MST})(f::MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}; check::Bool=true) where {BRT, BRET, RT, RET, MST<:MPolyComplementOfKPointIdeal}
  return L(numerator(f), denominator(f), check=check)
end

function (L::MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MST})(f::MPolyQuoElem{RET}) where {BRT, BRET, RT, RET, MST} 
  parent(f) == quotient_ring(L) || error("the given element does not belong to the correct ring") 
  return L(lift(f))
end

### additional functionality
@Markdown.doc """
    lift(f::MPolyQuoLocalizedRingElem)

For ``f = A//B ∈ (𝕜[x₁,…,xₙ]/I)[S⁻¹]`` this returns a representative 
``a//b ∈  𝕜[x₁,…,xₙ][S⁻¹]`` of the fraction. 
"""
lift(f::MPolyQuoLocalizedRingElem) = localized_ring(f)(lifted_numerator(f), lifted_denominator(f))

function isunit(f::MPolyQuoLocalizedRingElem) 
  lifted_numerator(f) in inverted_set(parent(f)) && return true
  return one(localized_ring(parent(f))) in localized_modulus(parent(f)) + ideal(localized_ring(parent(f)), lift(f))
end

function isunit(L::MPolyQuoLocalizedRing, f::MPolyLocalizedRingElem) 
  parent(f) == localized_ring(L) || error("element does not belong to the correct ring")
  numerator(f) in inverted_set(L) && return true
  one(localized_ring(L)) in localized_modulus(L) + ideal(localized_ring(L), f)
end

function isunit(L::MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MST}, f::RET) where {BRT, BRET, RT, RET, MST}
  parent(f) == base_ring(L) || error("element does not belong to the correct ring")
  f in inverted_set(L) && return true
  return one(localized_ring(L)) in localized_modulus(L) + ideal(localized_ring(L), localized_ring(L)(f))
end

function isunit(L::MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MST}, f::MPolyQuoElem{RET}) where {BRT, BRET, RT, RET, MST}
  parent(f) == quotient_ring(L) || error("element does not belong to the correct ring")
  lift(f) in inverted_set(L) && return true
  one(localized_ring(L)) in localized_modulus(L) + ideal(localized_ring(L), localized_ring(L)(f))
end

# WARNING: This routine runs forever if f is not a unit in L. 
# So this needs to be checked first!
function inv(L::MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}}, 
    f::MPolyQuoElem{RET}) where {BRT, BRET, RT, RET}
  Q = quotient_ring(L)
  parent(f) == quotient_ring(L) || error("element does not belong to the correct ring")
  W = localized_ring(L)
  R = base_ring(L)
  I = saturated_ideal(localized_modulus(L))
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
  Q = quotient_ring(L)
  parent(a) == base_ring(L) || error("element does not belong to the correct ring")
  W = localized_ring(L)
  R = base_ring(L)
  I = saturated_ideal(localized_modulus(L))
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
    return reduce_fraction((parent(a))(lifted_numerator(a) + lifted_numerator(b), lifted_denominator(a), check=false))
  end
  return reduce_fraction((parent(a))(lifted_numerator(a)*lifted_denominator(b) + lifted_numerator(b)*lifted_denominator(a), lifted_denominator(a)*lifted_denominator(b), check=false))
end

# TODO: improve this method.
function addeq!(a::T, b::T) where {T<:MPolyQuoLocalizedRingElem}
  a = a+b
  return a
end

function -(a::T, b::T) where {T<:MPolyQuoLocalizedRingElem}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  if lifted_denominator(a) == lifted_denominator(b) 
    return reduce_fraction((parent(a))(lifted_numerator(a) - lifted_numerator(b), lifted_denominator(a), check=false))
  end
  return reduce_fraction((parent(a))(lifted_numerator(a)*lifted_denominator(b) - lifted_numerator(b)*lifted_denominator(a), lifted_denominator(a)*lifted_denominator(b), check=false))
end

function *(a::T, b::T) where {T<:MPolyQuoLocalizedRingElem}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  return reduce_fraction((parent(a))(lifted_numerator(a)*lifted_numerator(b), lifted_denominator(a)*lifted_denominator(b), check=false))
end

function *(a::RET, b::MPolyQuoLocalizedRingElem{BRT, BRET, RT, RET, MST}) where {BRT<:Ring, BRET<:RingElem, RT<:Ring, RET <: RingElem, MST}
  return reduce_fraction((parent(b))(a*lifted_numerator(b), lifted_denominator(b), check=false))
end

function *(a::MPolyQuoLocalizedRingElem{BRT, BRET, RT, RET, MST}, b::RET) where {BRT<:Ring, BRET<:RingElem, RT<:Ring, RET <: RingElem, MST}
  return b*a
end

function *(a::BRET, b::MPolyQuoLocalizedRingElem{BRT, BRET, RT, RET, MST}) where {BRT<:Ring, BRET<:RingElem, RT<:Ring, RET <: RingElem, MST}
  return reduce_fraction((parent(b))(a*lifted_numerator(b), lifted_denominator(b), check=false))
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
  return lifted_numerator(a)*lifted_denominator(b) - lifted_numerator(b)*lifted_denominator(a) in localized_modulus(parent(a))
end

function ^(a::MPolyQuoLocalizedRingElem, i::fmpz)
  return parent(a)(lifted_numerator(a)^i, lifted_denominator(a)^i, check=false)
end

function ^(a::MPolyQuoLocalizedRingElem, i::Integer)
  return parent(a)(lifted_numerator(a)^i, lifted_denominator(a)^i, check=false)
end

function isone(a::MPolyQuoLocalizedRingElem) 
  return lifted_numerator(a) - lifted_denominator(a) in localized_modulus(parent(a))
end

function iszero(a::MPolyQuoLocalizedRingElem)
  iszero(lifted_numerator(a)) && return true
  return lifted_numerator(a) in localized_modulus(parent(a))
end

### enhancement of the arithmetic
function reduce_fraction(f::MPolyQuoLocalizedRingElem{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST<:MPolyPowersOfElement}
  return f
  h = lift(f)
  h = reduce(h, groebner_basis(localized_modulus(parent(f))))
  g = gcd(numerator(h), denominator(h))
  h = parent(h)(divexact(numerator(h), g), divexact(denominator(h), g), check=false)
  return parent(f)(h, check=false)
end

# for local orderings, reduction does not give the correct result.
function reduce_fraction(f::MPolyQuoLocalizedRingElem{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST<:MPolyComplementOfKPointIdeal}
  return f
  h = lift(f)
  g = gcd(numerator(h), denominator(h))
  h = parent(h)(divexact(numerator(h), g), divexact(denominator(h), g), check=false)
  return parent(f)(h, check=false)
end

### implementation of Oscar's general ring interface
one(W::MPolyQuoLocalizedRing) = W(one(base_ring(W)))
zero(W::MPolyQuoLocalizedRing)= W(zero(base_ring(W)))

elem_type(W::MPolyQuoLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyQuoLocalizedRingElem{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}
elem_type(T::Type{MPolyQuoLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyQuoLocalizedRingElem{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}

parent_type(W::MPolyQuoLocalizedRingElem{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyQuoLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}
parent_type(T::Type{MPolyQuoLocalizedRingElem{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyQuoLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}


### ideal constructors
# Note that by convention an ideal J in a localized algebra 
# L = (𝕜[x₁,…,xₙ]/I)[S⁻¹] is an ideal in 𝕜[x₁,…,xₙ][S⁻¹] 
# containing IS⁻¹. We provide the constructors here.

function ideal(L::MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MST}, 
    g::Vector{T}
  ) where {BRT, BRET, RT, RET, MST, T<:RingElement}
  gconv = L.(g)
  W = localized_ring(L)
  return ideal(W, lift.(gconv)) + localized_modulus(L)
end

ideal(L::MPolyQuoLocalizedRing, g::T) where {T<:RingElement} = ideal(L, [g])
ideal(L::MPolyQuoLocalizedRing, g::T) where {T<:MPolyQuoLocalizedRingElem} = ideal(L, [g])

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
    f::MPolyQuoLocalizedRingElem{BRT, BRET, RT, RET, MPolyComplementOfKPointIdeal{BRT, BRET, RT, RET}},
    g::Vector{MPolyQuoLocalizedRingElem{BRT, BRET, RT, RET, MPolyComplementOfKPointIdeal{BRT, BRET, RT, RET}}}
  ) where {BRT, BRET, RT, RET}
  n = length(g)
  L = parent(f)
  W = localized_ring(L)
  for a in g 
    parent(a) == L || error("elements do not belong to the same ring")
  end
  (d, a) = bring_to_common_denominator(vcat([f], g))
  h = [a[i+1]*lifted_numerator(g[i]) for i in 1:n]
  lbpa = LocalizedBiPolyArray(W.(h))
  p = a[1]*lifted_numerator(f)
  p_sing = to_singular_side(lbpa, p)
  S = singular_ring(lbpa)
  I = modulus(L)
  SI = to_singular_side(lbpa, gens(I))
  
  M, N, U = Singular.lift(
			  Singular.Module(S, vcat([Singular.vector(S, g) for g in gens(singular_gens(lbpa))], [Singular.vector(S, g) for g in SI])...),
			  Singular.Module(S, Singular.vector(S, p_sing)),
			  false, false, false)
  A = Singular.Matrix(M)
  iszero(N) || error("the first argument is not contained in the span of the second")
  u = L(one(base_ring(W)),to_oscar_side(lbpa, U[1,1]))
  lambda = [L(to_oscar_side(lbpa, A[i, 1]))*u for i in 1:n]
  return lambda
end

function write_as_linear_combination(
    f::MPolyQuoLocalizedRingElem{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}},
    g::Vector{MPolyQuoLocalizedRingElem{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}}}
  ) where {BRT, BRET, RT, RET}
  n = length(g)
  L = parent(f)
  for a in g 
    parent(a) == L || error("elements do not belong to the same ring")
  end
  (d, a) = bring_to_common_denominator(vcat([f], g))
  hg = [a[i+1]*lifted_numerator(g[i]) for i in 1:n]
  hf = lifted_numerator(f)*a[1]
  A, I, q, phi, theta = as_affine_algebra(L)
  SA, _ = Singular.PolynomialRing(Oscar.singular_ring(base_ring(A)), 
				  String.(symbols(A)),  
				  ordering=Singular.ordering_dp(1)
				  *Singular.ordering_dp(nvars(A)-1))
  Shg_ext = Singular.Ideal(SA, SA.(vcat(phi.(hg), gens(I))))
  M, N, U = Singular.lift(
                          Singular.Module(SA, [Singular.vector(SA, g) for g in gens(Shg_ext)]...),
			  Singular.Module(SA, Singular.vector(SA, SA(phi(hf)))),
			  false, false, false)
  iszero(N) || error("the first argument is not contained in the span of the second")
  W = localized_ring(L)
  evaluation_list = vcat([W(one(base_ring(L)), q)], gens(W))
  l = [Singular.Matrix(M)[i, 1] for i in 1:n]
  lambda = L.([evaluate(A(a), evaluation_list) for a in l])
  return lambda
end

write_as_linear_combination(f::MPolyQuoLocalizedRingElem, g::Vector) = write_as_linear_combination(f, parent(f).(g))

function write_as_linear_combination(f::T, g::Vector{T}) where {T<:MPolyQuoElem}
  Q = parent(f)
  R = base_ring(Q)
  I = modulus(Q)
  for b in g
    parent(b) == Q || error("elements do not belong to the same ring")
  end
  SR, _ = Singular.PolynomialRing(Oscar.singular_ring(base_ring(R)), 
				  String.(symbols(R)),  
				  ordering=Singular.ordering_dp(ngens(R)))
  Sg_ext = SR.(vcat(lift.(g), gens(I)))
  Sf = SR(lift(f))
  M, N, U = Singular.lift(
                          Singular.Module(SR, [Singular.vector(SR, a) for a in Sg_ext]...),
			  Singular.Module(SR, Singular.vector(SR, Sf)),
			  false, false, false)
  iszero(N) || error("the first argument is not contained in the span of the second")
  return Q.(R.([Singular.Matrix(M)[i, 1] for i in 1:length(g)]))
end


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
        isunit(S(res(f))) || error("map is not well defined")
      end
      for g in gens(modulus(L))
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

hom(L::MPolyQuoLocalizedRing, S::Ring, a::Vector{T}) where {T<:RingElem} = MPolyQuoLocalizedRingHom(L, S, a)

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
  if codomain(restricted_map(f)) === domain(g)
    return MPolyQuoLocalizedRingHom(domain(f), codomain(g), compose(restricted_map(f), g))
  elseif codomain(restricted_map(f)) === base_ring(domain(g)) 
    h = hom(base_ring(domain(g)), domain(g), domain(g).(gens(base_ring(domain(g)))))
    return MPolyQuoLocalizedRingHom(domain(f), codomain(g), compose(compose(restricted_map(f), h), g))
  end
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
  ) where {BRT, BRET, RT, RET}
  R = base_ring(L)
  A, phi, t = _add_variables_first(R, [inverse_name])
  theta = t[1]
  f = prod(denominators(inverted_set(L)))
  I = ideal(A, [phi(g) for g in gens(modulus(L))]) + ideal(A, [one(A)-theta*phi(f)])
  return A, I, f, phi, theta
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
  singC, _ = Singular.PolynomialRing(Oscar.singular_ring(base_ring(C)), 
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



function preimage(
    f::MPolyQuoLocalizedRingHom{
                                <:MPolyQuoLocalizedRing{<:Any, <:Any, <:Any, <:Any, 
                                                        <:MPolyPowersOfElement
                                                       },
                                <:MPolyQuoLocalizedRing{<:Any, <:Any, <:Any, <:Any, 
                                                        <:MPolyPowersOfElement
                                                       }
                               },
    I::MPolyLocalizedIdeal
  )
  base_ring(I) == localized_ring(codomain(f)) || error("the ideal does not belong to the codomain of the map")
  R = base_ring(domain(f))
  S = base_ring(codomain(f))
  Sc = helper_ring(f)
  lbpa = groebner_basis(I) # saturation takes place in this computation
  J = ideal(Sc, [helper_kappa(f)(g) for g in numerator.(oscar_gens(lbpa))]) + helper_ideal(f)
  return localized_ring(domain(f))(preimage(helper_eta(f), J))
end

function preimage(
    f::MPolyQuoLocalizedRingHom{
                                <:MPolyQuoLocalizedRing{<:Any, <:Any, <:Any, <:Any, 
                                                        <:MPolyPowersOfElement
                                                       },
                                <:MPolyQuoLocalizedRing{<:Any, <:Any, <:Any, <:Any, 
                                                        <:MPolyComplementOfKPointIdeal
                                                       }
                               },
    I::MPolyLocalizedIdeal
  )
  base_ring(I) == localized_ring(codomain(f)) || error("the ideal does not belong to the codomain of the map")
  J = ideal(helper_ring(f), helper_kappa(f).(gens(saturated_ideal(I)))) + helper_ideal(f)
  return localized_ring(domain(f))(preimage(helper_eta(f), J))
end
