import AbstractAlgebra: Ring, RingElem, Generic.Frac
import Base: issubset

export MPolyQuoLocalizedRing
export parent, inverted_set, base_ring, quotient_ring, localized_ring, modulus, localized_modulus
export Localization

export MPolyQuoLocalizedRingElem
export numerator, denominator, parent, lift

export MPolyQuoLocalizedRingHom
export domain, codomain, images
export helper_ring, helper_images, minimal_denominators, helper_eta, helper_phi, common_denominator, helper_ideal


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

mutable struct MPolyQuoLocalizedRing{
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

  # fields for caching
  J::MPolyLocalizedIdeal{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}

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
    L.J = W(I)
    return L
  end
end

### required getter functions 
base_ring(L::MPolyQuoLocalizedRing) = L.R
inverted_set(L::MPolyQuoLocalizedRing) = L.S

### additional getter functions
modulus(L::MPolyQuoLocalizedRing) = L.I
localized_modulus(L::MPolyQuoLocalizedRing) = L.J
quotient_ring(L::MPolyQuoLocalizedRing) = L.Q
localized_ring(L::MPolyQuoLocalizedRing) = L.W

### printing
function Base.show(io::IO, L::MPolyQuoLocalizedRing)
  print(io, "Localization of $(quotient_ring(L)) at the multiplicative set $(inverted_set(L))")
end

function Base.show(
    io::IO, 
    L::MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MPolyUnits{BRT, BRET, RT, RET}}
  ) where {BRT, BRET, RT, RET}
  print(io, quotient_ring(L))
end

### additional constructors
function quo(
    W::MPolyLocalizedRing{BRT, BRET, RT, RET, MST},
    I::MPolyLocalizedIdeal{BRT, BRET, RT, RET, MST}
  ) where {BRT, BRET, RT, RET, MST}
  R = base_ring(W)
  S = inverted_set(W)
  lbpa = groebner_basis(I) # In particular, this saturates the ideal
  J = ideal(R, numerator.(oscar_gens(lbpa))) # the preimage of I in R
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
  U = MPolyUnits(R)
  W = Localization(U)
  return MPolyQuoLocalizedRing(R, I, U, Q, W)
end

function MPolyQuoLocalizedRing(Q::RT) where {RT<:MPolyQuo}
  R = base_ring(Q)
  I = modulus(Q)
  U = MPolyUnits(R)
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


########################################################################
# Elements of localizations of polynomial algebras                     #
########################################################################

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
      b::RingElemType
    ) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}

    S = inverted_set(L)
    parent(a) == parent(b) == base_ring(L) || error("elements do not belong to the correct ring")
    b in inverted_set(L) || error("the given denominator is not admissible for this localization")
    return new{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}(L, a, b)
  end
end

### required getter functions 
parent(a::MPolyQuoLocalizedRingElem) = a.L
numerator(a::MPolyQuoLocalizedRingElem) = quotient_ring(parent(a))(a.numerator) 
denominator(a::MPolyQuoLocalizedRingElem) = quotient_ring(parent(a))(a.denominator) 

### additional getter functions
quotient_ring(a::MPolyQuoLocalizedRingElem) = quotient_ring(parent(a))
localized_ring(a::MPolyQuoLocalizedRingElem) = localized_ring(parent(a))
lifted_numerator(a::MPolyQuoLocalizedRingElem) = a.numerator
lifted_denominator(a::MPolyQuoLocalizedRingElem) = a.denominator

### required conversions
(L::MPolyQuoLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType})(f::RingElemType) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyQuoLocalizedRingElem(L, f, one(f))
(L::MPolyQuoLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType})(a::RingElemType, b::RingElemType) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyQuoLocalizedRingElem(L, a, b)

### additional conversions
(L::MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MST})(f::Frac{RET}) where {BRT, BRET, RT, RET, MST} = L(numerator(f), denominator(f))

function (L::MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MST})(f::MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST}
  parent(f) == localized_ring(L) || error("the given element does not belong to the correct localization")
  lbpa = groebner_basis(localized_modulus(L))
  f = reduce(f, lbpa)
  return MPolyQuoLocalizedRingElem(L, numerator(f), denominator(f))
end

function (L::MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MST})(f::MPolyQuoElem{RET}) where {BRT, BRET, RT, RET, MST} 
  parent(f) == quotient_ring(L) || error("the given element does not belong to the correct ring") 
  return L(lift(f))
end

### additional functionality
lift(f::MPolyQuoLocalizedRingElem) = localized_ring(f)(lifted_numerator(f), lifted_denominator(f))


### arithmetic #########################################################
function +(a::T, b::T) where {T<:MPolyQuoLocalizedRingElem}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  if denominator(a) == denominator(b) 
    return reduce_fraction((parent(a))(numerator(a) + numerator(b), denominator(a)))
  end
  return reduce_fraction((parent(a))(numerator(a)*denominator(b) + numerator(b)*denominator(a), denominator(a)*denominator(b)))
end

# TODO: improve this method.
function addeq!(a::T, b::T) where {T<:MPolyQuoLocalizedRingElem}
  a = a+b
  return a
end

function -(a::T, b::T) where {T<:MPolyQuoLocalizedRingElem}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  if denominator(a) == denominator(b) 
    return reduce_fraction((parent(a))(numerator(a) - numerator(b), denominator(a)))
  end
  return reduce_fraction((parent(a))(numerator(a)*denominator(b) - numerator(b)*denominator(a), denominator(a)*denominator(b)))
end

function *(a::T, b::T) where {T<:MPolyQuoLocalizedRingElem}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  return reduce_fraction((parent(a))(numerator(a)*numerator(b), denominator(a)*denominator(b)))
end

function *(a::RET, b::MPolyQuoLocalizedRingElem{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET <: RingElem, MST}
  return reduce_fraction((parent(b))(a*numerator(b), denominator(b)))
end

function *(a::MPolyQuoLocalizedRingElem{BRT, BRET, RT, RET, MST}, b::RET) where {BRT, BRET, RT, RET <: RingElem, MST}
  return b*a
end

function *(a::BRET, b::MPolyQuoLocalizedRingElem{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET <: RingElem, MST}
  return reduce_fraction((parent(b))(base_ring(b)(a)*numerator(b), denominator(b)))
end

function *(a::MPolyQuoLocalizedRingElem{BRT, BRET, RT, RET, MST}, b::BRET) where {BRT, BRET, RT, RET <: RingElem, MST}
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
  return numerator(a)*denominator(b) == numerator(b)*denominator(a)
end

function ^(a::MPolyQuoLocalizedRingElem, i::fmpz)
  return parent(a)(lifted_numerator(a)^i, lifted_denominator(a)^i)
end

function ^(a::MPolyQuoLocalizedRingElem, i::Integer)
  return parent(a)(lifted_numerator(a)^i, lifted_denominator(a)^i)
end

function isone(a::MPolyQuoLocalizedRingElem) 
  a = reduce_fraction(a)
  return (numerator(a) == denominator(a))
end

function iszero(a::MPolyQuoLocalizedRingElem)
  a = reduce_fraction(a)
  return iszero(numerator(a))
end

### enhancement of the arithmetic
function reduce_fraction(f::MPolyQuoLocalizedRingElem)
  h = lift(f)
  h = reduce(h, groebner_basis(localized_modulus(parent(f))))
  g = gcd(numerator(h), denominator(h))
  h = parent(h)(divexact(numerator(h), g), divexact(denominator(h), g))
  return parent(f)(h)
end


### implementation of Oscar's general ring interface
one(W::MPolyQuoLocalizedRing) = W(one(base_ring(W)))
zero(W::MPolyQuoLocalizedRing)= W(zero(base_ring(W)))

elem_type(W::MPolyQuoLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyQuoLocalizedRingElem{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}
elem_type(T::Type{MPolyQuoLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyQuoLocalizedRingElem{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}

parent_type(W::MPolyQuoLocalizedRingElem{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyQuoLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}
parent_type(T::Type{MPolyQuoLocalizedRingElem{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyQuoLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}


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
#             η
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

mutable struct MPolyQuoLocalizedRingHom{
    BaseRingType, 
    BaseRingElemType, 
    RingType, 
    RingElemType, 
    DomainMultSetType, 
    CodomainMultSetType
  } <: AbsLocalizedRingHom{
    RingType, RingElemType, DomainMultSetType, CodomainMultSetType
  }
  domain::MPolyQuoLocalizedRing
  codomain::MPolyQuoLocalizedRing
  images::Vector{MPolyLocalizedRingElem}

  # variables for caching
  helper_ring::RingType
  helper_images::Vector{RingElemType}
  minimal_denominators::Vector{RingElemType}
  eta::AlgHom{BaseRingElemType}
  phi::AlgHom{BaseRingElemType}

  function MPolyQuoLocalizedRingHom(
      L::MPolyQuoLocalizedRing{BRT, BRET, RT, RET, DMST}, 
      M::MPolyQuoLocalizedRing{BRT, BRET, RT, RET, CMST}, 
      a::Vector{MPolyLocalizedRingElem{BRT, BRET, RT, RET, CMST}}
    ) where {BRT, BRET, RT, RET, CMST, DMST}
    R = base_ring(L)
    S = base_ring(M)
    k = coefficient_ring(R) 
    k == coefficient_ring(S) || error("the two polynomial rings are not defined over the same coefficient ring")
    ngens(R) == length(a) || error("the number of images does not coincide with the number of variables")
    parent_check = true
    for x in a
      parent_check = parent_check && parent(x) == localized_ring(M)
    end
    parent_check || error("the images of the variables are not elements of the codomain")
    # Check whether this homomorphism is well defined
    # TODO: Implement that!
    return new{typeof(k), elem_type(k), typeof(R), elem_type(R), typeof(inverted_set(L)), typeof(inverted_set(M))}(L, M, a)
  end
end

### required getter functions
domain(f::MPolyQuoLocalizedRingHom) = f.domain
codomain(f::MPolyQuoLocalizedRingHom) = f.codomain
images(f::MPolyQuoLocalizedRingHom) = f.images

### required functionality
function (f::MPolyQuoLocalizedRingHom{BRT, BRET, RT, RET, DMST, CMST})(
    p::MPolyQuoLocalizedRingElem{BRT, BRET, RT, RET, DMST}
  ) where {BRT, BRET, RT, RET, DMST, CMST}
  parent(p) == domain(f) || error("the given element does not belong to the domain of the map")
  return codomain(f)(evaluate(lifted_numerator(p), images(f))//evaluate(lifted_denominator(p), images(f)))
end

### additional functionality 
function (f::MPolyQuoLocalizedRingHom{BRT, BRET, RT, RET, DMST, CMST})(
    p::MPolyLocalizedRingElem{BRT, BRET, RT, RET, DMST}
  ) where {BRT, BRET, RT, RET, DMST, CMST}
  parent(p) == localized_ring(domain(f)) || error("the given element does not belong to the domain of the map")
  return codomain(f)(evaluate(lifted_numerator(p), images(f))//evaluate(lifted_denominator(p), images(f)))
end

function (f::MPolyQuoLocalizedRingHom{BRT, BRET, RT, RET, DMST, CMST})(
    p::MPolyQuoElem{RET}
  ) where {BRT, BRET, RT, RET, DMST, CMST}
  parent(p) == quotient_ring(domain(f)) || error("the given element does not belong to the domain of the map")
  return codomain(f)(evaluate(lift(p), images(f)))
end

### overwriting of the generic method
function (f::MPolyQuoLocalizedRingHom{BRT, BRET, RT, RET, DMST, CMST})(
    p::RET
  ) where {BRT, BRET, RT, RET, DMST, CMST}
  parent(p) == base_ring(domain(f)) || error("the given element does not belong to the domain of the map")
  return codomain(f)(evaluate(p, images(f)))
end

### provide an extra method for elements of the base ring
function (f::MPolyQuoLocalizedRingHom{BRT, BRET, RT, RET, DMST, CMST})(p::BRET) where {BRT, BRET, RT, RET, DMST, CMST}
  parent(p) == coefficient_ring(base_ring(domain(f))) || error("the given element does not belong to the domain of the map")
  return codomain(f)(p)
end

### remove the ambiguity of methods in case the base ring is ZZ
function (f::MPolyQuoLocalizedRingHom)(p::fmpz) 
  return codomain(f)(p)
end

### implementing the Oscar map interface
identity_map(W::T) where {T<:MPolyQuoLocalizedRing} = MPolyQuoLocalizedRingHom(W, W, W.(gens(base_ring(W))))
function compose(
    f::MPolyQuoLocalizedRingHom{BRT, BRET, RT, RET, MST1, MST2}, 
    g::MPolyQuoLocalizedRingHom{BRT, BRET, RT, RET, MST2, MST3}
  ) where {BRT, BRET, RT, RET, MST1, MST2, MST3}
  codomain(f) == domain(g) || error("maps are not compatible")
  return MPolyQuoLocalizedRingHom(domain(f), codomain(g), g.(images(f)))
end

### helper_ring
# Sets up the ring S[c⁻¹] from the Lemma.
function helper_ring(f::MPolyQuoLocalizedRingHom{BRT, BRET, RT, RET, DMST, CMST}) where {BRT, BRET, RT, RET, DMST, CMST}
  if isdefined(f, :helper_ring)
    return f.helper_ring
  end
  f.minimal_denominators = Vector{RET}()
  R = base_ring(domain(f))
  S = base_ring(codomain(f))
  p = one(S)

  for d in [denominator(y) for y in images(f)]
    g = gcd(d, p)
    d_min = divexact(d, g)
    push!(f.minimal_denominators, d)
    p = p*d_min
  end
 
  help_ring, help_phi, theta = _add_variables(S, ["θ"])
  f.helper_ring = help_ring
  f.phi = help_phi
  c_inv = theta[1]
  f.helper_images = [f.phi(numerator(y))*c_inv*f.phi(divexact(p, denominator(y))) for y in images(f)]
  f.eta = AlgebraHomomorphism(R, help_ring, f.helper_images)
  return f.helper_ring
end

function helper_images(
    f::MPolyQuoLocalizedRingHom{BRT, BRET, RT, RET, DMST, CMST}
  ) where {BRT, BRET, RT, RET, DMST, CMST} 
  if !isdefined(f, :helper_images) 
    helper_ring(f)
  end
  return f.helper_images
end

function minimal_denominators(
    f::MPolyQuoLocalizedRingHom{BRT, BRET, RT, RET, DMST, CMST}
  ) where {BRT, BRET, RT, RET, DMST, CMST} 
  if !isdefined(f, :minimal_denominators) 
    helper_ring(f)
  end
  return f.minimal_denominators
end

function helper_eta(
    f::MPolyQuoLocalizedRingHom{BRT, BRET, RT, RET, DMST, CMST}
  ) where {BRT, BRET, RT, RET, DMST, CMST} 
  if !isdefined(f, :eta) 
    helper_ring(f)
  end
  return f.eta
end

function helper_phi(
    f::MPolyQuoLocalizedRingHom{BRT, BRET, RT, RET, DMST, CMST}
  ) where {BRT, BRET, RT, RET, DMST, CMST} 
  if !isdefined(f, :phi) 
    helper_ring(f)
  end
  return f.phi
end

function common_denominator(
    f::MPolyQuoLocalizedRingHom{BRT, BRET, RT, RET, DMST, CMST}
  ) where {BRT, BRET, RT, RET, DMST, CMST} 
  if !isdefined(f, :minimal_denominators) 
    helper_ring(f)
  end
  return (length(f.minimal_denominators) == 0 ? one(base_ring(codomain(f))) : prod(f.minimal_denominators))
end

function helper_ideal(
    f::MPolyQuoLocalizedRingHom{BRT, BRET, RT, RET, DMST, CMST}
  ) where {BRT, BRET, RT, RET, DMST, CMST} 
  if !isdefined(f, :helper_ring) 
    helper_ring(f)
  end
  Sc = helper_ring(f)
  return ideal(Sc, one(Sc)-last(gens(Sc))*helper_phi(f)(common_denominator(f)))
end



########################################################################
# Functionality for maps and ideals                                    #
########################################################################
# 
# The methods have to be adapted to the type of localization in the 
# target. It needs to be assured that all components which are invisible
# in the localization, are indeed discarded. 

### adds the variables with names specified in v to the polynomial 
# ring R and returns a triple consisting of the new ring, the embedding 
# of the original one, and a list of the new variables. 
function _add_variables(R::RingType, v::Vector{String}) where {RingType<:MPolyRing}
  ext_R, _ = PolynomialRing(coefficient_ring(R), vcat(symbols(R), Symbol.(v)))
  n = length(gens(R))
  phi = AlgebraHomomorphism(R, ext_R, gens(ext_R)[1:n])
  return ext_R, phi, gens(ext_R)[(length(gens(R))+1):length(gens(ext_R))]
end

function preimage(
    f::MPolyQuoLocalizedRingHom{BRT, BRET, RT, RET, DMST, CMST},
    I::MPolyLocalizedIdeal{BRT, BRET, RT, RET, CMST}
  ) where {BRT<:Ring, BRET<:RingElement, RT<:MPolyRing, RET<:MPolyElem, 
    DMST<:AbsMultSet{RT, RET}, CMST<:MPolyPowersOfElement{BRT, BRET, RT, RET}
  }
  base_ring(I) == localized_ring(codomain(f)) || error("the ideal does not belong to the codomain of the map")
  R = base_ring(domain(f))
  S = base_ring(codomain(f))
  Sc = helper_ring(f)
#  help_ring, phi, new_vars = _add_variables(S, ["θ"])
#  t = new_vars[1]
#  common_denom = one(R)
#  for y in [denominator(g) for g in images(f)]
#    g = gcd(common_denom, y)
#    common_denom = common_denom*divexact(y, g)
#  end
#  help_images = [numerator(y)*divexact(common_denom, denominator(y)) for y in images(f)]
#  help_hom = AlgebraHomomorphism(S, help_ring, [phi(x)*t for x in help_images])
  lbpa = groebner_basis(I) # saturation takes place in this computation
  #J = ideal(help_ring, [phi(g) for g in numerator.(oscar_gens(lbpa))]) + ideal(help_ring, one(help_ring)-t*phi(common_denom))
  J = ideal(Sc, [helper_phi(f)(g) for g in numerator.(oscar_gens(lbpa))]) + helper_ideal(f)
  return localized_ring(domain(f))(preimage(helper_eta(f), J))
end

function preimage(
    f::MPolyQuoLocalizedRingHom{BRT, BRET, RT, RET, DMST, CMST},
    I::MPolyLocalizedIdeal{BRT, BRET, RT, RET, CMST}
  ) where {BRT<:Ring, BRET<:RingElement, RT<:MPolyRing, RET<:MPolyElem, 
    DMST<:AbsMultSet{RT, RET}, CMST<:MPolyComplementOfKPointIdeal{BRT, BRET, RT, RET}
  }
  base_ring(I) == localized_ring(codomain(f)) || error("the ideal does not belong to the codomain of the map")
  R = base_ring(domain(f))
  S = base_ring(codomain(f))
  Sc = helper_ring(f)

  # throw away all components of the ideal I which do not touch the origin
  # TODO: This uses primary decomposition for now. Most probably, the 
  # computations can be carried out more efficiently.
  lbpa = LocalizedBiPolyArray(I)
  V = base_ring(I)
  decomp = [LocalizedBiPolyArray(V, K[1], shift=shift(lbpa)) for K in Singular.LibPrimdec.primdecGTZ(singular_ring(lbpa), singular_gens(lbpa))]
  if length(decomp) == 0
    return ideal(V, zero(V))
  end
  relevant_comp = decomp[1]
  for i in (2:length(decomp))
    relevant_comp = Singular.intersection(singular_gens(relevant_comp), singular_gens(decomp[i]))
  end

  # compute the preimage of the remaining ideal
  J = ideal(Sc, [helper_phi(f)(g) for g in numerator.(oscar_gens(relevant_comp))]) + helper_ideal(f)
  return localized_ring(domain(f))(preimage(helper_phi(f), J))
end
