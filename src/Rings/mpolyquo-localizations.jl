import AbstractAlgebra: Ring, RingElem, Generic.Frac

export MPolyQuoLocalizedRing
export parent, inverted_set, base_ring, loc_poly_ring, modulus
export Localization

export MPolyQuoLocalizedRingElem
export numerator, denominator, parent, lift


########################################################################
# Localizations of polynomial algebras                                 #
########################################################################
# 
# Let R = 𝕜[x₁,…,xₘ] be a polynomial ring, I ⊂ R some ideal 
# and Q = R/I its quotient. Then Q is naturally an R-module 
# and localization of Q as a ring coincides with localization 
# as an R-module in the sense that for every multiplicative 
# set S ⊂ R there is a commutative diagram 
#
#         R    →  Q = R/I
#         ↓       ↓
#   W = R[S⁻¹] → Q[S⁻¹].
#
# Observe that, moreover, for every multiplicative set 
# T ⊂ Q the preimage S of T in R is also a multiplicative set. 
#
# We may therefore treat localizations of polynomial algebras 
# as localizations of modules over free polynomial rings and 
# exclusively use the types of multiplicative sets which are 
# available for the latter.
#
# Note the following differences compared to the standard usage 
# of the localization interface:
#
#  * The `base_ring` returns neither Q, nor W, but R.
#  * The `BaseRingType` is the type of R and similar for 
#    the other ring-based type parameters.
#
# This is to make the data structure most accessible for 
# the computational backends.
#
#  * The type returned by `numerator` and `denominator` 
#    on an element of type `MPolyQuoLocalizedRingElem` is 
#    not `RingElemType`, but the type of `Q`. 
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
    modulus(Q) == I || error("the modulus of the quotient ring does not coincide with the ideal")
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

### additional constructors
function quo(
    W::MPolyLocalizedRing{BRT, BRET, RT, RET, MST},
    I::MPolyLocalizedIdeal{BRT, BRET, RT, RET, MST}
  ) where {BRT, BRET, RT, RET, MST}
  R = base_ring(W)
  S = inverted_set(W)
  lbpa = groebner_basis(I) # In particular, this saturates the ideal
  J = ideal(R, numerator.(oscar_gens(lbpa))) # the preimage of I in R
  return MPolyQuoLocalizedRing(R, J, S, quo(R, J), W)
end

function Localization(Q::MPolyQuo{RET}, S::MultSetType) where {RET <: RingElem, MultSetType <: AbsMultSet}
  return MPolyQuoLocalizedRing(base_ring(Q), modulus(Q), S, Q, Localization(S))
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


