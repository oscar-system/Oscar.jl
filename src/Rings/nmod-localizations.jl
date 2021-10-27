#################################################################
# Localizations of residue rings ℤ/nℤ                           #
#################################################################

import Nemo.NmodRing

export NmodComplementOfPrimeIdeal, NmodLocalizedRing, NmodLocalizedRingElem

export generator, ambient_ring, localize_at, parent, numerator, denominator


#######################################################################
# Types of multiplicatively closed sets in ℤ/nℤ                       #
#######################################################################

@Markdown.doc """
    NmodComplementOfPrimeIdeal <: AbsMultSet{NmodRing, nmod}

Complement of a prime ideal in a quotient ring `ℤ/nℤ`.
"""
mutable struct NmodComplementOfPrimeIdeal <: AbsMultSet{NmodRing, nmod}
  R::NmodRing # the ambient ring
  gen::nmod # a generator of the prime ideal

  function NmodComplementOfPrimeIdeal(gen::nmod)
    R = parent(gen)
    n = ZZ(modulus(R))
    a = lift(gen)
    r = gcd(n, a)
    isprime(r) || error("the given element does not generate a prime ideal")
    return new{}(R, R(r))
  end
end

### required functionality
function Base.in(b::nmod, S::NmodComplementOfPrimeIdeal) 
  return mod(lift(b),lift(generator(S))) != zero(b)
end

### required getter functions
ambient_ring(S::NmodComplementOfPrimeIdeal) = S.R

### additional constructors
NmodComplementOfPrimeIdeal(R::NmodRing, i::Oscar.IntegerUnion) = NmodComplementOfPrimeIdeal(R(i))

### additional functionality
generator(S::NmodComplementOfPrimeIdeal) = S.gen

#######################################################################
# Localizations of ℤ/nℤ                                               #
#######################################################################

@Markdown.doc """
    NmodLocalizedRing{MultSetType} <: AbsLocalizedRing{NmodRing, nmod, MultSetType}

Localization of a ring `ℤ/nℤ` at a multiplicatively closed set of type `MultSetType`.
"""
mutable struct NmodLocalizedRing{MultSetType} <: AbsLocalizedRing{NmodRing, nmod, MultSetType}
  R::NmodRing # the original ring before localization
  S::MultSetType # the set at which has been localized

  function NmodLocalizedRing(S::MultSetType) where {MultSetType}
    return new{MultSetType}(ambient_ring(S), S)
  end
end

### required getter functions
original_ring(W::NmodLocalizedRing{MultSetType}) where {MultSetType} = W.R
inverted_set(W::NmodLocalizedRing{MultSetType}) where {MultSetType} = W.S

### required extension of the localization function
localize_at(S::NmodComplementOfPrimeIdeal) = NmodLocalizedRing(S)


#######################################################################
# Elements in localizations of ℤ/nℤ                                   #
#######################################################################

@Markdown.doc """
    NmodLocalizedRingElem{MultSetType} <: AbsLocalizedRingElem{NmodRing, nmod, MultSetType}

Elements of localizations of quotient rings `ℤ/nℤ` at a 
multiplicatively closed set of type `MultSetType`.
"""
mutable struct NmodLocalizedRingElem{MultSetType} <: AbsLocalizedRingElem{NmodRing, nmod, MultSetType}
  numerator::nmod
  denominator::nmod
  W::NmodLocalizedRing{MultSetType} # the parent ring

  function NmodLocalizedRingElem(W::NmodLocalizedRing{MultSetType}, a::nmod, b::nmod) where {MultSetType} 
    original_ring(W) == parent(a) == parent(b) || error("elements do not belong to the original ring")
    b in inverted_set(W) || error("the given denominator is not an admissible unit in this ring")
    return new{MultSetType}(a, b, W)
  end
end

### required getter functions
parent(f::NmodLocalizedRingElem{MultSetType}) where {MultSetType} = f.W
numerator(f::NmodLocalizedRingElem{MultSetType}) where {MultSetType} = f.numerator
denominator(f::NmodLocalizedRingElem{MultSetType}) where {MultSetType} = f.denominator

### required conversions
(W::NmodLocalizedRing{MultSetType})(a::nmod, b::nmod) where {MultSetType} = NmodLocalizedRingElem(W, a, b)
(W::NmodLocalizedRing{MultSetType})(a::nmod) where {MultSetType} = NmodLocalizedRingElem(W, a, one(parent(a)))

### additional conversions
(W::NmodLocalizedRing{MultSetType})(a::Oscar.IntegerUnion, b::Oscar.IntegerUnion) where {MultSetType} = NmodLocalizedRingElem(W, original_ring(W)(a), original_ring(W)(b))
(W::NmodLocalizedRing{MultSetType})(a::Oscar.IntegerUnion) where {MultSetType} = NmodLocalizedRingElem(W, original_ring(W)(a), one(original_ring(W)))
