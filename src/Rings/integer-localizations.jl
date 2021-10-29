#################################################################################
# A minimal implementation of the localization interface for the ring ℤ         #
#################################################################################

export FmpzComplementOfPrimeIdeal, FmpzPowersOfElement, FmpzComplementOfZeroIdeal
export FmpzLocalizedRing, FmpzLocalizedRingElem
export generator, ambient_ring, localize_at, parent, numerator, denominator


#################################################################################
# Multiplicatively closed sets given by powers of a specific integer            #
#################################################################################

@Markdown.doc """
FmpzPowersOfElement

The multiplicative set given by the powers of some element.
"""
mutable struct FmpzPowersOfElement <: AbsMultSet{FlintIntegerRing, fmpz}
  d::Vector{fmpz} # a vector of admissible denominators. 

  function FmpzPowersOfElement(denominators::Vector{fmpz})
    return new(denominators)
  end
end

### required getter functions
ambient_ring(S::FmpzPowersOfElement) = ZZ

### required functionality
function Base.in(a::fmpz, S::FmpzPowersOfElement) 
  # check whether ∃ c ∈ ℤ, k ∈ ℕ₀: c ⋅ a = (∏ᵢ dᵢ)ᵏ, where dᵢ are the admissible denominators in S.
  d = prod(denominators(S))
  g = gcd(a, d)
  while(!isone(g))
    a = divexact(a, g)
    g = gcd(a, d)
  end
  return isone(a)
end

### additional constructors
FmpzPowersOfElement(d::fmpz) = FmpzPowersOfElement([d])
FmpzPowersOfElement(d::Oscar.IntegerUnion) = FmpzPowersOfElement(ZZ(d))
FmpzPowersOfElement(d::Vector{T}) where {T<:Oscar.IntegerUnion} = FmpzPowersOfElement(ZZ.(d))

### additional functionality
denominators(S::FmpzPowersOfElement) = S.d::Vector{fmpz}
Base.in(a::Oscar.IntegerUnion, S::FmpzPowersOfElement) = (ZZ(a) in S)


#################################################################################
# Complements of prime ideals in ℤ                                              #
#################################################################################

@Markdown.doc """
FmpzComplementOfPrimeIdeal

The multiplicatively closed set `S = ℤ ∖ ⟨p⟩` of integers outside a prime ideal `⟨p⟩`.
"""
mutable struct FmpzComplementOfPrimeIdeal <: AbsMultSet{FlintIntegerRing, fmpz}
  # essential fields
  p::fmpz

  function FmpzComplementOfPrimeIdeal(p::fmpz)
    isprime(p) || error("$(p) is not prime")
    return new(p)
  end
end

### required getter functions
ambient_ring(S::FmpzComplementOfPrimeIdeal) = ZZ

### required functionality
function Base.in(b::fmpz, S::FmpzComplementOfPrimeIdeal)
  return !(divides(b, S.p)[1])
end

### additional constructors
FmpzComplementOfPrimeIdeal(i::Integer) = FmpzComplementOfPrimeIdeal(ZZ(i))

### additional functionality
Base.in(a::Oscar.IntegerUnion, S::FmpzComplementOfPrimeIdeal) = (ZZ(a) in S)


#################################################################################
# The complement of the zero ideal in ℤ                                         #
#################################################################################

@Markdown.doc """
FmpzComplementOfZeroIdeal 

The complement of the zero ideal in `ℤ`.
"""
mutable struct FmpzComplementOfZeroIdeal <: AbsMultSet{FlintIntegerRing, fmpz}
  function FmpzComplementOfZeroIdeal() 
    return new{}()
  end
end

### required getter functions
ambient_ring(S::FmpzComplementOfZeroIdeal) = ZZ

### required functionality
function Base.in(b::fmpz, S::FmpzComplementOfZeroIdeal) 
  return !(b == zero(ZZ))
end

### additional functionality
Base.in(b::Oscar.IntegerUnion, S::FmpzComplementOfZeroIdeal) = (ZZ(b) in S)


#################################################################################
# Localizations of ℤ                                                            #
#################################################################################

@Markdown.doc """
    FmpzLocalizedRing{MultSetType} <: AbsLocalizedRing{FlintIntegerRing, fmpz, MultSetType} 

A minimal implementation for the localization `ℤ[S⁻¹]` of the ring of integers 
at a multiplicatively closed set `S` of type `MultSetType`.
"""
mutable struct FmpzLocalizedRing{MultSetType <: AbsMultSet{FlintIntegerRing, fmpz}} <: AbsLocalizedRing{FlintIntegerRing, fmpz, MultSetType} 
  S::MultSetType # The multiplicatively closed set S ⊂ R whose inverses are added to R

  function FmpzLocalizedRing(S::MultSetType) where {MultSetType}
    # Sanity check whether the multiplicatively closed set S is compatible with the 
    # given rings
    MultSetType <: AbsMultSet{FlintIntegerRing, fmpz} || error(
	"The type of the multiplicatively closed set is not compatible with the type of the ring"
	)
    return new{MultSetType}(S)
  end
end

### required getter functions
original_ring(W::FmpzLocalizedRing{MultSetType}) where {MultSetType} = ZZ::FlintIntegerRing
inverted_set(W::FmpzLocalizedRing{MultSetType}) where {MultSetType} = W.S::MultSetType

### required extensions of the localization function
localize_at(S::FmpzComplementOfPrimeIdeal) = FmpzLocalizedRing(S)
localize_at(S::FmpzComplementOfZeroIdeal) = FmpzLocalizedRing(S)
localize_at(S::FmpzPowersOfElement) = FmpzLocalizedRing(S)


#######################################################################
# Elements of localizations of ℤ                                      #
#######################################################################

@Markdown.doc """
    FmpzLocalizedRingElem{MultSetType}

Elements `a/b ∈ ℤ[S⁻¹]` of localizations of the ring of integers 
at a multiplicatively closed set `S` of type `MultSetType`.
"""
mutable struct FmpzLocalizedRingElem{MultSetType} <: AbsLocalizedRingElem{FlintIntegerRing, fmpz, MultSetType}
  numerator::fmpz
  denominator::fmpz
  R::FmpzLocalizedRing{MultSetType} # the parent ring

  function FmpzLocalizedRingElem(R::FmpzLocalizedRing{MultSetType}, a::fmpz, b::fmpz) where {MultSetType}
    # Perform some sanity checks
    b in inverted_set(R) || error("$b does not belong to the units of $R")
    return new{MultSetType}(a, b, R)
  end
end

### required getter functions
parent(f::FmpzLocalizedRingElem{MultSetType}) where {MultSetType} = f.R::FmpzLocalizedRing{MultSetType}
numerator(f::FmpzLocalizedRingElem{MultSetType}) where {MultSetType} = f.numerator::fmpz
denominator(f::FmpzLocalizedRingElem{MultSetType}) where {MultSetType} = f.denominator::fmpz

### required implementation of the arithmetic
function Base.:(//)(a::Oscar.IntegerUnion, b::FmpzLocalizedRingElem)
  c = ZZ(a)
  g = gcd(c, numerator(b))
  c = divexact(c, g)
  numerator(b) in inverted_set(parent(b)) || error("the second argument is not a unit in this local ring")
  return reduce_fraction((parent(b))(c*denominator(b), divexact(numerator(b), g)))
end

function Base.:(//)(a::T, b::T) where {T<:FmpzLocalizedRingElem}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  g = gcd(numerator(a), numerator(b))
  c = divexact(numerator(a), g)
  d = divexact(numerator(b), g)
  d in inverted_set(parent(b)) || error("the second argument is not a unit in this local ring")
  return reduce_fraction((parent(a))(c*denominator(b), d*denominator(a)))
end

### optional enhancement of the arithmetic
function reduce_fraction(f::FmpzLocalizedRingElem{MultSetType}) where {MultSetType}
  g = gcd(numerator(f), denominator(f))
  return FmpzLocalizedRingElem(parent(f), divexact(numerator(f), g), divexact(denominator(f), g))
end

### required conversions
(W::FmpzLocalizedRing{MultSetType})(f::fmpz) where {MultSetType} = FmpzLocalizedRingElem(W, f, one(ZZ))
(W::FmpzLocalizedRing{MultSetType})(a::fmpz, b::fmpz) where {MultSetType} = FmpzLocalizedRingElem(W, a, b)

### additional conversions
(W::FmpzLocalizedRing{MultSetType})(q::fmpq) where {MultSetType} = FmpzLocalizedRingElem(W, numerator(q), denominator(q))
(W::FmpzLocalizedRing{MultSetType})(q::Rational{Int}) where {MultSetType} = FmpzLocalizedRingElem(W, ZZ(numerator(q)), ZZ(denominator(q)))
@Markdown.doc """
    (W::FmpzLocalizedRing{MultSetType})(i::Oscar.IntegerUnion) where {MultSetType} = FmpzLocalizedRingElem(W, ZZ(i), one(ZZ))

Part of the minimal interface for localized rings. This routine returns the conversion of 
an integer i to an element i//1 ∈ W.
"""
(W::FmpzLocalizedRing{MultSetType})(i::Oscar.IntegerUnion) where {MultSetType} = FmpzLocalizedRingElem(W, ZZ(i), one(ZZ))
(W::FmpzLocalizedRing{MultSetType})(a::Oscar.IntegerUnion, b::Oscar.IntegerUnion) where {MultSetType} = FmpzLocalizedRingElem(W, ZZ(a), ZZ(b))

### implementation of Oscar's general ring interface
one(W::FmpzLocalizedRing{MultSetType}) where {MultSetType} = W(1)
zero(W::FmpzLocalizedRing{MultSetType}) where {MultSetType} = W(0)

elem_type(W::FmpzLocalizedRing{MultSetType}) where {MultSetType} = FmpzLocalizedRingElem{MultSetType}
elem_type(T::Type{FmpzLocalizedRing{MultSetType}}) where {MultSetType} = FmpzLocalizedRingElem{MultSetType}

parent_type(W::FmpzLocalizedRingElem{MultSetType}) where {MultSetType} = FmpzLocalizedRing{MultSetType}
parent_type(T::Type{FmpzLocalizedRingElem{MultSetType}}) where {MultSetType} = FmpzLocalizedRing{MultSetType}

