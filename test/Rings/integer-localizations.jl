#################################################################################
# A minimal implementation of the localization interface for the ring ℤ         #
#################################################################################

using Oscar

import Oscar: base_ring, inverted_set, ambient_ring, Localization, parent, numerator, denominator, one, zero, reduce_fraction
import Oscar.AbstractAlgebra: elem_type, parent_type

export FmpzComplementOfPrimeIdeal, FmpzPowersOfElement, FmpzComplementOfZeroIdeal
export FmpzLocalizedRing, FmpzLocalizedRingElem


#################################################################################
# Multiplicatively closed sets given by powers of a specific integer            #
#################################################################################

@doc raw"""
FmpzPowersOfElement

The multiplicative set given by the powers of some element.
"""
mutable struct FmpzPowersOfElement <: AbsMultSet{ZZRing, ZZRingElem}
  d::Vector{ZZRingElem} # a vector of admissible denominators. 

  function FmpzPowersOfElement(denominators::Vector{ZZRingElem})
    return new(denominators)
  end
end

### required getter functions
ambient_ring(S::FmpzPowersOfElement) = ZZ

### required functionality
function Base.in(a::ZZRingElem, S::FmpzPowersOfElement) 
  # check whether ∃ c ∈ ℤ, k ∈ ℕ₀: c ⋅ a = (∏ᵢ dᵢ)ᵏ, where dᵢ are the admissible denominators in S.
  d = prod(denominators(S))
  g = gcd(a, d)
  while !isone(g)
    a = divexact(a, g)
    g = gcd(a, d)
  end
  return isone(a)
end

### additional constructors
FmpzPowersOfElement(d::ZZRingElem) = FmpzPowersOfElement([d])
FmpzPowersOfElement(d::Oscar.IntegerUnion) = FmpzPowersOfElement(ZZ(d))
FmpzPowersOfElement(d::Vector{T}) where {T<:Oscar.IntegerUnion} = FmpzPowersOfElement(ZZ.(d))

### additional functionality
denominators(S::FmpzPowersOfElement) = S.d::Vector{ZZRingElem}
Base.in(a::Oscar.IntegerUnion, S::FmpzPowersOfElement) = (ZZ(a) in S)


#################################################################################
# Complements of prime ideals in ℤ                                              #
#################################################################################

@doc raw"""
FmpzComplementOfPrimeIdeal

The multiplicatively closed set `S = ℤ ∖ ⟨p⟩` of integers outside a prime ideal `⟨p⟩`.
"""
mutable struct FmpzComplementOfPrimeIdeal <: AbsMultSet{ZZRing, ZZRingElem}
  # essential fields
  p::ZZRingElem

  function FmpzComplementOfPrimeIdeal(p::ZZRingElem)
    is_prime(p) || error("$(p) is not prime")
    return new(p)
  end
end

### required getter functions
ambient_ring(S::FmpzComplementOfPrimeIdeal) = ZZ

### required functionality
function Base.in(b::ZZRingElem, S::FmpzComplementOfPrimeIdeal)
  return !(divides(b, S.p)[1])
end

### additional constructors
FmpzComplementOfPrimeIdeal(i::Integer) = FmpzComplementOfPrimeIdeal(ZZ(i))

### additional functionality
Base.in(a::Oscar.IntegerUnion, S::FmpzComplementOfPrimeIdeal) = (ZZ(a) in S)


#################################################################################
# The complement of the zero ideal in ℤ                                         #
#################################################################################

@doc raw"""
FmpzComplementOfZeroIdeal 

The complement of the zero ideal in `ℤ`.
"""
mutable struct FmpzComplementOfZeroIdeal <: AbsMultSet{ZZRing, ZZRingElem}
  function FmpzComplementOfZeroIdeal() 
    return new{}()
  end
end

### required getter functions
ambient_ring(S::FmpzComplementOfZeroIdeal) = ZZ

### required functionality
function Base.in(b::ZZRingElem, S::FmpzComplementOfZeroIdeal) 
  return !iszero(b)
end

### additional functionality
Base.in(b::Oscar.IntegerUnion, S::FmpzComplementOfZeroIdeal) = (ZZ(b) in S)


#################################################################################
# Localizations of ℤ                                                            #
#################################################################################

@doc raw"""
FmpzLocalizedRing{MultSetType <: AbsMultSet{ZZRing, ZZRingElem}} <: AbsLocalizedRing{ZZRing, ZZRingElem, MultSetType} 

A minimal implementation for the localization `ℤ[S⁻¹]` of the ring of integers 
at a multiplicatively closed set `S` of type `MultSetType`.
"""
mutable struct FmpzLocalizedRing{MultSetType <: AbsMultSet{ZZRing, ZZRingElem}} <: AbsLocalizedRing{ZZRing, ZZRingElem, MultSetType} 
  S::MultSetType # The multiplicatively closed set S ⊂ R whose inverses are added to R

  function FmpzLocalizedRing(S::MultSetType) where {MultSetType <: AbsMultSet{ZZRing, ZZRingElem}}
    # Sanity check whether the multiplicatively closed set S is compatible with the 
    # given rings
    MultSetType <: AbsMultSet{ZZRing, ZZRingElem} || error(
	"The type of the multiplicatively closed set is not compatible with the type of the ring"
	)
    return new{MultSetType}(S)
  end
end

### required getter functions
base_ring(W::FmpzLocalizedRing) = ZZ::ZZRing
inverted_set(W::FmpzLocalizedRing{MultSetType}) where {MultSetType} = W.S::MultSetType

### required extensions of the localization function
function Localization(S::FmpzComplementOfPrimeIdeal) 
  L = FmpzLocalizedRing(S)
  return L, MapFromFunc(x->L(x), base_ring(L), L)
end

function Localization(S::FmpzComplementOfZeroIdeal)
  L = FmpzLocalizedRing(S)
  return L, MapFromFunc(x->L(x), base_ring(L), L)
end

function Localization(S::FmpzPowersOfElement)
  L = FmpzLocalizedRing(S)
  return L, MapFromFunc(x->L(x), base_ring(L), L)
end


#######################################################################
# Elements of localizations of ℤ                                      #
#######################################################################

@doc raw"""
    FmpzLocalizedRingElem{MultSetType}

Elements `a/b ∈ ℤ[S⁻¹]` of localizations of the ring of integers 
at a multiplicatively closed set `S` of type `MultSetType`.
"""
mutable struct FmpzLocalizedRingElem{MultSetType} <: AbsLocalizedRingElem{ZZRing, ZZRingElem, MultSetType}
  numerator::ZZRingElem
  denominator::ZZRingElem
  R::FmpzLocalizedRing{MultSetType} # the parent ring

  function FmpzLocalizedRingElem(R::FmpzLocalizedRing{MultSetType}, a::ZZRingElem, b::ZZRingElem) where {MultSetType}
    # Perform some sanity checks
    b in inverted_set(R) || error("$b does not belong to the units of $R")
    return new{MultSetType}(a, b, R)
  end
end

### required getter functions
parent(f::FmpzLocalizedRingElem{MultSetType}) where {MultSetType} = f.R::FmpzLocalizedRing{MultSetType}
numerator(f::FmpzLocalizedRingElem) = f.numerator::ZZRingElem
denominator(f::FmpzLocalizedRingElem) = f.denominator::ZZRingElem

### required implementation of the arithmetic
function Base.:(//)(a::ZZRingElem, b::FmpzLocalizedRingElem)
  c = ZZ(a)
  g = gcd(c, numerator(b))
  c = divexact(c, g)
  numerator(b) in inverted_set(parent(b)) || error("the second argument is not a unit in this local ring")
  return reduce_fraction((parent(b))(c*denominator(b), divexact(numerator(b), g)))
end

Base.:(//)(a::Integer, b::FmpzLocalizedRingElem) = ZZ(a)//b

function Base.:(//)(a::T, b::T) where {T<:FmpzLocalizedRingElem}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  g = gcd(numerator(a), numerator(b))
  c = divexact(numerator(a), g)
  d = divexact(numerator(b), g)
  d in inverted_set(parent(b)) || error("the second argument is not a unit in this local ring")
  return reduce_fraction((parent(a))(c*denominator(b), d*denominator(a)))
end

### optional enhancement of the arithmetic
function reduce_fraction(f::FmpzLocalizedRingElem) 
  g = gcd(numerator(f), denominator(f))
  return FmpzLocalizedRingElem(parent(f), divexact(numerator(f), g), divexact(denominator(f), g))
end

### required conversions
(W::FmpzLocalizedRing)(f::ZZRingElem) = FmpzLocalizedRingElem(W, f, ZZ(1))
(W::FmpzLocalizedRing)(a::ZZRingElem, b::ZZRingElem) = FmpzLocalizedRingElem(W, a, b)

### additional conversions
(W::FmpzLocalizedRing)(q::QQFieldElem) = FmpzLocalizedRingElem(W, numerator(q), denominator(q))
(W::FmpzLocalizedRing)(q::Rational{T}) where {T<:Oscar.IntegerUnion} = FmpzLocalizedRingElem(W, ZZ(numerator(q)), ZZ(denominator(q)))
(W::FmpzLocalizedRing)(a::Oscar.IntegerUnion, b::Oscar.IntegerUnion) = FmpzLocalizedRingElem(W, ZZ(a), ZZ(b))

### implementation of Oscar's general ring interface
one(W::FmpzLocalizedRing) = W(1)
zero(W::FmpzLocalizedRing) = W(0)

elem_type(W::FmpzLocalizedRing{MultSetType}) where {MultSetType} = FmpzLocalizedRingElem{MultSetType}
elem_type(T::Type{FmpzLocalizedRing{MultSetType}}) where {MultSetType} = FmpzLocalizedRingElem{MultSetType}

parent_type(W::FmpzLocalizedRingElem{MultSetType}) where {MultSetType} = FmpzLocalizedRing{MultSetType}
parent_type(T::Type{FmpzLocalizedRingElem{MultSetType}}) where {MultSetType} = FmpzLocalizedRing{MultSetType}


########################################################################
# The actual tests for the above implementation                        #
########################################################################

@testset "integer-localizations" begin
  S = FmpzPowersOfElement([5,14])
  @test 25 in S
  @test !(33 in S)
  @test 5*5*7 in S
  @test ambient_ring(S) == ZZ
  W, _ = Localization(S)
  @test base_ring(W) == ambient_ring(S)
  @test inverted_set(W) == S
  a = W(3)
  b = W(7, 5)
  @test a+b == W(3*5+7, 5)
  @test a-b == W(3*5-7, 5)
  @test 1//b == W(5,7)
  @test b^2 == W(49, 25)
  c = W(4, 5)
  @test c//c == one(W)

  U = FmpzComplementOfPrimeIdeal(13)
  @test !(13^5 in U)
  @test !(13*4289729837 in U)
  @test 5783790198374098 in U
  @test ambient_ring(U) == ZZ
  W, _ = Localization(U)
  @test base_ring(W) == ambient_ring(U)
  @test inverted_set(W) == U
  a = W(4, 17)
  b = W(4*17)
  b = b//W(19)
  @test a//b == W(19//17^2)
  @test a - b == W( 4//17 - 4*17//19 )
  @test a + b == W( 4//17 + 4*17//19 )
  @test a * b == W( 4//17 * 4*17//19 )

  O = FmpzComplementOfZeroIdeal()
  @test 234890 in O
  @test !(0 in O)
  @test ambient_ring(O) == ZZ
  W, _ = Localization(O)
  @test base_ring(W) == ambient_ring(O)
  @test inverted_set(W) == O
  a = W(4, 17)
  b = W(4*17)
  b = b//W(19)
  @test a//b == W(19//17^2)
  @test a - b == W( 4//17 - 4*17//19 )
  @test a + b == W( 4//17 + 4*17//19 )
  @test a * b == W( 4//17 * 4*17//19 )
end
