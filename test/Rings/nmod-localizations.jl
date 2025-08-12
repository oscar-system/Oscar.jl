#################################################################
# Localizations of residue rings ℤ/nℤ                           #
#################################################################

using Oscar

import Oscar: base_ring, base_ring_type, inverted_set, ring, localization, parent, numerator, denominator, one, zero
import Oscar.AbstractAlgebra: elem_type, parent_type

export NmodComplementOfPrimeIdeal, NmodLocalizedRing, NmodLocalizedRingElem

export generator, ring, localization, parent, numerator, denominator


#######################################################################
# Types of multiplicatively closed sets in ℤ/nℤ                       #
#######################################################################

@doc raw"""
    NmodComplementOfPrimeIdeal <: AbsMultSet{zzModRing, zzModRingElem}

Complement of a prime ideal in a quotient ring `ℤ/nℤ`.
"""
mutable struct NmodComplementOfPrimeIdeal <: AbsMultSet{zzModRing, zzModRingElem}
  R::zzModRing # the ambient ring
  gen::zzModRingElem # a generator of the prime ideal

  function NmodComplementOfPrimeIdeal(gen::zzModRingElem)
    R = parent(gen)
    n = ZZ(modulus(R))
    a = lift(gen)
    r = gcd(n, a)
    is_prime(r) || error("the given element does not generate a prime ideal")
    return new{}(R, R(r))
  end
end

### required functionality
function Base.in(b::zzModRingElem, S::NmodComplementOfPrimeIdeal) 
  return mod(lift(b), lift(generator(S))) != zero(b)
end

### required getter functions
ring(S::NmodComplementOfPrimeIdeal) = S.R

### additional constructors
NmodComplementOfPrimeIdeal(R::zzModRing, i::Oscar.IntegerUnion) = NmodComplementOfPrimeIdeal(R(i))

### additional functionality
generator(S::NmodComplementOfPrimeIdeal) = S.gen
Base.in(b::Oscar.IntegerUnion, S::NmodComplementOfPrimeIdeal) = (ring(S)(b) in S)

#######################################################################
# Localizations of ℤ/nℤ                                               #
#######################################################################

@doc raw"""
NmodLocalizedRing{MultSetType <: AbsMultSet{zzModRing, zzModRingElem}} <: AbsLocalizedRing{zzModRing, zzModRingElem, MultSetType}

Localization of a ring `ℤ/nℤ` at a multiplicatively closed set of type `MultSetType`.
"""
mutable struct NmodLocalizedRing{MultSetType <: AbsMultSet{zzModRing, zzModRingElem}} <: AbsLocalizedRing{zzModRing, zzModRingElem, MultSetType}
  R::zzModRing # the original ring before localization
  S::MultSetType # the set at which has been localized

  function NmodLocalizedRing(S::MultSetType) where {MultSetType <: AbsMultSet{zzModRing, zzModRingElem}}
    return new{MultSetType}(ring(S), S)
  end
end

### required getter functions
base_ring(W::NmodLocalizedRing) = W.R::zzModRing
base_ring_type(::Type{<:NmodLocalizedRing}) = zzModRing
inverted_set(W::NmodLocalizedRing{MultSetType}) where {MultSetType} = W.S::MultSetType

### required extension of the localization function
function localization(S::NmodComplementOfPrimeIdeal)
  L = NmodLocalizedRing(S)
  return L, MapFromFunc(base_ring(L), L, L)
end


#######################################################################
# Elements in localizations of ℤ/nℤ                                   #
#######################################################################

@doc raw"""
    NmodLocalizedRingElem{MultSetType} <: AbsLocalizedRingElem{zzModRing, zzModRingElem, MultSetType}

Elements of localizations of quotient rings `ℤ/nℤ` at a 
multiplicatively closed set of type `MultSetType`.
"""
mutable struct NmodLocalizedRingElem{MultSetType} <: AbsLocalizedRingElem{zzModRing, zzModRingElem, MultSetType}
  numerator::zzModRingElem
  denominator::zzModRingElem
  W::NmodLocalizedRing{MultSetType} # the parent ring

  function NmodLocalizedRingElem(W::NmodLocalizedRing{MultSetType}, a::zzModRingElem, b::zzModRingElem; check::Bool=true) where {MultSetType} 
    base_ring(W) == parent(a) == parent(b) || error("elements do not belong to the original ring")
    @check b in inverted_set(W) "the given denominator is not an admissible unit in this ring"
    return new{MultSetType}(a, b, W)
  end
end

### required getter functions
parent(f::NmodLocalizedRingElem) = f.W
numerator(f::NmodLocalizedRingElem) = f.numerator
denominator(f::NmodLocalizedRingElem) = f.denominator

### required conversions
(W::NmodLocalizedRing)(a::zzModRingElem, b::zzModRingElem; check::Bool=true) = NmodLocalizedRingElem(W, a, b, check=check)
(W::NmodLocalizedRing)(a::zzModRingElem) = NmodLocalizedRingElem(W, a, one(parent(a)), check=false)

### required implementation of the arithmetic
Base.:(//)(a::Oscar.IntegerUnion, b::NmodLocalizedRingElem) = ((parent(b)(a))//b)

function Base.:(//)(a::T, b::T) where {T<:NmodLocalizedRingElem}
  parent(a) == parent(b) || error("the arguments do not belong to the same ring")
  g = gcd(numerator(a), numerator(b))
  c = divexact(numerator(b), g)
  c in inverted_set(parent(b)) || error("the second argument is not a unit in this local ring")
  return reduce_fraction((parent(a))(divexact(numerator(a),g)*denominator(b), c*denominator(a)))
end

### additional conversions
(W::NmodLocalizedRing)(a::T, b::T; check::Bool=true) where {T<:Oscar.IntegerUnion} = W(base_ring(W)(a), base_ring(W)(b), check=check)
(W::NmodLocalizedRing)(a::Oscar.IntegerUnion) = W(base_ring(W)(a), one(base_ring(W)), check=false)
(W::NmodLocalizedRing)(q::QQFieldElem; check::Bool=true) = W(numerator(q), denominator(q), check=check)
(W::NmodLocalizedRing)(i::Int64) = W(base_ring(W)(i), one(base_ring(W)), check=false)
(W::NmodLocalizedRing)(q::Rational{T}; check::Bool=true) where {T<:Oscar.IntegerUnion} = W(numerator(q), denominator(q), check=check)

### implementation of Oscar's general ring interface
one(W::NmodLocalizedRing) = W(1)
zero(W::NmodLocalizedRing) = W(0)

elem_type(W::NmodLocalizedRing{MultSetType}) where {MultSetType} = NmodLocalizedRingElem{MultSetType}
elem_type(T::Type{NmodLocalizedRing{MultSetType}}) where {MultSetType} = NmodLocalizedRingElem{MultSetType}

parent_type(W::NmodLocalizedRingElem{MultSetType}) where {MultSetType} = NmodLocalizedRing{MultSetType}
parent_type(T::Type{NmodLocalizedRingElem{MultSetType}}) where {MultSetType} = NmodLocalizedRing{MultSetType}


########################################################################
# The actual tests for the above implementation                        #
########################################################################

@testset "zzModRingElem-localizations" begin
  R, p = quo(ZZ, 101*13)

  U = NmodComplementOfPrimeIdeal(R(13))
  @test !(13^5 in U)
  @test !(13*4289729837 in U)
  @test 5783790198374098 in U
  @test ring(U) == R
  W, _ = localization(U)
  @test base_ring(W) == ring(U)
  @test inverted_set(W) == U
  a = W(4, 17)
  b = W(4*17)
  b = b//W(19)
  @test a//b == W(19//17^2)
  @test a - b == W( 4//17 - 4*17//19 )
  @test a + b == W( 4//17 + 4*17//19 )
  @test a * b == W( 4//17 * 4*17//19 )
  b = W(13//19)
  @test b//b == one(W)
end
