export MPolyUnits, MPolyComplementOfPrimeIdeal, MPolyComplementOfKPointIdeal, MPolyPowersOfElement, MPolyUnionOfMultSets
export rand
export sets

export MPolyLocalizedRing
export ambient_ring, point_coordinates, inverted_set, denominators

export MPolyLocalizedRingElem
export numerator, denominator, fraction, parent
export reduce_fraction

export MPolyLocalizedIdeal
export gens, base_ring, groebner_bases, default_ordering, dim 

export Localization, ideal

export LocalizedBiPolyArray
export oscar_gens, oscar_ring, singular_ring, singular_gens, ordering, shift
export groebner_basis, groebner_assure

export MPolyLocalizedRingHom
export domain, codomain, images

import AbstractAlgebra: Ring, RingElem

########################################################################
# General framework for localizations of multivariate polynomial rings #
########################################################################

########################################################################
# Units in polynomial rings; localization does nothing in this case    #
########################################################################

mutable struct MPolyUnits{
    BaseRingType, 
    BaseRingElemType,
    RingType,
    RingElemType
  } <: AbsMultSet{
    RingType,
    RingElemType
  }

  R::RingType

  function MPolyUnits(R::MPolyRing)
    return new{typeof(coefficient_ring(R)), elem_type(coefficient_ring(R)), typeof(R), elem_type(R)}(R)
  end
end

### required getter functions
ambient_ring(S::MPolyUnits) = S.R

### required functionality
function Base.in(
    f::RingElemType, 
    S::MPolyUnits{BaseRingType, BaseRingElemType, RingType, RingElemType}
  ) where {BaseRingType, BaseRingElemType, RingType, RingElemType}
  return divides(one(ambient_ring(S)), f)[1]
end

### printing
function Base.show(io::IO, S::MPolyUnits)
  print(io, "units of ")
  print(io, ambient_ring(S))
end

### generation of random elements 
function rand(S::MPolyUnits, v1::UnitRange{Int}, v2::UnitRange{Int}, v3::UnitRange{Int})
  return one(ambient_ring(S))
end


########################################################################
# Complements of prime ideals                                          #
########################################################################

@Markdown.doc """
MPolyComplementOfPrimeIdeal{
    BaseRingType, 
    BaseRingElemType,
    RingType,
    RingElemType
  } <: AbsMultSet{
    RingType,
    RingElemType
  }

The complement of a prime ideal `P ⊂ 𝕜[x₁,…,xₙ]` in a multivariate polynomial ring 
with elements of type `RingElemType` over a base ring `𝕜` of type `BaseRingType`.
"""
mutable struct MPolyComplementOfPrimeIdeal{
    BaseRingType, 
    BaseRingElemType,
    RingType,
    RingElemType
  } <: AbsMultSet{
    RingType,
    RingElemType
  }

  # The parent polynomial ring 𝕜[x₁,…,xₙ]
  R::RingType
  # The prime ideal whose complement this is
  P::MPolyIdeal{RingElemType}

  function MPolyComplementOfPrimeIdeal(
      P::MPolyIdeal{RingElemType}; 
      check::Bool=false
    ) where {RingElemType}
    R = base_ring(P)
    check && (isprime(P) || error("the ideal $P is not prime"))
    return new{typeof(coefficient_ring(R)), elem_type(coefficient_ring(R)), typeof(R), elem_type(R)}(R, P)
  end
end

### required getter functions
ambient_ring(
    S::MPolyComplementOfPrimeIdeal) = S.R

### required functionality
function Base.in(
    f::RingElemType, 
    S::MPolyComplementOfPrimeIdeal{BaseRingType, BaseRingElemType, RingType, RingElemType}
  ) where {BaseRingType, BaseRingElemType, RingType, RingElemType}
  return !(f in prime_ideal(S))
end

### additional functionality
prime_ideal(S::MPolyComplementOfPrimeIdeal) = S.P

### printing
function Base.show(io::IO, S::MPolyComplementOfPrimeIdeal)
  print(io, "complement of ")
  print(io, prime_ideal(S))
end

### generation of random elements 
function rand(S::MPolyComplementOfPrimeIdeal, v1::UnitRange{Int}, v2::UnitRange{Int}, v3::UnitRange{Int})
  f = rand(ambient_ring(S), v1, v2, v3)
  if f in prime_ideal(S)
    return rand(S, v1, v2, v3)
  end
  return f
end

########################################################################
# Complements of maximal ideals corresponding to 𝕜-points              #
########################################################################

@Markdown.doc """
MPolyComplementOfKPointIdeal{
    BaseRingType,
    BaseRingElemType, 
    RingType,
    RingElemType
  } <: AbsMultSet{
    RingType, 
    RingElemType
  }

Complement of a maximal ideal ``𝔪 = ⟨x₁-a₁,…,xₙ-aₙ⟩⊂ 𝕜[x₁,…xₙ]`` with ``aᵢ∈ 𝕜``.
"""
mutable struct MPolyComplementOfKPointIdeal{
    BaseRingType,
    BaseRingElemType, 
    RingType,
    RingElemType
  } <: AbsMultSet{
    RingType, 
    RingElemType
  }

  # The parent polynomial ring 𝕜[x₁,…,xₙ]
  R::RingType
  # The coordinates aᵢ of the point in 𝕜ⁿ corresponding to the maximal ideal
  a::Vector{BaseRingElemType}

  function MPolyComplementOfKPointIdeal(R::RingType, a::Vector{BaseRingElemType}) where {RingType<:MPolyRing, BaseRingElemType}
    length(a) == length(gens(R)) || error("the number of variables in the ring does not coincide with the number of coordinates")
    n = length(a)
    k = coefficient_ring(R)
    if n > 0 
      base_ring(R) == parent(a[1]) || error("the coordinates are not elements of the base ring")
    else
      elem_type(k) == BaseRingElemType || error("the type of the coordinates does not match the elem_type of the base ring")
    end
    S = new{typeof(k), BaseRingElemType, RingType, elem_type(R)}(R, a)
    return S
  end
end

### required getter functions
ambient_ring(S::MPolyComplementOfKPointIdeal) = S.R

### additional getter functions 
point_coordinates(S::MPolyComplementOfKPointIdeal) = S.a

### required functionality
function Base.in(
    f::RingElemType, 
    S::MPolyComplementOfKPointIdeal{BaseRingType, BaseRingElemType, RingType, RingElemType}
  ) where {BaseRingType, BaseRingElemType, RingType, RingElemType}
  parent(f) == ambient_ring(S) || return false
  return !(evaluate(f, point_coordinates(S)) == zero(ambient_ring(S)))
end

### printing
function Base.show(io::IO, S::MPolyComplementOfKPointIdeal)
  print(io, "point with coordinates ")
  print(io, point_coordinates(S))
end

### generation of random elements 
function rand(S::MPolyComplementOfKPointIdeal, v1::UnitRange{Int}, v2::UnitRange{Int}, v3::UnitRange{Int})
  f = rand(ambient_ring(S), v1, v2, v3)
  if !(f in S)
    return rand(S, v1, v2, v3)
  end
  return f
end

########################################################################
# Powers of elements                                                   #
########################################################################

@Markdown.doc """
MPolyPowersOfElement{
    BaseRingType,
    BaseRingElemType, 
    RingType,
    RingElemType
  } <: AbsMultSet{
    RingType, 
    RingElemType
  }

The set `S = { aᵏ : k ∈ ℕ₀ }` for some ``a ∈ R`` with ``R`` of type `BaseRingType`.
"""
mutable struct MPolyPowersOfElement{
    BaseRingType,
    BaseRingElemType, 
    RingType,
    RingElemType
  } <: AbsMultSet{
    RingType, 
    RingElemType
  }

  R::RingType # the parent ring
  a::Vector{RingElemType} # the list of elements whose powers belong to this set

  function MPolyPowersOfElement(R::RingType, a::Vector{RingElemType}) where {RingType<:MPolyRing, RingElemType<:MPolyElem}
    for f in a 
      parent(f) == R || error("element does not belong to the given ring")
      !iszero(f) || error("can not localize at the zero element")
    end
    k = coefficient_ring(R)
    return new{typeof(k), elem_type(k), RingType, RingElemType}(R, a)
  end
end

### required getter functions
ambient_ring(S::MPolyPowersOfElement) = S.R

### additional constructors
MPolyPowersOfElement(f::RET) where {RET<:MPolyElem} = MPolyPowersOfElement([f])

### additional functionality
denominators(S::MPolyPowersOfElement) = copy(S.a)

### required functionality
function Base.in(
    f::RingElemType, 
    S::MPolyPowersOfElement{BaseRingType, BaseRingElemType, RingType, RingElemType}
  ) where {BaseRingType, BaseRingElemType, RingType, RingElemType}
  R = parent(f)
  R == ambient_ring(S) || return false
  if iszero(f) 
    return false
  end
  d = (length(denominators(S)) == 0 ? one(R) : prod(denominators(S)))
  # We need to check whether for some a ∈ R and k ∈ ℕ we have 
  #   a⋅f = dᵏ.
  (i, o) = ppio(f, d)
  return divides(one(R), o)[1]
end

### printing
function Base.show(io::IO, S::MPolyPowersOfElement)
  print(io, "powers of ")
  print(io, denominators(S))
end

### generation of random elements 
function rand(S::MPolyPowersOfElement, v1::UnitRange{Int}, v2::UnitRange{Int}, v3::UnitRange{Int})
  return prod([f^(abs(rand(Int))%10) for f in denominators(S)])::elem_type(ambient_ring(S))
end

### Taking the product of such sets
function product(
    T::MPolyPowersOfElement{BRT, BRET, RT, RET}, 
    U::MPolyPowersOfElement{BRT, BRET, RT, RET}
  ) where {BRT, BRET, RT, RET}
  a = denominators(T)
  for f in denominators(U)
    for d in a
      (_, f) = ppio(f, d) 
    end
    if !(divides(one(parent(f)), f)[1])
      push!(a, f)
    end
  end
  return MPolyPowersOfElement(a)
end

*(T::MPolyPowersOfElement{BRT, BRET, RT, RET}, 
  U::MPolyPowersOfElement{BRT, BRET, RT, RET}
 ) where {BRT, BRET, RT, RET} = product(T,U)



@Markdown.doc """
MPolyUnionOfMultSets{
    BaseRingType,
    BaseRingElemType, 
    RingType,
    RingElemType
  } <: AbsMultSet{
    RingType, 
    RingElemType
  }

A finite union of arbitrary other multiplicative sets.
"""
mutable struct MPolyUnionOfMultSets{
    BaseRingType,
    BaseRingElemType, 
    RingType,
    RingElemType
  } <: AbsMultSet{
    RingType, 
    RingElemType
  }
  R::RingType
  U::Vector{<:AbsMultSet{RingType, RingElemType}}

  function MPolyUnionOfMultSets(R::RT, U::Vector{<:AbsMultSet{RT, RET}}) where {RT<:MPolyRing, RET<:MPolyElem}
    for s in U
      ambient_ring(s) == R || error("multiplicative set does not live in the given ring")
    end
    return new{typeof(coefficient_ring(R)), elem_type(coefficient_ring(R)), typeof(R), elem_type(R)}(R, U)
  end
end

### required getter functions
ambient_ring(S::MPolyUnionOfMultSets) = S.R

### additional functionality
get_index(S::MPolyUnionOfMultSets, i::Integer) = S.U[i]
sets(S::MPolyUnionOfMultSets) = S.U

### required functionality
function Base.in(
    f::RingElemType, 
    S::MPolyUnionOfMultSets{BaseRingType, BaseRingElemType, RingType, RingElemType}
  ) where {BaseRingType, BaseRingElemType, RingType, RingElemType}
  R = ambient_ring(S)
  divides(one(R), f)[1] && return true
  U = sets(S)
  for s in U
    f in s && return true
  end
  # From here onwards, computations might not be cheap anymore
  a = factor(f)
  for fac in a
    # check for each factor whether it belongs to one of the admissible sets
    tmp_result = false
    for s in U
      if fac[1] in s 
	tmp_result = true
	break
      end
    end
    tmp_result || return false
  end
  return true
end

### printing
function Base.show(io::IO, S::MPolyUnionOfMultSets)
  print(io, "union of the multiplicative sets ")
  for s in sets(S)
    print(io, s)
    print(io, " ")
  end
end

### generation of random elements 
function rand(S::MPolyUnionOfMultSets, v1::UnitRange{Int}, v2::UnitRange{Int}, v3::UnitRange{Int})
  return prod([rand(s, v1, v2, v3) for s in sets(S)])::elem_type(ambient_ring(S))
end


########################################################################
# Localizations of polynomial rings over admissible fields             #
########################################################################

@Markdown.doc """
MPolyLocalizedRing{
    BaseRingType,
    BaseRingElemType,
    RingType,
    RingElemType,
    MultSetType
  } <: AbsLocalizedRing{
    RingType,
    RingType,
    MultSetType
  }

The localization of a multivariate polynomial ring ``R = 𝕜[x₁,…,xₙ]`` over a 
base field ``𝕜`` of type `BaseRingType` and with elements of type `RingElemType` 
at a multiplicative set ``S ⊂ R`` of type `MultSetType`.
"""
mutable struct MPolyLocalizedRing{
    BaseRingType,
    BaseRingElemType,
    RingType,
    RingElemType,
    MultSetType <: AbsMultSet{RingType, RingElemType}
  } <: AbsLocalizedRing{
    RingType,
    RingType,
    MultSetType
  }
  R::RingType # The parent ring which is being localized
  S::MultSetType # The multiplicatively closed set that has been inverted 

  function MPolyLocalizedRing(
      R::RingType, 
      S::MultSetType
    ) where {RingType<:MPolyRing, MultSetType<:AbsMultSet}
    # TODO: Add some sanity checks here?
    ambient_ring(S) == R || error("the multiplicative set is not contained in the given ring")
    k = coefficient_ring(R)
    R_loc = new{typeof(k), elem_type(k), RingType, elem_type(R), MultSetType}(R, S)
    return R_loc
  end
end

### required getter functions 
base_ring(W::MPolyLocalizedRing) = W.R

inverted_set(W::MPolyLocalizedRing) = W.S

### required extension of the localization function
Localization(S::MPolyUnits) = MPolyLocalizedRing(ambient_ring(S), S)

Localization(S::MPolyComplementOfPrimeIdeal) = MPolyLocalizedRing(ambient_ring(S), S)

Localization(S::MPolyComplementOfKPointIdeal) = MPolyLocalizedRing(ambient_ring(S), S)

Localization(S::MPolyPowersOfElement) = MPolyLocalizedRing(ambient_ring(S), S)

Localization(S::MPolyUnionOfMultSets) = MPolyLocalizedRing(ambient_ring(S), S)

function Localization(
    W::MPolyLocalizedRing{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}}, 
    S::MPolyPowersOfElement{BRT, BRET, RT, RET}
  ) where {BRT, BRET, RT, RET}
  return Localization(S*inverted_set(W))
end

### additional constructors
MPolyLocalizedRing(R::RingType, P::MPolyIdeal{RingElemType}) where {RingType, RingElemType} = MPolyLocalizedRing(R, MPolyComplementOfPrimeIdeal(P))

Localization(R::MPolyRing, f::MPolyElem) = Localization(MPolyPowersOfElement(R, [f]))
Localization(R::MPolyRing, v::Vector{T}) where {T<:MPolyElem} = Localization(MPolyPowersOfElement(R, v))

function Localization(
    W::MPolyLocalizedRing{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}}, 
    f::RET
  ) where {BRT, BRET, RT, RET<:RingElement}
  R = base_ring(W)
  parent(f) == R || error("the given element does not belong to the correct ring")
  S = inverted_set(W)
  if f in S
    return W
  end
  g = denominators(S)
  h = gcd(prod(g), f)
  return MPolyLocalizedRing(R, MPolyPowersOfElement(R, vcat(g, divexact(f, h))))
end

function Localization(
    W::MPolyLocalizedRing{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}}, 
    v::Vector{RET}
  ) where {BRT, BRET, RT, RET}
  V = W
  for f in v
    V = Localization(V, f)
  end
  return V
end

### generation of random elements 
function rand(W::MPolyLocalizedRing, v1::UnitRange{Int}, v2::UnitRange{Int}, v3::UnitRange{Int})
  return W(rand(base_ring(W), v1, v2, v3), rand(inverted_set(W), v1, v2, v3))
end


########################################################################
# Elements of localized polynomial rings                               #
########################################################################

@Markdown.doc """
MPolyLocalizedRingElem{
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

Elements of localizations of polynomial rings.
"""
mutable struct MPolyLocalizedRingElem{
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

  frac::AbstractAlgebra.Generic.Frac{RingElemType}
  R_loc::MPolyLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}

  function MPolyLocalizedRingElem(
      f::AbstractAlgebra.Generic.Frac{RingElemType}, 
      R_loc::MPolyLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}
    ) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}

    base_ring(parent(f)) == base_ring(R_loc) || error(
	"the numerator and denominator of the given fraction do not belong to the original ring before localization"
      )
    denominator(f) in inverted_set(R_loc) || error(
	"the given denominator is not admissible for this localization"
      )
    return new{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}(f, R_loc)
  end
end

### required getter functions 
numerator(a::MPolyLocalizedRingElem) = numerator(a.frac)

denominator(a::MPolyLocalizedRingElem) = denominator(a.frac)

parent(a::MPolyLocalizedRingElem) = a.R_loc

### additional getter functions
fraction(a::MPolyLocalizedRingElem) = a.frac

### required conversions
(W::MPolyLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType})(f::RingElemType) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyLocalizedRingElem((FractionField(base_ring(W)))(f), W)
function (W::MPolyLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType})(a::RingElemType, b::RingElemType) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} 
  return MPolyLocalizedRingElem(a//b, W)
end

### additional conversions
(W::MPolyLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType})(f::AbstractAlgebra.Generic.Frac{RingElemType}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyLocalizedRingElem(f, W)

### additional promotions 
AbstractAlgebra.promote_rule(::Type{RET}, ::Type{MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT<:Ring, RET<:RingElement, MST} = MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}

AbstractAlgebra.promote_rule(::Type{BRET}, ::Type{MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}}) where {BRT<:Ring, BRET<:RingElement, RT<:Ring, RET<:RingElement, MST} = MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}

AbstractAlgebra.promote_rule(::Type{Integer}, ::Type{MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}}) where {BRT<:Ring, BRET<:RingElement, RT<:Ring, RET<:RingElement, MST} = MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}

AbstractAlgebra.promote_rule(::Type{fmpz}, ::Type{MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}}) where {BRT<:Ring, BRET<:RingElement, RT<:Ring, RET<:RingElement, MST} = MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}

### overwriting the arithmetic using the fractions from AbstractAlgebra
function +(a::T, b::T) where {T<:MPolyLocalizedRingElem}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  return (parent(a))(fraction(a) + fraction(b))
end

function -(a::T, b::T) where {T<:MPolyLocalizedRingElem}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  return (parent(a))(fraction(a) - fraction(b))
end

function -(a::T) where {T<:MPolyLocalizedRingElem}
  return (parent(a))((-1)*fraction(a))
end

function *(a::T, b::T) where {T<:MPolyLocalizedRingElem}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  return (parent(a))(fraction(a) * fraction(b))
end

function *(a::RET, b::MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET <: RingElem, MST}
  return (parent(b))(a*fraction(b))
end

function *(a::MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}, b::RET) where {BRT, BRET, RT, RET <: RingElem, MST}
  return b*a
end

function *(a::BRET, b::MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}) where {BRT, BRET <: RingElem, RT, RET, MST}
  return (parent(b))(a*fraction(b))
end

function *(a::MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}, b::BRET) where {BRT, BRET <: RingElem, RT, RET, MST}
  return b*a
end

function Base.:(//)(a::Integer, b::T) where {T<:MPolyLocalizedRingElem}
  return (parent(b))(a//fraction(b))
end

function Base.:(//)(a::fmpz, b::T) where {T<:MPolyLocalizedRingElem}
  return (parent(b))(a//fraction(b))
end

function Base.:(//)(a::T, b::T) where {T<:MPolyLocalizedRingElem}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  g = gcd(numerator(a), numerator(b))
  c = divexact(numerator(a), g)
  d = divexact(numerator(b), g)
  numerator(fraction(b)) in inverted_set(parent(b)) || error("the second argument is not a unit in this local ring")
  return (parent(a))(fraction(a) // fraction(b))
end

function ==(a::T, b::T) where {T<:MPolyLocalizedRingElem}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  return fraction(a) == fraction(b)
end

# We need to manually split this into three methods, because 
# otherwise it seems that Julia can not dispatch this function.
function ^(a::MPolyLocalizedRingElem, i::Int64)
  return parent(a)(fraction(a)^i)
end
function ^(a::MPolyLocalizedRingElem, i::Integer)
  return parent(a)(fraction(a)^i)
end
function ^(a::MPolyLocalizedRingElem, i::fmpz)
  return parent(a)(fraction(a)^i)
end

function divexact(p::T, q::T; check::Bool=false) where {T<:MPolyLocalizedRingElem} 
  W = parent(p)
  S = inverted_set(W)
  parent(q) == W || error("incompatible rings")
  a = numerator(p)
  b = denominator(p)
  c = numerator(q)
  d = denominator(q)
  ### a/b = (m/n)*(c/d) <=> a*d*n = c*b*m.
  # Using the factoriality of polynomial rings, 
  # we can cancel common factors of a*d and c*b 
  # and read off the minimal m and n from that.
  a = a*d
  c = c*b
  g = gcd(a, c)
  m = divexact(a, g)
  n = divexact(c, g)
  if !(n in S) 
    error("not an exact division")
  end
  return W(m, n)
end

########################################################################
# implementation of Oscar's general ring interface                     #
########################################################################

one(W::MPolyLocalizedRing) = W(one(base_ring(W)))
zero(W::MPolyLocalizedRing) = W(zero(base_ring(W)))

elem_type(W::MPolyLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyLocalizedRingElem{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}
elem_type(T::Type{MPolyLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyLocalizedRingElem{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}

parent_type(f::MPolyLocalizedRingElem{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}
parent_type(T::Type{MPolyLocalizedRingElem{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}

(W::MPolyLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType})(f::MPolyLocalizedRingElem{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyLocalizedRingElem(fraction(f), W)

(W::MPolyLocalizedRing)() = zero(W)
(W::MPolyLocalizedRing)(a::Integer) = W(base_ring(W)(a))

isdomain_type(T::Type{MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = true 
isexact_type(T::Type{MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = true

### promotion rules
AbstractAlgebra.promote_rule(::Type{MPolyLocalizedRingElem{RT, RET, MST}}, ::Type{MPolyLocalizedRingElem{RT, RET, MST}}) where {RT<:Ring, RET<:RingElement, MST} = MPolyLocalizedRingElem{RT, RET, MST}

function AbstractAlgebra.promote_rule(::Type{MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}}, ::Type{T}) where {BRT<:Ring, BRET<:RingElement, RT<:Ring, RET<:RingElement, MST, T<:RingElement} 
  AbstractAlgebra.promote_rule(RET, T) == RET && return MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}
  return AbstractAlgebra.promote_rule(BRET, T) == BRET ? MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST} : Union{}
end

########################################################################
# Singular functionality                                               #
########################################################################
@Markdown.doc """
LocalizedBiPolyArray{BRT, BRET, RT, RET, MST}

Main workhorse for binding of ideals in localizations ``R[S⁻¹]`` of 
multivariate polynomial rings ``R = 𝕜[x₁,…,xₙ]`` to Singular. 
"""
mutable struct LocalizedBiPolyArray{BRT, BRET, RT, RET, MST}
  # The generators on the Oscar side
  oscar_gens::Vector{MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}}
  # The numerators of the above fractions as elements in the singular version 
  # of the original ring before localization.
  singular_gens::Singular.sideal
  # The localized ring on the Oscar side.
  oscar_ring::MPolyLocalizedRing{BRT, BRET, RT, RET, MST}
  # The polynomial ring on the singular side
  singular_ring::Singular.PolyRing
  # The ordering used for the above singular ring.
  ordering::Symbol
  # An optional shift vector applied to the polynomials in Oscar when 
  # translating them to the Singular ring. 
  # This is to make local orderings useful for localizations at 𝕜-points.
  shift::Vector{BRET}
  # Flag for caching
  is_groebner_basis::Bool

  function LocalizedBiPolyArray(
      oscar_ring::MPolyLocalizedRing{BRT, BRET, RT, RET, MST},
      oscar_gens::Vector{MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}};
      ordering::Symbol=:degrevlex, 
      shift::Vector{BRET}=Vector{BRET}()
    ) where {BRT, BRET, RT, RET, MST}
    lbpa = new{BRT, BRET, RT, RET, MST}()
    # TODO: Add some sanity checks here
    lbpa.oscar_ring = oscar_ring
    lbpa.oscar_gens = oscar_gens
    lbpa.ordering = ordering
    # fill up the shift vector with zeroes if it is not provided in full length
    for i in (length(shift)+1:nvars(base_ring(oscar_ring)))
      push!(shift, zero(coefficient_ring(base_ring(oscar_ring))))
    end
    lbpa.shift = shift
    lbpa.is_groebner_basis=false
    return lbpa
  end
  
  function LocalizedBiPolyArray(
      oscar_ring::MPolyLocalizedRing{BRT, BRET, RT, RET, MST},
      singular_gens::Singular.sideal; 
      shift::Vector{BRET}=Vector{BRET}(), 
      is_groebner_basis::Bool=false,
      ordering::Symbol=:degrevlex
    ) where {BRT, BRET, RT, RET, MST}
    lbpa = new{BRT, BRET, RT, RET, MST}()
    # TODO: Add some sanity checks here
    lbpa.oscar_ring = oscar_ring
    lbpa.singular_gens = singular_gens
    lbpa.singular_ring = base_ring(singular_gens)
    R = base_ring(oscar_ring)
    lbpa.ordering = Singular.ordering_as_symbol(lbpa.singular_ring)
    k = coefficient_ring(R)
    # fill up the shift vector with zeroes if it is not provided in full length
    for i in (length(shift)+1:nvars(R))
      push!(shift, zero(k))
    end
    lbpa.shift = shift
    inv_shift_hom = AlgebraHomomorphism(R, R, [gen(R, i) - R(shift[i]) for i in (1:nvars(R))])
    lbpa.oscar_gens = [ oscar_ring(y) for y in (inv_shift_hom.(R.([x for x in gens(singular_gens)])))]
    lbpa.is_groebner_basis=is_groebner_basis
    return lbpa
  end
end

oscar_gens(lbpa::LocalizedBiPolyArray) = lbpa.oscar_gens
oscar_ring(lbpa::LocalizedBiPolyArray) = lbpa.oscar_ring
ordering(lbpa::LocalizedBiPolyArray) = lbpa.ordering
shift(lbpa::LocalizedBiPolyArray) = lbpa.shift

function singular_gens(lbpa::LocalizedBiPolyArray)
  singular_assure(lbpa)
  return lbpa.singular_gens
end

function singular_ring(lbpa::LocalizedBiPolyArray)
  singular_assure(lbpa)
  return lbpa.singular_ring
end

function singular_assure(lbpa::LocalizedBiPolyArray)
  if !isdefined(lbpa, :singular_ring)
    lbpa.singular_ring = Singular.PolynomialRing(
	Oscar.singular_ring(base_ring(base_ring(oscar_ring(lbpa)))), 
        [string(a) for a = Nemo.symbols(base_ring(oscar_ring(lbpa)))], 
        ordering = ordering(lbpa), 
        cached = false
      )[1]
  end
  if !isdefined(lbpa, :singular_gens)
    shift_hom = hom(base_ring(oscar_ring(lbpa)), base_ring(oscar_ring(lbpa)), 
        [gen(base_ring(oscar_ring(lbpa)), i) + lbpa.shift[i] for i in (1:nvars(base_ring(oscar_ring(lbpa))))])
    lbpa.singular_gens = Singular.Ideal(lbpa.singular_ring,
	[lbpa.singular_ring(shift_hom(numerator(x))) for x in oscar_gens(lbpa)])
  end
end

function std(lbpa::LocalizedBiPolyArray)
  i = Singular.std(singular_gens(lbpa))
  return LocalizedBiPolyArray(
             oscar_ring(lbpa), 
	     i, 
	     shift=shift(lbpa), 
	     ordering=ordering(lbpa), 
	     is_groebner_basis=true
	   )
end


########################################################################
# Ideals in localizations of multivariate polynomial rings             #
########################################################################

@Markdown.doc """
MPolyLocalizedIdeal{BRT, BRET, RT, RET, MST} <: AbsLocalizedIdeal{RT, RET, MST}

Ideals in localizations of polynomial rings.
"""
mutable struct MPolyLocalizedIdeal{BRT, BRET, RT, RET, MST} <: AbsLocalizedIdeal{RT, RET, MST}
  # the initial set of generators, not to be changed ever!
  gens::Vector{MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}}
  # the ambient ring for this ideal
  W::MPolyLocalizedRing{BRT, BRET, RT, RET, MST} 
  
  # Fields for caching
  groebner_bases::Dict{Symbol, LocalizedBiPolyArray{BRT, BRET, RT, RET, MST}}
  dimension::Int
 
  function MPolyLocalizedIdeal(
      W::MPolyLocalizedRing{BRT, BRET, RT, RET, MST}, 
      gens::Vector{MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}};
      check::Bool=false
    ) where {BRT, BRET, RT, RET, MST}
    R = base_ring(W)
    k = base_ring(R)
    S = inverted_set(W)
    for f in gens
      parent(f) == W || error("generator is not an element of the given ring")
      check && (denominator(f) in S || error("fraction is not an element of the localization"))
    end
    I = new{BRT, BRET, RT, RET, MST}()
    I.gens = gens
    I.W = W
    I.groebner_bases = Dict{Symbol, LocalizedBiPolyArray{BRT, BRET, RT, RET, MST}}()
    return I
  end
end
 
### required getter functions
gens(I::MPolyLocalizedIdeal) = I.gens
base_ring(I::MPolyLocalizedIdeal) = I.W

### additional getter functions
groebner_bases(I::MPolyLocalizedIdeal) = I.groebner_bases

# the default ordering; probably mathematically useless
default_ordering(I::MPolyLocalizedIdeal) = :degrevlex

# specific default orderings for other cases
default_ordering(
    I::MPolyLocalizedIdeal{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}}
  ) where {BRT, BRET, RT, RET} = :degrevlex

default_ordering(
    I::MPolyLocalizedIdeal{BRT, BRET, RT, RET, MPolyComplementOfKPointIdeal{BRT, BRET, RT, RET}}
  ) where {BRT, BRET, RT, RET} = :negdegrevlex

function dim(I::MPolyLocalizedIdeal)
  if isdefined(I,:dimension)
    return I.dimension
  end
  error("not implemented")
end

### Conversion of ideals in the original ring to localized ideals
function (W::MPolyLocalizedRing{BRT, BRET, RT, RET, MST})(I::MPolyIdeal{RET}) where {BRT, BRET, RT, RET, MST}
  return MPolyLocalizedIdeal(W, W.(gens(I)))
end

### required constructors 
function ideal(
    W::MPolyLocalizedRing{BRT, BRET, RT, RET, MST}, 
    f::MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}
  ) where {BRT, BRET, RT, RET, MST}
  return MPolyLocalizedIdeal(W, [f])
end

function ideal(
    W::MPolyLocalizedRing{BRT, BRET, RT, RET, MST}, 
    gens::Vector{MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}}
  ) where {BRT, BRET, RT, RET, MST}
  return MPolyLocalizedIdeal(W, gens)
end

function ideal(
    W::MPolyLocalizedRing{BRT, BRET, RT, RET, MST}, 
    f::RET
  ) where {BRT, BRET, RT, RET, MST}
  return MPolyLocalizedIdeal(W, [W(f)])
end

function ideal(
    W::MPolyLocalizedRing{BRT, BRET, RT, RET, MST}, 
    gens::Vector{RET}
  ) where {BRT, BRET, RT, RET, MST}
  return MPolyLocalizedIdeal(W, W.(gens))
end

### required functionality
function Base.in(
    f::MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}, 
    I::MPolyLocalizedIdeal{BRT, BRET, RT, RET, MST} 
  ) where {BRT, BRET, RT, RET, MST}
  parent(f) == base_ring(I) || return false
  lbpa = groebner_basis(I)
  return iszero(reduce(f, lbpa))
end

function Base.in(
    f::RET, 
    I::MPolyLocalizedIdeal{BRT, BRET, RT, RET, MST} 
  ) where {BRT, BRET, RT, RET, MST}
  return base_ring(I)(f) in I
end

### Default constructors 
# The ordering is determined from the type of the multiplicative set
LocalizedBiPolyArray(I::MPolyLocalizedIdeal{BRT, BRET, RT, RET, MPolyComplementOfKPointIdeal{BRT, BRET, RT, RET}}) where {BRT, BRET, RT, RET} = LocalizedBiPolyArray(base_ring(I), gens(I), ordering=:negdegrevlex, shift=point_coordinates(inverted_set(base_ring(I)))) 

LocalizedBiPolyArray(I::MPolyLocalizedIdeal{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}}) where {BRT, BRET, RT, RET} = LocalizedBiPolyArray(base_ring(I), gens(I), ordering=:degrevlex) 

function Base.show(io::IO, I::MPolyLocalizedIdeal) 
  print(io, "ideal in $(base_ring(I)) generated by the elements ") 
  n = length(gens(I))
  for i in (1:n-1)
    print(io, "$(gens(I)[i]), ")
  end
  print(io, last(gens(I)))
end

########################################################################
# Groebner and standard bases                                          #
########################################################################

### the catchall implementation; most probably mathematically useless!
function groebner_basis(
    I::MPolyLocalizedIdeal,
    ordering::Symbol
  )
  D = groebner_bases(I)
  # check whether a standard basis has already been computed for this ordering
  if haskey(D, ordering)
    return D[ordering]
  end
  # if not, set up a LocalizedBiPolyArray
  W = base_ring(I)
  R = base_ring(W)
  S = inverted_set(W)
  lbpa = LocalizedBiPolyArray(W, gens(I), ordering=ordering)
  # compute the standard basis and cache the result
  D[ordering] = std(lbpa)
  return D[ordering]
end

function groebner_basis(I::MPolyLocalizedIdeal)
  return groebner_basis(I, default_ordering(I))
end

function groebner_assure(I::MPolyLocalizedIdeal)
  D = groebner_bases(I)
  if length(D) > 0 
    return
  end
  W = base_ring(I)
  R = base_ring(W)
  S = inverted_set(W)
  lbpa = LocalizedBiPolyArray(I)
  D[default_ordering(I)] = std(lbpa)
  return
end

### special routines for localizations at 𝕜-points
function groebner_basis(
    I::MPolyLocalizedIdeal{BRT, BRET, RT, RET, MPolyComplementOfKPointIdeal{BRT, BRET, RT, RET}}; 
    ordering::Symbol=:negdegrevlex
  ) where {BRT, BRET, RT, RET}
  D = groebner_bases(I)
  # check whether a standard basis has already been computed for this ordering
  if haskey(D, ordering)
    return D[ordering]
  end
  # if not, set up a LocalizedBiPolyArray
  W = base_ring(I)
  R = base_ring(W)
  S = inverted_set(W)::MPolyComplementOfKPointIdeal{BRT, BRET, RT, RET}
  a = point_coordinates(S)
  lbpa = LocalizedBiPolyArray(W, gens(I), ordering=ordering, shift=a)
  # Check whether this ordering is admissible
  !Singular.has_local_ordering(singular_ring(lbpa)) && error("The ordering has to be a local ordering.")
  # compute the standard basis and cache the result
  D[ordering] = std(lbpa)
  return D[ordering]
end

function groebner_assure(
    I::MPolyLocalizedIdeal{BRT, BRET, RT, RET, MPolyComplementOfKPointIdeal{BRT, BRET, RT, RET}}
  ) where {BRT, BRET, RT, RET}
  D = groebner_bases(I)
  if length(D) > 0 
    return
  end
  W = base_ring(I)
  R = base_ring(W)
  S = inverted_set(W)::MPolyComplementOfKPointIdeal{BRT, BRET, RT, RET}
  a = point_coordinates(S)
  # Choose negdegrevlex as a default local ordering
  lbpa = LocalizedBiPolyArray(W, gens(I), ordering=:negdegrevlex, shift=a)
  D[:negdegrevlex] = std(lbpa)
  return
end
  
### special routines for localizations at powers of elements
function groebner_basis(
    I::MPolyLocalizedIdeal{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}}; 
    ordering::Symbol=:degrevlex
  ) where {BRT, BRET, RT, RET}
  D = groebner_bases(I)
  # check whether a standard basis has already been computed for this ordering
  if haskey(D, ordering)
    return D[ordering]
  end
  # if not, set up a LocalizedBiPolyArray
  W = base_ring(I)
  R = base_ring(W)
  S = inverted_set(W)::MPolyPowersOfElement{BRT, BRET, RT, RET}
  lbpa = LocalizedBiPolyArray(W, gens(I), ordering=ordering)
  a = denominators(S)
  sing_ring = singular_ring(lbpa)
  sing_a = [sing_ring(x) for x in a]
  sing_ideal = singular_gens(lbpa)
  for h in sing_a
    sing_ideal = Singular.saturation(sing_ideal, Singular.Ideal(sing_ring, [h]))[1]
  end
  sing_ideal = Singular.std(sing_ideal)
  lbpa = LocalizedBiPolyArray(W, sing_ideal, ordering=ordering)
  return lbpa
end


### reduction to normal form with respect to a list of elements
function Base.reduce(
    f::MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}, 
    lbpa::LocalizedBiPolyArray{BRT, BRET, RT, RET, MST}
  ) where {BRT, BRET, RT, RET, MST}

  W = parent(f)
  W == oscar_ring(lbpa) || error("element does not belong to the Oscar ring of the biPolyArray")
  R = base_ring(W)
  shift_hom = hom(R, R, [gen(R, i) + lbpa.shift[i] for i in (1:nvars(R))])
  singular_n = singular_ring(lbpa)(shift_hom(numerator(f)))
  singular_n = Singular.reduce(singular_n, singular_gens(lbpa))
  if iszero(singular_n) 
    return zero(W)
  end
  #singular_d = singular_ring(lbpa)(shift_hom(denominator(f)))
  #singular_d = Singular.reduce(singular_d, singular_gens(lbpa))
  inv_shift_hom = hom(R,R, [gen(R, i) - lbpa.shift[i] for i in (1:nvars(R))])
  #return W(inv_shift_hom(R(singular_n)), inv_shift_hom(R(singular_d)))
  return W(inv_shift_hom(R(singular_n)), denominator(f))
end


########################################################################
# Homomorphisms of localized polynomial rings                          #
########################################################################
#
# Let P = 𝕜[x₁,…,xₘ] and Q = 𝕜[y₁,…,yₙ] be polynomial rings 
# Any homomorphism ϕ : P[U⁻¹] → Q[V⁻¹] is completely determined 
# by the images of the variables 
#
#     ϕ(xᵢ) = aᵢ(y)/ bᵢ(y)
#
# where bᵢ(y) ∈ V for all i. Given such a list of images
# one can always define the associated homomorphism ϕ' : P → Q[V⁻¹]. 
# This extends to a well defined homomorphism ϕ as above iff
#
#     for all f ∈ U : ϕ'(f) is a unit in Q[V⁻¹]            (*)
#
# and this is easily seen to be the case iff there exists an 
# element c ∈ V such that c⋅ϕ'(f) ∈ V

mutable struct MPolyLocalizedRingHom{BaseRingType, BaseRingElemType, 
    RingType, RingElemType, DomainMultSetType, CodomainMultSetType
  } <: AbsLocalizedRingHom{
    RingType, RingElemType, DomainMultSetType, CodomainMultSetType
  }
  # the domain of definition
  V::MPolyLocalizedRing{BaseRingType, BaseRingElemType, 
			RingType, RingElemType, DomainMultSetType}
  # the codomain 
  W::MPolyLocalizedRing{BaseRingType, BaseRingElemType, 
			RingType, RingElemType, CodomainMultSetType}
  # the images of the variables
  images::Vector{MPolyLocalizedRingElem{
    BaseRingType, BaseRingElemType, RingType, RingElemType, 
    CodomainMultSetType}}

  function MPolyLocalizedRingHom(
      V::MPolyLocalizedRing{BRT, BRET, RT, RET, DMST}, 
      W::MPolyLocalizedRing{BRT, BRET, RT, RET, CMST}, 
      a::Vector{MPolyLocalizedRingElem{BRT, BRET, RT, RET, CMST}}
    ) where {BRT, BRET, RT, RET, CMST, DMST}
    R = base_ring(V)
    S = base_ring(W)
    k = coefficient_ring(R) 
    k == coefficient_ring(S) || error("the two polynomial rings are not defined over the same coefficient ring")
    ngens(R) == length(a) || error("the number of images does not coincide with the number of variables")
    parent_check = true
    for x in a
      parent_check = parent_check && parent(x) == W
      denominator(x) in inverted_set(W) || error("the homomorphism is not well defined")
    end
    parent_check || error("the images of the variables are not elements of the codomain")
    # Check whether this homomorphism is well defined
    # TODO: Implement that!
    return new{typeof(k), elem_type(k), typeof(R), elem_type(R), typeof(inverted_set(V)), typeof(inverted_set(W))}(V, W, a)
  end
end
  
### additional constructors
function MPolyLocalizedRingHom(
      R::RT,
      W::MPolyLocalizedRing{BRT, BRET, RT, RET, CMST}, 
      a::Vector{MPolyLocalizedRingElem{BRT, BRET, RT, RET, CMST}}
    ) where {BRT, BRET, RT, RET, CMST}
  return MPolyLocalizedRingHom(Localization(MPolyUnits(R)), W, a)
end

function MPolyLocalizedRingHom(
      V::MPolyLocalizedRing{BRT, BRET, RT, RET, DMST}, 
      S::RT,
      a::Vector{RET}
    ) where {BRT, BRET, RT, RET, DMST}
  W = Localization(MPolyUnits(S))
  return MPolyLocalizedRingHom(V, W, W.(a))
end

### required getter functions
domain(f::MPolyLocalizedRingHom) = f.V
codomain(f::MPolyLocalizedRingHom) = f.W
images(f::MPolyLocalizedRingHom) = f.images

### required functionality
function (f::MPolyLocalizedRingHom{BRT, BRET, RT, RET, DMST, CMST})(p::MPolyLocalizedRingElem{BRT, BRET, RT, RET, DMST}) where {BRT, BRET, RT, RET, DMST, CMST}
  parent(p) == domain(f) || error("the given element does not belong to the domain of the map")
  return evaluate(numerator(p), images(f))//evaluate(denominator(p), images(f))
end

### overwriting of the generic method
function (f::MPolyLocalizedRingHom{BRT, BRET, RT, RET, DMST, CMST})(p::RET) where {BRT, BRET, RT, RET, DMST, CMST}
  parent(p) == base_ring(domain(f)) || error("the given element does not belong to the domain of the map")
  return evaluate(p, images(f))
end

### provide an extra method for elements of the base ring
function (f::MPolyLocalizedRingHom{BRT, BRET, RT, RET, DMST, CMST})(p::BRET) where {BRT, BRET, RT, RET, DMST, CMST}
  parent(p) == coefficient_ring(base_ring(domain(f))) || error("the given element does not belong to the domain of the map")
  return codomain(f)(p)
end

### remove the ambiguity of methods in case the base ring is ZZ
function (f::MPolyLocalizedRingHom)(p::fmpz) 
  return codomain(f)(p)
end

### implementing the Oscar map interface
identity_map(W::T) where {T<:MPolyLocalizedRing} = MPolyLocalizedRingHom(W, W, W.(gens(base_ring(W))))
function compose(
    f::MPolyLocalizedRingHom{BRT, BRET, RT, RET, MST1, MST2}, 
    g::MPolyLocalizedRingHom{BRT, BRET, RT, RET, MST2, MST3}
  ) where {BRT, BRET, RT, RET, MST1, MST2, MST3}
  codomain(f) == domain(g) || error("maps are not compatible")
  return MPolyLocalizedRingHom(domain(f), codomain(g), g.(images(f)))
end

function preimage(f::MPolyLocalizedRingHom, I::MPolyLocalizedIdeal)
  base_ring(I) == codomain(f) || error("the ideal does not belong to the codomain of the map")
  R = base_ring(domain(f))
  error("not implemented")
end
