export MPolyComplementOfPrimeIdeal, MPolyComplementOfKPointIdeal, MPolyPowersOfElement, MPolyProductOfMultSets, MPolyLeadingMonOne
export prime_ideal
export rand, sets, issubset, units_of, simplify!, is_trivial
import Base: issubset

export MPolyLocalizedRing
export ambient_ring, point_coordinates, inverted_set, denominators, gens

export MPolyLocalizedRingElem
export numerator, denominator, fraction, parent, is_unit, divexact
export reduce_fraction

export MPolyLocalizedIdeal
export gens, base_ring, saturated_ideal, intersect, quotient

export Localization, ideal
export bring_to_common_denominator, write_as_linear_combination

export MPolyLocalizedRingHom
export domain, codomain, images

import AbstractAlgebra: Ring, RingElem

export complement_of_ideal, powers_of_element

########################################################################
# General framework for localizations of multivariate polynomial rings #
########################################################################

########################################################################
# Multiplicatively closed sets in multivariate polynomial rings        #
########################################################################

abstract type AbsMPolyMultSet{BRT, BRET, RT, RET} <: AbsMultSet{RT, RET} end

########################################################################
# Powers of elements                                                   #
########################################################################

@Markdown.doc """
    MPolyPowersOfElement{
        BaseRingType,
        BaseRingElemType, 
        RingType,
        RingElemType
      } <: AbsMPolyMultSet{
        BaseRingType,
        BaseRingElemType, 
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
  } <: AbsMPolyMultSet{
    BaseRingType,
    BaseRingElemType, 
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
MPolyPowersOfElement(f::RET) where {RET<:MPolyElem} = MPolyPowersOfElement(parent(f), [f])
units_of(R::RT) where {RT<:MPolyRing} = MPolyPowersOfElement(R, [one(R)])

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
  #return divides(one(R), o)[1]
  return is_unit(o)
end

### iteration 
Base.iterate(U::MPolyPowersOfElement) = (length(U.a)>0 ? (U.a[1], 1) : nothing)
Base.iterate(U::MPolyPowersOfElement, a::Tuple{<:MPolyElem, Int}) = (a[2] < length(U.a) ? (U.a[a[2]+1], a[2]+1) : nothing)
Base.iterate(U::MPolyPowersOfElement, i::Int) = (i < length(U.a) ? (U.a[i+1], i+1) : nothing)

is_trivial(U::MPolyPowersOfElement) = (U == units_of(ambient_ring(U)))

### printing
function Base.show(io::IO, S::MPolyPowersOfElement)
  print(io, "powers of ")
  print(io, denominators(S))
end

### generation of random elements 
function rand(S::MPolyPowersOfElement, v1::UnitRange{Int}, v2::UnitRange{Int}, v3::UnitRange{Int})
  R = ambient_ring(S)
  return prod(f^rand(0:5) for f in denominators(S); init = one(R))::elem_type(R)
end

### simplification. 
# Replaces each element d by its list of square free prime divisors.
function simplify!(S::MPolyPowersOfElement)
  R = ambient_ring(S)
  new_denom = Vector{elem_type(R)}()
  for d in denominators(S)
    for a in factor(d)
      push!(new_denom, a[1])
    end
  end
  S.a = new_denom
  return S
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
      } <: AbsMPolyMultSet{
        BaseRingType, 
        BaseRingElemType,
        RingType,
        RingElemType
      }

The complement of a prime ideal ``P ⊂ 𝕜[x₁,…,xₙ]`` in a multivariate polynomial ring 
with elements of type `RingElemType` over a base ring ``𝕜`` of type `BaseRingType`.
"""
mutable struct MPolyComplementOfPrimeIdeal{
    BaseRingType, 
    BaseRingElemType,
    RingType,
    RingElemType
  } <: AbsMPolyMultSet{
    BaseRingType, 
    BaseRingElemType,
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
    check && (is_prime(P) || error("the ideal $P is not prime"))
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
      } <: AbsMPolyMultSet{
        BaseRingType,
        BaseRingElemType, 
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
  } <: AbsMPolyMultSet{
    BaseRingType,
    BaseRingElemType, 
    RingType, 
    RingElemType
  }

  # The parent polynomial ring 𝕜[x₁,…,xₙ]
  R::RingType
  # The coordinates aᵢ of the point in 𝕜ⁿ corresponding to the maximal ideal
  a::Vector{BaseRingElemType}

  function MPolyComplementOfKPointIdeal(R::RingType, a::Vector{T}) where {RingType<:MPolyRing, T<:RingElement}
    length(a) == length(gens(R)) || error("the number of variables in the ring does not coincide with the number of coordinates")
    n = length(a)
    kk = coefficient_ring(R)
    b = kk.(a) # fails if the input is not compatible
    S = new{typeof(kk), elem_type(kk), RingType, elem_type(R)}(R, b)
    return S
  end
end

@doc Markdown.doc"""  
    complement_of_ideal(R::MPolyRing, a::Vector)

Given a polynomial ring ``R``, say ``R = K[x_1,\dots, x_n]``, and given a vector 
``a = (a_1, \dots, a_n)`` of ``n`` elements of ``K``, return the multiplicatively 
closed subset ``R\setminus m``, where ``m`` is the maximal ideal 

$$m = \langle x_1-a_1,\dots, x_n-a_n\rangle \subset R.$$

    complement_of_ideal(P::MPolyIdeal; check::Bool=false)

Given a prime ideal ``P`` of a polynomial ring ``R``, say,
return the multiplicatively closed subset ``R\setminus P.``

!!! note
    If  `check` is set to `true`, the function checks whether ``P`` is indeed a prime ideal. 

    This may take some time.

# Examples
```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"]);

julia> U = complement_of_ideal(R, [0, 0 ,0])
complement of maximal ideal corresponding to point with coordinates fmpq[0, 0, 0]

julia> P = ideal(R, [x])
ideal(x)

julia> U = complement_of_ideal(P)
complement of ideal(x)
```
"""
complement_of_ideal(R::MPolyRing, a::Vector) = MPolyComplementOfKPointIdeal(R, a)
complement_of_ideal(P::MPolyIdeal; check::Bool=false) = MPolyComplementOfPrimeIdeal(P; check)

@doc Markdown.doc"""  
    powers_of_element(f::MPolyElem)

Given an element `f` of a polynomial ring, return the multiplicatively 
closed subset of the polynomial ring which is formed by the powers of `f`.

# Examples
```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"]);

julia> f = x
x

julia> U = powers_of_element(f)
powers of fmpq_mpoly[x]
```
"""
powers_of_element(f::MPolyElem) = MPolyPowersOfElement(f)

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
  print(io, "complement of maximal ideal corresponding to point with coordinates ")
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


@Markdown.doc """
    MPolyProductOfMultSets{
        BaseRingType,
        BaseRingElemType, 
        RingType,
        RingElemType
      } <: AbsMPolyMultSet{
        BaseRingType,
        BaseRingElemType, 
        RingType, 
        RingElemType
      }

A finite product `T⋅U = { a⋅b : a ∈ T, b∈ U}` of arbitrary other 
multiplicative sets in a multivariate polynomial ring.
"""
mutable struct MPolyProductOfMultSets{
    BaseRingType,
    BaseRingElemType, 
    RingType,
    RingElemType
  } <: AbsMPolyMultSet{
    BaseRingType,
    BaseRingElemType, 
    RingType, 
    RingElemType
  }
  R::RingType
  U::Vector{<:AbsMPolyMultSet{BaseRingType, BaseRingElemType, RingType, RingElemType}}

  function MPolyProductOfMultSets(R::RT, U::Vector{<:AbsMPolyMultSet{BRT, BRET, RT, RET}}) where {BRT<:Ring, BRET<:RingElement, RT<:MPolyRing, RET<:MPolyElem}
    for s in U
      ambient_ring(s) == R || error("multiplicative set does not live in the given ring")
    end
    return new{typeof(coefficient_ring(R)), elem_type(coefficient_ring(R)), typeof(R), elem_type(R)}(R, U)
  end
end

### required getter functions
ambient_ring(S::MPolyProductOfMultSets) = S.R

### additional functionality
getindex(S::MPolyProductOfMultSets, i::Integer) = S.U[i]
sets(S::MPolyProductOfMultSets) = copy(S.U)

### required functionality
function Base.in(
    f::RingElemType, 
    S::MPolyProductOfMultSets{BaseRingType, BaseRingElemType, RingType, RingElemType}
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
function Base.show(io::IO, S::MPolyProductOfMultSets)
  print(io, "product of the multiplicative sets [")
  for s in sets(S)
    print(io, s)
    s != last(sets(S)) && print(io, ", ")
  end
  print(io, "]")
end

### generation of random elements 
function rand(S::MPolyProductOfMultSets, v1::UnitRange{Int}, v2::UnitRange{Int}, v3::UnitRange{Int})
  return prod([rand(s, v1, v2, v3) for s in sets(S)])::elem_type(ambient_ring(S))
end

########################################################################
# Localization associated to a monomial ordering                       #
########################################################################

@Markdown.doc """
    MPolyLeadingMonOne{
        BaseRingType,
        BaseRingElemType,
        RingType,
        RingElemType
      } <: AbsMPolyMultSet{
        BaseRingType,
        BaseRingElemType,
        RingType,
        RingElemType
      }

The set `S = { a in R : leading_monomial(a, ord) = 1 }` for a fixed
monomial ordering `ord`.
"""
mutable struct MPolyLeadingMonOne{
    BaseRingType,
    BaseRingElemType,
    RingType,
    RingElemType
  } <: AbsMPolyMultSet{
    BaseRingType,
    BaseRingElemType,
    RingType,
    RingElemType
  }

  R::RingType # the parent ring
  ord::MonomialOrdering

  function MPolyLeadingMonOne(R::RingType, ord::MonomialOrdering) where {RingType <: MPolyRing}
    @assert R === ord.R "Ordering does not belong to the given ring"
    k = coefficient_ring(R)
    return new{typeof(k), elem_type(k), RingType, elem_type(R)}(R, ord)
  end
end

### required getter functions
ambient_ring(S::MPolyLeadingMonOne) = S.R

### additional constructors
MPolyLeadingMonOne(ord::MonomialOrdering) = MPolyLeadingMonOne(ord.R, ord)

ordering(S::MPolyLeadingMonOne) = S.ord

### required functionality
function Base.in(
    f::RingElemType,
    S::MPolyLeadingMonOne{BaseRingType, BaseRingElemType, RingType, RingElemType}
  ) where {BaseRingType, BaseRingElemType, RingType, RingElemType}
  R = parent(f)
  R == ambient_ring(S) || return false
  if iszero(f)
    return false
  end
  return isone(leading_monomial(f; ordering = ordering(S)))
end

### printing
function Base.show(io::IO, S::MPolyLeadingMonOne)
  print(io, "elements of ")
  print(io, ambient_ring(S))
  print(io, " with leading monomial 1 w.r.t. ")
  print(io, ordering(S))
end

########################################################################
# Arithmetic of multiplicative sets                                    #
########################################################################

### containment ########################################################
⊂(T::AbsMPolyMultSet, U::AbsMPolyMultSet) = issubset(T, U)

function issubset(T::AbsMPolyMultSet, U::AbsMPolyMultSet) 
  ambient_ring(T) == ambient_ring(U) || return false
  error("comparison of multiplicative sets of type $(typeof(T)) and $(typeof(U)) is not implemented")
end

function ==(T::AbsMPolyMultSet, U::AbsMPolyMultSet) 
  return (issubset(T, U) && issubset(U, T))
end

function issubset(
    T::MPolyComplementOfPrimeIdeal{BRT, BRET, RT, RET},
    U::MPolyComplementOfPrimeIdeal{BRT, BRET, RT, RET}
  ) where {BRT, BRET, RT, RET}
  return issubset(prime_ideal(U), prime_ideal(T))
end

function issubset(
    T::MPolyComplementOfKPointIdeal{BRT, BRET, RT, RET},
    U::MPolyComplementOfKPointIdeal{BRT, BRET, RT, RET}
  ) where {BRT, BRET, RT, RET}
  R = ambient_ring(T)
  R == ambient_ring(U) || error("multiplicative sets do not belong to the same ring")
  a = point_coordinates(U)
  b = point_coordinates(T)
  for i in 1:length(a)
    a[i] == b[i] || return false
  end
  return true
end

function issubset(
    T::MPolyComplementOfKPointIdeal{BRT, BRET, RT, RET},
    U::MPolyComplementOfPrimeIdeal{BRT, BRET, RT, RET}
  ) where {BRT, BRET, RT, RET}
  R = ambient_ring(T)
  R == ambient_ring(U) || error("multiplicative sets do not belong to the same ring")
  a = point_coordinates(T)
  for i in 1:length(a)
    (gens(R)[i]- R(a[i])) in prime_ideal(U) || return false
  end
  return true
end

function issubset(
    T::MPolyComplementOfPrimeIdeal{BRT, BRET, RT, RET},
    U::MPolyComplementOfKPointIdeal{BRT, BRET, RT, RET}
  ) where {BRT, BRET, RT, RET}
  R = ambient_ring(T)
  R == ambient_ring(U) || error("multiplicative sets do not belong to the same ring")
  a = point_coordinates(U)
  for f in gens(prime_ideal(T))
    iszero(evaluate(f, a)) || return false
  end
  return true
end

function issubset(
    T::MPolyPowersOfElement{BRT, BRET, RT, RET},
    U::AbsMPolyMultSet{BRT, BRET, RT, RET}
  ) where {BRT, BRET, RT, RET}
  for a in denominators(T)
    a in U || return false
  end
  return true
end

function issubset(
    T::MPolyComplementOfPrimeIdeal{BRT, BRET, RT, RET},
    U::MPolyPowersOfElement{BRT, BRET, RT, RET},
  ) where {BRT, BRET, RT, RET}
  return false
end

function issubset(
    T::MPolyComplementOfKPointIdeal{BRT, BRET, RT, RET},
    U::MPolyPowersOfElement{BRT, BRET, RT, RET},
  ) where {BRT, BRET, RT, RET}
  return false
end

function issubset(
    T::MPolyProductOfMultSets{BRT, BRET, RT, RET},
    U::MPolyProductOfMultSets{BRT, BRET, RT, RET}
  ) where {BRT, BRET, RT, RET}
  for V in sets(T)
    issubset(V, U) || return false
  end
  return true
end

function issubset(
    T::MPolyProductOfMultSets{BRT, BRET, RT, RET},
    U::AbsMPolyMultSet{BRT, BRET, RT, RET}
  ) where {BRT, BRET, RT, RET}
  for V in sets(T)
    issubset(V, U) || return false
  end
  return true
end

function issubset(
    T::AbsMPolyMultSet{BRT, BRET, RT, RET},
    U::MPolyProductOfMultSets{BRT, BRET, RT, RET}
  ) where {BRT, BRET, RT, RET}
  for V in sets(U)
    issubset(T, V) && return true
  end
  error("containment can not be checked")
end

function issubset(
    T::MPolyPowersOfElement{BRT, BRET, RT, RET},
    U::MPolyProductOfMultSets{BRT, BRET, RT, RET}
  ) where {BRT, BRET, RT, RET}
  for d in denominators(T)
    d in U || return false
  end
  return true
end

### intersections ######################################################
function intersect(T::AbsMPolyMultSet, U::AbsMPolyMultSet) 
  error("intersection of multiplicative sets of type $(typeof(T)) and $(typeof(U)) is not implemented")
end

# TODO: Implement this if necessary!

### functionality for taking products
#
# Definition.
# Let T and U be multiplicative sets in a commutative ring R. The product 
# of T and U is defined as 
#
#   T⋅U = { f⋅g : f ∈ T and g ∈ U }.
#
# A product of multiplicative sets U = U₁⋅…⋅Uₙ is called interreduced 
# if neither one of the factors Uᵢ is contained in one of the others Uⱼ, j≠i.
#
# Lemma. 
# Any product of multiplicative sets U = U₁⋅…⋅Uₙ may be replaced by 
# an interreduced one. 
#
# Remark. 
# An interreduced factorization of a product of multiplicative sets may 
# not be unique: Consider the ring ℤ[x] and the multiplicative sets 
#   T  = {powers of 5x}
#   T' = {powers of x}
#   S  = {constant polynomials outside 7ℤ }.
# Then 
#   T⋅S = { a⋅xᵏ : a ∉ 7ℤ, k ∈ ℕ₀ } = T'⋅S.
#
# Upshot: Whenever a product is taken, some interreduced form of the 
# entire product is returned. Besides the obvious simplification in 
# case all factors are contained in a single one, it is difficult to 
# determine which interreduction is the best one. 

*(T::AbsMPolyMultSet, U::AbsMPolyMultSet) = product(T, U)

@Markdown.doc """
    product(T::AbsMPolyMultSet, U::AbsMPolyMultSet)

Return the product of the multiplicative subsets `T` and `U`. 

Alternatively, write `T*U`.

# Examples
```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"]);

julia> T = complement_of_ideal(R, [0, 0 ,0])
complement of maximal ideal corresponding to point with coordinates fmpq[0, 0, 0]

julia> f = x
x

julia> U = powers_of_element(f)
powers of fmpq_mpoly[x]

julia> S = product(T, U)
product of the multiplicative sets [complement of maximal ideal corresponding to point with coordinates fmpq[0, 0, 0], powers of fmpq_mpoly[x]]
```
"""
function product(T::AbsMPolyMultSet, U::AbsMPolyMultSet)
  R = ambient_ring(T)
  R == ambient_ring(U) || error("multiplicative sets do not belong to the same ring")
  issubset(T, U) && return U
  issubset(U, T) && return T
  return MPolyProductOfMultSets(R, [T, U])
end

function product(T::MST, U::MST) where {MST<:MPolyProductOfMultSets} 
  R = ambient_ring(T)
  R == ambient_ring(U) || error("multiplicative sets do not belong to the same ring")
  new_sets = Vector()
  for S in sets(T)
    push!(new_sets, S)
    for V in sets(U)
      if issubset(S, V) 
	pop!(new_sets)
	break
      end
    end
  end
  n = length(new_sets)
  for V in sets(U)
    push!(new_sets, V)
    for U in new_sets[1:n]
      if issubset(V, U) 
	pop!(new_sets)
	break
      end
    end
  end
  return MPolyProductOfMultSets(R, [x for x in new_sets])
end

function product(T::MPolyProductOfMultSets{BRT, BRET, RT, RET}, U::MST) where {BRT, BRET, RT, RET, MST<:AbsMPolyMultSet{BRT, BRET, RT, RET}}
  R = ambient_ring(T)
  R == ambient_ring(U) || error("multiplicative sets do not belong to the same ring")
  for V in sets(T)
    issubset(U, T) && return T
  end
  new_sets = U
  for V in sets(T)
    issubset(V, U) || push!(new_sets, V)
  end
  return MPolyProductOfMultSets(R, new_sets)
end

product(U::MST, T::MPolyProductOfMultSets{BRT, BRET, RT, RET}) where {BRT, BRET, RT, RET, MST<:AbsMPolyMultSet{BRT, BRET, RT, RET}} = product(T, U)

function product(T::MPolyComplementOfPrimeIdeal{BRT, BRET, RT, RET}, U::MPolyComplementOfKPointIdeal{BRT, BRET, RT, RET}) where {BRT, BRET, RT, RET}
  R = ambient_ring(T)
  R == ambient_ring(U) || error("multiplicative sets do not belong to the same ring")
  P = prime_ideal(T)
  for f in gens(P)
    if iszero(evaluate(f, point_coordinates(U)))
      return MPolyProductOfMultSets(R, [U, T])
    end
  end
  return U
end

function product(
    U::MPolyComplementOfKPointIdeal{BRT, BRET, RT, RET},
    T::MPolyComplementOfPrimeIdeal{BRT, BRET, RT, RET}
  ) where {BRT, BRET, RT, RET}
  return product(T, U)
end

function product(T::MST, U::MST) where {MST<:MPolyComplementOfKPointIdeal}
  R = ambient_ring(T)
  R == ambient_ring(U) || error("multiplicative sets do not belong to the same ring")
  a = point_coordinates(U)
  b = point_coordinates(T)
  for i in 1:length(a)
    a[i] == b[i] || return MPolyProductOfMultSets(R, [U, T])
  end
  return T
end

function product(T::MST, U::MST) where {MST<:MPolyComplementOfPrimeIdeal}
  R = ambient_ring(T)
  R == ambient_ring(U) || error("multiplicative sets do not belong to the same ring")
  if issubset(prime_ideal(T), prime_ideal(U))
    return T
  end
  if issubset(prime_ideal(U), prime_ideal(T))
    return U
  end
  return MPolyProductOfMultSets(R, [U, T])
end

function product(T::MST, U::MST) where {MST<:MPolyPowersOfElement}
  R = ambient_ring(T) 
  R == ambient_ring(U) || error("multiplicative sets do not belong to the same ring")
  new_denoms = Vector{elem_type(R)}()
  for f in denominators(T)
    for g in denominators(U)
      (_, f) = ppio(f, g)
    end
    if !(divides(one(parent(f)), f)[1])
      push!(new_denoms, f)
    end
  end
  n = length(new_denoms)
  for g in denominators(U)
    for f in new_denoms[1:n]
      (_, g) = ppio(g, f)
    end
    if !(divides(one(parent(g)), g)[1])
      push!(new_denoms, g)
    end
  end
  return (length(new_denoms) == 0 ? units_of(R) : MPolyPowersOfElement(R, new_denoms))
end

function product(
    U::AbsMPolyMultSet{BRT, BRET, RT, RET},
    T::MPolyPowersOfElement{BRT, BRET, RT, RET}
  ) where {BRT, BRET, RT, RET}
  R = ambient_ring(T)
  R == ambient_ring(U) || error("multiplicative sets do not belong to the same ring")
  for a in denominators(T)
    a in U || return MPolyProductOfMultSets(R, [U, T])
  end
  return U
end

function product(
    T::MPolyPowersOfElement{BRT, BRET, RT, RET},
    U::AbsMPolyMultSet{BRT, BRET, RT, RET}
  ) where {BRT, BRET, RT, RET}
  R = ambient_ring(T)
  R == ambient_ring(U) || error("multiplicative sets do not belong to the same ring")
  for a in denominators(T)
    a in U || return MPolyProductOfMultSets(R, [U, T])
  end
  return U
end

function product(
    T::MPolyPowersOfElement{BRT, BRET, RT, RET},
    U::MPolyProductOfMultSets{BRT, BRET, RT, RET}
  ) where {BRT, BRET, RT, RET}
  R = ambient_ring(T)
  R == ambient_ring(U) || error("multiplicative sets do not belong to the same ring")
  keep_denom = Vector{RET}()
  for a in denominators(T)
    a in U || (push!(keep_denom, a))
  end
  length(keep_denom) == 0 && return U
  return MPolyProductOfMultSets(vcat(sets(U), MPolyPowersOfElement(keep_denom)))
end
function product(
    U::MPolyProductOfMultSets{BRT, BRET, RT, RET},
    T::MPolyPowersOfElement{BRT, BRET, RT, RET}
  ) where {BRT, BRET, RT, RET} 
  return T*U
end

### Preimages of multiplicative sets.
# In general for a ring homomorphism f : R → S and a multiplicative 
# set U ⊂ S the preimage V = f⁻¹(U) is a multiplicative set in R. 
# Membership in V can easily be tested, but introducing a new type 
# for preimages makes it necessary to extend all dispatch routines. 
# It is not clear what is the best strategy for all this. 

function preimage(f::Oscar.AffAlgHom, U::MST) where {MST<:AbsMPolyMultSet}
  error("not implemented")
end

### Transfer of multiplicative sets along ring homomorphisms 
function (phi::MPolyAnyMap{<:MPolyRing, <:MPolyRing, Nothing})(U::MPolyPowersOfElement;
                                                               check::Bool=true
                                                              )
  ambient_ring(U) === domain(phi) || error("multiplicative set does not lay in the domain of the morphism")
  S = codomain(phi) 
  SU = MPolyPowersOfElement(S, phi.(denominators(U)))
  return SU
end

function (phi::MPolyAnyMap{<:MPolyRing, <:MPolyRing, Nothing})(U::MPolyComplementOfPrimeIdeal;
                                                               check::Bool=true
                                                              )
  ambient_ring(U) === domain(phi) || error("multiplicative set does not lay in the domain of the morphism")
  S = codomain(phi) 
  Q = ideal(S, phi.(gens(prime_ideal(U))))
  SU = MPolyComplementOfPrimeIdeal(S, Q, check=check)
  return SU
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
@attributes mutable struct MPolyLocalizedRing{
    BaseRingType,
    BaseRingElemType,
    RingType,
    RingElemType,
    MultSetType <: AbsMPolyMultSet{BaseRingType, BaseRingElemType, RingType, RingElemType}
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
    ) where {RingType<:MPolyRing, MultSetType<:AbsMPolyMultSet}
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

### additional getter functions
gens(W::MPolyLocalizedRing) = W.(gens(base_ring(W)))
ngens(W::MPolyLocalizedRing) = ngens(base_ring(W))

### required extension of the localization function
@Markdown.doc """

    localization(R::MPolyRing, U::AbsMPolyMultSet)   

Return the localization of `R` at `U`, together with the localization map.

# Examples
```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"]);

julia> P = ideal(R, [x])
ideal(x)

julia> U = complement_of_ideal(P)
complement of ideal(x)

julia> Rloc, iota = localization(R, U);

julia> Rloc
localization of Multivariate Polynomial Ring in x, y, z over Rational Field at the complement of ideal(x)

julia> iota
Map with following data
Domain:
=======
Multivariate Polynomial Ring in x, y, z over Rational Field
Codomain:
=========
localization of Multivariate Polynomial Ring in x, y, z over Rational Field at the complement of ideal(x)
```
""" localization(R::MPolyRing, U::AbsMPolyMultSet)

###localization is an Abstract Algebra alias for Localization

function Localization(S::AbsMPolyMultSet)
    R = ambient_ring(S)
    Rloc = MPolyLocalizedRing(R, S)
    #iota = MapFromFunc(x -> Rloc(x), R, Rloc)
    iota = hom(R, Rloc, Rloc.(gens(R)))
    return Rloc, iota
end

function Localization(R::MPolyRing, ord::MonomialOrdering)
  @assert R === ord.R
  return Localization(MPolyLeadingMonOne(ord))
end

### Successive localizations are handled by the dispatch for products
function Localization(
    W::MPolyLocalizedRing{BRT, BRET, RT, RET, MST}, 
    S::AbsMPolyMultSet{BRT, BRET, RT, RET}
  ) where {BRT, BRET, RT, RET, MST}
  issubset(S, inverted_set(W)) && return W, identity_map(W)
  U = S*inverted_set(W)
  L, _ = Localization(U)
  #return L, MapFromFunc((x->(L(numerator(x), denominator(x), check=false))), W, L)
  return L, MPolyLocalizedRingHom(W, L, hom(base_ring(W), L, L.(gens(base_ring(W)))), check=false)
end

### additional constructors
MPolyLocalizedRing(R::RingType, P::MPolyIdeal{RingElemType}) where {RingType, RingElemType} = MPolyLocalizedRing(R, MPolyComplementOfPrimeIdeal(P))

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
  L = MPolyLocalizedRing(R, MPolyPowersOfElement(R, vcat(g, divexact(f, h))))
  return L, MPolyLocalizedRingHom(W, L, hom(base_ring(W), L, L.(gens(base_ring(W)))), check=false)
  #return L, MapFromFunc((x->L(numerator(x), denominator(x), check=false)), W, L)
end

function Localization(
    W::MPolyLocalizedRing{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}}, 
    v::Vector{RET}
  ) where {BRT, BRET, RT, RET}
  V = W
  for f in v
    V = Localization(V, f)
  end
  return V, MPolyLocalizedRingHom(W, V, hom(base_ring(W), V, V.(gens(base_ring(W)))), check=false)
  #return V, MapFromFunc((x->V(numerator(x), denominator(x), check=false)), W, V)
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

  W::MPolyLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}
  frac::AbstractAlgebra.Generic.Frac{RingElemType}

  function MPolyLocalizedRingElem(
      W::MPolyLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType},
      f::AbstractAlgebra.Generic.Frac{RingElemType};
      check::Bool=true
    ) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}
    base_ring(parent(f)) == base_ring(W) || error(
	"the numerator and denominator of the given fraction do not belong to the original ring before localization"
      )
    if check && !iszero(f) && !is_unit(denominator(f))
      denominator(f) in inverted_set(W) || error("the given denominator is not admissible for this localization")
    end
    return new{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}(W, f)
  end
end

### required getter functions 
numerator(a::MPolyLocalizedRingElem) = numerator(a.frac)

denominator(a::MPolyLocalizedRingElem) = denominator(a.frac)

parent(a::MPolyLocalizedRingElem) = a.W

### additional getter functions
fraction(a::MPolyLocalizedRingElem) = a.frac
# to assure compatibility with generic code for MPolyQuoLocalizedRings:
lifted_numerator(a::MPolyLocalizedRingElem) = numerator(a)
lifted_denominator(a::MPolyLocalizedRingElem) = denominator(a)

### required conversions
function (W::MPolyLocalizedRing{
    BaseRingType, 
    BaseRingElemType, 
    RingType, 
    RingElemType, 
    MultSetType
  })(f::RingElemType; check::Bool=true) where {
    BaseRingType, 
    BaseRingElemType, 
    RingType, 
    RingElemType<:RingElem, 
    MultSetType
  } 
  return MPolyLocalizedRingElem(W, FractionField(base_ring(W))(f), check=false)
end

(W::MPolyLocalizedRing{
    BaseRingType, 
    BaseRingElemType, 
    RingType, 
    RingElemType, 
    MultSetType
  })(f::BaseRingElemType) where {
    BaseRingType<:Ring, 
    BaseRingElemType<:RingElem, 
    RingType<:Ring, 
    RingElemType<:RingElem, 
    MultSetType<:AbsMultSet
  } = MPolyLocalizedRingElem(W, FractionField(base_ring(W))(base_ring(W)(f)), check=false)

# Remove ambiguities
(W::MPolyLocalizedRing{
    BaseRingType, 
    BaseRingElemType, 
    RingType, 
    RingElemType, 
    MultSetType
  })(f::BaseRingElemType) where {
    BaseRingType<:Ring, 
    BaseRingElemType<:fmpz, 
    RingType<:Ring, 
    RingElemType<:RingElem, 
    MultSetType<:AbsMultSet
  } = MPolyLocalizedRingElem(W, FractionField(base_ring(W))(base_ring(W)(f)), check=false)

function (W::MPolyLocalizedRing{
    BaseRingType, 
    BaseRingElemType, 
    RingType, 
    RingElemType, 
    MultSetType
  })(a::RingElemType, b::RingElemType; check::Bool=true) where {
    BaseRingType, 
    BaseRingElemType, 
    RingType, 
    RingElemType, 
    MultSetType
  } 
  return MPolyLocalizedRingElem(W, a//b, check=check)
end

### additional conversions
(W::MPolyLocalizedRing{
    BaseRingType, 
    BaseRingElemType, 
    RingType, 
    RingElemType, 
    MultSetType
  })(f::AbstractAlgebra.Generic.Frac{RingElemType}; check::Bool=true) where {
    BaseRingType, 
    BaseRingElemType, 
    RingType, 
    RingElemType, 
    MultSetType
  } = MPolyLocalizedRingElem(W, f, check=check)

### additional promotions 
AbstractAlgebra.promote_rule(::Type{RET}, ::Type{MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT<:Ring, RET<:RingElement, MST} = MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}

AbstractAlgebra.promote_rule(::Type{BRET}, ::Type{MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}}) where {BRT<:Ring, BRET<:RingElement, RT<:Ring, RET<:RingElement, MST} = MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}

AbstractAlgebra.promote_rule(::Type{Integer}, ::Type{MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}}) where {BRT<:Ring, BRET<:RingElement, RT<:Ring, RET<:RingElement, MST} = MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}

AbstractAlgebra.promote_rule(::Type{fmpz}, ::Type{MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}}) where {BRT<:Ring, BRET<:RingElement, RT<:Ring, RET<:RingElement, MST} = MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}

### overwriting the arithmetic using the fractions from AbstractAlgebra
function +(a::T, b::T) where {T<:MPolyLocalizedRingElem}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  return (parent(a))(fraction(a) + fraction(b), check=false)
end

function -(a::T, b::T) where {T<:MPolyLocalizedRingElem}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  return (parent(a))(fraction(a) - fraction(b), check=false)
end

function -(a::T) where {T<:MPolyLocalizedRingElem}
  return (parent(a))((-1)*fraction(a), check=false)
end

function *(a::T, b::T) where {T<:MPolyLocalizedRingElem}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  return (parent(a))(fraction(a) * fraction(b), check=false)
end

function *(a::RET, b::MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET <: RingElem, MST}
  return (parent(b))(a*fraction(b), check=false)
end

function *(a::MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}, b::RET) where {BRT, BRET, RT, RET <: RingElem, MST}
  return b*a
end

function *(a::BRET, b::MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}) where {BRT, BRET <: RingElem, RT, RET, MST}
  return (parent(b))(a*numerator(b), denominator(b), check=false)
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
  return (parent(a))(fraction(a) // fraction(b), check=false)
end

function ==(a::T, b::T) where {T<:MPolyLocalizedRingElem}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  return fraction(a) == fraction(b)
end

# We need to manually split this into three methods, because 
# otherwise it seems that Julia can not dispatch this function.
function ^(a::MPolyLocalizedRingElem, i::Int64)
  return parent(a)(fraction(a)^i, check=false)
end
function ^(a::MPolyLocalizedRingElem, i::Integer)
  return parent(a)(fraction(a)^i, check=false)
end
function ^(a::MPolyLocalizedRingElem, i::fmpz)
  return parent(a)(fraction(a)^i, check=false)
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
  return W(m, n, check=false)
end

@Markdown.doc """
    is_unit(f::MPolyLocalizedRingElem)

Return `true`, if `f` is a unit of `parent(f)`, `false` otherwise.

# Examples

```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"]);

julia> P = ideal(R, [x])
ideal(x)

julia> U = complement_of_ideal(P)
complement of ideal(x)

julia> Rloc, iota = localization(U);

julia> is_unit(iota(x))
false

julia> is_unit(iota(y))
true
```
""" 
is_unit(f::MPolyLocalizedRingElem) = numerator(f) in inverted_set(parent(f))


########################################################################
# implementation of Oscar's general ring interface                     #
########################################################################

one(W::MPolyLocalizedRing) = W(one(base_ring(W)))
zero(W::MPolyLocalizedRing) = W(zero(base_ring(W)))

elem_type(W::MPolyLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyLocalizedRingElem{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}
elem_type(T::Type{MPolyLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyLocalizedRingElem{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}

parent_type(f::MPolyLocalizedRingElem{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}
parent_type(T::Type{MPolyLocalizedRingElem{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}

function (W::MPolyLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType})(f::MPolyLocalizedRingElem{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}; check::Bool=true) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} 
  parent(f) === W && return f
  return MPolyLocalizedRingElem(W, fraction(f), check=check)
end

(W::MPolyLocalizedRing)() = zero(W)
(W::MPolyLocalizedRing)(a::Integer) = W(base_ring(W)(a), check=false)
(W::MPolyLocalizedRing)(a::Int64) = W(base_ring(W)(a), check=false)
(W::MPolyLocalizedRing)(a::fmpz) = W(base_ring(W)(a), check=false)

is_domain_type(T::Type{MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = true 
is_exact_type(T::Type{MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = true

### promotion rules
AbstractAlgebra.promote_rule(::Type{MPolyLocalizedRingElem{RT, RET, MST}}, ::Type{MPolyLocalizedRingElem{RT, RET, MST}}) where {RT<:Ring, RET<:RingElement, MST} = MPolyLocalizedRingElem{RT, RET, MST}

function AbstractAlgebra.promote_rule(::Type{MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}}, ::Type{T}) where {BRT<:Ring, BRET<:RingElement, RT<:Ring, RET<:RingElement, MST, T<:RingElement} 
  AbstractAlgebra.promote_rule(RET, T) == RET && return MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}
  return AbstractAlgebra.promote_rule(BRET, T) == BRET ? MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST} : Union{}
end

@attr function base_ring_shifts(L::MPolyLocalizedRing{<:Any, <:Any, <:Any, <:Any, <:MPolyComplementOfKPointIdeal}) 
  a = point_coordinates(inverted_set(L))
  R = base_ring(L)
  shift = hom(R, R, gens(R)+R.(a))
  back_shift = hom(R, R, gens(R)-R.(a))
  return shift, back_shift
end

########################################################################
# Ideals in localizations of multivariate polynomial rings             #
########################################################################
# 
# If 𝕜 is a Noetherian ring, any localization W = R[U⁻¹] of a multi-
# variate polynomial ring R = 𝕜[x₁,…,xₙ] is again Noetherian and 
# any ideal I ⊂ W is of the form I = I'⋅ W for some ideal I' ⊂ R. 
# This correspondence is not 1:1 but for any ideal I ⊂ W we always 
# have that 
#
#   J = { x∈ R | ∃ u ∈ U : u⋅x ∈ I }
#
# is the unique element which is maximal among all ideals I' in R for 
# which I = I'⋅W. We call this the `saturated ideal` I': U of the localization.
#
# Since computation of the saturated ideal might be expensive, 
# we usually only cache a 'pre-saturated ideal' I' ⊂ J ⊂ I': U.

@Markdown.doc """
    MPolyLocalizedIdeal{
        LocRingType<:MPolyLocalizedRing, 
        LocRingElemType<:MPolyLocalizedRingElem
      } <: AbsLocalizedIdeal{LocRingElemType}

Ideals in localizations of polynomial rings.
"""
@attributes mutable struct MPolyLocalizedIdeal{
    LocRingType<:MPolyLocalizedRing, 
    LocRingElemType<:MPolyLocalizedRingElem
  } <: AbsLocalizedIdeal{LocRingElemType}
  # the initial set of generators, not to be changed ever!
  gens::Vector{LocRingElemType}
  # the ambient ring for this ideal
  W::LocRingType

  # fields for caching 
  map_from_base_ring::Hecke.Map

  # The pre_saturated_ideal can be any ideal I in the `base_ring` of W
  # such that W(I) is this ideal
  pre_saturated_ideal::MPolyIdeal
  # The following field contains a matrix A such that the *transpose* 
  # of A can be used to compute coordinates v of elements a w.r.t the 
  # generators of the pre_saturated_ideal to the coordinates w 
  # of W(a) w.r.t the generators of this ideal: 
  #   
  #     v*A^T = w, or equivalently (A*v^T)^T = w.
  #
  # In the code below we will often use the latter, because multiplication 
  # with sparse matrices is only implemented from the left. 
  pre_saturation_data::SMat

  is_saturated::Bool
  saturated_ideal::MPolyIdeal
 
  function MPolyLocalizedIdeal(
      W::MPolyLocalizedRing, 
      gens::Vector{LocRingElemType};
      map_from_base_ring::Hecke.Map = MapFromFunc(
          x->W(x),
          y->(isone(denominator(y)) ? numerator(y) : divexact(numerator(y), denominator(y))),
          base_ring(W), 
          W
        )
    ) where {LocRingElemType<:AbsLocalizedRingElem}
    for f in gens
      parent(f) == W || error("generator is not an element of the given ring")
    end

    I = new{typeof(W), LocRingElemType}()
    I.gens = gens
    I.W = W
    I.map_from_base_ring = map_from_base_ring
    I.is_saturated=false
    return I
  end
end
 
### required getter functions
gens(I::MPolyLocalizedIdeal) = copy(I.gens)
base_ring(I::MPolyLocalizedIdeal) = I.W

### type getters
ideal_type(::Type{MPolyLocalizedRingType}) where {MPolyLocalizedRingType<:MPolyLocalizedRing} = MPolyLocalizedIdeal{MPolyLocalizedRingType, elem_type(MPolyLocalizedRingType)}
ideal_type(L::MPolyLocalizedRing) = ideal_type(typeof(L))

### additional getter functions 
map_from_base_ring(I::MPolyLocalizedIdeal) = I.map_from_base_ring
is_saturated(I::MPolyLocalizedIdeal) = I.is_saturated
ngens(I::MPolyLocalizedIdeal) = length(I.gens)
getindex(I::MPolyLocalizedIdeal, k::Int) = copy(I.gens[k])

function ideal_membership(a::RingElem, I::MPolyLocalizedIdeal)
  L = base_ring(I)
  parent(a) == L || return L(a) in I
  b = numerator(a)
  b in pre_saturated_ideal(I) && return true
  is_saturated(I) && return false
  R = base_ring(L)
  J = pre_saturated_ideal(I)
  (success, x, u) = has_solution(generator_matrix(J), MatrixSpace(R, 1, 1)([b]), inverted_set(L))
  !success && return false
  # cache the intermediate result
  extend_pre_saturated_ideal!(I, b, x, u, check=false)
  return true
end

# TODO: Also add a special dispatch for localizations at 𝕜-points
@attr function is_prime(I::MPolyLocalizedIdeal)
  return is_prime(saturated_ideal(I))
end

### Additional constructors
function intersect(I::MPolyLocalizedIdeal, J::MPolyLocalizedIdeal)
  L = base_ring(I)
  L == base_ring(J) || error("ideals must be defined in the same ring")
  preI = Oscar.pre_saturated_ideal(I)
  preJ = Oscar.pre_saturated_ideal(J) 
  R = base_ring(L)
  K = intersect(I, J)
  return L(K)
end

### Further functionality
function coordinates(a::RingElem, I::MPolyLocalizedIdeal; check::Bool=true)
  L = base_ring(I)
  parent(a) === L || return coordinates(L(a), I, check=check)
  if check 
    a in I || error("the given element is not in the ideal")
  end
  J = pre_saturated_ideal(I)
  R = base_ring(J)
  p = numerator(a)
  if p in J
    q = denominator(a)
    # caching has been done during the call of `in`, so the following will work
    x = coordinates(p, pre_saturated_ideal(I))
    # multiplications sparse*dense have to be carried out this way round.
    return transpose(mul(pre_saturation_data(I), transpose(L(one(q), q, check=false)*change_base_ring(L, x))))
  else
    (success, x, u) = has_solution(generator_matrix(J), MatrixSpace(R, 1, 1)([p]), inverted_set(L), check=false)
    !success && error("check for membership was disabled, but element is not in the ideal")
    # cache the intermediate result
    #result = L(one(R), u*denominator(a), check=false)*change_base_ring(L, x)*pre_saturation_data(I)
    result = transpose(mul(pre_saturation_data(I), transpose(L(one(R), u*denominator(a), check=false)*change_base_ring(L, x))))
    extend_pre_saturated_ideal!(I, p, x, u, check=false)
    return result
  end
end

function coordinates(
    a::RingElem, I::MPolyLocalizedIdeal{LRT}; check::Bool=true
  ) where {LRT<:MPolyLocalizedRing{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement}}
  L = base_ring(I)
  parent(a) === L || return coordinates(L(a), I, check=check)
  if check 
    a in I || error("the given element is not in the ideal")
  end
  saturated_ideal(I, with_generator_transition=true) # Computing the saturation first is cheaper than the generic Posur method
  # Note that a call to saturated_ideal overwrites the cache for pre_saturated_ideal
  J = pre_saturated_ideal(I)
  R = base_ring(J)
  p = numerator(a)
  x = coordinates(p, J)
  q = denominator(a)
  return transpose(mul(pre_saturation_data(I), transpose(L(one(q), q, check=false)*change_base_ring(L, x))))
end

generator_matrix(J::MPolyIdeal) = MatrixSpace(base_ring(J), ngens(J), 1)(gens(J))

@Markdown.doc """
    saturated_ideal(I::MPolyLocalizedIdeal)

Given an ideal `I` of a localization, say, `Rloc` of a multivariate polynomial ring, say, `R`,
return the saturation of `I` over `R`. That is, return the largest ideal of `R` whose extension to 
`Rloc` is `I`. This is the preimage of `I` under the localization map.

    saturated_ideal(I::MPolyQuoLocalizedIdeal)

Given an ideal `I` of a localization, say, `RQL` of a quotient, say, `RQ` of a multivariate 
polynomial ring, say, `R`, return the preimage of the saturation of `I` over `RQ` under the 
projection map `R -> RQ`.

# Examples
```jldoctest
julia> R, (x,) = PolynomialRing(QQ, ["x"]);

julia> U = powers_of_element(x)
powers of fmpq_mpoly[x]

julia> Rloc, iota = localization(R, U);

julia> I = ideal(Rloc, [x+x^2])
ideal in localization of Multivariate Polynomial Ring in x over Rational Field at the powers of fmpq_mpoly[x] generated by the elements (x^2 + x)//1

julia> SI = saturated_ideal(I)
ideal(x + 1)

julia> base_ring(SI)
Multivariate Polynomial Ring in x over Rational Field

julia> U = complement_of_ideal(R, [0])
complement of maximal ideal corresponding to point with coordinates fmpq[0]

julia> Rloc, iota = localization(R, U);

julia> I = ideal(Rloc, [x+x^2])
ideal in localization of Multivariate Polynomial Ring in x over Rational Field at the complement of maximal ideal corresponding to point with coordinates fmpq[0] generated by the elements (x^2 + x)//1

julia> saturated_ideal(I)
ideal(x)
```
"""
function saturated_ideal(I::MPolyLocalizedIdeal)
  if !isdefined(I, :saturated_ideal)
    error("method `saturated_ideal` is not implemented for ideals of type $(typeof(I))")
  end
  return I.saturated_ideal
end

function pre_saturated_ideal(I::MPolyLocalizedIdeal)
  if !isdefined(I, :pre_saturated_ideal)
    W = base_ring(I)
    I.pre_saturated_ideal = ideal(base_ring(W), numerator.(gens(I)))
    r = length(gens(I))
    A = zero_matrix(SMat, W, 0, r)
    #A = zero(MatrixSpace(W, r, r))
    for i in 1:r
      push!(A, sparse_row(W, [(i, W(denominator(gens(I)[i])))]))
      #A[i, i] = denominator(gens(I)[i])
    end
    I.pre_saturation_data = A
  end
  return I.pre_saturated_ideal
end

function pre_saturation_data(I::MPolyLocalizedIdeal)
  if !isdefined(I, :pre_saturation_data)
    pre_saturated_ideal(I)
  end
  return I.pre_saturation_data
end

function extend_pre_saturated_ideal!(
    I::MPolyLocalizedIdeal, f::PT, x::MatrixElem{PT}, u::PT;
    check::Bool=true
  ) where {PT <: MPolyElem}
  nrows(x) == 1 || error("matrix must be a row vector")
  L = base_ring(I)
  R = base_ring(L)
  J = pre_saturated_ideal(I)
  if check
    u*f == dot(x, gens(J)) || error("input is not coherent")
  end
  J_ext = ideal(R, vcat(gens(J), [f]))
  T = pre_saturation_data(I)
  y = mul(T, transpose(L(one(u), u, check=false)*change_base_ring(L, x)))
  T_ext = zero_matrix(SMat, L, 0, ncols(T)+1)
  for i in 1:length(y)
    push!(T_ext, T[i] + sparse_row(L, [(ncols(T) + 1, y[i, 1])]))
  end
  #T_ext = vcat(T, L(one(u), u, check=false)*change_base_ring(L, x)*T)
  I.pre_saturated_ideal = J_ext
  I.pre_saturation_data = T_ext
  return I
end

function extend_pre_saturated_ideal!(
    I::MPolyLocalizedIdeal, f::Vector{PT}, x::MatrixElem{PT}, u::Vector{PT};
    check::Bool=true
  ) where {PT <: MPolyElem}
  L = base_ring(I)
  R = base_ring(L)
  J = pre_saturated_ideal(I)
  n = length(f)
  n == length(u) == nrows(x) || error("input dimensions are not compatible")
  if check
    for i in 1:n
      u[i]*f[i] == dot(x[i, :], gens(J)) || error("input is not coherent")
    end
  end
  J_ext = ideal(R, vcat(gens(J), f))
  T = pre_saturation_data(I)
  #T_ext = vcat(T, 
  #             diagonal_matrix([L(one(v), v, check=false) for v in u])*
  #             change_base_ring(L, x)*T
  #            )
  #y = T * transpose(L(one(u), u, check=false)*change_base_ring(L, x))
  y = mul(T, transpose(change_base_ring(L, x)))
  for i in 1:ncols(y)
    for j in 1:n
      y[i, j] = y[i, j]*L(one(u[i]), u[i], check=false) 
    end
  end
  T_ext = zero_matrix(SMat, L, 0, ncols(T)+n)
  for i in 1:length(y)
    push!(T_ext, T[i] + sparse_row(L, [(ncols(T) + j, y[i, j]) for j in 1:n]))
  end
  I.pre_saturated_ideal = J_ext
  I.pre_saturation_data = T_ext
  return I
end

function _diagonal_sparse_matrix(L::Ring, v::Vector{T}) where {T<:RingElem}
  A = sparse_matrix(L)
  for i in 1:length(v)
    push!(A, sparse_row(L, [(i, v[i])]))
  end
  return A
end

function coordinate_shift(
    L::MPolyLocalizedRing{<:Any, <:Any, <:Any, <:Any, <:MPolyComplementOfKPointIdeal}
  )
  if !has_attribute(L, :coordinate_shift)
    U = inverted_set(L)
    a = point_coordinates(U)
    R = base_ring(L)
    Ls = MPolyLocalizedRing(R, MPolyComplementOfKPointIdeal(R, [0 for i in 1:ngens(R)]))
    xs = [ x + a for (x, a) in zip(gens(base_ring(L)), a) ]
    xs_inv = [ x - a for (x, a) in zip(gens(base_ring(L)), a) ]
    shift = MapFromFunc(
                f -> Ls(evaluate(numerator(f), xs), evaluate(denominator(f), xs), check=false),
                g -> L(evaluate(numerator(g), xs_inv), evaluate(denominator(g), xs_inv), check=false),
                L, Ls
              )
    set_attribute!(L, :coordinate_shift, shift)
  end
  return get_attribute(L, :coordinate_shift)::Hecke.Map
end


function saturated_ideal(
    I::MPolyLocalizedIdeal{LRT};
    with_generator_transition::Bool=false
  ) where {LRT<:MPolyLocalizedRing{<:Any, <:Any, <:Any, <:Any, <:MPolyComplementOfKPointIdeal}}
  if !isdefined(I, :saturated_ideal)
    is_saturated(I) && return pre_saturated_ideal(I)
    L = base_ring(I)
    R = base_ring(base_ring(I))
    result = ideal(R, [one(R)])
    J = pre_saturated_ideal(I)
    if !all(x->iszero(evaluate(x, point_coordinates(inverted_set(L)))), gens(J)) 
      I.saturated_ideal = result
      return result
    end
    pdec = primary_decomposition(J)
    for (Q, P) in pdec
      if all(x->iszero(evaluate(x, point_coordinates(inverted_set(L)))), gens(P))
        result = intersect(result, Q)
      end
    end
    I.saturated_ideal = result
    if with_generator_transition
      error("computation of the transition matrix for the generators is not supposed to happen because of using local orderings")
    end
  end
  return I.saturated_ideal
end

function saturated_ideal(
    I::MPolyLocalizedIdeal{LRT};
    with_generator_transition::Bool=false
  ) where {LRT<:MPolyLocalizedRing{<:Any, <:Any, <:Any, <:Any, <:MPolyComplementOfPrimeIdeal}}
  if !isdefined(I, :saturated_ideal)
    is_saturated(I) && return pre_saturated_ideal(I)
    L = base_ring(I)
    R = base_ring(L)
    result = ideal(R, [one(R)])
    U = inverted_set(base_ring(I))
    if !issubset(I, L(prime_ideal(U)))
      I.saturated_ideal = result
      return result
    end
    J = pre_saturated_ideal(I)
    pdec = primary_decomposition(J)
    for (Q, P) in pdec
      if issubset(P,prime_ideal(U))
        result = intersect(result, Q)
      end
    end
    I.saturated_ideal = result
    if with_generator_transition
      error("computation of the transition matrix for the generators is not supposed to happen for localizations at complements of prime ideals")
    end
  end
  return I.saturated_ideal
end

function saturated_ideal(
    I::MPolyLocalizedIdeal{LRT};
    strategy::Symbol=:iterative_saturation,
    with_generator_transition::Bool=true
  ) where {LRT<:MPolyLocalizedRing{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement}}
  if !isdefined(I, :saturated_ideal)
    is_saturated(I) && return pre_saturated_ideal(I)
    L = base_ring(I)
    R = base_ring(L)
    if strategy==:iterative_saturation
      Jsat = pre_saturated_ideal(I)
      for d in denominators(inverted_set(L))
        if !is_unit(d) && !iszero(Jsat)
          Jsat = saturation(Jsat, ideal(R, d))
        end
        if with_generator_transition
          # We completely overwrite the pre_saturated_ideal with the generators 
          # of the saturated ideal
          cache = Vector()
          for g in gens(Jsat)
            (k, dttk) = Oscar._minimal_power_such_that(d, p->(p*g in pre_saturated_ideal(I)))
            if k > 0
              push!(cache, (g, coordinates(dttk*g, pre_saturated_ideal(I)), dttk))
            else
              # We hit the unit ideal. Now we could simply return that; but we don't for 
              # the moment.
              push!(cache, (g, coordinates(dttk*g, pre_saturated_ideal(I)), dttk))
            end
          end
          if length(cache) > 0
            #A = zero(MatrixSpace(L, ngens(Jsat), ngens(I)))

            #for i in 1:length(cache)
            #  (g, a, dttk) = cache[i]
            #  A[i, :] = L(one(dttk), dttk, check=false)*change_base_ring(L, a)*pre_saturation_data(I)
            #end
            cols = Vector()
            for i in 1:length(cache)
              (g, a, dttk) = cache[i]
              push!(cols, mul(pre_saturation_data(I), transpose(L(one(dttk), dttk, check=false)*change_base_ring(L, a))))
            end
            A = zero_matrix(SMat, L, 0, length(cache))
            for i in 1:ngens(I)
              v = sparse_row(L, [(j, cols[j][i, 1]) for j in 1:length(cols)])
              push!(A, v)
            end
            I.pre_saturated_ideal = ideal(R, gens(Jsat))
            I.pre_saturation_data = A
#            extend_pre_saturated_ideal!(I, 
#                                        elem_type(R)[g for (g, x, u) in cache],
#                                        vcat(dense_matrix_type(R)[x for (g, x, u) in cache]),
#                                        [u for (g, x, u) in cache], 
#                                        check=false
#                                       )
          end
          I.is_saturated=true
        end
      end
      I.saturated_ideal = Jsat
    elseif strategy==:single_saturation
      d = prod(denominators(inverted_set(L)))
      Jsat = pre_saturated_ideal(I)
      if !is_unit(d) && !iszero(pre_saturated_ideal(I))
        Jsat = saturation(Jsat, ideal(R, d))
      end
      if with_generator_transition
        cache = Vector()
        for g in gens(Jsat)
          (k, dttk) = Oscar._minimal_power_such_that(d, p->(p*g in pre_saturated_ideal(I)))
          #if k > 0
            push!(cache, (g, coordinates(dttk*g, pre_saturated_ideal(I)), dttk))
          #end
        end
        if length(cache) > 0
          # We completely overwrite the pre_saturated_ideal with the generators 
          # of the saturated ideal
          #A = zero(MatrixSpace(L, ngens(Jsat), ngens(I)))
          #for i in 1:length(cache)
          #  (g, a, dttk) = cache[i]
          #  A[i, :] = L(one(dttk), dttk, check=false)*change_base_ring(L, a)*pre_saturation_data(I)
          #end
          cols = Vector()
          for i in 1:length(cache)
            (g, a, dttk) = cache[i]
            push!(cols, mul(pre_saturation_data(I), transpose(L(one(dttk), dttk, check=false)*change_base_ring(L, a))))
          end
          A = zero_matrix(SMat, L, 0, length(cache))
          for i in 1:ngens(I)
            push!(A, sparse_row(L, [(j, cols[j][i, 1]) for j in 1:length(cols)]))
          end
          I.pre_saturated_ideal = Jsat
          I.pre_saturation_data = A
         #extend_pre_saturated_ideal!(I, 
         #                            elem_type(R)[g for (g, x, u) in cache],
         #                            vcat(dense_matrix_type(R)[x for (g, x, u) in cache]),
         #                            [u for (g, x, u) in cache], 
         #                            check=false
         #                           )
        end
        I.is_saturated=true
      end
      I.saturated_ideal = Jsat
    else
      error("strategy $strategy not implemented")
    end
  end
  return I.saturated_ideal
end

@doc Markdown.doc"""
    saturated_ideal(I::MPolyLocalizedIdeal)
    saturated_ideal(I::MPolyQuoLocalizedIdeal)
    saturated_ideal(I::MPolyIdeal)
    saturated_ideal(I::MPolyQuoIdeal)

Returns the largest ideal in $R$ mapping to $I$ under the canonical map $R \longrightarrow S$ for an ideal $I \in S$, where $S$ may be any of the
types 'MPolyLocalizedRing', 'MPolyQuoLocalizedRing', 'MPolyRing' and 'MPolyQuoRing' and $R$ is the underlying ring of type 'MPolyRing'.
Note that the last two variants are only provided to allow a coherent usage.  

# Examples
``` jldoctest
julia> R, (x, y) = QQ["x", "y"]
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> I = ideal(R, [x, y^2+1])
ideal(x, y^2 + 1)

julia> U = MPolyComplementOfPrimeIdeal(I)
complement of ideal(x, y^2 + 1)

julia> L = MPolyLocalizedRing(R, U)
localization of Multivariate Polynomial Ring in x, y over Rational Field at the complement of ideal(x, y^2 + 1)

julia> J = ideal(L,[y*(x^2+(y^2+1)^2)])
ideal in localization of Multivariate Polynomial Ring in x, y over Rational Field at the complement of ideal(x, y^2 + 1) generated by the elements (x^2*y + y^5 + 2*y^3 + y)//1

julia> saturated_ideal(J)
ideal(x^2 + y^4 + 2*y^2 + 1)

julia> JJ = ideal(R,[y*(x^2+(y^2+1)^2)])
ideal(x^2*y + y^5 + 2*y^3 + y)

julia> saturated_ideal(JJ)
ideal(x^2*y + y^5 + 2*y^3 + y)

```
"""
function saturated_ideal(
    I::MPolyLocalizedIdeal{LRT};
    with_generator_transition::Bool=false
  ) where {LRT<:MPolyLocalizedRing{<:Any, <:Any, <:Any, <:Any, <:MPolyProductOfMultSets}}
  if !is_saturated(I)
    W = base_ring(I)
    R = base_ring(W)
    J = ideal(R, numerator.(gens(I)))
    for U in sets(inverted_set(W))
      L, _ = Localization(U)
      J = saturated_ideal(L(J))
    end
    if with_generator_transition
      for g in gens(J)
        g in I
      end
    end
    I.is_saturated = true
  end
  return pre_saturated_ideal(I)
end

function cache_transitions_for_saturation(I::MPolyLocalizedIdeal) 
  for g in saturated_ideal(I)
    g in I
  end
  return I
end

# for convenience of scripting user allow I::MPolyIdeal as input to
# saturated_ideal and return the ideal itself 
saturated_ideal(I::MPolyIdeal) = I

# the following overwrites the membership test 
# assuming that direct computation of the saturation 
# is cheaper when localizing at powers of elements.
function ideal_membership(
    a::RingElem, 
    I::MPolyLocalizedIdeal{LocRingType}
  ) where {
    LocRingType<:MPolyLocalizedRing{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement}
  }
  L = base_ring(I)
  parent(a) == L || return L(a) in I
  b = numerator(a)
  return b in saturated_ideal(I)
end


### Conversion of ideals in the original ring to localized ideals
function (W::MPolyLocalizedRing{BRT, BRET, RT, RET, MST})(I::MPolyIdeal{RET}) where {BRT, BRET, RT, RET, MST}
  return MPolyLocalizedIdeal(W, W.(gens(I)))
end

### required constructors 
function ideal(
    W::MPolyLocalizedRing, f
  )
  return MPolyLocalizedIdeal(W, [W(f)])
end

@doc Markdown.doc"""    
    ideal(Rloc::MPolyLocalizedRing, V::Vector)

Given a localization `Rloc` of a multivariate polynomial ring, and given a vector `V` of 
elements of `Rloc`, create the ideal of `Rloc` which is generated by the entries of `V`.
"""
function ideal(
    W::MPolyLocalizedRing,
    gens::Vector
  )
  length(gens) == 0 && return MPolyLocalizedIdeal(W, elem_type(W)[])
  return MPolyLocalizedIdeal(W, W.(gens))
end

function ideal(
    W::MPolyLocalizedRing,
    I::Ideal
  )
  return MPolyLocalizedIdeal(W, W.(gens(I)))
end

### additional functionality
function issubset(I::IdealType, J::IdealType) where {IdealType<:MPolyLocalizedIdeal}
  base_ring(I) == base_ring(J) || error("ideals do not belong to the same ring")
  for g in gens(I)
    g in J || return false
  end
  return true
end

==(I::IdealType, J::IdealType) where {IdealType<:MPolyLocalizedIdeal} = (issubset(I, J) && issubset(J, I))

function +(I::IdealType, J::IdealType) where {IdealType<:MPolyLocalizedIdeal}
  return ideal(base_ring(I), vcat(gens(I), gens(J)))
end

# TODO: The following method can probably be fine tuned for specific localizations.
function intersect(I::IdealType, J::IdealType) where {IdealType<:MPolyLocalizedIdeal}
  return base_ring(I)(intersect(pre_saturated_ideal(I), pre_saturated_ideal(J)))
end

# TODO: The following method can probably be fine tuned for specific localizations.
function quotient(I::IdealType, J::IdealType) where {IdealType<:MPolyLocalizedIdeal}
  #return base_ring(I)(quotient(saturated_ideal(I), saturated_ideal(J)))
  return base_ring(I)(quotient(pre_saturated_ideal(I), pre_saturated_ideal(J)))
end

Base.:(:)(I::IdealType, J::IdealType) where {IdealType<:MPolyLocalizedIdeal} = quotient(I, J)


function Base.show(io::IO, I::MPolyLocalizedIdeal) 
  print(io, "ideal in $(base_ring(I)) generated by the elements ") 
  n = length(gens(I))
  for i in (1:n-1)
    print(io, "$(gens(I)[i]), ")
  end
  print(io, last(gens(I)))
end

@attr function shifted_ideal(
    I::MPolyLocalizedIdeal{LRT, LRET}
  ) where {LRT<:MPolyLocalizedRing{<:Any, <:Any, <:Any, <:Any, <:MPolyComplementOfKPointIdeal}, LRET}
  L = base_ring(I)
  R = base_ring(L)
  shift, tfihs = base_ring_shifts(L)
  return ideal(R, shift.(gens(pre_saturated_ideal(I))))
end

function ideal_membership(
    a::RingElem, 
    I::MPolyLocalizedIdeal{LocRingType}
  ) where {
    LocRingType<:MPolyLocalizedRing{<:Any, <:Any, <:Any, <:Any, <:MPolyComplementOfKPointIdeal}
  }
  L = base_ring(I)
  parent(a) == L || return L(a) in I
  L = base_ring(I)
  R = base_ring(L)
  shift, tfihs = base_ring_shifts(L)
  Is = shifted_ideal(I)
  # We have to call for that groebner basis once manually. 
  # Otherwise the ideal membership will complain about the ordering not being global.
  o = negdegrevlex(gens(R))
  standard_basis(Is, ordering=o)
  return ideal_membership(shift(numerator(a)), Is, ordering=o)
end

function coordinates(
    a::RingElem,
    I::MPolyLocalizedIdeal{LRT} 
  ) where {LRT<:MPolyLocalizedRing{<:Any, <:Any, <:Any, <:Any, <:MPolyComplementOfKPointIdeal}}
  L = base_ring(I)
  L == parent(a) || return coordinates(L(a), I)
  R = base_ring(L)
  J = shifted_ideal(I)
  shift, tfihs = base_ring_shifts(L)
  p = shift(numerator(a))
  o = negdegrevlex(gens(R))
  x, u = Oscar.lift(p, J, o)
  T = pre_saturation_data(I)
  return transpose(mul(T, transpose(L(one(base_ring(L)), tfihs(u)*denominator(a), check=false)*change_base_ring(L, map_entries(tfihs,x)))))
  #return L(one(base_ring(L)), tfihs(u)*denominator(a), check=false)*change_base_ring(L, map_entries(tfihs,x))*T
end

########################################################################
# The next method is based on the following observation: 
#
#     p//q ∈ I ⋅ U⁻¹ ⋅ V⁻¹ ⊂ R[U⁻¹⋅V⁻¹] 
#   ⇔ ∃ u ∈ U, v∈ V : u ⋅ v ⋅ p ∈ I
#   ⇔ ∃ v ∈ V : v ⋅ p ∈ I : U 
#
# So, to compute the coordinates of p//q in I, we can proceed 
# inductively, by caching J = I : U and computing the coordinates 
# of p in J over R[V⁻¹]. 
#
# In particular for the case where U is the localization at a 
# hypersurface and V the localization at the complement of a prime 
# ideal, computing and working with the saturation in the first case 
# is quicker, while using local orderings/ the Posur method is favourable 
# in the second case. Since this situation will be of particular 
# interest, we provide a special path for that.
function coordinates(
    a::RingElem,
    I::MPolyLocalizedIdeal{LRT} 
  ) where {LRT<:MPolyLocalizedRing{<:Any, <:Any, <:Any, <:Any, <:MPolyProductOfMultSets}}
  L = base_ring(I)
  parent(a) == L || return coordinates(L(a), I)

  R = base_ring(L)
  U = sets(inverted_set(base_ring(I)))

  if length(U) == 1 
    if !has_attribute(I, :popped_ideal)
      W, _ = Localization(R, U[1])
      popped_ideal = W(pre_saturated_ideal(I))
      set_attribute!(I, :popped_ideal, popped_ideal)
    end
    popped_ideal = get_attribute(I, :popped_ideal)
    return transpose(mul(pre_saturation_data(I), transpose(L(one(R), denominator(a), check=false)*map_entries(x->L(x, check=false), coordinates(numerator(a), popped_ideal)))))
    #return L(one(R), denominator(a), check=false)*map_entries(x->L(x, check=false), coordinates(numerator(a), popped_ideal))*pre_saturation_data(I)
  end

  if numerator(a) in pre_saturated_ideal(I) 
    return transpose(mul(T, transpose(L(one(R), denominator(a), check=false)*map_entries(L, coordinates(numerator(a), pre_saturated_ideal(I))))))
  end

  i = findfirst(x->(typeof(x)<:MPolyPowersOfElement), U)
  if !isnothing(i)
    if !has_attribute(I, :popped_ideal)
      S = popat!(U, i)
      W, _ = Localization(base_ring(L), S)
      popped_ideal = ideal(W, pre_saturated_ideal(I))
      saturated_ideal(popped_ideal, with_generator_transition=true)
      set_attribute!(I, :popped_ideal, popped_ideal)
      Wnext, _ = Localization(R, MPolyProductOfMultSets(R, U))
      next_ideal = Wnext(pre_saturated_ideal(popped_ideal))
      set_attribute!(I, :next_ideal, next_ideal)
    end
    popped_ideal = get_attribute(I, :popped_ideal)
    next_ideal = get_attribute(I, :next_ideal)
    y = coordinates(numerator(a), next_ideal)
    x = transpose(mul(map_entries(x->L(x, check=false), pre_saturation_data(popped_ideal)), 
                      transpose(map_entries(x->L(x, check=false), y))))
    #x = map_entries(x->L(x, check=false), y)*map_entries(x->L(x, check=false), pre_saturation_data(popped_ideal))
    return L(one(R), denominator(a), check=false)*x
  else
    if !has_attribute(I, :popped_ideal)
      S = pop!(U)
      W, _ = Localization(base_ring(L), S)
      popped_ideal = ideal(W, pre_saturated_ideal(I))
      saturated_ideal(popped_ideal, with_generator_transition=true)
      set_attribute!(I, :popped_ideal, popped_ideal)
      Wnext, _ = Localization(R, MPolyProductOfMultSets(R, U))
      next_ideal = Wnext(pre_saturated_ideal(popped_ideal))
      set_attribute!(I, :next_ideal, next_ideal)
    end
    popped_ideal = get_attribute(I, :popped_ideal)
    next_ideal = get_attribute(I, :next_ideal)
    y = coordinates(numerator(a), next_ideal)
    x = transpose(mul(map_entries(x->L(x, check=false), pre_saturation_data(popped_ideal)), 
                      transpose(map_entries(x->L(x, check=false), y))))
    #x = map_entries(x->L(x, check=false), y)*map_entries(x->L(x, check=false), pre_saturation_data(popped_ideal))
    return L(one(R), denominator(a), check=false)*x
  end
end

function ideal_membership(
    a::RingElem,
    I::MPolyLocalizedIdeal{LRT} 
  ) where {LRT<:MPolyLocalizedRing{<:Any, <:Any, <:Any, <:Any, <:MPolyProductOfMultSets}}
  L = base_ring(I)
  parent(a) == L || return L(a) in I

  R = base_ring(L)
  U = sets(inverted_set(base_ring(I)))

  if length(U) == 1 
    if !has_attribute(I, :popped_ideal)
      W, _ = Localization(R, U[1])
      popped_ideal = W(pre_saturated_ideal(I))
      set_attribute!(I, :popped_ideal, popped_ideal)
    end
    popped_ideal = get_attribute(I, :popped_ideal)
    return numerator(a) in popped_ideal
  end

  numerator(a) in pre_saturated_ideal(I) && return true

  i = findfirst(x->(typeof(x)<:MPolyPowersOfElement), U)
  if !isnothing(i)
    if !has_attribute(I, :popped_ideal)
      S = popat!(U, i)
      W, _ = Localization(base_ring(L), S)
      popped_ideal = ideal(W, pre_saturated_ideal(I))
      saturated_ideal(popped_ideal, with_generator_transition=true)
      set_attribute!(I, :popped_ideal, popped_ideal)
      Wnext, _ = Localization(R, MPolyProductOfMultSets(R, U))
      next_ideal = Wnext(pre_saturated_ideal(popped_ideal))
      set_attribute!(I, :next_ideal, next_ideal)
    end
    next_ideal = get_attribute(I, :next_ideal)
    return numerator(a) in next_ideal
  else
    if !has_attribute(I, :popped_ideal)
      S = pop!(U)
      W, _ = Localization(base_ring(L), S)
      popped_ideal = ideal(W, pre_saturated_ideal(I))
      saturated_ideal(popped_ideal, with_generator_transition=true)
      set_attribute!(I, :popped_ideal, popped_ideal)
      Wnext, _ = Localization(R, MPolyProductOfMultSets(R, U))
      next_ideal = Wnext(pre_saturated_ideal(popped_ideal))
      set_attribute!(I, :next_ideal, next_ideal)
    end
    next_ideal = get_attribute(I, :next_ideal)
    return numerator(a) in next_ideal
  end
end


########################################################################
# special treatment of localization at orderings                       #
########################################################################

function ideal_membership(
    a::RingElem,
    I::MPolyLocalizedIdeal{LRT} 
  ) where {LRT<:MPolyLocalizedRing{<:Any, <:Any, <:Any, <:Any, <:MPolyLeadingMonOne}}
  parent(a) == base_ring(I) || return base_ring(I)(a) in I
  J = pre_saturated_ideal(I)
  p = numerator(a)
  o = ordering(inverted_set(parent(a)))
  # We have to call for that groebner basis once manually. 
  # Otherwise the ideal membership will complain about the ordering not being global.
  standard_basis(J, ordering=o)
  return ideal_membership(p, J, ordering=o)
end

function coordinates(
    a::RingElem,
    I::MPolyLocalizedIdeal{LRT} 
  ) where {LRT<:MPolyLocalizedRing{<:Any, <:Any, <:Any, <:Any, <:MPolyLeadingMonOne}}
  L = base_ring(I)
  L == parent(a) || return coordinates(L(a), I)
  J = pre_saturated_ideal(I)
  p = numerator(a)
  o = ordering(inverted_set(parent(a)))
  x, u = Oscar.lift(p, J, o)
  T = pre_saturation_data(I)
  return transpose(mul(T, transpose(L(one(base_ring(L)), u*denominator(a), check=false)*
                                    change_base_ring(L, x))))
end

@Markdown.doc """
    bring_to_common_denominator(f::Vector{T}) where {T<:MPolyLocalizedRingElem}

Given a vector of fractions [a₁//b₁,…,aₙ//bₙ] return a pair 
(d, λ) consisting of a common denominator d and a vector 
λ = [λ₁,…,λₙ] such that aᵢ//bᵢ = λᵢ⋅aᵢ//d
"""
function bring_to_common_denominator(f::Vector{T}) where {T<:MPolyLocalizedRingElem}
  length(f) == 0 && error("need at least one argument to determine the return type")
  R = base_ring(parent(f[1]))
  for a in f
    R == base_ring(parent(a)) || error("elements do not belong to the same ring")
  end
  d = one(R)
  a = Vector{elem_type(R)}()
  for den in denominator.(f)
    b = gcd(d, den)
    c = divexact(den, b)
    e = divexact(d, b)
    d = d*c
    a = [c*k for k in a]
    push!(a, e)
  end
  return d, a
end

write_as_linear_combination(f::MPolyLocalizedRingElem, g::Vector) = write_as_linear_combination(f, parent(f).(g))

# return the localized ring as a quotient of a polynomial ring using Rabinowitsch's trick.
@Markdown.doc """
    as_affine_algebra(
      L::MPolyLocalizedRing{BRT, BRET, RT, RET, 
      MPolyPowersOfElement{BRT, BRET, RT, RET}}; 
      inverse_name::String="θ"
    ) where {BRT, BRET, RT, RET}

For a localized polynomial ring ``L = 𝕜[x₁,…,xₙ][f⁻¹]`` this returns a 
quintuple ``(A, I, d, ϕ, θ)`` consisting of 

  * an `AffineAlgebra` ``A = 𝕜[x₁,…,xₙ,θ]/⟨1 - θ⋅d⟩`` isomorphic to ``L``
  * the ideal ``⟨1 - θ⋅d⟩``
  * an element ``d ∈ 𝕜[x₁,…,xₙ]`` at which has been localized
  * the natural inclusion ``ϕ : 𝕜[x₁,…,xₙ] ↪ A``
  * the localization variable ``θ`` corresponding to ``d⁻¹``.
"""
function as_affine_algebra(
    L::MPolyLocalizedRing{BRT, BRET, RT, RET, 
			     MPolyPowersOfElement{BRT, BRET, RT, RET}}; 
    inverse_name::String="θ"
  ) where {BRT, BRET, RT, RET}
  R = base_ring(L)
  A, phi, t = _add_variables_first(R, [inverse_name])
  theta = t[1]
  f = prod(denominators(inverted_set(L)))
  I = ideal(A, [one(A)-theta*phi(f)])
  return A, I, f, phi, theta
end


########################################################################
# Homomorphisms of localized polynomial rings                          #
########################################################################

mutable struct MPolyLocalizedRingHom{
                                     DomainType<:MPolyLocalizedRing, 
                                     CodomainType<:Ring, 
                                     RestrictedMapType<:Map
                                    } <: AbsLocalizedRingHom{
                                                             DomainType, 
                                                             CodomainType, 
                                                             RestrictedMapType
                                                            }
  # the domain of definition
  W::DomainType

  ### the codomain
  # Why do we need to store the codomain explicitly and not extract it from 
  # the restricted map? 
  # Because the restricted map res probably takes values in a strictly 
  # smaller or even completely different ring than S. If C is the codomain 
  # of res, then we impliticly assume (and check) that coercion from elements 
  # of C to S is possible. This is strictly weaker than requiring res to implement 
  # the full map interface. In that case, one would -- for example -- need to 
  # implement a new type just for the inclusion map R → R[U⁻¹] of a ring into 
  # its localization.
  S::CodomainType

  # the restriction of the map to the base_ring of the domain
  res::RestrictedMapType

  function MPolyLocalizedRingHom(
      W::DomainType,
      S::CodomainType,
      res::RestrictedMapType;
      check::Bool=true
    ) where {DomainType<:MPolyLocalizedRing, CodomainType<:Ring, RestrictedMapType<:Map}
    R = base_ring(W)
    U = inverted_set(W)
    domain(res) === R || error("the domain of the restricted map does not coincide with the base ring of the localization")
    for f in U
      is_unit(S(res(f))) || error("image of $f is not a unit in the codomain")
    end
    return new{DomainType, CodomainType, RestrictedMapType}(W, S, res) 
  end
end
  
### additional constructors
function MPolyLocalizedRingHom(
      R::MPolyRing,
      S::Ring,
      a::Vector{RingElemType}
  ) where {RingElemType<:RingElem}
  res = hom(R, S, a)
  W, _ = Localization(units_of(R))
  return MPolyLocalizedRingHom(W, S, res)
end

function MPolyLocalizedRingHom(
      W::MPolyLocalizedRing,
      S::Ring,
      a::Vector{RingElemType};
      check::Bool=true
  ) where {RingElemType<:RingElem}
  R = base_ring(W)
  res = hom(R, S, a)
  return MPolyLocalizedRingHom(W, S, res, check=check)
end

@doc Markdown.doc"""
    hom(Rloc::MPolyLocalizedRing, S::Ring, phi::Map)

Given a localized ring `Rloc`of type `MPolyLocalizedRing`, say `Rloc` is the localization 
of a multivariate polynomial ring `R` at the multiplicatively closed subset `U` of `R`, and 
given a homomorphism `phi` from `R` to `S` sending elements of `U` to units in `S`, return 
the homomorphism from `Rloc` to `S` whose composition with the localization map is `phi`.

    hom(Rloc::MPolyLocalizedRing, S::Ring, V::Vector)

Given a localized ring `Rloc` as above, and given a vector `V` of `nvars` elements of `S`, let `phi` 
be the homomorphism from `R` to `S` which is determined by the entries of `V` as the images of
the generators of `R`, and proceed as above.

    hom(RQL::MPolyQuoLocalizedRing, S::Ring, phi::Map)

Given a localized ring `RQL`of type `MPolyQuoLocalizedRing`, say `RQL` is the localization 
of a quotient ring `RQ` of a multivariate polynomial ring `R` at the multiplicatively closed subset `U` of `R`, and 
given a homomorphism `phi` from `R` to `S` sending elements of `U` to units in `S` and elements of the modulus
of `RQ` to zero, return the homomorphism from `Rloc` to `S` whose composition with the localization map `RQ -> RQL`
and the projection map `R -> RQ` is `phi`.

    hom(RQL::MPolyQuoLocalizedRing, S::Ring, V::Vector)

Given a localized ring `RQL`as above, and given a vector `V` of `nvars` elements of `S`, let `phi` 
be the homomorphism from `R` to `S` which is determined by the entries of `V` as the images of
the generators of `R`, and proceed as above.

!!! warning
    Except from the case where the type of `U` is `<: MPolyPowersOfElement`, the condition on `phi` requiring that elements of `U` are send to units in `S` is not checked by the `hom` constructor.
 
# Examples
```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"]);

julia> I = ideal(R, [y-x^2, z-x^3]);

julia> RQ, p = quo(R, I);

julia> UR = complement_of_ideal(R, [0, 0, 0]);

julia> RQL, _ = localization(RQ, UR);

julia> T, (t,) =  PolynomialRing(QQ, ["t"]);

julia> UT = complement_of_ideal(T, [0]);

julia> TL, _ =  localization(T, UT);

julia> PHI = hom(RQL, TL, TL.([t, t^2, t^3]));

julia> PSI = hom(TL, RQL, RQL.([x]));

julia> PHI(RQL(z))
t^3//1

julia> PSI(TL(t))
x//1
```
"""
hom(W::MPolyLocalizedRing, S::Ring, res::Map; check::Bool=true) = MPolyLocalizedRingHom(W, S, res, check=check)
hom(W::MPolyLocalizedRing, S::Ring, a::Vector{RET}; check::Bool=true) where {RET<:RingElem} = MPolyLocalizedRingHom(W, S, a, check=check)

### required getter functions
domain(PHI::MPolyLocalizedRingHom) = PHI.W
codomain(PHI::MPolyLocalizedRingHom) = PHI.S

@doc Markdown.doc"""
    restricted_map(PHI::MPolyLocalizedRingHom)

    restricted_map(PHI::MPolyQuoLocalizedRingHom)

Given a ring homomorphism `PHI` starting from a localized multivariate polynomial ring
(a localized quotient of a multivariate polynomial ring), return the composition of `PHI` 
with the localization map (with the composition of the localization map and the projection map).

# Examples
```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"]);

julia> I = ideal(R, [y-x^2, z-x^3]);

julia> RQ, p = quo(R, I);

julia> UR = complement_of_ideal(R, [0, 0, 0]);

julia> RQL, _ = localization(RQ, UR);

julia> T, (t,) =  PolynomialRing(QQ, ["t"]);

julia> UT = complement_of_ideal(T, [0]);

julia> TL, _ =  localization(T, UT);

julia> PHI = hom(RQL, TL, TL.([t, t^2, t^3]));

julia> PSI = hom(TL, RQL, RQL.([x]));

julia> phi = restricted_map(PHI)
Map with following data
Domain:
=======
Multivariate Polynomial Ring in x, y, z over Rational Field
Codomain:
=========
localization of Multivariate Polynomial Ring in t over Rational Field at the complement of maximal ideal corresponding to point with coordinates fmpq[0]

julia> psi = restricted_map(PSI)
Map with following data
Domain:
=======
Multivariate Polynomial Ring in t over Rational Field
Codomain:
=========
Localization of Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(-x^2 + y, -x^3 + z) at the multiplicative set complement of maximal ideal corresponding to point with coordinates fmpq[0, 0, 0]
```
"""
restricted_map(PHI::MPolyLocalizedRingHom) = PHI.res

### implementing the Oscar map interface
function identity_map(W::T) where {T<:MPolyLocalizedRing} 
  MPolyLocalizedRingHom(W, W, identity_map(base_ring(W)))
end

function compose(
    f::MPolyLocalizedRingHom, 
    g::Hecke.Map{<:Ring, <:Ring}
  )
  codomain(f) === domain(g) || error("maps are not compatible")
  R = base_ring(domain(f))
  return MPolyLocalizedRingHom(domain(f), codomain(g), hom(R, codomain(g), [g(f(x)) for x in gens(R)]))
end

(f::MPolyLocalizedRingHom)(I::Ideal) = ideal(codomain(f), domain(f).(gens(I)))

function ==(f::MPolyLocalizedRingHom, g::MPolyLocalizedRingHom) 
  domain(f) === domain(g) || return false
  codomain(f) === codomain(g) || return false
  for x in gens(base_ring(domain(f)))
    f(x) == g(x) || return false
  end
  return true
end

function divides(a::MPolyLocalizedRingElem, b::MPolyLocalizedRingElem)
  W = parent(a)
  W == parent(b) || error("elements do not belong to the same ring")
  F = FreeMod(W, 1)
  A = MatrixSpace(W, 1, 1)([b])
  M, _ = sub(F, A)
  represents_element(a*F[1], M) || return (false, zero(W))
  x = coordinates(a*F[1], M)
  return true, W(x[1])
end

# This had to be moved after the definition of the elements.
function Localization(R::MPolyRing, f::MPolyElem)
  U = MPolyPowersOfElement(R, [f])
  L = MPolyLocalizedRing(R, U)
  function func(a::MPolyElem)
    parent(a) == R || error("element does not belong to the correct ring")
    return L(a)
  end
  function func_inv(a::MPolyLocalizedRingElem{<:Any, <:Any, <:Any, <:Any, 
                                              <:MPolyPowersOfElement}
    )
    L == parent(a) || error("element does not belong to the correct ring")
    isone(denominator(a)) && return numerator(a)
    return divexact(numerator(a), denominator(a))
  end
  return L, MapFromFunc(func, func_inv, R, L)
end

#############################################################################
# Further functionality for elements of localized rings
#############################################################################

function derivative(f::MPolyLocalizedRingElem, i::Int)
  num = derivative(numerator(f), i)*denominator(f) - derivative(denominator(f), i)*numerator(f)
  den = denominator(f)^2
  g = gcd(num, den)
  return parent(f)(divexact(num, g), divexact(den, g), check=false)
end

function jacobi_matrix(f::MPolyLocalizedRingElem)
  L = parent(f)
  n = nvars(base_ring(L))
  return matrix(L, n, 1, [derivative(f, i) for i=1:n])
end

function jacobi_matrix(g::Vector{<:MPolyLocalizedRingElem})
  R = parent(g[1])
  n = nvars(base_ring(R))
  @assert all(x->parent(x) == R, g)
  return matrix(R, n, length(g), [derivative(x, i) for i=1:n for x = g])
end

function _is_integral_domain(R::MPolyLocalizedRing)
  return true
end

### Method compatible with `localized_ring`: This is the trivial lift to itself. 
function lift(a::MPolyLocalizedRingElem)
  return a
end

### Allow computation of kernels of maps to MPolyLocalizedRings
function kernel(f::MPolyAnyMap{<:MPolyRing, <:MPolyLocalizedRing})
  P = domain(f)
  L = codomain(f)
  I = ideal(L, zero(L))
  R = base_ring(L)
  J = saturated_ideal(I)
  d = [lifted_denominator(g) for g in f.(gens(domain(f)))]
  W = MPolyQuoLocalizedRing(R, ideal(R, zero(R)), MPolyPowersOfElement(R, d))
  id =  _as_affine_algebra(W)
  A = codomain(id)
  h = hom(P, A, id.(f.(gens(P))))
  return preimage(h, ideal(A, id.(W.(gens(J)))))
end

### For introducing in Function to docu
########################################

@Markdown.doc """

    in(f::MPolyElem, U::AbsMPolyMultSet)   

Return `true` if `f` is contained in `U`, `false` otherwise.

# Examples
```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"]);

julia> S = complement_of_ideal(R, [0, 0 ,0])
complement of maximal ideal corresponding to point with coordinates fmpq[0, 0, 0]

julia> y in S
false

julia> P = ideal(R, [x])
ideal(x)

julia> T = complement_of_ideal(P)
complement of ideal(x)

julia> y in T
true

julia> U = powers_of_element(x)
powers of fmpq_mpoly[x]

julia> x^3 in U
true

julia> (1+y)*x^2 in product(S, U)
true
```
""" Base.in(f::MPolyElem, U::AbsMPolyMultSet) 

