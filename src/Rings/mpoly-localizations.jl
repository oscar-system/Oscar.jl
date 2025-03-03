import Base: issubset

import AbstractAlgebra: Ring, RingElem

########################################################################
# General framework for localizations of multivariate polynomial rings #
########################################################################

########################################################################
# Powers of elements                                                   #
########################################################################

### required getter functions
ring(S::MPolyPowersOfElement) = S.R

### additional constructors
MPolyPowersOfElement(f::RET) where {RET<:MPolyRingElem} = MPolyPowersOfElement(parent(f), [f])
units_of(R::RT) where {RT<:MPolyRing} = MPolyPowersOfElement(R, [one(R)])

### additional functionality
denominators(S::MPolyPowersOfElement) = copy(S.a)

### required functionality
function Base.in(
    f::RingElemType, 
    S::MPolyPowersOfElement{BaseRingType, BaseRingElemType, RingType, RingElemType}
  ) where {BaseRingType, BaseRingElemType, RingType, RingElemType}
  R = parent(f)
  R === ring(S) || return false
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
Base.iterate(U::MPolyPowersOfElement, a::Tuple{<:MPolyRingElem, Int}) = (a[2] < length(U.a) ? (U.a[a[2]+1], a[2]+1) : nothing)
Base.iterate(U::MPolyPowersOfElement, i::Int) = (i < length(U.a) ? (U.a[i+1], i+1) : nothing)

is_trivial(U::MPolyPowersOfElement) = (U == units_of(ring(U)))

### printing

function Base.show(io::IO, ::MIME"text/plain", S::MPolyPowersOfElement)
  io = pretty(io)
  println(io, "Multiplicative subset")
  print(io, Indent())
  println(io, "of ", Lowercase(), ring(S))
  print(io, "given by the products of ")
  print(io, "[")
  join(io, denominators(S), ", ")
  print(io, "]")
  print(io, Dedent())
end

function Base.show(io::IO, S::MPolyPowersOfElement)
  if is_terse(io)
    print(io, "Products of ")
    print(io, ItemQuantity(length(denominators(S)), "element"))
  else
    print(io, "Products of (")
    join(io, denominators(S), ", ")
    print(io, ")")
  end
end

### generation of random elements 
function rand(S::MPolyPowersOfElement, v1::AbstractUnitRange{Int}, v2::AbstractUnitRange{Int}, v3::AbstractUnitRange{Int})
  R = ring(S)
  return prod(f^rand(0:2) for f in denominators(S); init = one(R))::elem_type(R)
end

function rand(rng::Random.AbstractRNG, S::MPolyPowersOfElement, v1::AbstractUnitRange{Int}, v2::AbstractUnitRange{Int}, v3::AbstractUnitRange{Int})
  R = ring(S)
  return prod(f^rand(rng, 0:2) for f in denominators(S); init = one(R))::elem_type(R)
end

### simplification. 
# Replaces each element d by its list of square free prime divisors.
function simplify!(S::MPolyPowersOfElement)
  S.is_simplified && return S
  S.a = denominators(simplify(S))
  S.is_simplified = true
  return S
end

function simplify(S::MPolyPowersOfElement)
  if S.is_simplified 
    result = MPolyPowersOfElement(ring(S), denominators(S))
    result.is_simplified = true
    return result
  end
  R = ring(S)
  new_denom = Vector{elem_type(R)}()
  for d in denominators(S)
    for (a, _) in factor(d)
      a in new_denom && continue
      push!(new_denom, a)
    end
  end
  return MPolyPowersOfElement(R, new_denom)
end

function simplify_light!(S::MPolyPowersOfElement)
  (S.is_simplified || S.is_lightly_simplified) && return S
  S.a = denominators(simplify_light(S))
  return S
end

function simplify_light(S::MPolyPowersOfElement)
  if S.is_simplified || S.is_lightly_simplified
    result = MPolyPowersOfElement(ring(S), denominators(S))
    result.is_lightly_simplified = true
    return result
  end
    
  # TODO: replace this by coprime_base once this has a method for multivariate polynomial rings.
  R = ring(S)
  new_denom = Vector{elem_type(R)}()

  for d in denominators(S)
    # get rid of the square free part
    for j in 1:ngens(R)
      dd = derivative(d, j)
      is_zero(dd) && continue
      c = gcd(d, dd)
      if !isone(c) 
        push!(new_denom, divexact(d, c))
        push!(new_denom, c)
      end
    end
    push!(new_denom, d)
  end

  new_denom = elem_type(R)[a for a in unique!(new_denom) if !is_unit(a)]
  fac_denom = elem_type(R)[]

  # Try to factor the denominators if this is deemed possible
  for a in new_denom
    if total_degree(a) > 50 || length(a) > 10000 || coefficient_ring(R) isa Hecke.RelSimpleNumField{AbsSimpleNumFieldElem}
      push!(fac_denom, a)
    else
      for (b, _) in factor(a)
        push!(fac_denom, b)
      end
    end
  end

  new_denom = fac_denom

  # Otherwise, it is still possible to play one denominator against another 
  # with gcd to get more factors in a relatively efficient way.
  simplified = true
  while simplified
    new_denom = elem_type(R)[a for a in unique!(new_denom) if !is_unit(a)]
    sort!(new_denom, by=total_degree)
    simplified = false
    for i in 1:length(new_denom)-1
      for j in i+1:length(new_denom)
        success, q = divides(new_denom[j], new_denom[i])
        if success
          new_denom[j] = q
          simplified = true
        end
      end
    end
  end

  result = MPolyPowersOfElement(R, new_denom)
  result.is_lightly_simplified = true
  return result
end


########################################################################
# Complements of prime ideals                                          #
########################################################################

### required getter functions
ring(S::MPolyComplementOfPrimeIdeal) = S.R

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
function Base.show(io::IO, ::MIME"text/plain", S::MPolyComplementOfPrimeIdeal)
  io = pretty(io)
  println(io, "Complement")
  println(io, Indent(), "of prime ", Lowercase(), prime_ideal(S))
  print(io, "in ", Lowercase(), ring(S))
  print(io, Dedent())
end

function Base.show(io::IO, S::MPolyComplementOfPrimeIdeal)
  if is_terse(io)
    print(io, "Complement of prime ideal")
  else
    io = pretty(io)
    print(io, "Complement of prime ", Lowercase(), prime_ideal(S))
  end
end

### generation of random elements 
function rand(S::MPolyComplementOfPrimeIdeal, v1::AbstractUnitRange{Int}, v2::AbstractUnitRange{Int}, v3::AbstractUnitRange{Int})
  f = rand(ring(S), v1, v2, v3)
  if f in prime_ideal(S)
    return rand(S, v1, v2, v3)
  end
  return f
end

function rand(rng::Random.AbstractRNG, S::MPolyComplementOfPrimeIdeal, v1::AbstractUnitRange{Int}, v2::AbstractUnitRange{Int}, v3::AbstractUnitRange{Int})
  f = rand(rng, ring(S), v1, v2, v3)
  if f in prime_ideal(S)
    return rand(rng, S, v1, v2, v3)
  end
  return f
end


########################################################################
# Complements of maximal ideals corresponding to 𝕜-points              #
########################################################################

# A function with this name exists with a method for 
# MPolyComplementOfPrimeIdeal. In order to treat them in 
# parallel, we introduce the analogous method here. 
function prime_ideal(S::MPolyComplementOfKPointIdeal) 
  kk = coefficient_ring(ring(S))
  # TODO: Test whether kk is an integral domain
  if !isdefined(S, :m)
    R = ring(S)
    x = gens(R)
    n = ngens(R)
    m = ideal(R, [x[i]-a for (i, a) in enumerate(point_coordinates(S))])
    set_attribute!(m, :is_prime=>true)
    kk isa Field && set_attribute!(m, :is_maximal=>true)
    set_attribute!(m, :dim=>dim(coefficient_ring(ring(S))))
    S.m = m
  end
  return S.m
end

dim(kk::Field) = 0

maximal_ideal(S::MPolyComplementOfKPointIdeal{<:Field}) = prime_ideal(S)

@doc raw"""  
    complement_of_point_ideal(R::MPolyRing, a::Vector)

Given a polynomial ring ``R``, say ``R = K[x_1,\dots, x_n]``, and given a vector 
``a = (a_1, \dots, a_n)`` of ``n`` elements of ``K``, return the multiplicatively 
closed subset ``R\setminus m``, where ``m`` is the maximal ideal 

$$m = \langle x_1-a_1,\dots, x_n-a_n\rangle \subset R.$$

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> U = complement_of_point_ideal(R, [0, 0 ,0])
Complement
  of maximal ideal corresponding to rational point with coordinates (0, 0, 0)
  in multivariate polynomial ring in 3 variables over QQ

```
"""
complement_of_point_ideal(R::MPolyRing, a::Vector) = MPolyComplementOfKPointIdeal(R, a)

@doc raw"""
    complement_of_prime_ideal(P::MPolyIdeal; check::Bool=true)

Given a prime ideal ``P`` of a polynomial ring ``R``, say,
return the multiplicatively closed subset ``R\setminus P.``

!!! note
    Since  `check` is set to `true`, the function checks whether ``P`` is indeed a prime ideal. 

    This may take some time.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> P = ideal(R, [x])
Ideal generated by
  x

julia> U = complement_of_prime_ideal(P)
Complement
  of prime ideal (x)
  in multivariate polynomial ring in 3 variables over QQ
```
"""
complement_of_prime_ideal(P::MPolyIdeal; check::Bool=true) = MPolyComplementOfPrimeIdeal(P; check)

@doc raw"""  
    powers_of_element(f::MPolyRingElem)

Given an element `f` of a polynomial ring, return the multiplicatively 
closed subset of the polynomial ring which is formed by the powers of `f`.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> f = x
x

julia> U = powers_of_element(f)
Multiplicative subset
  of multivariate polynomial ring in 3 variables over QQ
  given by the products of [x]
```
"""
powers_of_element(f::MPolyRingElem) = MPolyPowersOfElement(f)

### required getter functions
ring(S::MPolyComplementOfKPointIdeal) = S.R

### additional getter functions 
point_coordinates(S::MPolyComplementOfKPointIdeal) = S.a

### required functionality
function Base.in(
    f::RingElemType, 
    S::MPolyComplementOfKPointIdeal{BaseRingType, BaseRingElemType, RingType, RingElemType}
  ) where {BaseRingType, BaseRingElemType, RingType, RingElemType}
  parent(f) == ring(S) || return false
  return !(evaluate(f, point_coordinates(S)) == zero(ring(S)))
end

### printing
function Base.show(io::IO, ::MIME"text/plain", S::MPolyComplementOfKPointIdeal)
  io = pretty(io)
  println(io, "Complement")
  print(io, Indent(), "of maximal ideal corresponding to rational point ")
  print(io, Indent(), "with coordinates ")
  print(io, "(")
  join(io, point_coordinates(S), ", ")
  println(io, ")")
  print(io, Dedent())
  print(io, "in ", Lowercase(), ring(S))
  print(io, Dedent())
end

function Base.show(io::IO, S::MPolyComplementOfKPointIdeal)
  if is_terse(io)
    print(io, "Complement of maximal ideal")
  else
    print(io, "Complement of maximal ideal of point ")
    print(io, "(")
    join(io, point_coordinates(S), ", ")
    print(io, ")")
  end
end

### generation of random elements 
function rand(S::MPolyComplementOfKPointIdeal, v1::AbstractUnitRange{Int}, v2::AbstractUnitRange{Int}, v3::AbstractUnitRange{Int})
  f = rand(ring(S), v1, v2, v3)
  if !(f in S)
    return rand(S, v1, v2, v3)
  end
  return f
end
function rand(rng::Random.AbstractRNG, S::MPolyComplementOfKPointIdeal, v1::AbstractUnitRange{Int}, v2::AbstractUnitRange{Int}, v3::AbstractUnitRange{Int})
  f = rand(rng, ring(S), v1, v2, v3)
  if !(f in S)
    return rand(rng, S, v1, v2, v3)
  end
  return f
end


### required getter functions
ring(S::MPolyProductOfMultSets) = S.R

### additional functionality
getindex(S::MPolyProductOfMultSets, i::Integer) = S.U[i]
sets(S::MPolyProductOfMultSets) = copy(S.U)

### required functionality
function Base.in(
    f::RingElemType, 
    S::MPolyProductOfMultSets{BaseRingType, BaseRingElemType, RingType, RingElemType}
  ) where {BaseRingType, BaseRingElemType, RingType, RingElemType}
  R = ring(S)
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
function Base.show(io::IO, ::MIME"text/plain", S::MPolyProductOfMultSets)
  io = pretty(io)
  println(io, "Product of the multiplicative sets")
  print(io, Indent())
  M = sets(S)
  for i in 1:length(M)-1
    println(io, Lowercase(), M[i])
  end
  print(io, Lowercase(), M[end])
  print(io, Dedent())
end

function Base.show(io::IO, S::MPolyProductOfMultSets)
  io = pretty(io)
  if is_terse(io)
    print(io, "Product of multiplicatively closed subsets")
  else
    print(io, "Product of the multiplicative subsets [")
    for s in sets(S)
      print(terse(io), Lowercase(), s)
      s !== last(sets(S)) && print(io, ", ")
    end
    print(io, "]")
  end
end

### generation of random elements 
function rand(S::MPolyProductOfMultSets, v1::AbstractUnitRange{Int}, v2::AbstractUnitRange{Int}, v3::AbstractUnitRange{Int})
  return prod([rand(s, v1, v2, v3) for s in sets(S)])::elem_type(ring(S))
end
function rand(rng::Random.AbstractRNG, S::MPolyProductOfMultSets, v1::AbstractUnitRange{Int}, v2::AbstractUnitRange{Int}, v3::AbstractUnitRange{Int})
  return prod([rand(rng, s, v1, v2, v3) for s in sets(S)])::elem_type(ring(S))
end

########################################################################
# Localization associated to a monomial ordering                       #
########################################################################

### required getter functions
ring(S::MPolyLeadingMonOne) = S.R

### additional constructors
MPolyLeadingMonOne(ord::MonomialOrdering) = MPolyLeadingMonOne(ord.R, ord)

ordering(S::MPolyLeadingMonOne) = S.ord

### required functionality
function Base.in(
    f::RingElemType,
    S::MPolyLeadingMonOne{BaseRingType, BaseRingElemType, RingType, RingElemType}
  ) where {BaseRingType, BaseRingElemType, RingType, RingElemType}
  R = parent(f)
  R == ring(S) || return false
  if iszero(f)
    return false
  end
  return isone(leading_monomial(f; ordering = ordering(S)))
end

### printing
function Base.show(io::IO,::MIME"text/plain", S::MPolyLeadingMonOne)
  io = pretty(io)
  println(io, "Elements")
  print(io, Indent())
  println(io, "of ", ring(S))
  print(io, Dedent())
  println(io, "with leading monomial 1 w.r.t. ")
  print(io, Indent())
  print(io, ordering(S))
  print(io, Dedent())
end

function Base.show(io::IO, S::MPolyLeadingMonOne)
  io = pretty(io)
  if is_terse(io)
    print(io, "Elements with leading monomial 1")
  else
    print(io, "Elements with leading monomial 1 w.r.t. ")
    print(terse(io), ordering(S))
  end
end

########################################################################
# Arithmetic of multiplicative sets                                    #
########################################################################

### containment ########################################################
⊂(T::AbsMPolyMultSet, U::AbsMPolyMultSet) = issubset(T, U)

function issubset(T::AbsMPolyMultSet, U::AbsMPolyMultSet) 
  ring(T) == ring(U) || return false
  error("comparison of multiplicative sets of type $(typeof(T)) and $(typeof(U)) is not implemented")
end

function ==(T::AbsMPolyMultSet, U::AbsMPolyMultSet) 
  T === U && return true
  return (issubset(T, U) && issubset(U, T))
end

function Base.hash(T::AbsMPolyMultSet{BRT, BRET, RT, RET}, h::UInt) where {BRT,BRET,RT,RET}
  error("not implemented")
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
  R = ring(T)
  R == ring(U) || error("multiplicative sets do not belong to the same ring")
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
  R = ring(T)
  R == ring(U) || error("multiplicative sets do not belong to the same ring")
  a = point_coordinates(T)
  for i in 1:length(a)
    (gen(R, i)- R(a[i])) in prime_ideal(U) || return false
  end
  return true
end

function issubset(
    T::MPolyComplementOfPrimeIdeal{BRT, BRET, RT, RET},
    U::MPolyComplementOfKPointIdeal{BRT, BRET, RT, RET}
  ) where {BRT, BRET, RT, RET}
  R = ring(T)
  R == ring(U) || error("multiplicative sets do not belong to the same ring")
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

@doc raw"""
    product(T::AbsMPolyMultSet, U::AbsMPolyMultSet)

Return the product of the multiplicative subsets `T` and `U`. 

Alternatively, write `T*U`.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> T = complement_of_point_ideal(R, [0, 0 ,0])
Complement
  of maximal ideal corresponding to rational point with coordinates (0, 0, 0)
  in multivariate polynomial ring in 3 variables over QQ

julia> f = x
x

julia> U = powers_of_element(f)
Multiplicative subset
  of multivariate polynomial ring in 3 variables over QQ
  given by the products of [x]

julia> S = product(T, U)
Product of the multiplicative sets
  complement of maximal ideal of point (0, 0, 0)
  products of (x)
```
"""
function product(T::AbsMPolyMultSet, U::AbsMPolyMultSet)
  R = ring(T)
  R == ring(U) || error("multiplicative sets do not belong to the same ring")
  issubset(T, U) && return U
  issubset(U, T) && return T
  return MPolyProductOfMultSets(R, [T, U])
end

function product(T::MST, U::MST) where {MST<:MPolyProductOfMultSets} 
  R = ring(T)
  R == ring(U) || error("multiplicative sets do not belong to the same ring")
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
  R = ring(T)
  R == ring(U) || error("multiplicative sets do not belong to the same ring")
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
  R = ring(T)
  R == ring(U) || error("multiplicative sets do not belong to the same ring")
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
  R = ring(T)
  R == ring(U) || error("multiplicative sets do not belong to the same ring")
  a = point_coordinates(U)
  b = point_coordinates(T)
  for i in 1:length(a)
    a[i] == b[i] || return MPolyProductOfMultSets(R, [U, T])
  end
  return T
end

function product(T::MST, U::MST) where {MST<:MPolyComplementOfPrimeIdeal}
  R = ring(T)
  R == ring(U) || error("multiplicative sets do not belong to the same ring")
  if issubset(prime_ideal(T), prime_ideal(U))
    return T
  end
  if issubset(prime_ideal(U), prime_ideal(T))
    return U
  end
  return MPolyProductOfMultSets(R, [U, T])
end

function product(T::MST, U::MST) where {MST<:MPolyPowersOfElement}
  R = ring(T) 
  R == ring(U) || error("multiplicative sets do not belong to the same ring")
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
  R = ring(T)
  R == ring(U) || error("multiplicative sets do not belong to the same ring")
  for a in denominators(T)
    a in U || return MPolyProductOfMultSets(R, [U, T])
  end
  return U
end

function product(
    T::MPolyPowersOfElement{BRT, BRET, RT, RET},
    U::AbsMPolyMultSet{BRT, BRET, RT, RET}
  ) where {BRT, BRET, RT, RET}
  R = ring(T)
  R == ring(U) || error("multiplicative sets do not belong to the same ring")
  for a in denominators(T)
    a in U || return MPolyProductOfMultSets(R, [U, T])
  end
  return U
end

function product(
    T::MPolyPowersOfElement{BRT, BRET, RT, RET},
    U::MPolyProductOfMultSets{BRT, BRET, RT, RET}
  ) where {BRT, BRET, RT, RET}
  R = ring(T)
  R == ring(U) || error("multiplicative sets do not belong to the same ring")
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
  ring(U) === domain(phi) || error("multiplicative set does not lay in the domain of the morphism")
  S = codomain(phi) 
  SU = MPolyPowersOfElement(S, phi.(denominators(U)))
  return SU
end

function (phi::MPolyAnyMap{<:MPolyRing, <:MPolyRing, Nothing})(U::MPolyComplementOfPrimeIdeal;
                                                               check::Bool=true
                                                              )
  ring(U) === domain(phi) || error("multiplicative set does not lay in the domain of the morphism")
  S = codomain(phi) 
  Q = ideal(S, phi.(gens(prime_ideal(U))))
  SU = MPolyComplementOfPrimeIdeal(S, Q, check=check)
  return SU
end

########################################################################
# Localizations of polynomial rings over admissible fields             #
########################################################################

### required getter functions 
base_ring(W::MPolyLocRing) = W.R
base_ring_type(::Type{MPolyLocRing{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = BRT
inverted_set(W::MPolyLocRing) = W.S

### additional getter functions
gens(W::MPolyLocRing) = W.(gens(base_ring(W)))
number_of_generators(W::MPolyLocRing) = number_of_generators(base_ring(W))

### required extension of the localization function
@doc raw"""

    localization(R::MPolyRing, U::AbsMPolyMultSet)   

Return the localization of `R` at `U`, together with the localization map.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> P = ideal(R, [x])
Ideal generated by
  x

julia> U = complement_of_prime_ideal(P)
Complement
  of prime ideal (x)
  in multivariate polynomial ring in 3 variables over QQ

julia> Rloc, iota = localization(R, U);

julia> Rloc
Localization
  of multivariate polynomial ring in 3 variables x, y, z
    over rational field
  at complement of prime ideal (x)

julia> iota
Ring homomorphism
  from multivariate polynomial ring in 3 variables over QQ
  to localization of multivariate polynomial ring in 3 variables over QQ at complement of prime ideal (x)
defined by
  x -> x
  y -> y
  z -> z
```
""" localization(R::MPolyRing, U::AbsMPolyMultSet)

function localization(S::AbsMPolyMultSet)
    R = ring(S)
    Rloc = MPolyLocRing(R, S)
    #iota = MapFromFunc(R, Rloc, Rloc)
    iota = hom(R, Rloc, Rloc.(gens(R)), check=false)
    return Rloc, iota
end

function localization(R::MPolyRing, ord::MonomialOrdering)
  @assert R === ord.R
  return localization(MPolyLeadingMonOne(ord))
end

### Successive localizations are handled by the dispatch for products
function localization(
    W::MPolyLocRing{BRT, BRET, RT, RET, MST}, 
    S::AbsMPolyMultSet{BRT, BRET, RT, RET}
  ) where {BRT, BRET, RT, RET, MST}
  issubset(S, inverted_set(W)) && return W, identity_map(W)
  U = S*inverted_set(W)
  L, _ = localization(U)
  #return L, MapFromFunc(W, L, (x->(L(numerator(x), denominator(x), check=false))))
  return L, MPolyLocalizedRingHom(W, L, hom(base_ring(W), L, L.(gens(base_ring(W)))), check=false)
end

### additional constructors
MPolyLocRing(R::RingType, P::MPolyIdeal{RingElemType}; check::Bool=true) where {RingType, RingElemType} = MPolyLocRing(R, MPolyComplementOfPrimeIdeal(P); check)

localization(R::MPolyRing, v::Vector{T}) where {T<:MPolyRingElem} = localization(MPolyPowersOfElement(R, v))

function localization(
    W::MPolyLocRing{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}}, 
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
  L = MPolyLocRing(R, MPolyPowersOfElement(R, vcat(g, divexact(f, h))))
  return L, MPolyLocalizedRingHom(W, L, hom(base_ring(W), L, L.(gens(base_ring(W)))), check=false)
  #return L, MapFromFunc(W, L, (x->L(numerator(x), denominator(x), check=false)))
end

function localization(
    W::MPolyLocRing{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}}, 
    v::Vector{RET}
  ) where {BRT, BRET, RT, RET}
  V = W
  for f in v
    V = localization(V, f)
  end
  return V, MPolyLocalizedRingHom(W, V, hom(base_ring(W), V, V.(gens(base_ring(W)))), check=false)
  #return V, MapFromFunc(W, V, (x->V(numerator(x), denominator(x), check=false)))
end

### generation of random elements 
function rand(W::MPolyLocRing, v1::AbstractUnitRange{Int}, v2::AbstractUnitRange{Int}, v3::AbstractUnitRange{Int})
  return W(rand(base_ring(W), v1, v2, v3), rand(inverted_set(W), v1, v2, v3))
end

function rand(rng::Random.AbstractRNG, W::MPolyLocRing, v1::AbstractUnitRange{Int}, v2::AbstractUnitRange{Int}, v3::AbstractUnitRange{Int})
  return W(rand(rng, base_ring(W), v1, v2, v3), rand(rng, inverted_set(W), v1, v2, v3))
end

### Conformance test element generation
function ConformanceTests.generate_element(W::Oscar.MPolyLocRing)
  return rand(W, 0:3, 0:4, 0:3)
end


########################################################################
# Elements of localized polynomial rings                               #
########################################################################

### required getter functions 
numerator(a::MPolyLocRingElem) = numerator(a.frac)

denominator(a::MPolyLocRingElem) = denominator(a.frac)

parent(a::MPolyLocRingElem) = a.W

### additional getter functions
fraction(a::MPolyLocRingElem) = a.frac
# to assure compatibility with generic code for MPolyQuoLocalizedRings:
lifted_numerator(a::MPolyLocRingElem) = numerator(a)
lifted_denominator(a::MPolyLocRingElem) = denominator(a)

### required conversions
function (W::MPolyLocRing{
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
  return MPolyLocRingElem(W, fraction_field(base_ring(W))(f), check=false)
end

(W::MPolyLocRing{
    BaseRingType, 
    BaseRingElemType, 
    RingType, 
    RingElemType, 
    MultSetType
  })(f::BaseRingElemType; check::Bool=true) where {
    BaseRingType<:Ring, 
    BaseRingElemType<:RingElem, 
    RingType<:Ring, 
    RingElemType<:RingElem, 
    MultSetType<:AbsMultSet
  } = MPolyLocRingElem(W, fraction_field(base_ring(W))(base_ring(W)(f)), check=check)

# Remove ambiguities
(W::MPolyLocRing{
    BaseRingType, 
    BaseRingElemType, 
    RingType, 
    RingElemType, 
    MultSetType
  })(f::BaseRingElemType) where {
    BaseRingType<:Ring, 
    BaseRingElemType<:ZZRingElem, 
    RingType<:Ring, 
    RingElemType<:RingElem, 
    MultSetType<:AbsMultSet
  } = MPolyLocRingElem(W, fraction_field(base_ring(W))(base_ring(W)(f)), check=false)

function (W::MPolyLocRing{
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
  return MPolyLocRingElem(W, a//b, check=check)
end

### additional conversions
(W::MPolyLocRing{
    BaseRingType, 
    BaseRingElemType, 
    RingType, 
    RingElemType, 
    MultSetType
  })(f::AbstractAlgebra.Generic.FracFieldElem{RingElemType}; check::Bool=true) where {
    BaseRingType, 
    BaseRingElemType, 
    RingType, 
    RingElemType, 
    MultSetType
  } = MPolyLocRingElem(W, f, check=check)

### additional promotions 
AbstractAlgebra.promote_rule(::Type{RET}, ::Type{MPolyLocRingElem{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT<:Ring, RET<:RingElement, MST} = MPolyLocRingElem{BRT, BRET, RT, RET, MST}

AbstractAlgebra.promote_rule(::Type{BRET}, ::Type{MPolyLocRingElem{BRT, BRET, RT, RET, MST}}) where {BRT<:Ring, BRET<:RingElement, RT<:Ring, RET<:RingElement, MST} = MPolyLocRingElem{BRT, BRET, RT, RET, MST}

AbstractAlgebra.promote_rule(::Type{Integer}, ::Type{MPolyLocRingElem{BRT, BRET, RT, RET, MST}}) where {BRT<:Ring, BRET<:RingElement, RT<:Ring, RET<:RingElement, MST} = MPolyLocRingElem{BRT, BRET, RT, RET, MST}

AbstractAlgebra.promote_rule(::Type{ZZRingElem}, ::Type{MPolyLocRingElem{BRT, BRET, RT, RET, MST}}) where {BRT<:Ring, BRET<:RingElement, RT<:Ring, RET<:RingElement, MST} = MPolyLocRingElem{BRT, BRET, RT, RET, MST}

### overwriting the arithmetic using the fractions from AbstractAlgebra
function +(a::T, b::T) where {T<:MPolyLocRingElem}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  return (parent(a))(fraction(a) + fraction(b), check=false)
end

function -(a::T, b::T) where {T<:MPolyLocRingElem}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  return (parent(a))(fraction(a) - fraction(b), check=false)
end

function -(a::T) where {T<:MPolyLocRingElem}
  return (parent(a))((-1)*fraction(a), check=false)
end

function *(a::T, b::T) where {T<:MPolyLocRingElem}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  return (parent(a))(fraction(a) * fraction(b), check=false)
end

function *(a::RET, b::MPolyLocRingElem{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET <: RingElem, MST}
  return (parent(b))(a*fraction(b), check=false)
end

function *(a::MPolyLocRingElem{BRT, BRET, RT, RET, MST}, b::RET) where {BRT, BRET, RT, RET <: RingElem, MST}
  return b*a
end

function *(a::BRET, b::MPolyLocRingElem{BRT, BRET, RT, RET, MST}) where {BRT, BRET <: RingElem, RT, RET, MST}
  return (parent(b))(a*fraction(b), check=false)
end

function *(a::MPolyLocRingElem{BRT, BRET, RT, RET, MST}, b::BRET) where {BRT, BRET <: RingElem, RT, RET, MST}
  return b*a
end

function Base.:(/)(a::Integer, b::T) where {T<:MPolyLocRingElem}
  return (parent(b))(a//fraction(b))
end

function Base.:(/)(a::ZZRingElem, b::T) where {T<:MPolyLocRingElem}
  return (parent(b))(a//fraction(b))
end

function ==(a::T, b::T) where {T<:MPolyLocRingElem}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  return fraction(a) == fraction(b)
end

# We need to manually split this into three methods, because 
# otherwise it seems that Julia can not dispatch this function.
function ^(a::MPolyLocRingElem, i::Int64)
  return parent(a)(fraction(a)^i, check=false)
end
function ^(a::MPolyLocRingElem, i::Integer)
  return parent(a)(fraction(a)^i, check=false)
end
function ^(a::MPolyLocRingElem, i::ZZRingElem)
  return parent(a)(fraction(a)^i, check=false)
end

function divexact(p::T, q::T; check::Bool=false) where {T<:MPolyLocRingElem} 
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

@doc raw"""
    is_unit(f::MPolyLocRingElem)

Return `true`, if `f` is a unit of `parent(f)`, `false` otherwise.

# Examples

```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> P = ideal(R, [x])
Ideal generated by
  x

julia> U = complement_of_prime_ideal(P)
Complement
  of prime ideal (x)
  in multivariate polynomial ring in 3 variables over QQ

julia> Rloc, iota = localization(U);

julia> is_unit(iota(x))
false

julia> is_unit(iota(y))
true
```
""" 
is_unit(f::MPolyLocRingElem) = numerator(f) in inverted_set(parent(f))

function is_zero_divisor(f::MPolyLocRingElem)
  iszero(f) && return true
  if is_constant(f)
    c = first(AbstractAlgebra.coefficients(numerator(f)))
    return is_zero_divisor(c)
  end
  return is_zero_divisor(numerator(f))
end


function is_constant(a::MPolyLocRingElem)
  reduce_fraction(a)
  return is_constant(numerator(a)) && is_constant(denominator(a))
end

########################################################################
# implementation of Oscar's general ring interface                     #
########################################################################

one(W::MPolyLocRing) = W(one(base_ring(W)))
zero(W::MPolyLocRing) = W(zero(base_ring(W)))

elem_type(::Type{MPolyLocRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyLocRingElem{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}
parent_type(::Type{MPolyLocRingElem{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyLocRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}

function (W::MPolyLocRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType})(f::MPolyLocRingElem{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}; check::Bool=true) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} 
  parent(f) === W && return f
  return MPolyLocRingElem(W, fraction(f), check=check)
end

(W::MPolyLocRing)() = zero(W)
(W::MPolyLocRing)(a::Integer) = W(base_ring(W)(a), check=false)
(W::MPolyLocRing)(a::Int64) = W(base_ring(W)(a), check=false)
(W::MPolyLocRing)(a::ZZRingElem) = W(base_ring(W)(a), check=false)

is_domain_type(T::Type{MPolyLocRingElem{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = true 
is_exact_type(T::Type{MPolyLocRingElem{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = true

@attr Any function base_ring_shifts(L::MPolyLocRing{<:Any, <:Any, <:Any, <:Any, <:MPolyComplementOfKPointIdeal}) 
  a = point_coordinates(inverted_set(L))
  R = base_ring(L)
  shift = hom(R, R, gens(R)+R.(a), check=false)
  back_shift = hom(R, R, gens(R)-R.(a), check=false)
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

### required getter functions
gens(I::MPolyLocalizedIdeal) = copy(I.gens)
gen(I::MPolyLocalizedIdeal, i::Int) = I.gens[i]
base_ring(I::MPolyLocalizedIdeal) = I.W
base_ring_type(::Type{MPolyLocalizedIdeal{LRT, LRET}}) where {LRT, LRET} = LRT

### type getters
ideal_type(::Type{MPolyLocalizedRingType}) where {MPolyLocalizedRingType<:MPolyLocRing} = MPolyLocalizedIdeal{MPolyLocalizedRingType, elem_type(MPolyLocalizedRingType)}
ideal_type(L::MPolyLocRing) = ideal_type(typeof(L))

### additional getter functions 
map_from_base_ring(I::MPolyLocalizedIdeal) = I.map_from_base_ring
is_saturated(I::MPolyLocalizedIdeal) = I.is_saturated
number_of_generators(I::MPolyLocalizedIdeal) = length(I.gens)

function ideal_membership(a::RingElem, I::MPolyLocalizedIdeal)
  L = base_ring(I)
  parent(a) == L || return L(a) in I
  b = numerator(a)
  b in pre_saturated_ideal(I) && return true
  is_saturated(I) && return false
  R = base_ring(L)
  J = pre_saturated_ideal(I)
  (success, x, u) = has_solution(generator_matrix(J), matrix(R, 1, 1, [b]), inverted_set(L))
  !success && return false
  # cache the intermediate result
  extend_pre_saturated_ideal!(I, b, x, u, check=false)
  return true
end

# TODO: Also add a special dispatch for localizations at 𝕜-points
@attr Bool function is_prime(I::MPolyLocalizedIdeal)
  return is_prime(saturated_ideal(I))
end

### Additional constructors
@doc raw"""
    intersect(I::MPolyLocalizedIdeal, J::MPolyLocalizedIdeal)
    intersect(I::MPolyLocalizedIdeal, J::MPolyLocalizedIdeal...)
    intersect(VI::Vector{<:MPolyLocalizedIdeal{T}}) where T
Return the intersection of two or more ideals.

# Examples
```jldoctest
julia> R,(x,y,z,w) = QQ[:x, :y, :z, :w];

julia> f = x+y+z+w-1;

julia> T = powers_of_element(f);

julia> RL,phiL = localization(R,T);

julia> I=ideal(RL,RL.([x+y+z,w-1]))
Ideal generated by
  x + y + z
  w - 1

julia> is_one(I)
true

julia> J=ideal(RL,RL.([x^2,y]))
Ideal generated by
  x^2
  y

julia> K=ideal(RL,RL.([x]))
Ideal generated by
  x

julia> intersect(J,K)
Ideal generated by
  x*y
  x^2

julia> intersect(I,J,K)
Ideal generated by
  x*y*w - x*y
  x^2*w - x^2
  x^2*y + x*y^2 + x*y*z
  x^3 + x^2*z - x*y^2 - x*y*z

julia> intersect([I,J,K])
Ideal generated by
  x*y*w - x*y
  x^2*w - x^2 + x*y*w - x*y
  x^2*y + x*y^2 + x*y*z
  x^3 + x^2*z - x*y^2 - x*y*z
```
"""
function intersect(I::MPolyLocalizedIdeal, J::MPolyLocalizedIdeal)
  L = base_ring(I)
  L == base_ring(J) || error("ideals must be defined in the same ring")
  preI = pre_saturated_ideal(I)
  preJ = pre_saturated_ideal(J)
  K = intersect(I, J)
  return L(K)
end

function intersect(I::MPolyLocalizedIdeal,J::MPolyLocalizedIdeal...) 
  L = base_ring(I)
  erg = pre_saturated_ideal(I)
  for K in J
    base_ring(K) == L || error("base rings must match")
    erg = intersect(erg,pre_saturated_ideal(K))
  end
  return L(erg)
end

function intersect(VI::Vector{<:MPolyLocalizedIdeal{T}}) where T
  @assert length(VI)!=0
  L = base_ring(VI[1])
  all(J -> base_ring(J) == L,VI) || error("base rings must match")
  VIpre = [pre_saturated_ideal(J) for J in VI]
  erg = Base.intersect(VIpre)
  return L(erg)
end

### Further functionality
function coordinates(a::RingElem, I::MPolyLocalizedIdeal; check::Bool=true)
  L = base_ring(I)
  parent(a) === L || return coordinates(L(a), I, check=check)
  @check a in I "the given element is not in the ideal"
  J = pre_saturated_ideal(I)
  R = base_ring(J)
  p = numerator(a)
  if p in J
    q = denominator(a)
    # caching has been done during the call of `in`, so the following will work
    x = coordinates(p, pre_saturated_ideal(I))
    # multiplications sparse*dense have to be carried out this way round.
    return transpose(pre_saturation_data(I) * transpose(L(one(q), q, check=false)*change_base_ring(L, x)))
  else
    (success, x, u) = has_solution(generator_matrix(J), matrix(R, 1, 1, [p]), inverted_set(L), check=false)
    !success && error("check for membership was disabled, but element is not in the ideal")
    # cache the intermediate result
    #result = L(one(R), u*denominator(a), check=false)*change_base_ring(L, x)*pre_saturation_data(I)
    result = transpose(pre_saturation_data(I) * transpose(L(one(R), u*denominator(a), check=false)*change_base_ring(L, x)))
    extend_pre_saturated_ideal!(I, p, x, u, check=false)
    return result
  end
end

function coordinates(
    a::RingElem, I::MPolyLocalizedIdeal{LRT}; check::Bool=true
  ) where {LRT<:MPolyLocRing{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement}}
  L = base_ring(I)
  parent(a) === L || return coordinates(L(a), I, check=check)

  b = numerator(a)
  if b in pre_saturated_ideal(I)
    x = coordinates(b, pre_saturated_ideal(I))
    q = denominator(a)
    return transpose(pre_saturation_data(I) * transpose(L(one(q), q, check=false)*change_base_ring(L, x)))
  end

  @check a in I "the given element is not in the ideal"
  !is_saturated(I) && _replace_pre_saturated_ideal(I, saturated_ideal(I), prod(denominators(inverted_set(L)); init=one(base_ring(L)))) # Computing the saturation first is cheaper than the generic Posur method
  J = pre_saturated_ideal(I)
  R = base_ring(J)
  p = numerator(a)
  x = coordinates(p, J)
  q = denominator(a)
  return transpose(pre_saturation_data(I) * transpose(L(one(q), q, check=false)*change_base_ring(L, x)))
end

generator_matrix(J::MPolyIdeal) = matrix(base_ring(J), ngens(J), 1, gens(J))

@doc raw"""
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
julia> R, (x,) = polynomial_ring(QQ, [:x]);

julia> U = powers_of_element(x)
Multiplicative subset
  of multivariate polynomial ring in 1 variable over QQ
  given by the products of [x]

julia> Rloc, iota = localization(R, U);

julia> I = ideal(Rloc, [x+x^2])
Ideal generated by
  x^2 + x

julia> SI = saturated_ideal(I)
Ideal generated by
  x + 1

julia> base_ring(SI)
Multivariate polynomial ring in 1 variable x
  over rational field

julia> U = complement_of_point_ideal(R, [0])
Complement
  of maximal ideal corresponding to rational point with coordinates (0)
  in multivariate polynomial ring in 1 variable over QQ

julia> Rloc, iota = localization(R, U);

julia> I = ideal(Rloc, [x+x^2])
Ideal generated by
  x^2 + x

julia> saturated_ideal(I)
Ideal generated by
  x
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
    r = ngens(I)
    A = zero_matrix(SMat, W, 0, r)
    #A = zero_matrix(W, r, r)
    for i in 1:r
      push!(A, sparse_row(W, [(i, W(denominator(gen(I, i))))]))
      #A[i, i] = denominator(gen(I, i))
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
  ) where {PT <: MPolyRingElem}
  nrows(x) == 1 || error("matrix must be a row vector")
  L = base_ring(I)
  R = base_ring(L)
  J = pre_saturated_ideal(I)
  @check u*f == dot(x, gens(J)) "input is not coherent"
  J_ext = ideal(R, vcat(gens(J), [f]))
  T = pre_saturation_data(I)
  y = T * transpose(L(one(u), u, check=false)*change_base_ring(L, x))
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
  ) where {PT <: MPolyRingElem}
  L = base_ring(I)
  R = base_ring(L)
  J = pre_saturated_ideal(I)
  n = length(f)
  n == length(u) == nrows(x) || error("input dimensions are not compatible")
  @check begin
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
  y = T * transpose(change_base_ring(L, x))
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

@attr Map function coordinate_shift(
    L::MPolyLocRing{<:Any, <:Any, <:Any, <:Any, <:MPolyComplementOfKPointIdeal}
  )
  U = inverted_set(L)
  a = point_coordinates(U)
  R = base_ring(L)
  Ls = MPolyLocRing(R, MPolyComplementOfKPointIdeal(R, [0 for i in 1:ngens(R)]))
  xs = [ x + a for (x, a) in zip(gens(base_ring(L)), a) ]
  xs_inv = [ x - a for (x, a) in zip(gens(base_ring(L)), a) ]
  shift = MapFromFunc(
              L, Ls,
              f -> Ls(evaluate(numerator(f), xs), evaluate(denominator(f), xs), check=false),
              g -> L(evaluate(numerator(g), xs_inv), evaluate(denominator(g), xs_inv), check=false),
            )
  return shift
end


function saturated_ideal(
    I::MPolyLocalizedIdeal{LRT};
    with_generator_transition::Bool=false
  ) where {LRT<:MPolyLocRing{<:Any, <:Any, <:Any, <:Any, <:MPolyComplementOfKPointIdeal}}
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
  ) where {LRT<:MPolyLocRing{<:Any, <:Any, <:Any, <:Any, <:MPolyComplementOfPrimeIdeal}}
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
    with_generator_transition::Bool=false
  ) where {LRT<:MPolyLocRing{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement}}
  if !isdefined(I, :saturated_ideal)
    is_saturated(I) && return pre_saturated_ideal(I)
    L = base_ring(I)
    R = base_ring(L)
    if iszero(I)
      I.saturated_ideal = ideal(R, elem_type(R)[])
      return I.saturated_ideal
    end
    if strategy==:iterative_saturation
      Jsat = pre_saturated_ideal(I)
      U = inverted_set(L)
      for d in denominators(U)
        if !is_unit(d) && !iszero(Jsat)
          Jsat = saturation(Jsat, ideal(R, d))
          if with_generator_transition
            _replace_pre_saturated_ideal(I, Jsat, d)
          end
        end
      end
      I.saturated_ideal = Jsat
    elseif strategy==:single_saturation
      d = prod(denominators(inverted_set(L)))
      Jsat = pre_saturated_ideal(I)
      if !is_unit(d) && !iszero(pre_saturated_ideal(I))
        Jsat = saturation(Jsat, ideal(R, d))
        if with_generator_transition
          _replace_pre_saturated_ideal(I, Jsat, d)
        end
      end
      I.saturated_ideal = Jsat
    else
      error("strategy $strategy not implemented")
    end
  end
  return I.saturated_ideal
end

@doc raw"""
    _replace_pre_saturated_ideal(
        I::MPolyLocalizedIdeal{LRT}, J::MPolyIdeal, d::MPolyRingElem
      ) where {LRT<:MPolyLocRing{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement}}

For ``I ⊂ R[S⁻¹]`` with `saturated_ideal` ``J ⊂ R`` this replaces the 
`pre_saturated_ideal` of ``I`` by ``J`` and computes the transition data for all
the generators using only powers of `d` as denominators.

"""
function _replace_pre_saturated_ideal(
    I::MPolyLocalizedIdeal{LRT}, Jsat::MPolyIdeal, d::MPolyRingElem
  ) where {LRT<:MPolyLocRing{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement}}
  L = base_ring(I)
  R = base_ring(L)
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
    #A = zero_matrix(L, ngens(Jsat), ngens(I))

    #for i in 1:length(cache)
    #  (g, a, dttk) = cache[i]
    #  A[i, :] = L(one(dttk), dttk, check=false)*change_base_ring(L, a)*pre_saturation_data(I)
    #end
    cols = Vector()
    for i in 1:length(cache)
      (g, a, dttk) = cache[i]
      push!(cols, pre_saturation_data(I) * transpose(L(one(dttk), dttk, check=false)*change_base_ring(L, a)))
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

@doc raw"""
    saturated_ideal(I::MPolyLocalizedIdeal)
    saturated_ideal(I::MPolyQuoLocalizedIdeal)
    saturated_ideal(I::MPolyIdeal)
    saturated_ideal(I::MPolyQuoIdeal)

Return the largest ideal in $R$ mapping to $I$ under the canonical map $R \longrightarrow S$ for an ideal $I \in S$, where $S$ may be any of the
types 'MPolyLocRing', 'MPolyQuoLocRing', 'MPolyRing' and 'MPolyQuoRing' and $R$ is the underlying ring of type 'MPolyRing'.
Note that the last two variants are only provided to allow a coherent usage.  

# Examples
```jldoctest
julia> R, (x, y) = QQ[:x, :y]
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

julia> I = ideal(R, [x, y^2+1])
Ideal generated by
  x
  y^2 + 1

julia> U = complement_of_prime_ideal(I)
Complement
  of prime ideal (x, y^2 + 1)
  in multivariate polynomial ring in 2 variables over QQ

julia> L, _ = localization(R, U)
(Localization of multivariate polynomial ring in 2 variables over QQ at complement of prime ideal (x, y^2 + 1), Hom: R -> localized ring)

julia> J = ideal(L,[y*(x^2+(y^2+1)^2)])
Ideal generated by
  x^2*y + y^5 + 2*y^3 + y

julia> saturated_ideal(J)
Ideal generated by
  x^2 + y^4 + 2*y^2 + 1

julia> JJ = ideal(R,[y*(x^2+(y^2+1)^2)])
Ideal generated by
  x^2*y + y^5 + 2*y^3 + y

julia> saturated_ideal(JJ)
Ideal generated by
  x^2*y + y^5 + 2*y^3 + y

```
"""
function saturated_ideal(
    I::MPolyLocalizedIdeal{LRT};
    with_generator_transition::Bool=false
  ) where {LRT<:MPolyLocRing{<:Any, <:Any, <:Any, <:Any, <:MPolyProductOfMultSets}}
  if !is_saturated(I)
    W = base_ring(I)
    R = base_ring(W)
    J = ideal(R, numerator.(gens(I)))
    for U in sets(inverted_set(W))
      L, _ = localization(U)
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
    LocRingType<:MPolyLocRing{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement}
  }
  L = base_ring(I)
  parent(a) == L || return L(a) in I
  b = numerator(a)
  b in pre_saturated_ideal(I) && return true
  return b in saturated_ideal(I)
end


### Conversion of ideals in the original ring to localized ideals
function (W::MPolyLocRing{BRT, BRET, RT, RET, MST})(I::MPolyIdeal{RET}) where {BRT, BRET, RT, RET, MST}
  result = MPolyLocalizedIdeal(W, W.(gens(I)))
  # Make sure we keep previously computed groebner and standard bases.
  result.pre_saturated_ideal = I
  r = ngens(I)
  A = zero_matrix(SMat, W, 0, r)
  for i in 1:r
    push!(A, sparse_row(W, [(i, one(W))]))
  end
  result.pre_saturation_data = A
  return result
end

### required constructors 
function ideal(
    W::MPolyLocRing, f
  )
  return MPolyLocalizedIdeal(W, [W(f)])
end

@doc raw"""    
    ideal(Rloc::MPolyLocRing, V::Vector)

Given a localization `Rloc` of a multivariate polynomial ring, and given a vector `V` of 
elements of `Rloc`, create the ideal of `Rloc` which is generated by the entries of `V`.
"""
function ideal(
    W::MPolyLocRing,
    gens::Vector
  )
  length(gens) == 0 && return MPolyLocalizedIdeal(W, elem_type(W)[])
  return MPolyLocalizedIdeal(W, W.(gens))
end

function ideal(
    W::MPolyLocRing,
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

function ==(I::IdealType, J::IdealType) where {IdealType<:MPolyLocalizedIdeal}
  I === J && return true
  (issubset(I, J) && issubset(J, I))
end

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

@attr Any function shifted_ideal(
    I::MPolyLocalizedIdeal{LRT, LRET}
  ) where {LRT<:MPolyLocRing{<:Any, <:Any, <:Any, <:Any, <:MPolyComplementOfKPointIdeal}, LRET}
  L = base_ring(I)
  R = base_ring(L)
  shift, tfihs = base_ring_shifts(L)
  return ideal(R, shift.(gens(pre_saturated_ideal(I))))
end

function ideal_membership(
    a::RingElem, 
    I::MPolyLocalizedIdeal{LocRingType}
  ) where {
    LocRingType<:MPolyLocRing{<:Any, <:Any, <:Any, <:Any, <:MPolyComplementOfKPointIdeal}
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
  ) where {LRT<:MPolyLocRing{<:Any, <:Any, <:Any, <:Any, <:MPolyComplementOfKPointIdeal}}
  L = base_ring(I)
  L == parent(a) || return coordinates(L(a), I)
  R = base_ring(L)
  J = shifted_ideal(I)
  shift, tfihs = base_ring_shifts(L)
  p = shift(numerator(a))
  o = negdegrevlex(gens(R))
  x, u = Oscar.lift(p, J, o)
  T = pre_saturation_data(I)
  return transpose(T * transpose(L(one(base_ring(L)), tfihs(u)*denominator(a), check=false)*change_base_ring(L, map_entries(tfihs,x))))
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
  ) where {LRT<:MPolyLocRing{<:Any, <:Any, <:Any, <:Any, <:MPolyProductOfMultSets}}
  L = base_ring(I)
  parent(a) == L || return coordinates(L(a), I)

  R = base_ring(L)
  U = sets(inverted_set(base_ring(I)))

  if length(U) == 1 
    if !has_attribute(I, :popped_ideal)
      W, _ = localization(R, U[1])
      popped_ideal = W(pre_saturated_ideal(I))
      set_attribute!(I, :popped_ideal, popped_ideal)
    end
    popped_ideal = get_attribute(I, :popped_ideal)
    return transpose(pre_saturation_data(I) * transpose(L(one(R), denominator(a), check=false)*map_entries(x->L(x, check=false), coordinates(numerator(a), popped_ideal))))
    #return L(one(R), denominator(a), check=false)*map_entries(x->L(x, check=false), coordinates(numerator(a), popped_ideal))*pre_saturation_data(I)
  end

  if numerator(a) in pre_saturated_ideal(I) 
    return transpose(T * transpose(L(one(R), denominator(a), check=false)*map_entries(L, coordinates(numerator(a), pre_saturated_ideal(I)))))
  end

  i = findfirst(x->(x isa MPolyPowersOfElement), U)
  if !isnothing(i)
    if !has_attribute(I, :popped_ideal)
      S = popat!(U, i)
      W, _ = localization(base_ring(L), S)
      popped_ideal = ideal(W, pre_saturated_ideal(I))
      saturated_ideal(popped_ideal, with_generator_transition=true)
      set_attribute!(I, :popped_ideal, popped_ideal)
      Wnext, _ = localization(R, MPolyProductOfMultSets(R, U))
      next_ideal = Wnext(pre_saturated_ideal(popped_ideal))
      set_attribute!(I, :next_ideal, next_ideal)
    end
    popped_ideal = get_attribute(I, :popped_ideal)
    next_ideal = get_attribute(I, :next_ideal)
    y = coordinates(numerator(a), next_ideal)
    x = transpose(map_entries(x->L(x, check=false), pre_saturation_data(popped_ideal)) *
                      transpose(map_entries(x->L(x, check=false), y)))
    #x = map_entries(x->L(x, check=false), y)*map_entries(x->L(x, check=false), pre_saturation_data(popped_ideal))
    return L(one(R), denominator(a), check=false)*x
  else
    if !has_attribute(I, :popped_ideal)
      S = pop!(U)
      W, _ = localization(base_ring(L), S)
      popped_ideal = ideal(W, pre_saturated_ideal(I))
      saturated_ideal(popped_ideal, with_generator_transition=true)
      set_attribute!(I, :popped_ideal, popped_ideal)
      Wnext, _ = localization(R, MPolyProductOfMultSets(R, U))
      next_ideal = Wnext(pre_saturated_ideal(popped_ideal))
      set_attribute!(I, :next_ideal, next_ideal)
    end
    popped_ideal = get_attribute(I, :popped_ideal)
    next_ideal = get_attribute(I, :next_ideal)
    y = coordinates(numerator(a), next_ideal)
    x = transpose(map_entries(x->L(x, check=false), pre_saturation_data(popped_ideal)) *
                      transpose(map_entries(x->L(x, check=false), y)))
    #x = map_entries(x->L(x, check=false), y)*map_entries(x->L(x, check=false), pre_saturation_data(popped_ideal))
    return L(one(R), denominator(a), check=false)*x
  end
end

function ideal_membership(
    a::RingElem,
    I::MPolyLocalizedIdeal{LRT} 
  ) where {LRT<:MPolyLocRing{<:Any, <:Any, <:Any, <:Any, <:MPolyProductOfMultSets}}
  L = base_ring(I)
  parent(a) == L || return L(a) in I

  R = base_ring(L)
  U = sets(inverted_set(base_ring(I)))

  if length(U) == 1 
    if !has_attribute(I, :popped_ideal)
      W, _ = localization(R, U[1])
      popped_ideal = W(pre_saturated_ideal(I))
      set_attribute!(I, :popped_ideal, popped_ideal)
    end
    popped_ideal = get_attribute(I, :popped_ideal)
    return numerator(a) in popped_ideal
  end

  numerator(a) in pre_saturated_ideal(I) && return true

  i = findfirst(x->(x isa MPolyPowersOfElement), U)
  if !isnothing(i)
    if !has_attribute(I, :popped_ideal)
      S = popat!(U, i)
      W, _ = localization(base_ring(L), S)
      popped_ideal = ideal(W, pre_saturated_ideal(I))
      saturated_ideal(popped_ideal, with_generator_transition=true)
      set_attribute!(I, :popped_ideal, popped_ideal)
      Wnext, _ = localization(R, MPolyProductOfMultSets(R, U))
      next_ideal = Wnext(pre_saturated_ideal(popped_ideal))
      set_attribute!(I, :next_ideal, next_ideal)
    end
    next_ideal = get_attribute(I, :next_ideal)
    return numerator(a) in next_ideal
  else
    if !has_attribute(I, :popped_ideal)
      S = pop!(U)
      W, _ = localization(base_ring(L), S)
      popped_ideal = ideal(W, pre_saturated_ideal(I))
      saturated_ideal(popped_ideal, with_generator_transition=true)
      set_attribute!(I, :popped_ideal, popped_ideal)
      Wnext, _ = localization(R, MPolyProductOfMultSets(R, U))
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
  ) where {LRT<:MPolyLocRing{<:Any, <:Any, <:Any, <:Any, <:MPolyLeadingMonOne}}
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
  ) where {LRT<:MPolyLocRing{<:Any, <:Any, <:Any, <:Any, <:MPolyLeadingMonOne}}
  L = base_ring(I)
  L == parent(a) || return coordinates(L(a), I)
  J = pre_saturated_ideal(I)
  p = numerator(a)
  o = ordering(inverted_set(parent(a)))
  x, u = Oscar.lift(p, J, o)
  T = pre_saturation_data(I)
  return transpose(T * transpose(L(one(base_ring(L)), u*denominator(a), check=false)*
                                    change_base_ring(L, x)))
end

@doc raw"""
    bring_to_common_denominator(f::Vector{T}) where {T<:MPolyLocRingElem}

Given a vector of fractions [a₁//b₁,…,aₙ//bₙ] return a pair 
(d, λ) consisting of a common denominator d and a vector 
λ = [λ₁,…,λₙ] such that aᵢ//bᵢ = λᵢ⋅aᵢ//d
"""
function bring_to_common_denominator(f::Vector{T}) where {T<:MPolyLocRingElem}
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

write_as_linear_combination(f::MPolyLocRingElem, g::Vector) = write_as_linear_combination(f, parent(f).(g))

# return the localized ring as a quotient of a polynomial ring using Rabinowitsch's trick.
@doc raw"""
    as_affine_algebra(
      L::MPolyLocRing{BRT, BRET, RT, RET, 
      MPolyPowersOfElement{BRT, BRET, RT, RET}}; 
      inverse_name::VarName=:0
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
    L::MPolyLocRing{BRT, BRET, RT, RET, 
               MPolyPowersOfElement{BRT, BRET, RT, RET}}; 
    inverse_name::VarName=:_0
  ) where {BRT, BRET, RT, RET}
  R = base_ring(L)
  A, phi, t = _add_variables_first(R, [Symbol(inverse_name)])
  theta = t[1]
  f = prod(denominators(inverted_set(L)))
  I = ideal(A, [one(A)-theta*phi(f)])
  return A, I, f, phi, theta
end


########################################################################
# Homomorphisms of localized polynomial rings                          #
########################################################################

### Mapping of elements
function (f::MPolyLocalizedRingHom)(a::AbsLocalizedRingElem)
  parent(a) === domain(f) || return f(domain(f)(a))
  if total_degree(denominator(a)) > 10
    res = restricted_map(f)
    img_num = res(numerator(a))
    den = denominator(a)
    img_den = one(img_num)
    fac_den = factor(den)
    for (a, k) in fac_den
      img_den = img_den * inv(res(a))^k
    end
    img_den = img_den * inv(res(unit(fac_den)))
    return img_num * img_den
  end
  return codomain(f)(restricted_map(f)(numerator(a)))*inv(codomain(f)(restricted_map(f)(denominator(a))))
end

### additional constructors
function MPolyLocalizedRingHom(
      R::MPolyRing,
      S::Ring,
      a::Vector{RingElemType}
  ) where {RingElemType<:RingElem}
  res = hom(R, S, a, check=false)
  W, _ = localization(units_of(R))
  return MPolyLocalizedRingHom(W, S, res)
end

function MPolyLocalizedRingHom(
      W::MPolyLocRing,
      S::Ring,
      a::Vector{RingElemType};
      check::Bool=true
  ) where {RingElemType<:RingElem}
  R = base_ring(W)
  res = hom(R, S, a, check=false)
  return MPolyLocalizedRingHom(W, S, res, check=check)
end

### printing
"""
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> I = ideal(R, [y-x^2, z-x^3]);

julia> RQ, p = quo(R, I);

julia> UR = complement_of_point_ideal(R, [0, 0, 0]);

julia> RQL, _ = localization(RQ, UR);

julia> T, (t,) =  polynomial_ring(QQ, [:t]);

julia> UT = complement_of_point_ideal(T, [0]);

julia> TL, _ =  localization(T, UT);

julia> PSI = hom(TL, RQL, RQL.([x]))
Ring homomorphism
  from localization of multivariate polynomial ring in 1 variable over QQ at complement of maximal ideal of point (0)
  to localization of RQ at complement of maximal ideal
defined by
  t -> x

julia> (PSI, )
(hom: Localized ring -> Localized quotient of multivariate polynomial ring,)

```
"""
function Base.show(io::IO, ::MIME"text/plain", phi::MPolyLocalizedRingHom)
  io = pretty(io)
  println(terse(io), phi)
  print(io, Indent())
  println(io, "from ", Lowercase(), domain(phi))
  println(io, "to ", Lowercase(), codomain(phi))
  println(io, Dedent(), "defined by", Indent())
  R = base_ring(domain(phi))
  psi = restricted_map(phi)
  to = is_unicode_allowed() ? " ↦ " : " -> "
  for i in 1:ngens(R)-1
    println(io, R[i], to, psi(R[i]))
  end
  n = ngens(R)
  print(io, R[n], to, psi(R[n]))
  print(io, Dedent())
end

function Base.show(io::IO, phi::MPolyLocalizedRingHom)
  if is_terse(io)
    print(io, "Ring homomorphism")
  else
    R = base_ring(domain(phi))
    psi = restricted_map(phi)
    io = pretty(io)
    io = terse(io)
    print(io, "hom: ", domain(phi))
    if is_unicode_allowed()
      print(io, " → ")
    else
      print(io, " -> ")
    end
    print(io, codomain(phi))
  end
end


@doc raw"""
    hom(Rloc::MPolyLocRing, S::Ring, phi::Map)

Given a localized ring `Rloc`of type `MPolyLocRing`, say `Rloc` is the localization 
of a multivariate polynomial ring `R` at the multiplicatively closed subset `U` of `R`, and 
given a homomorphism `phi` from `R` to `S` sending elements of `U` to units in `S`, return 
the homomorphism from `Rloc` to `S` whose composition with the localization map is `phi`.

    hom(Rloc::MPolyLocRing, S::Ring, V::Vector)

Given a localized ring `Rloc` as above, and given a vector `V` of `nvars` elements of `S`, let `phi` 
be the homomorphism from `R` to `S` which is determined by the entries of `V` as the images of
the generators of `R`, and proceed as above.

    hom(RQL::MPolyQuoLocRing, S::Ring, phi::Map)

Given a localized ring `RQL`of type `MPolyQuoLocRing`, say `RQL` is the localization 
of a quotient ring `RQ` of a multivariate polynomial ring `R` at the multiplicatively closed subset `U` of `R`, and 
given a homomorphism `phi` from `R` to `S` sending elements of `U` to units in `S` and elements of the modulus
of `RQ` to zero, return the homomorphism from `Rloc` to `S` whose composition with the localization map `RQ -> RQL`
and the projection map `R -> RQ` is `phi`.

    hom(RQL::MPolyQuoLocRing, S::Ring, V::Vector)

Given a localized ring `RQL`as above, and given a vector `V` of `nvars` elements of `S`, let `phi` 
be the homomorphism from `R` to `S` which is determined by the entries of `V` as the images of
the generators of `R`, and proceed as above.

!!! warning
    Except from the case where the type of `U` is `<: MPolyPowersOfElement`, the condition on `phi` requiring that elements of `U` are send to units in `S` is not checked by the `hom` constructor.
 
# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> I = ideal(R, [y-x^2, z-x^3]);

julia> RQ, p = quo(R, I);

julia> UR = complement_of_point_ideal(R, [0, 0, 0]);

julia> RQL, _ = localization(RQ, UR);

julia> T, (t,) =  polynomial_ring(QQ, [:t]);

julia> UT = complement_of_point_ideal(T, [0]);

julia> TL, _ =  localization(T, UT);

julia> PHI = hom(RQL, TL, TL.([t, t^2, t^3]))
Ring homomorphism
  from localization of RQ at complement of maximal ideal
  to localization of multivariate polynomial ring in 1 variable over QQ at complement of maximal ideal of point (0)
defined by
  x -> t
  y -> t^2
  z -> t^3

julia> PSI = hom(TL, RQL, RQL.([x]))
Ring homomorphism
  from localization of multivariate polynomial ring in 1 variable over QQ at complement of maximal ideal of point (0)
  to localization of RQ at complement of maximal ideal
defined by
  t -> x

julia> PHI(RQL(z))
t^3

julia> PSI(TL(t))
x
```
"""
hom(W::MPolyLocRing, S::Ring, res::Map; check::Bool=true) = MPolyLocalizedRingHom(W, S, res, check=check)
hom(W::MPolyLocRing, S::Ring, a::Vector{RET}; check::Bool=true) where {RET<:RingElem} = MPolyLocalizedRingHom(W, S, a, check=check)

### required getter functions
domain(PHI::MPolyLocalizedRingHom) = PHI.W
codomain(PHI::MPolyLocalizedRingHom) = PHI.S

@doc raw"""
    restricted_map(PHI::MPolyLocalizedRingHom)

    restricted_map(PHI::MPolyQuoLocalizedRingHom)

Given a ring homomorphism `PHI` starting from a localized multivariate polynomial ring
(a localized quotient of a multivariate polynomial ring), return the composition of `PHI` 
with the localization map (with the composition of the localization map and the projection map).

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> I = ideal(R, [y-x^2, z-x^3]);

julia> RQ, p = quo(R, I);

julia> UR = complement_of_point_ideal(R, [0, 0, 0]);

julia> RQL, _ = localization(RQ, UR);

julia> T, (t,) =  polynomial_ring(QQ, [:t]);

julia> UT = complement_of_point_ideal(T, [0]);

julia> TL, _ =  localization(T, UT);

julia> PHI = hom(RQL, TL, TL.([t, t^2, t^3]));

julia> PSI = hom(TL, RQL, RQL.([x]));

julia> phi = restricted_map(PHI)
Ring homomorphism
  from multivariate polynomial ring in 3 variables over QQ
  to localization of multivariate polynomial ring in 1 variable over QQ at complement of maximal ideal of point (0)
defined by
  x -> t
  y -> t^2
  z -> t^3

julia> psi = restricted_map(PSI)
Ring homomorphism
  from multivariate polynomial ring in 1 variable over QQ
  to localization of RQ at complement of maximal ideal
defined by
  t -> x
```
"""
restricted_map(PHI::MPolyLocalizedRingHom) = PHI.res

### implementing the Oscar map interface
function identity_map(W::T) where {T<:MPolyLocRing} 
  MPolyLocalizedRingHom(W, W, identity_map(base_ring(W)), check=false)
end

function compose(
    f::MPolyLocalizedRingHom, 
    g::Map{<:Ring, <:Ring}
  )
  codomain(f) === domain(g) || error("maps are not compatible")
  R = base_ring(domain(f))
  res = restricted_map(f)
  if codomain(res) === domain(g)
    return MPolyLocalizedRingHom(domain(f), codomain(g), compose(res, g), check=false)
  else
    new_res = MapFromFunc(domain(res), codomain(g), x->g(domain(g)(res(x))))
    return MPolyLocalizedRingHom(domain(f), codomain(g), new_res, check=false)
  end
end

(f::MPolyLocalizedRingHom)(I::Ideal) = ideal(codomain(f), f.(domain(f).(gens(I))))

function ==(f::MPolyLocalizedRingHom, g::MPolyLocalizedRingHom) 
  f === g && return true
  domain(f) === domain(g) || return false
  codomain(f) === codomain(g) || return false
  for x in gens(base_ring(domain(f)))
    f(x) == g(x) || return false
  end
  return true
end

function Base.hash(f::MPolyLocalizedRingHom, h::UInt)
  b = 0xe29cdc1ecb68b7a0 % UInt
  h = hash(domain(f), h)
  h = hash(codomain(f), h)
  for x in gens(base_ring(domain(f)))
    h = hash(f(x), h)
  end
  return xor(h, b)
end

function divides(a::MPolyLocRingElem, b::MPolyLocRingElem)
  W = parent(a)
  W == parent(b) || error("elements do not belong to the same ring")
  F = FreeMod(W, 1)
  A = matrix(W, 1, 1, [b])
  M, _ = sub(F, A)
  represents_element(a*F[1], M) || return (false, zero(W))
  x = coordinates(a*F[1], M)
  return true, W(x[1])
end

# This had to be moved after the definition of the elements.
function localization(R::MPolyRing, f::MPolyRingElem)
  U = MPolyPowersOfElement(R, [f])
  L = MPolyLocRing(R, U)
  function func(a::MPolyRingElem)
    parent(a) == R || error("element does not belong to the correct ring")
    return L(a)
  end
  function func_inv(a::MPolyLocRingElem{<:Any, <:Any, <:Any, <:Any, 
                                              <:MPolyPowersOfElement}
    )
    L == parent(a) || error("element does not belong to the correct ring")
    isone(denominator(a)) && return numerator(a)
    return divexact(numerator(a), denominator(a))
  end
  return L, MapFromFunc(R, L, func, func_inv)
end

#############################################################################
# Further functionality for elements of localized rings
#############################################################################

function derivative(f::MPolyLocRingElem, i::Int)
  num = derivative(numerator(f), i)*denominator(f) - derivative(denominator(f), i)*numerator(f)
  den = denominator(f)^2
  g = gcd(num, den)
  return parent(f)(divexact(num, g), divexact(den, g), check=false)
end

function jacobian_matrix(f::MPolyLocRingElem)
  L = parent(f)
  n = nvars(base_ring(L))
  return matrix(L, n, 1, [derivative(f, i) for i=1:n])
end

function jacobian_matrix(g::Vector{<:MPolyLocRingElem})
  R = parent(g[1])
  n = nvars(base_ring(R))
  @assert all(x->parent(x) == R, g)
  return matrix(R, n, length(g), [derivative(x, i) for i=1:n for x = g])
end

function _is_integral_domain(R::MPolyLocRing)
  return true
end

### Method compatible with `localized_ring`: This is the trivial lift to itself. 
function lift(a::MPolyLocRingElem)
  return a
end

### Allow computation of kernels of maps to MPolyLocalizedRings
function kernel(f::MPolyAnyMap{<:MPolyRing, <:MPolyLocRing})
  P = domain(f)
  L = codomain(f)
  I = ideal(L, zero(L))
  R = base_ring(L)
  J = saturated_ideal(I)
  d = [lifted_denominator(g) for g in f.(gens(domain(f)))]
  U = simplify!(MPolyPowersOfElement(R, d))
  W = MPolyQuoLocRing(R, ideal(R, zero(R)), U)
  id =  _as_affine_algebra(W)
  A = codomain(id)
  h = hom(P, A, id.(f.(gens(P))))
  # Even though `id` seems to be type stable, the type is not inferred for empty lists for some reason.
  # Thus, we need to to that manually. 
  gg = Vector{elem_type(A)}(id.(W.(gens(J))))
  return preimage(h, ideal(A, gg))
end

### For introducing in Function to docu
########################################

@doc raw"""

    in(f::MPolyRingElem, U::AbsMPolyMultSet)   

Return `true` if `f` is contained in `U`, `false` otherwise.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> S = complement_of_point_ideal(R, [0, 0 ,0])
Complement
  of maximal ideal corresponding to rational point with coordinates (0, 0, 0)
  in multivariate polynomial ring in 3 variables over QQ

julia> y in S
false

julia> P = ideal(R, [x])
Ideal generated by
  x

julia> T = complement_of_prime_ideal(P)
Complement
  of prime ideal (x)
  in multivariate polynomial ring in 3 variables over QQ

julia> y in T
true

julia> U = powers_of_element(x)
Multiplicative subset
  of multivariate polynomial ring in 3 variables over QQ
  given by the products of [x]

julia> x^3 in U
true

julia> (1+y)*x^2 in product(S, U)
true
```
""" Base.in(f::MPolyRingElem, U::AbsMPolyMultSet) 

### Some auxiliary functions

@attr MPolyLocalizedIdeal function radical(I::MPolyLocalizedIdeal)
  get_attribute(I, :is_prime, false) && return I
  # heuristics have shown that the radical of the saturated ideal 
  # is often quicker 
  #J = pre_saturated_ideal(I)
  J = saturated_ideal(I)
  R = base_ring(I)
  return ideal(base_ring(I), [g for g in R.(gens(radical(J))) if !is_zero(g)])
end

@attr Union{Int, NegInf} function dim(I::MPolyLocalizedIdeal{RT}) where {RT<:MPolyLocRing{<:Ring, <:RingElem, <:MPolyRing, <:MPolyRingElem, <:MPolyPowersOfElement}}
  return dim(saturated_ideal(I))
end

@attr Union{Int, NegInf} function dim(I::MPolyLocalizedIdeal{RT}) where {RT<:MPolyLocRing{<:Ring, <:RingElem, <:MPolyRing, <:MPolyRingElem, <:MPolyComplementOfPrimeIdeal}}
  return dim(saturated_ideal(I)) - dim(prime_ideal(inverted_set(base_ring(I))))
end

@attr Union{Int, NegInf} function dim(I::MPolyLocalizedIdeal{RT}) where {RT<:MPolyLocRing{<:Field, <:FieldElem, <:MPolyRing, <:MPolyRingElem, <:MPolyComplementOfKPointIdeal}}
  J = shifted_ideal(I)
  # TODO: Is there a more conceptual way to do this???
  oo = negdegrevlex(gens(base_ring(J)))
  return dim(leading_ideal(J, ordering=oo))
end


################################################################################
#
# Minimal generating set
#
################################################################################

@doc raw"""
    minimal_generating_set(I::MPolyLocalizedIdeal)

Given an ideal `I` in the localization of a multivariate polynomial ring
over a field at a point, return an array containing a minimal set of
generators of `I`. If `I` is the zero ideal an empty list is returned.

# Examples
```jldoctest
julia> R, (x, y) = QQ[:x, :y];

julia> L, phi = localization(R, complement_of_point_ideal(R, [1, 2]));

julia> I = ideal(L, [x-1, y-2])^2
Ideal generated by
  x^2 - 2*x + 1
  x*y - 2*x - y + 2
  x*y - 2*x - y + 2
  y^2 - 4*y + 4

julia> minimal_generating_set(I)
3-element Vector{Oscar.MPolyLocRingElem{QQField, QQFieldElem, QQMPolyRing, QQMPolyRingElem, Oscar.MPolyComplementOfKPointIdeal{QQField, QQFieldElem, QQMPolyRing, QQMPolyRingElem}}}:
 y^2 - 4*y + 4
 x*y - 2*x - y + 2
 x^2 - 2*x + 1
```
"""
@attr Vector{<:MPolyLocRingElem} function minimal_generating_set(
    I::MPolyLocalizedIdeal{<:MPolyLocRing{<:Field, <:FieldElem,
                                          <:MPolyRing, <:MPolyRingElem,
                                          <:MPolyComplementOfKPointIdeal}
                          }
  )
  L = base_ring(I)
  R = base_ring(L)

  @req coefficient_ring(R) isa AbstractAlgebra.Field "The coefficient ring must be a field"

  shift, back_shift = base_ring_shifts(L)
  I_shift = shifted_ideal(I)
  o = negdegrevlex(R) # the default local ordering

  if !isempty(I_shift.gb)

    # use available gb for shifted ideal is available
    G = first(values(I_shift.gb))
    is_local(G.ord) || error("inconsistent data: local ring, but global ordering for I_shift")
    G.gens.S.isGB = true
    _, I_shift_min = Singular.mstd(singular_generators(G, G.ord))
  else

    # no gb for shifted ideal available
    gensI_shift = singular_generators(I_shift, o)
    I_shift_gb, I_shift_min = Singular.mstd(gensI_shift)

    # store gb of I_shift
    computed_gb = IdealGens(base_ring(I_shift),I_shift_gb, false)
    computed_gb.isGB = true
    I_shift.gb[o] = computed_gb
  end

  # return result after shifting back and moving to the correct ring
  computed_min = IdealGens(base_ring(I_shift), I_shift_min, false)
  I_min = L.(back_shift.(R.(gens(computed_min))))
  return filter(!iszero, I_min)

end

@doc raw"""
    small_generating_set(I::MPolyLocalizedIdeal)

Given an ideal `I` in a localization of a multivariate polynomial ring
over a field, return an array containing a set of generators of `I`,
which is usually smaller than the original one. 

If `I` is the zero ideal an empty list is returned.

If the localization is at a point, a minimal set of generators is returned.

# Note: not available for localizations at prime ideals.

# Examples
```jldoctest
julia> R, (x, y) = QQ[:x, :y];

julia> L, phi = localization(R, powers_of_element(x));

julia> I = ideal(L, [x*(x-1), y])^2
Ideal generated by
  x^4 - 2*x^3 + x^2
  x^2*y - x*y
  x^2*y - x*y
  y^2

julia> small_generating_set(I)
3-element Vector{Oscar.MPolyLocRingElem{QQField, QQFieldElem, QQMPolyRing, QQMPolyRingElem, Oscar.MPolyPowersOfElement{QQField, QQFieldElem, QQMPolyRing, QQMPolyRingElem}}}:
 y^2
 x^2*y - x*y
 x^4 - 2*x^3 + x^2

```
"""
small_generating_set(I::MPolyLocalizedIdeal{<:MPolyLocRing{<:Field, <:FieldElem,
                        <:MPolyRing, <:MPolyRingElem,
                        <:MPolyComplementOfKPointIdeal}})  = minimal_generating_set(I)

@attr Vector{elem_type(base_ring(I))} function small_generating_set(
    I::MPolyLocalizedIdeal{<:MPolyLocRing{<:Field, <:FieldElem,
                                          <:MPolyRing, <:MPolyRingElem,
                                          <:MPolyPowersOfElement}};
    algorithm::Symbol=:simple
  )
  L = base_ring(I)
  R = base_ring(L)
  #I_min = L.(small_generating_set(saturated_ideal(I)))
  I_min = L.(small_generating_set(ideal(R, numerator.(gens(I)))))
  filter(!iszero, I_min)
end

dim(R::MPolyLocRing{<:Field, <:FieldElem, <:MPolyRing, <:MPolyRingElem, <:MPolyComplementOfPrimeIdeal}) = nvars(base_ring(R)) - dim(prime_ideal(inverted_set(R)))

########################################################################
# Localizations of graded rings                                        #
########################################################################
function is_graded(L::MPolyLocRing{<:Ring, <:RingElem, <:MPolyDecRing})
  return true
end

function grading_group(L::MPolyLocRing{<:Ring, <:RingElem, <:MPolyDecRing})
  return grading_group(base_ring(L))
end

function degree(a::MPolyLocRingElem{<:Ring, <:RingElem, <:MPolyDecRing}; check::Bool=true)
  return degree(numerator(a); check) - degree(denominator(a); check)
end

function degree(::Type{Int}, a::MPolyLocRingElem{<:Ring, <:RingElem, <:MPolyDecRing}; check::Bool=true)
  return degree(Int, numerator(a); check) - degree(Int, denominator(a); check)
end

function degree(::Type{Vector{Int}}, a::MPolyLocRingElem{<:Ring, <:RingElem, <:MPolyDecRing}; check::Bool=true)
  return degree(Vector{Int}, numerator(a); check) - degree(Vector{Int}, denominator(a); check)
end

function _degree_fast(a::MPolyLocRingElem{<:Ring, <:RingElem, <:MPolyDecRing})
  return _degree_fast(numerator(a)) - _degree_fast(denominator(a))
end

function is_homogeneous(a::MPolyLocRingElem{<:Ring, <:RingElem, <:MPolyDecRing})
  return is_homogeneous(numerator(a)) && is_homogeneous(denominator(a))
end

coefficient_ring(L::MPolyLocRing) = coefficient_ring(base_ring(L))

########################################################################
# Ported from old implementation in src/Rings/mpoly-local.jl
########################################################################
function Oscar.localization(R::MPolyRing{S}, m::MPolyIdeal) where S
  return MPolyLocRing(R, complement_of_point_ideal(m))
end

complement_of_point_ideal(m::MPolyIdeal) = complement_of_point_ideal(base_ring(m), rational_point_coordinates(m))

