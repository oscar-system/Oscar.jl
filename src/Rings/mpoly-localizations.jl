export MPolyComplementOfPrimeIdeal, MPolyComplementOfKPointIdeal, MPolyPowersOfElement

export ambient_ring, point_coordinates, inverted_set

export reduce_fraction

export fraction, parent

export MPolyLocalizedRing

export MPolyLocalizedIdeal
export gens, base_ring, groebner_bases, default_ordering, dim 

export localize_at, ideal

export LocalizedBiPolyArray
export oscar_gens, oscar_ring, singular_ring, singular_gens, ordering, shift
export groebner_basis, groebner_assure

import AbstractAlgebra.Ring

########################################################################
# General framework for localizations of multivariate polynomial rings #
########################################################################


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
    end
    k = coefficient_ring(R)
    return new{typeof(k), elem_type(k), RingType, RingElemType}(R, a)
  end
end

### required getter functions
ambient_ring(S::MPolyPowersOfElement) = S.R

### additional functionality
denominators(S::MPolyPowersOfElement) = S.a

### required functionality
function Base.in(
    f::RingElemType, 
    S::MPolyPowersOfElement{BaseRingType, BaseRingElemType, RingType, RingElemType}
  ) where {BaseRingType, BaseRingElemType, RingType, RingElemType}
  parent(f) == ambient_ring(S) || return false

  # TODO: This algorithm can most probably be improved!
  R = ambient_ring(S)
  d = prod(denominators(S))
  g = gcd(f, d)
  while !isone(g)
    f = divexact(f, g)
    g = gcd(f, d)
  end
  # return true iff the remaining f is a unit in R
  return divides(one(R), f)[1]
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
localize_at(S::MPolyComplementOfPrimeIdeal) = MPolyLocalizedRing(ambient_ring(S), S)

localize_at(S::MPolyComplementOfKPointIdeal) = MPolyLocalizedRing(ambient_ring(S), S)

localize_at(S::MPolyPowersOfElement) = MPolyLocalizedRing(ambient_ring(S), S)


### additional constructors
MPolyLocalizedRing(R::RingType, P::MPolyIdeal{RingElemType}) where {RingType, RingElemType} = MPolyLocalizedRing(R, MPolyComplementOfPrimeIdeal(P))

localize_at(R::MPolyRing, f::MPolyElem) = localize_at(MPolyPowersOfElement(R, [f]))
localize_at(R::MPolyRing, v::Vector{T}) where {T<:MPolyElem} = localize_at(MPolyPowersOfElement(R, v))

function localize_at(
    W::MPolyLocalizedRing{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}}, 
    f::RET
  ) where {BRT, BRET, RT, RET}
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

function localize_at(
    W::MPolyLocalizedRing{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}}, 
    v::Vector{RET}
  ) where {BRT, BRET, RT, RET}
  V = W
  for f in v
    V = localize_at(V, f)
  end
  return V
end

########################################################################
# Elements of local polynomial rings                                   #
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
(W::MPolyLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType})(a::RingElemType, b::RingElemType) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyLocalizedRingElem(a//b, W)

### additional conversions
(W::MPolyLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType})(f::AbstractAlgebra.Generic.Frac{RingElemType}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyLocalizedRingElem(f, W)

(W::MPolyLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType})(a::Oscar.IntegerUnion) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = W(base_ring(W)(a))

### overwriting the arithmetic using the fractions from AbstractAlgebra
function +(a::T, b::T) where {T<:MPolyLocalizedRingElem}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  return (parent(a))(fraction(a) + fraction(b))
end

function -(a::T, b::T) where {T<:MPolyLocalizedRingElem}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  return (parent(a))(fraction(a) - fraction(b))
end

function *(a::T, b::T) where {T<:MPolyLocalizedRingElem}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  return (parent(a))(fraction(a) * fraction(b))
end

function Base.:(//)(a::Oscar.IntegerUnion, b::T) where {T<:MPolyLocalizedRingElem}
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

function ^(a::MPolyLocalizedRingElem, i::Oscar.IntegerUnion)
  return parent(a)(fraction(a)^i)
end

### implementation of Oscar's general ring interface
one(W::MPolyLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = W(one(base_ring(W)))
zero(W::MPolyLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = W(zero(base_ring(W)))

elem_type(W::MPolyLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyLocalizedRingElem{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}
elem_type(T::Type{MPolyLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyLocalizedRingElem{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}

parent_type(W::MPolyLocalizedRingElem{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}
parent_type(T::Type{MPolyLocalizedRingElem{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}


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
  default_ordering::Symbol
  groebner_bases::Dict{Symbol, LocalizedBiPolyArray{BRT, BRET, RT, RET, MST}}
  dimension::Int
 
  function MPolyLocalizedIdeal(
      W::MPolyLocalizedRing{BRT, BRET, RT, RET, MST}, 
      gens::Vector{MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}};
      default_ordering::Symbol=:degrevlex, 
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
    I.default_ordering = default_ordering
    I.groebner_bases = Dict{Symbol, LocalizedBiPolyArray{BRT, BRET, RT, RET, MST}}()
    return I
  end
end
 
### required getter functions
gens(I::MPolyLocalizedIdeal) = I.gens
base_ring(I::MPolyLocalizedIdeal) = I.W

### additional getter functions
groebner_bases(I::MPolyLocalizedIdeal) = I.groebner_bases

default_ordering(I::MPolyLocalizedIdeal) = I.default_ordering

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
  return reduce(f, lbpa) == base_ring(I)(0)
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
  inv_shift_hom = hom(R,R, [gen(R, i) - lbpa.shift[i] for i in (1:nvars(R))])
  return W(inv_shift_hom(R(singular_n)), denominator(f))
end



