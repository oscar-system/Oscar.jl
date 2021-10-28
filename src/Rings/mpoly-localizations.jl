export MPolyComplementOfPrimeIdeal, MPolyComplementOfKPointIdeal

export ambient_ring, point_coordinates

export original_ring, inverted_set
export reduce_fraction

export fraction, parent

export MPolyLocalizedRing

export MPolyLocalizedIdeal
export gens, base_ring, groebner_bases, default_ordering, dim 

export localize_at, ideal

export LocalizedBiPolyArray
export oscar_gens, oscar_ring, singular_ring, singular_gens, ordering, shift
export groebner_basis, groebner_assure, reduce

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
    S::MPolyComplementOfPrimeIdeal{BaseRingType, BaseRingElemType, RingType, RingElemType}
  ) where {BaseRingType, BaseRingElemType, RingType, RingElemType} = S.R

### required functionality
function Base.in(
    f::RingElemType, 
    S::MPolyComplementOfPrimeIdeal{BaseRingType, BaseRingElemType, RingType, RingElemType}
  ) where {BaseRingType, BaseRingElemType, RingType, RingElemType}
  return !(f in prime_ideal(S))
end

### additional functionality
prime_ideal(
    S::MPolyComplementOfPrimeIdeal{BaseRingType, BaseRingElemType, RingType, RingElemType}
  ) where {BaseRingType, BaseRingElemType, RingType, RingElemType} = S.P


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
ambient_ring(
    S::MPolyComplementOfKPointIdeal{BaseRingType, BaseRingElemType, RingType, RingElemType}
  ) where {BaseRingType, BaseRingElemType, RingType, RingElemType} = S.R

### additional getter functions 
function point_coordinates(
    S::MPolyComplementOfKPointIdeal{BaseRingType, BaseRingElemType, RingType, RingElemType}
  ) where {BaseRingType, BaseRingElemType, RingType, RingElemType}
  return S.a
end

### required functionality
function Base.in(
    f::RingElemType, 
    S::MPolyComplementOfKPointIdeal{BaseRingType, BaseRingElemType, RingType, RingElemType}
  ) where {BaseRingType, BaseRingElemType, RingType, RingElemType}
  parent(f) == ambient_ring(S) || error("the given element does not belong to the same ring as the multiplicative set")
  return !(evaluate(f, point_coordinates(S)) == zero(ambient_ring(S)))
end


########################################################################
# Localizations of polynomial rings over admissible fields             #
########################################################################

@Markdown.doc """
    MPolyLocalizedRing{
        BaseRingType, 
	BaseRingElemType,
        RingType,
        RingElemType
      } <: AbsLocalizedRing{
        RingType,
        MPolyElem{BaseRingType}, 
        MPolyComplementOfPrimeIdeal{BaseRingType, RingType, RingElemType}
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
    MultSetType
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
original_ring(
    W::MPolyLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}
  ) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = W.R

inverted_set(
    W::MPolyLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}
  ) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = W.S

### required extension of the localization function
localize_at(
    S::MPolyComplementOfPrimeIdeal{
            BaseRingType, 
	    BaseRingElemType,
    	    RingType, 
	    RingElemType
        }
    ) where {
	BaseRingType,
	BaseRingElemType,
	RingType, 
	RingElemType
    } = MPolyLocalizedRing(ambient_ring(S), S)

localize_at(
        S::MPolyComplementOfKPointIdeal{
	    BaseRingType,
    	    BaseRingElemType, 
    	    RingType, 
    	    RingElemType
        }
    ) where {
	BaseRingType,
	BaseRingElemType, 
	RingType, 
	RingElemType
    } = MPolyLocalizedRing(ambient_ring(S), S)

### additional constructors
MPolyLocalizedRing(R::RingType, P::MPolyIdeal{RingElemType}) where {RingType, RingElemType} = MPolyLocalizedRing(R, MPolyComplementOfPrimeIdeal(P))


########################################################################
# Elements of local polynomial rings                                   #
########################################################################

@Markdown.doc """
    MPolyLocalizedRingElem{BaseRingType, RingElemType} <: AbsLocalizedRingElem{MPolyRing{BaseRingType}, MPolyElem{BaseRingType}, MPolyComplementOfPrimeIdeal{BaseRingType, RingElemType}} 


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

    base_ring(parent(f)) == original_ring(R_loc) || error(
	"the numerator and denominator of the given fraction do not belong to the original ring before localization"
      )
    denominator(f) in inverted_set(R_loc) || error(
	"the given denominator is not admissible for this localization"
      )
    return new{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}(f, R_loc)
  end
end

### required getter functions 
numerator(
    a::MPolyLocalizedRingElem{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}
  ) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = numerator(a.frac)

denominator(
    a::MPolyLocalizedRingElem{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}
  ) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = denominator(a.frac)

parent(
    a::MPolyLocalizedRingElem{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}
  ) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = a.R_loc

### additional getter functions
fraction(
    a::MPolyLocalizedRingElem{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}
  ) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = a.frac

### required conversions
(W::MPolyLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType})(f::RingElemType) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyLocalizedRingElem((FractionField(original_ring(W)))(f), W)
(W::MPolyLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType})(a::RingElemType, b::RingElemType) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyLocalizedRingElem(a//b, W)

### additional conversions
(W::MPolyLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType})(f::AbstractAlgebra.Generic.Frac{RingElemType}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyLocalizedRingElem(f, W)

(W::MPolyLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType})(a::Oscar.IntegerUnion) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = W(original_ring(W)(a))

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
  numerator(fraction(b)) in inverted_set(parent(b)) || error("the second argument is not a unit in this local ring")
  return (parent(a))(fraction(a) // fraction(b))
end

function ==(a::T, b::T) where {T<:MPolyLocalizedRingElem}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  return fraction(a) == fraction(b)
end

function ^(a::T, i::Oscar.IntegerUnion) where {T<:MPolyLocalizedRingElem}
  return parent(a)(fraction(a)^i)
end

### implementation of Oscar's general ring interface
one(W::MPolyLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = W(one(original_ring(W)))
zero(W::MPolyLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = W(zero(original_ring(W)))

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
    for i in (length(shift)+1:length(gens(original_ring(oscar_ring))))
      append!(shift, zero(oscar_ring))
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
    R = original_ring(oscar_ring)
    lbpa.ordering = Singular.ordering_as_symbol(lbpa.singular_ring)
    k = coefficient_ring(R)
    # fill up the shift vector with zeroes if it is not provided in full length
    for i in (length(shift)+1:nvars(R))
	      append!(shift, zero(k))
    end
    lbpa.shift = shift
    inv_shift_hom = AlgebraHomomorphism(R, R, [gen(R, i) - R(shift[i]) for i in (1:nvars(R))])
    lbpa.oscar_gens = [ oscar_ring(y) for y in (inv_shift_hom.(R.([x for x in gens(singular_gens)])))]
    lbpa.is_groebner_basis=is_groebner_basis
    return lbpa
  end
end

oscar_gens(lbpa::LocalizedBiPolyArray{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST} = lbpa.oscar_gens
oscar_ring(lbpa::LocalizedBiPolyArray{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST} = lbpa.oscar_ring
ordering(lbpa::LocalizedBiPolyArray{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST} = lbpa.ordering
shift(lbpa::LocalizedBiPolyArray{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST} = lbpa.shift

function singular_gens(lbpa::LocalizedBiPolyArray{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST}
  singular_assure(lbpa)
  return lbpa.singular_gens
end

function singular_ring(lbpa::LocalizedBiPolyArray{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST}
  singular_assure(lbpa)
  return lbpa.singular_ring
end

function singular_assure(lbpa::LocalizedBiPolyArray{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST}
  if !isdefined(lbpa, :singular_ring)
    lbpa.singular_ring = Singular.PolynomialRing(
	Oscar.singular_ring(base_ring(original_ring(oscar_ring(lbpa)))), 
        [string(a) for a = Nemo.symbols(original_ring(oscar_ring(lbpa)))], 
        ordering = ordering(lbpa), 
        cached = false
      )[1]
  end
  if !isdefined(lbpa, :singular_gens)
    shift_hom = hom(original_ring(oscar_ring(lbpa)), original_ring(oscar_ring(lbpa)), 
        [gen(original_ring(oscar_ring(lbpa)), i) + lbpa.shift[i] for i in (1:nvars(original_ring(oscar_ring(lbpa))))])
    lbpa.singular_gens = Singular.Ideal(lbpa.singular_ring,
	[lbpa.singular_ring(shift_hom(numerator(x))) for x in oscar_gens(lbpa)])
  end
end

function std(lbpa::LocalizedBiPolyArray{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST}
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
    R = original_ring(W)
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
gens(I::MPolyLocalizedIdeal{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST} = I.gens
base_ring(I::MPolyLocalizedIdeal{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST} = I.W

### additional getter functions
groebner_bases(I::MPolyLocalizedIdeal{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST} = I.groebner_bases

default_ordering(I::MPolyLocalizedIdeal{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST} = I.default_ordering

function dim(I::MPolyLocalizedIdeal{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST}
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


########################################################################
# Groebner and standard bases                                          #
########################################################################

### the catchall implementation; most probably mathematically useless!
function groebner_basis(
    I::MPolyLocalizedIdeal{BRT, BRET, RT, RET, MST}, 
    ordering::Symbol
  ) where {BRT, BRET, RT, RET, MST} 
  D = groebner_bases(I)
  # check whether a standard basis has already been computed for this ordering
  if haskey(D, ordering)
    return D[ordering]
  end
  # if not, set up a LocalizedBiPolyArray
  W = parent(I)
  R = original_ring(W)
  S = inverted_set(W)::MPolyComplementOfKPointIdeal{BRT, BRET, RT, RET}
  a = point_coordinates(S)
  lbpa = LocalizedBiPolyArray(W, gens(I), ordering=ordering, shift=a)
  # compute the standard basis and cache the result
  D[ordering] = std(lbpa)
  return D[ordering]
end

function groebner_basis(
    I::MPolyLocalizedIdeal{BRT, BRET, RT, RET, MST}, 
  ) where {BRT, BRET, RT, RET, MST} 
  return groebner_basis(I, default_ordering(I))
end

function groebner_assure(
    I::MPolyLocalizedIdeal{BRT, BRET, RT, RET, MST}
  ) where {BRT, BRET, RT, RET, MST}
  D = groebner_bases(I)
  if length(D) > 0 
    return
  end
  W = parent(I)
  R = original_ring(W)
  S = inverted_set(W)
  lbpa = LocalizedBiPolyArray(W, gens(I), ordering=default_ordering(I))
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
  R = original_ring(W)
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
  W = parent(I)
  R = original_ring(W)
  S = inverted_set(W)::MPolyComplementOfKPointIdeal{BRT, BRET, RT, RET}
  a = point_coordinates(S)
  # Choose negdegrevlex as a default local ordering
  lbpa = LocalizedBiPolyArray(W, gens(I), ordering=:negdegrevlex, shift=a)
  D[:negdegrevlex] = std(lbpa)
  return
end
  
function reduce(
    f::MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}, 
    lbpa::LocalizedBiPolyArray{BRT, BRET, RT, RET, MST}
  ) where {BRT, BRET, RT, RET, MST}

  W = parent(f)
  W == oscar_ring(lbpa) || error("element does not belong to the Oscar ring of the biPolyArray")
  R = original_ring(W)
  shift_hom = hom(R, R, [gen(R, i) + lbpa.shift[i] for i in (1:nvars(R))])
  singular_n = singular_ring(lbpa)(shift_hom(numerator(f)))
  singular_n = Singular.reduce(singular_n, singular_gens(lbpa))
  inv_shift_hom = hom(R,R, [gen(R, i) - lbpa.shift[i] for i in (1:nvars(R))])
  return W(inv_shift_hom(R(singular_n)), denominator(f))
end



