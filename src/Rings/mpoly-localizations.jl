export MPolyComplementOfPrimeIdeal, MPolyComplementOfMaxIdeal

export original_ring, inverted_set
export reduce_fraction

export fraction, parent

export MPolyLocalRing, MPolyMaxLocalRing

export MPolyLocalizedIdeal, MPolyLocalRingIdeal, MPolyMaxLocalIdeal

export localize_at, ideal

import AbstractAlgebra.Ring

########################################################################
# General framework for localizations of multivariate polynomial rings #
########################################################################


########################################################################
# Complements of prime ideals                                          #
########################################################################

@Markdown.doc """
    MPolyComplementOfPrimeIdeal{BaseRingType, RingElemType} <: AbsComplementOfPrimeIdeal{MPolyRing{BaseRingType}, MPolyIdeal{RingElemType}, RingElemType}

The complement of a prime ideal `P ⊂ 𝕜[x₁,…,xₙ]` in a multivariate polynomial ring 
with elements of type `RingElemType` over a base ring `𝕜` of type `BaseRingType`.
"""
mutable struct MPolyComplementOfPrimeIdeal{BaseRingType, RingElemType} <: AbsComplementOfPrimeIdeal{MPolyRing{BaseRingType}, MPolyIdeal{RingElemType}, RingElemType}
  # The parent polynomial ring 𝕜[x₁,…,xₙ]
  R::MPolyRing{BaseRingType}
  # The prime ideal whose complement this is
  P::MPolyIdeal{RingElemType}

  function MPolyComplementOfPrimeIdeal(
      P::MPolyIdeal{RingElemType}; 
      check::Bool=false
    ) where {RingElemType}
    R = base_ring(P)
    check && (isprime(P) || error("the ideal $P is not prime"))
    return new{typeof(coefficient_ring(R)), elem_type(R)}(R, P)
  end
end

### required getter functions
ambient_ring(
    W::MPolyComplementOfPrimeIdeal{BaseRingType, RingElemType}
  ) where {BaseRingType, RingElemType} = S.R

### required functionality
function Base.in(f::RingElemType, S::MPolyComplementOfPrimeIdeal{BaseRingType, RingElemType}) where {BaseRingType, RingElemType}
  return !(f in prime_ideal(S))
end

### additional functionality
prime_ideal(
    S::MPolyComplementOfPrimeIdeal{BaseRingType, RingElemType}
  ) where {BaseRingType, RingElemType} = S.P


########################################################################
# Complements of certain maximal ideals                                #
########################################################################

@Markdown.doc """
    MPolyComplementOfMaxIdeal{BaseRingType, BaseRingElemType, RingElemType} <: AbsComplementOfPrimeIdeal{MPolyRing{BaseRingType}, MPolyIdeal{RingElemType}, RingElemType}

Complement of a maximal ideal ``𝔪 = ⟨x₁-a₁,…,xₙ-aₙ⟩⊂ 𝕜[x₁,…xₙ]`` with ``aᵢ∈ 𝕜``.
"""
mutable struct MPolyComplementOfMaxIdeal{BaseRingType, BaseRingElemType, RingElemType} <: AbsComplementOfPrimeIdeal{MPolyRing{BaseRingType}, MPolyIdeal{RingElemType}, RingElemType}
  # The parent polynomial ring 𝕜[x₁,…,xₙ]
  R::MPolyRing{BaseRingType}
  # The coordinates aᵢ of the point in 𝕜ⁿ corresponding to the maximal ideal
  a::Vector{BaseRingElemType}

  function MPolyComplementOfMaxIdeal(R::MPolyRing{BaseRingType}, a::Vector{BaseRingElemType}) where {BaseRingType, BaseRingElemType}
    length(a) == length(gens(R)) || error("the number of variables in the ring does not coincide with the number of coordinates")
    n = length(a)
    if n > 0 
      base_ring(R) == parent(a[1]) || error("the coordinates are not elements of the base ring")
    else
      elem_type(base_ring(R)) == BaseRingElemType || error("the type of the coordinates does not match the elem_type of the base ring")
    end
    S = new{BaseRingType, BaseRingElemType, elem_type(R)}(R, a)
    return S
  end
end

### required getter functions
ambient_ring(
    S::MPolyComplementOfMaxIdeal{BaseRingType, BaseRingElemType, RingElemType}
  ) where {BaseRingType, BaseRingElemType, RingElemType} = S.R

### required functionality
function Base.in(f::RingElemType, S::MPolyComplementOfMaxIdeal{BaseRingType, BaseRingElemType, RingElemType}) where {BaseRingType, BaseRingElemType, RingElemType}
  return !(evaluate(f, S.a) == zero(S.R))
end

### additional functionality
function point_coordinates(S::MPolyComplementOfMaxIdeal{BaseRingType, BaseRingElemType, RingElemType}) where {BaseRingType, BaseRingElemType, RingElemType}
  return S.a
end


#################################################################
# General localizations of polynomial rings                     #
#################################################################

@Markdown.doc """
    MPolyLocalizedRing{BaseRingType, RingElemType, MultSetType} <: AbsLocalizedRing{MPolyRing{BaseRingType}, MPolyElem{BaseRingType}}

Localizations ``R[S⁻¹]`` of free polynomial rings ``R = 𝕜[x₁,…,xₙ]`` with elements of type 
`RingElemType` over some base ring ``𝕜`` of type `BaseRingType` and with multiplicatively 
closed set ``S`` of type `MultSetType`.
"""
mutable struct MPolyLocalizedRing{BaseRingType, RingElemType, MultSetType} <: AbsLocalizedRing{MPolyRing{BaseRingType}, MPolyElem{BaseRingType}, MultSetType}
  R::MPolyRing{BaseRingType} # The parent ring which is being localized
  S::MultSetType

  function MPolyLocalizedRing(R::MPolyRing{BaseRingType}, 
      S::MultSetType) where {BaseRingType, MultSetType}
    # TODO: Add some sanity checks here?
    R_loc = new{BaseRingType, elem_type(R), MultSetType}(R,S)
    return R_loc
  end
end

### required getter functions
original_ring(W::MPolyLocalizedRing{BaseRingType, RingElemType, MultSetType}) where {BaseRingType, RingElemType, MultSetType} = W.R
inverted_set(W::MPolyLocalizedRing{BaseRingType, RingElemType, MultSetType}) where {BaseRingType, RingElemType, MultSetType} = W.S

### required extension of the localization function
localize_at(S::AbsMultSet{MPolyRing{BaseRingType}, MPolyElem{BaseRingType}}) where {BaseRingType<:Ring} = MPolyLocalizedRing(ambient_ring(S), S)


############################################################################
# Localizations of polynomial rings over admissible fields at prime ideals #
############################################################################
@Markdown.doc """
    MPolyLocalRing{BaseRingType, RingElemType} <: AbsLocalizedRing{MPolyRing{BaseRingType}, MPolyElem{BaseRingType}, MPolyComplementOfPrimeIdeal{BaseRingType, RingElemType}}

The localization of a multivariate polynomial ring R = 𝕜[x₁,…,xₙ] over a 
base field 𝕜 of type `BaseRingType` and with elements of type `RingElemType` 
at a prime ideal P ⊂ R.
"""
mutable struct MPolyLocalRing{
    BaseRingType, 
    RingElemType
  } <: AbsLocalizedRing{
    MPolyRing{BaseRingType}, 
    MPolyElem{BaseRingType}, 
    MPolyComplementOfPrimeIdeal{BaseRingType, RingElemType}
  }
  R::MPolyRing{BaseRingType} # The parent ring which is being localized
  S::MPolyComplementOfPrimeIdeal{BaseRingType, RingElemType} 

  function MPolyLocalRing(R::MPolyRing{BaseRingType}, 
      S::MPolyComplementOfPrimeIdeal{BaseRingType, RingElemType}
    ) where {BaseRingType, RingElemType}
    # TODO: Add some sanity checks here?
    R_loc = new{BaseRingType, RingElemType}(R,S)
    return R_loc
  end
end

### required getter functions 
original_ring(W::MPolyLocalRing{BaseRingType, RingElemType}) where {BaseRingType, RingElemType} = W.R
inverted_set(W::MPolyLocalRing{BaseRingType, RingElemType}) where {BaseRingType, RingElemType} = W.S

### required extension of the localization function
localize_at(S::MPolyComplementOfPrimeIdeal{BaseRingType, RingElemType}) where {BaseRingType, RingElemType} = MPolyLocalRing(ambient_ring(S), S)

### additional constructors
MPolyLocalRing(R::MPolyRing{BaseRingType}, P::MPolyIdeal{RingElemType}) where {BaseRingType, RingElemType} = MPolyLocalRing(R, MPolyComplementOfPrimeIdeal(P))


########################################################################
# Elements of local polynomial rings                                   #
########################################################################

mutable struct MPolyLocalRingElem{BaseRingType, RingElemType} <: AbsLocalizedRingElem{MPolyRing{BaseRingType}, MPolyElem{BaseRingType}, MPolyComplementOfPrimeIdeal{BaseRingType, RingElemType}} 
  frac::AbstractAlgebra.Generic.Frac{RingElemType}
  R_loc::MPolyLocalRing{BaseRingType, RingElemType}

  function MPolyLocalRingElem(f::AbstractAlgebra.Generic.Frac{RingElemType}, R_loc::MPolyLocalRing{BaseRingType, RingElemType}) where {BaseRingType, RingElemType}

    base_ring(parent(f)) == original_ring(R_loc) || error("the numerator and denominator of the given fraction do not belong to the original ring before localization")
    denominator(f) in inverted_set(R_loc) || error("the given denominator is not admissible for this localization")
    return new{BaseRingType, RingElemType}(f, R_loc)
  end
end

### required getter functions 
numerator(a::MPolyLocalRingElem{BaseRingType, RingElemType}) where {BaseRingType, RingElemType} = numerator(a.frac)
denominator(a::MPolyLocalRingElem{BaseRingType, RingElemType}) where {BaseRingType, RingElemType} = denominator(a.frac)
parent(a::MPolyLocalRingElem{BaseRingType, RingElemType}) where {BaseRingType, RingElemType} = a.R_loc

### additional getter functions
fraction(a::MPolyLocalRingElem{BaseRingType, RingElemType}) where {BaseRingType, RingElemType} = a.frac

### required conversions
(W::MPolyLocalRing{BaseRingType, RingElemType})(f::RingElemType) where {BaseRingType, RingElemType} = MPolyLocalRingElem((FractionField(original_ring(W)))(f), W)
(W::MPolyLocalRing{BaseRingType, RingElemType})(a::RingElemType, b::RingElemType) where {BaseRingType, RingElemType} = MPolyLocalRingElem(a//b, W)

### additional conversions
(W::MPolyLocalRing{BaseRingType, RingElemType})(f::AbstractAlgebra.Generic.Frac{RingElemType}) where {BaseRingType, RingElemType} = MPolyLocalRingElem(f, W)

### overwriting the arithmetic using the fractions from AbstractAlgebra
function +(a::MPolyLocalRingElem{BaseRingType, RingElemType}, b::MPolyLocalRingElem{BaseRingType, RingElemType}) where {BaseRingType, RingElemType}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  return (parent(a))(fraction(a) + fraction(b))
end

function -(a::MPolyLocalRingElem{BaseRingType, RingElemType}, b::MPolyLocalRingElem{BaseRingType, RingElemType}) where {BaseRingType, RingElemType}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  return (parent(a))(fraction(a) - fraction(b))
end

function *(a::MPolyLocalRingElem{BaseRingType, RingElemType}, b::MPolyLocalRingElem{BaseRingType, RingElemType}) where {BaseRingType, RingElemType}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  return (parent(a))(fraction(a) * fraction(b))
end

function Base.:(//)(a::MPolyLocalRingElem{BaseRingType, RingElemType}, b::MPolyLocalRingElem{BaseRingType, RingElemType}) where {BaseRingType, RingElemType}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  numerator(fraction(b)) in inverted_set(parent(b)) || error("the second argument is not a unit in this local ring")
  return (parent(a))(fraction(a) // fraction(b))
end

function ==(a::T, b::T) where {T<:MPolyLocalRingElem}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  return fraction(a) == fraction(b)
end

function ^(a::T, i::Oscar.IntegerUnion) where {T<:MPolyLocalRingElem}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  return parent(a)(fraction(a)^i)
end


############################################################################
# Localizations of polynomial rings at certain maximal ideals              #
############################################################################

@Markdown.doc """
    MPolyMaxLocalRing{BaseRingType, BaseRingElemType, RingElemType} <: AbsLocalizedRing{MPolyRing{BaseRingType}, MPolyElem{BaseRingType}, MPolyComplementOfPrimeIdeal{BaseRingType, RingElemType}

The localization of a multivariate polynomial ring R = 𝕜[x₁,…,xₙ]
with elements of type `RingElemType` over a base field 𝕜 of type
`BaseRingType` and with elements of type `BaseRingElemType` at a
maximal ideal ``𝔪 = ⟨x₁-a₁,…,xₙ-aₙ⟩`` with coefficients
``aᵢ∈ 𝕜``.
"""
mutable struct MPolyMaxLocalRing{BaseRingType, BaseRingElemType, RingElemType} <: AbsLocalizedRing{MPolyRing{BaseRingType}, MPolyElem{BaseRingType}, MPolyComplementOfPrimeIdeal{BaseRingType, RingElemType}}
  R::MPolyRing{BaseRingType} # The parent ring which is being localized
  S::MPolyComplementOfMaxIdeal{BaseRingType, BaseRingElemType, RingElemType} 

  function MPolyMaxLocalRing(
      R::MPolyRing{BaseRingType}, 
      S::MPolyComplementOfMaxIdeal{BaseRingType, BaseRingElemType, RingElemType}
    ) where {BaseRingType, BaseRingElemType, RingElemType}
    # TODO: Add some sanity checks here?
    R_loc = new{BaseRingType, BaseRingElemType, RingElemType}(R, S)
    return R_loc
  end
end

function localize_at(S::MPolyComplementOfMaxIdeal{BaseRingType, BaseRingElemType, RingElemType}) where {BaseRingType, BaseRingElemType, RingElemType}
  return MPolyMaxLocalRing(ambient_ring(S), S)
end

function MPolyMaxLocalRing( R::MPolyRing{BaseRingType}, P::Ideal{RingElemType} ) where {BaseRingType, RingElemType}
  return MPolyLocalRing(R, ComplementOfMaxIdeal(P))
end

function original_ring(W::MPolyMaxLocalRing{BaseRingType, BaseRingElemType, RingElemType}) where {BaseRingType, BaseRingElemType, RingElemType}
  return W.R
end

function inverted_set(W::MPolyMaxLocalRing{BaseRingType, BaseRingElemType, RingElemType}) where {BaseRingType, BaseRingElemType, RingElemType}
  return W.S
end


########################################################################
# Elements of localizations of polynomial rings at maximal ideals      #
########################################################################

mutable struct MPolyMaxLocalRingElem{BaseRingType, BaseRingElemType, RingElemType} <: AbsLocalizedRingElem{MPolyRing{BaseRingType}, MPolyElem{BaseRingType}, MPolyComplementOfPrimeIdeal{BaseRingType, RingElemType}} 
  frac::AbstractAlgebra.Generic.Frac{RingElemType}
  R_loc::MPolyMaxLocalRing{BaseRingType, BaseRingElemType, RingElemType}

  function MPolyMaxLocalRingElem(f::AbstractAlgebra.Generic.Frac{RingElemType}, R_loc::MPolyMaxLocalRing{BaseRingType, BaseRingElemType, RingElemType}) where {BaseRingType, BaseRingElemType, RingElemType}
    base_ring(parent(f)) == original_ring(R_loc) || error("the numerator and denominator of the given fraction do not belong to the original ring before localization")
    denominator(f) in inverted_set(R_loc) || error("the given denominator is not admissible for this localization")
    return new{BaseRingType, BaseRingElemType, RingElemType}(f, R_loc)
  end
end

fraction(a::MPolyMaxLocalRingElem{BaseRingType, BaseRingElemType, RingElemType}) where {BaseRingType, BaseRingElemType, RingElemType} = a.frac
numerator(a::MPolyMaxLocalRingElem{BaseRingType, BaseRingElemType, RingElemType}) where {BaseRingType, BaseRingElemType, RingElemType} = numerator(a.frac)
denominator(a::MPolyMaxLocalRingElem{BaseRingType, BaseRingElemType, RingElemType}) where {BaseRingType, BaseRingElemType, RingElemType} = denominator(a.frac)

parent(a::MPolyMaxLocalRingElem{BaseRingType, BaseRingElemType, RingElemType}) where {BaseRingType, BaseRingElemType, RingElemType} = a.R_loc

(W::MPolyMaxLocalRing{BaseRingType, BaseRingElemType, RingElemType})(f::RingElemType) where {BaseRingType, BaseRingElemType, RingElemType} = MPolyLocalRingElem((FractionField(original_ring(W)))(f), W)

(W::MPolyMaxLocalRing{BaseRingType, BaseRingElemType, RingElemType})(f::AbstractAlgebra.Generic.Frac{RingElemType}) where {BaseRingType, BaseRingElemType, RingElemType} = MPolyLocalRingElem(f, W)


### The main workhorse for binding to the Singular backend
@Markdown.doc """
    LocalizedBiPolyArray{LocRingType, RingElemType}

Main workhorse for binding of ideals in localizations ``R[S⁻¹]`` of type `LocRingType` 
of polynomial rings ``R`` with elements of type `RingElemType` to the singular backend. 

This is supposed to be used for all types of localizations, but with different 
orderings for the polynomial ring on the singular side, depending on the 
algorithmic implementation.
"""
mutable struct LocalizedBiPolyArray{LocRingType, RingElemType}
  # The generators on the Oscar side, given simply as fractions
  oscar_gens::Vector{AbstractAlgebra.Generic.Frac{RingElemType}}
  # The numerators of the above fractions as elements in the singular version 
  # of the original ring before localization.
  singular_gens::Singular.sideal
  # The localized ring on the Oscar side.
  oscar_ring::LocRingType
  # A polynomial ring on the singular side, mapping surjectively to the 
  # original ring before localization. 
  singular_ring::Singular.PolyRing
  # The ordering used for the above singular ring.
  ordering::Symbol
  # An optional shift vector applied to the polynomials in Oscar when 
  # translating them to the Singular ring. 
  shift::Vector{RingElemType}
  # Flag for caching
  is_groebner_basis::Bool

  function LocalizedBiPolyArray(oscar_ring::LocRingType, 
      oscar_gens::Vector{AbstractAlgebra.Generic.Frac{RingElemType}};
      ordering=:degrevlex, shift=Vector{RingElemType}()
    ) where {LocRingType <: AbsLocalizedRing, RingElemType}
    BPA = new{LocRingType, RingElemType}()
    # TODO: Add some sanity checks here
    BPA.oscar_ring = oscar_ring
    BPA.oscar_gens = oscar_gens
    BPA.ordering = ordering
    # fill up the shift vector with zeroes if it is not provided in full length
    for i in (length(shift)+1:length(gens(original_ring(oscar_ring))))
      append!(shift, zero(oscar_ring))
    end
    BPA.shift = shift
    BPA.is_groebner_basis=false
    return BPA
  end
  
  function LocalizedBiPolyArray(oscar_ring::LocRingType, 
      singular_gens::Singular.sideal, shift::Vector{RingElemType}
    ) where {LocRingType <: AbsLocalizedRing, RingElemType}
    BPA = new{LocRingType, RingElemType}()
    # TODO: Add some sanity checks here
    BPA.oscar_ring = oscar_ring
    BPA.singular_gens = singular_gens
    BPA.singular_ring = base_ring(singular_gens)
    R = original_ring(oscar_ring)
    @show BPA.singular_ring
    @show typeof(BPA.singular_ring)
    BPA.ordering = Singular.ordering_as_symbol(BPA.singular_ring)
    # fill up the shift vector with zeroes if it is not provided in full length
    for i in (length(shift)+1:length(gens(R)))
	      append!(shift, zero(R))
    end
    BPA.shift = shift
    nvars(R) == length(shift) || error("the number of variables does not coincide with the number of coordinates")
    inv_shift_hom = AlgebraHomomorphism(R, R, [gen(R, i) + shift[i] for i in (1:nvars(R))])
    BPA.oscar_gens = [ y//one(y) for y in (inv_shift_hom.(R.([x for x in gens(singular_gens)])))]
    BPA.is_groebner_basis=false
    return BPA
  end
end

oscar_gens(I::LocalizedBiPolyArray{LocRingType, RingElemType}) where {LocRingType, RingElemType} = I.oscar_gens
oscar_ring(I::LocalizedBiPolyArray{LocRingType, RingElemType}) where {LocRingType, RingElemType} = I.oscar_ring

function singular_gens(lbpa::LocalizedBiPolyArray{LocRingType, RingElemType}) where {LocRingType, RingElemType}
  singular_assure!(lbpa)
  return lbpa.singular_gens
end

function singular_ring(lbpa::LocalizedBiPolyArray{LocRingType, RingElemType}) where {LocRingType, RingElemType} 
  singular_assure!(lbpa)
  return lbpa.singular_ring
end
ordering(lbpa::LocalizedBiPolyArray{LocRingType, RingElemType}) where {LocRingType, RingElemType} = lbpa.ordering
shift(lbpa::LocalizedBiPolyArray{LocRingType, RingElemType}) where {LocRingType, RingElemType} = lbpa.shift


function _singular_ring(oscar_ring::MPolyMaxLocalRing{BaseRingType, RingElemType, MultSetType}; ord::Symbol = :degrevlex) where {BaseRingType, RingElemType, MultSetType}
  return Singular.PolynomialRing(Oscar.singular_ring(base_ring(original_ring(oscar_ring))), 
				 [string(a) for a = Nemo.symbols(original_ring(oscar_ring))], 
				 ordering = ord, 
				 cached = false)[1]
end

function singular_assure!(lbpa::LocalizedBiPolyArray{RingElemType, MultSetType}) where {RingElemType, MultSetType}
  if !isdefined(lbpa, :singular_ring)
    lbpa.singular_ring = _singular_ring(oscar_ring(lbpa), ord=ordering(lbpa))
  end
  if !isdefined(lbpa, :singular_gens)
    shift_hom = hom(original_ring(oscar_ring(lbpa)), original_ring(oscar_ring(lbpa)), 
        [gen(original_ring(oscar_ring(lbpa)), i) - lbpa.shift[i] for i in (1:nvars(original_ring(oscar_ring(lbpa))))])
    lbpa.singular_gens = Singular.Ideal(lbpa.singular_ring,
	[lbpa.singular_ring(shift_hom(numerator(x))) for x in oscar_gens(lbpa)])
  end
end


@Markdown.doc """
    MPolyMaxLocalIdeal{BaseRingType, BaseRingElemType, RingElemType} <: AbsMaxLocalRingIdeal{MPolyRing{BaseRingType}, MPolyMaxLocalRing{BaseRingType, BaseRingElemType, RingElemType}, RingElemType} 

Ideals in localizations of polynomial rings ``R = 𝕜[x₁,…,xₙ]`` at maximal ideals 
``𝔪 = ⟨x₁-a₁,…,xₙ-aₙ⟩`` with coefficients ``aᵢ∈ 𝕜``.
"""
mutable struct MPolyMaxLocalIdeal{
    BaseRingType, 
    BaseRingElemType, 
    RingElemType
  } <: AbsLocalizedRingIdeal{
    MPolyRing{BaseRingType}, 
    MPolyMaxLocalRingElem{BaseRingType, BaseRingElemType, MPolyElem{BaseRingType}}, 
    MPolyComplementOfMaxIdeal{BaseRingType, BaseRingElemType, RingElemType}
  }
  # the initial set of generators, not to be changed ever!
  gens::LocalizedBiPolyArray{MPolyMaxLocalRing{BaseRingType, BaseRingElemType, RingElemType}, RingElemType}
  # the ambient ring for this ideal
  R_loc::MPolyMaxLocalRing{BaseRingType, BaseRingElemType, RingElemType}
  
  # Fields for caching
  groebner_basis::LocalizedBiPolyArray{MPolyMaxLocalRing{BaseRingType, BaseRingElemType, RingElemType}, RingElemType}
  dimension::Int
 

  function MPolyMaxLocalIdeal(R_loc::MPolyMaxLocalRing{BaseRingType, BaseRingElemType, RingElemType}, f::Vector{AbstractAlgebra.Generic.Frac{RingElemType}}) where {BaseRingType, BaseRingElemType, RingElemType}
    R = original_ring(R_loc)
    k = base_ring(R)
    S = inverted_set(R_loc)
    a = point_coordinates(S)
    for x in f
      denominator(x) in S || error("fraction is not an element of the localization")
    end
    I = new{BaseRingType, BaseRingElemType, RingElemType}()
    I.gens = LocalizedBiPolyArray(R_loc, f, ordering=:negdegrevlex, shift=R.(a))
    I.R_loc = R_loc
    return I
  end
end

gens(I::MPolyMaxLocalIdeal{BaseRingType, BaseRingElemType, RingElemType}) where {BaseRingType, BaseRingElemType, RingElemType} = I.gens
ambient_ring(I::MPolyMaxLocalIdeal{BaseRingType, BaseRingElemType, RingElemType}) where {BaseRingType, BaseRingElemType, RingElemType} = I.R_loc

function shift(I::MPolyMaxLocalIdeal{BaseRingType, BaseRingElemType, RingElemType}) where {BaseRingType, BaseRingElemType, RingElemType}
  lbpa = I.gens
  return shift(lbpa)
end

function dimension(I::MPolyMaxLocalIdeal{BaseRingType, BaseRingElemType, RingElemType}) where {BaseRingType, BaseRingElemType, RingElemType} 
  error("not implemented")
end

function MPolyMaxLocalIdeal(R_loc::MPolyMaxLocalRing{BaseRingType, BaseRingElemType, RingElemType}, I::MPolyIdeal{RingElemType}) where {BaseRingType, BaseRingElemType, RingElemType} 
  if length(gens(I))==0 
    return MPolyMaxLocalIdeal(R_loc, Vector{AbstractAlgebra.Generic.Frac{RingElemType}}())
  end
  R = original_ring(R_loc)
  Q = AbstractAlgebra.Generic.FractionField(R)
  base_ring(I) == R || error("ideal does not belong to the original ring before localization")
  return MPolyMaxLocalIdeal(R_loc, Q.(gens(I)))
end

function ideal(R_loc::MPolyMaxLocalRing{BaseRingType, BaseRingElemType, RingElemType}, I::MPolyIdeal{RingElemType}) where {BaseRingType, BaseRingElemType, RingElemType} 
  return MPolyMaxLocalIdeal(R_loc, I)
end


###############################################################################
# Groebner bases                                                              #
###############################################################################

function groebner_assure(I::MPolyMaxLocalIdeal{BaseRingType, BaseRingElemType, RingElemType}) where {BaseRingType, BaseRingElemType, RingElemType}
  if isdefined(I, :groebner_basis) 
    return I.groebner_basis
  end
  gb = Singular.std(singular_gens(gens(I)))
end

function groebner_basis(I::MPolyMaxLocalIdeal{BaseRingType, BaseRingElemType, RingElemType}; ord::Symbol = :negdegrevlex) where {BaseRingType, BaseRingElemType, RingElemType}
  if ord != ordering(gens(I))
    B = LocalizedBiPolyArray(ambient_ring(I), oscar_gens(gens(I)), ordering = ord)
    singular_assure!(B)
    R = singular_ring(B)
    !Oscar.Singular.has_local_ordering(R) && error("The ordering has to be a local ordering.")
    gb = Singular.std(singular_gens(B))
    return LocalizedBiPolyArray(ambient_ring(I), gb, shift(I))
  end
  if !isdefined(I, :groebner_basis)
    B = LocalizedBiPolyArray(ambient_ring(I), oscar_gens(gens(I)); ordering = ord, shift = shift(I))
    singular_assure!(B)
    R = singular_ring(B)
    !Oscar.Singular.has_local_ordering(R) && error("The ordering has to be a local ordering.")
    gb = Singular.std(singular_gens(B))
    I.groebner_basis = LocalizedBiPolyArray(ambient_ring(I), gb, shift(I))
  end
  return oscar_gens(I.groebner_basis)
end

