export MPolyComplementOfPrimeIdeal, MPolyComplementOfKPointIdeal

export original_ring, inverted_set
export reduce_fraction

export fraction, parent

export MPolyLocalizedRing

export MPolyLocalizedIdeal

export localize_at, ideal

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
        RingElemType
      } <: AbsMultSet{
        MPolyRing{BaseRingType}, 
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


############################################################################
# Localizations of polynomial rings over admissible fields at prime ideals #
############################################################################
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
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  return parent(a)(fraction(a)^i)
end

### implementation of Oscar's general ring interface
one(W::MPolyLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = W(one(original_ring(W)))
zero(W::MPolyLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = W(zero(original_ring(W)))

elem_type(W::MPolyLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyLocalizedRingElem{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}
elem_type(T::Type{MPolyLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyLocalizedRingElem{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}

parent_type(W::MPolyLocalizedRingElem{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}
parent_type(T::Type{MPolyLocalizedRingElem{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}}) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType} = MPolyLocalizedRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}


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
    MPolyComplementOfKPointIdeal{BaseRingType, BaseRingElemType, RingElemType}
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

