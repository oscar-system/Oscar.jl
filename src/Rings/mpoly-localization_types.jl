########################################################################
# Multiplicatively closed sets in multivariate polynomial rings        #
########################################################################

abstract type AbsMPolyMultSet{BRT, BRET, RT, RET} <: AbsMultSet{RT, RET} end


########################################################################
# Powers of elements                                                   #
########################################################################

@doc raw"""
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

```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> powers_of_element(x)
Multiplicative subset
  of multivariate polynomial ring in 2 variables over QQ
  given by the products of [x]

julia> (powers_of_element(x),)
(Products of (x),)

```
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
  is_simplified::Bool
  is_lightly_simplified::Bool

  function MPolyPowersOfElement(R::RingType, a::Vector{RingElemType}) where {RingType<:MPolyRing, RingElemType<:MPolyRingElem}
    for f in a 
      parent(f) == R || error("element does not belong to the given ring")
      !iszero(f) || error("can not localize at the zero element")
    end
    k = coefficient_ring(R)
    return new{typeof(k), elem_type(k), RingType, RingElemType}(R, a, false, false)
  end
end


########################################################################
# Complements of prime ideals                                          #
########################################################################
@doc raw"""
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

```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> S = complement_of_prime_ideal(ideal([x]))
Complement
  of prime ideal (x)
  in multivariate polynomial ring in 2 variables over QQ

julia> (S,)
(Complement of prime ideal (x),)

```
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
      check::Bool=true
    ) where {RingElemType}
    R = base_ring(P)
    @check is_prime(P) "the ideal $P is not prime"
    return new{typeof(coefficient_ring(R)), elem_type(coefficient_ring(R)), typeof(R), elem_type(R)}(R, P)
  end
end

########################################################################
# Complements of maximal ideals corresponding to 𝕜-points              #
########################################################################

@doc raw"""
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

!!! note 
    The coefficient ring is required to be a field.

```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> S = complement_of_point_ideal(R,[1,2])
Complement
  of maximal ideal corresponding to rational point with coordinates (1, 2)
  in multivariate polynomial ring in 2 variables over QQ

julia> (S,)
(Complement of maximal ideal of point (1, 2),)

```
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
  m::MPolyIdeal # Field for caching the associated maximal ideal

  function MPolyComplementOfKPointIdeal(R::RingType, a::Vector{T}) where {RingType<:MPolyRing, T<:RingElement}
    length(a) == ngens(R) || error("the number of variables in the ring does not coincide with the number of coordinates")
    n = length(a)
    kk = coefficient_ring(R)
    @req kk isa Field "This localization is only available over fields"
    b = kk.(a) # fails if the input is not compatible
    S = new{typeof(kk), elem_type(kk), RingType, elem_type(R)}(R, b)
    return S
  end
end

@doc raw"""
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

```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> S = complement_of_point_ideal(R,[1,2]);

julia> T = powers_of_element(x);

julia> ST = Oscar.MPolyProductOfMultSets(R, [S, T])
Product of the multiplicative sets
  complement of maximal ideal of point (1, 2)
  products of (x)

julia> (ST,)
(Product of the multiplicative subsets [complement of maximal ideal, products of 1 element],)

```
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

  function MPolyProductOfMultSets(R::RT, U::Vector{<:AbsMPolyMultSet{BRT, BRET, RT, RET}}) where {BRT<:Ring, BRET<:RingElement, RT<:MPolyRing, RET<:MPolyRingElem}
    for s in U
      ring(s) == R || error("multiplicative set does not live in the given ring")
    end
    return new{typeof(coefficient_ring(R)), elem_type(coefficient_ring(R)), typeof(R), elem_type(R)}(R, U)
  end
end

########################################################################
# Localization associated to a monomial ordering                       #
########################################################################

@doc raw"""
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

########################################################################
# Localizations of polynomial rings over admissible fields             #
########################################################################

@doc raw"""
    MPolyLocRing{
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
@attributes mutable struct MPolyLocRing{
    BaseRingType,
    BaseRingElemType,
    RingType,
    RingElemType,
    MultSetType <: AbsMPolyMultSet{BaseRingType, BaseRingElemType, RingType, RingElemType}
  } <: AbsLocalizedRing{
    RingType,
    RingElemType,
    MultSetType
  }
  R::RingType # The parent ring which is being localized
  S::MultSetType # The multiplicatively closed set that has been inverted 

  function MPolyLocRing(
      R::RingType, 
      S::MultSetType
    ) where {RingType<:MPolyRing, MultSetType<:AbsMPolyMultSet}
    # TODO: Add some sanity checks here?
    ring(S) == R || error("the multiplicative set is not contained in the given ring")
    k = coefficient_ring(R)
    R_loc = new{typeof(k), elem_type(k), RingType, elem_type(R), MultSetType}(R, S)
    return R_loc
  end
end

########################################################################
# Elements of localized polynomial rings                               #
########################################################################

@doc raw"""
    MPolyLocRingElem{
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
mutable struct MPolyLocRingElem{
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

  W::MPolyLocRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}
  frac::AbstractAlgebra.Generic.FracFieldElem{RingElemType}

  function MPolyLocRingElem(
      W::MPolyLocRing{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType},
      f::AbstractAlgebra.Generic.FracFieldElem{RingElemType};
      check::Bool=true
    ) where {BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}
    base_ring(parent(f)) === base_ring(W) || error(
      "the numerator and denominator of the given fraction do not belong to the original ring before localization"
      )
    @check begin
      if !iszero(f) && !is_unit(denominator(f))
        denominator(f) in inverted_set(W) || error("the given denominator is not admissible for this localization")
      end
    end
    return new{BaseRingType, BaseRingElemType, RingType, RingElemType, MultSetType}(W, f)
  end
end

########################################################################
# Ideals in localizations of multivariate polynomial rings             #
########################################################################

@doc raw"""
    MPolyLocalizedIdeal{
        LocRingType<:MPolyLocRing, 
        LocRingElemType<:MPolyLocRingElem
      } <: AbsLocalizedIdeal{LocRingElemType}

Ideals in localizations of polynomial rings.
"""
@attributes mutable struct MPolyLocalizedIdeal{
    LocRingType<:MPolyLocRing, 
    LocRingElemType<:MPolyLocRingElem
  } <: AbsLocalizedIdeal{LocRingElemType}
  # the initial set of generators, not to be changed ever!
  gens::Vector{LocRingElemType}
  # the ambient ring for this ideal
  W::LocRingType

  # fields for caching 
  map_from_base_ring::Map

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
      W::MPolyLocRing, 
      gens::Vector{LocRingElemType};
      map_from_base_ring::Map = MapFromFunc(
          base_ring(W), 
          W,
          W,
          y->(isone(denominator(y)) ? numerator(y) : divexact(numerator(y), denominator(y))),
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

########################################################################
# Homomorphisms of localized polynomial rings                          #
########################################################################

@attributes mutable struct MPolyLocalizedRingHom{
                                     DomainType<:MPolyLocRing, 
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
    ) where {DomainType<:MPolyLocRing, CodomainType<:Ring, RestrictedMapType<:Map}
    R = base_ring(W)
    U = inverted_set(W)
    domain(res) === R || error("the domain of the restricted map does not coincide with the base ring of the localization")
    @check begin
      for f in U
        is_unit(S(res(f))) || error("image of $f is not a unit in the codomain")
      end
    end
    return new{DomainType, CodomainType, RestrictedMapType}(W, S, res) 
  end
end
  
