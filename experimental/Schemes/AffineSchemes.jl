import AbstractAlgebra.Ring
import Base: intersect

export Scheme
export Spec, OO, defining_ideal
export spec_type, ring_type
export base_ring_type, base_ring_elem_type, poly_type, poly_ring_type, mult_set_type, ring_type
export affine_space, empty_spec
export EmptyScheme

export is_open_embedding, is_closed_embedding, canonically_isomorphic, hypersurface_complement, subscheme, name_of, set_name!
export closure, product

export SpecMor, morphism_type
export pullback, domain, codomain, preimage, restrict, graph, identity_map, inclusion_map, is_isomorphism, is_inverse_of
# TODO for Tommy: Find out why the following are necessary
AbstractAlgebra.promote_rule(::Type{gfp_mpoly}, ::Type{fmpz}) = gfp_mpoly
AbstractAlgebra.promote_rule(::Type{gfp_elem}, ::Type{fmpz}) = gfp_elem
AbstractAlgebra.promote_rule(::Type{gfp_elem}, ::Type{AbstractAlgebra.Generic.Frac{gfp_mpoly}}) = AbstractAlgebra.Generic.Frac{gfp_mpoly}

@Markdown.doc """
    Scheme{BaseRingType<:Ring, BaseRingElemType<:RingElement} 

A scheme over a ring ``𝕜`` of type `BaseRingType` with elements 
of type `BaseRingElemType`.
"""
abstract type Scheme{BaseRingType<:Ring, BaseRingElemType<:RingElement} end

struct EmptyScheme{BaseRingType, BaseRingElemType}<:Scheme{BaseRingType, BaseRingElemType} 
  k::BaseRingType
  function EmptyScheme(k::BaseRingType) where {BaseRingType<:Ring}
    return new{BaseRingType, elem_type(k)}(k)
  end
end

@Markdown.doc """
    Spec{BRT, BRET, RT, RET, MST} <: Scheme{BRT, BRET}

An affine scheme ``X = Spec ((R/I)[S⁻¹])`` with ``R = k[x₁,…,xₙ]`` a free 
polynomial algebra of type `RT` over a base ring ``k`` of type 
`BRT`, ``I ⊂ R`` a finitely generated ideal 
with elements of type `RET`, and ``S`` a multiplicative set in ``R`` of 
type `MST`.
"""
mutable struct Spec{BRT, BRET, RT, RET, MST} <: Scheme{BRT, BRET}
  # the basic fields 
  OO::MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MST}
  # fields for caching
  name::String # the name of this scheme for printing

  function Spec(OO::MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST}
    return new{BRT, BRET, RT, RET, MST}(OO)
  end
end

### Type getters

ring_type(::Type{Spec{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MST}
ring_type(X::Spec) = ring_type(typeof(X))

base_ring_type(X::Spec{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST} = BRT
base_ring_elem_type(X::Spec{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST} = BRET
mult_set_type(X::Spec{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST} = MST
poly_ring_type(X::Spec{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST} = RT
poly_type(X::Spec{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST} = RET

base_ring_type(::Type{Spec{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = BRT
base_ring_elem_type(::Type{Spec{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = BRET
mult_set_type(::Type{Spec{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = MST
poly_ring_type(::Type{Spec{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = RT
poly_type(::Type{Spec{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = RET

### type constructors

# this defaults to Specs over localizations of polynomial rings at hypersurfaces
# and does not cover localizations for germs!
spec_type(R::T) where {T<:AbstractAlgebra.Ring} = Spec{T, elem_type(T), mpoly_ring_type(T), mpoly_type(T), MPolyPowersOfElement{T, elem_type(T), mpoly_ring_type(T), mpoly_type(T)}}
spec_type(::Type{T}) where {T<:AbstractAlgebra.Ring} = Spec{T, elem_type(T), mpoly_ring_type(T), mpoly_type(T), MPolyPowersOfElement{T, elem_type(T), mpoly_ring_type(T), mpoly_type(T)}}
spec_type(L::MPolyQuoLocalizedRing{S, T, U, V, W}) where {S, T, U, V, W} = Spec{S, T, U, V, W}
spec_type(::Type{MPolyQuoLocalizedRing{S, T, U, V, W}}) where {S, T, U, V, W} = Spec{S, T, U, V, W}



### Getter functions

@Markdown.doc """
    OO(X::Spec)

For ``X = Spec ((𝕜[x₁,…,xₙ]/I)[S⁻¹])`` this returns ``(𝕜[x₁,…,xₙ]/I)[S⁻¹]``.
"""
OO(X::Spec) = X.OO

function name_of(X::Spec) 
  if isdefined(X, :name)
    return X.name
  end
  return "unnamed affine variety"
end

function set_name!(X::Spec, name::String) 
  X.name = name
end

function Base.show(io::IO, X::Spec) 
  if isdefined(X, :name)
    print(io, name_of(X))
    return
  end
  print(io, "Spec of $(OO(X))")
end

base_ring(X::Spec) = coefficient_ring(base_ring(OO(X)))
defining_ideal(X::Spec) = modulus(OO(X))

### Copy constructor
Spec(X::Spec) = Spec(OO(X))

### internal copy routines
Base.deepcopy_internal(X::Spec, dict::IdDict) = Spec(deepcopy_internal(OO(X), dict))

### additional constructors 
Spec(R::MPolyRing) = Spec(MPolyQuoLocalizedRing(R))

Spec(Q::MPolyQuo) = Spec(MPolyQuoLocalizedRing(Q))

Spec(W::MPolyLocalizedRing) = Spec(MPolyQuoLocalizedRing(W))

Spec(R::MPolyRing, I::MPolyIdeal, U::AbsMPolyMultSet) = Spec(MPolyQuoLocalizedRing(R, I, U))
Spec(R::MPolyRing, U::AbsMPolyMultSet) = Spec(MPolyQuoLocalizedRing(R, ideal(R, [zero(R)]), U))
Spec(R::MPolyRing, I::MPolyIdeal) = Spec(MPolyQuoLocalizedRing(R, I, units_of(R)))

# Hack for the construction of the empty scheme over kk 
# as an instance of Spec
function empty_spec(kk::BRT) where {BRT<:AbstractAlgebra.Ring} 
  R, (x,) = PolynomialRing(kk, ["x"])
  return Spec(R, ideal(R, [x]), MPolyPowersOfElement(x))
end

### closed subschemes defined by ideals
function subscheme(X::Spec{BRT, BRET, RT, RET, MST}, I::MPolyLocalizedIdeal{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST}
  localized_ring(OO(X)) == base_ring(I) || error("ideal does not live in the correct ring")
  return Spec(quo(localized_ring(OO(X)), I + localized_modulus(OO(X))))
end

function subscheme(X::Spec{BRT, BRET, RT, RET, MST}, I::MPolyIdeal{RET}) where {BRT, BRET, RT, RET, MST}
  base_ring(OO(X)) == base_ring(I) || error("ideal does not live in the correct ring")
  return Spec(quo(OO(X), I))
end
  
@Markdown.doc """
    subscheme(X::Spec{BRT, BRET, RT, RET, MST}, f::RET) where {BRT, BRET, RT, RET, MST}

For a scheme ``X = Spec ((𝕜[x₁,…,xₙ]/I)[S⁻¹])`` and an element ``f ∈ 𝕜[x₁,…,xₙ]`` 
this returns the closed subscheme defined by the ideal ``I' = I + ⟨f⟩``.
"""
function subscheme(X::Spec{BRT, BRET, RT, RET, MST}, f::RET) where {BRT, BRET, RT, RET, MST}
  R = base_ring(OO(X))
  I = ideal(R, f)
  return subscheme(X, I)
end

function subscheme(X::Spec{BRT, BRET, RT, RET, MST}, f::Vector{RET}) where {BRT, BRET, RT, RET, MST}
  R = base_ring(OO(X))
  I = ideal(R, f)
  return subscheme(X, I)
end

function subscheme(X::Spec, f::RET) where {RET<:MPolyQuoLocalizedRingElem}
  I = ideal(OO(X), [f])
  return subscheme(X, I)
end

function subscheme(X::Spec, f::Vector{RET}) where {RET<:MPolyQuoLocalizedRingElem}
  I = ideal(OO(X), f)
  return subscheme(X, I)
end

function subscheme(X::Spec, f::RET) where {RET<:MPolyLocalizedRingElem}
  I = ideal(OO(X), [f])
  return subscheme(X, I)
end

function subscheme(X::Spec, f::Vector{RET}) where {RET<:MPolyLocalizedRingElem}
  I = ideal(OO(X), f)
  return subscheme(X, I)
end

function subscheme(X::Spec{BRT, BRET, RT, RET, MST}, f::BRET) where {BRT, BRET, RT, RET, MST}
  R = base_ring(OO(X))
  I = ideal(R, R(f))
  return subscheme(X, I)
end

function subscheme(X::Spec{BRT, BRET, RT, RET, MST}, f::Vector{BRET}) where {BRT, BRET, RT, RET, MST}
  R = base_ring(OO(X))
  I = ideal(R, R.(f))
  return subscheme(X, I)
end

### open subschemes defined by complements of hypersurfaces
@Markdown.doc """
    hypersurface_complement(X::Spec{BRT, BRET, RT, RET, MST}, f::RET) where {BRT, BRET, RT, RET, MST<:MPolyPowersOfElement{BRT, BRET, RT, RET}}

For a scheme ``X = Spec ((𝕜[x₁,…,xₙ]/I)[S⁻¹])`` and an element ``f ∈ 𝕜[x₁,…,xₙ]`` 
this returns the open subscheme ``U = X ∖ V(f)`` defined by the complement of the vanishing 
locus of ``f``.
"""
function hypersurface_complement(X::Spec{BRT, BRET, RT, RET, MST}, f::RET) where {BRT, BRET, RT, RET, MST<:MPolyPowersOfElement{BRT, BRET, RT, RET}}
  R = base_ring(OO(X))
  parent(f) == R || error("the element does not belong to the correct ring")
  iszero(f) && return subscheme(X, [one(R)])
  return Spec(Localization(OO(X), MPolyPowersOfElement(f)))
end

function hypersurface_complement(
    X::Spec{BRT, BRET, RT, RET, MST}, 
    f::MPolyLocalizedRingElem{BRT, BRET, RT, RET, MST}
  ) where {BRT, BRET, RT, RET, MST<:MPolyPowersOfElement{BRT, BRET, RT, RET}}
  return Spec(Localization(OO(X), MPolyPowersOfElement(numerator(f))))
end

function hypersurface_complement(
    X::Spec{BRT, BRET, RT, RET, MST}, 
    f::MPolyQuoLocalizedRingElem{BRT, BRET, RT, RET, MST}
  ) where {BRT, BRET, RT, RET, MST<:MPolyPowersOfElement{BRT, BRET, RT, RET}}
  parent(f) == OO(X) || error("the element does not belong to the correct ring")
  return Spec(Localization(OO(X), MPolyPowersOfElement(lifted_numerator(f))))
end


### testing containment
issubset(X::EmptyScheme{BRT, BRET}, Y::Scheme{BRT, BRET}) where {BRT, BRET} = true

function issubset(Y::Spec{BRT, BRET, RT, RET, MST1}, X::EmptyScheme{BRT, BRET}) where {BRT, BRET, RT, RET, MST1} 
  return iszero(one(OO(Y)))
end

@Markdown.doc """
    issubset(
      X::Spec{BRT, BRET, RT, RET, MST1}, 
      Y::Spec{BRT, BRET, RT, RET, MST2}
    ) where {BRT, BRET, RT, RET, MST1<:MPolyPowersOfElement{BRT, BRET, RT, RET}, MST2<:MPolyPowersOfElement{BRT, BRET, RT, RET}}

Checks whether ``X`` is a subset of ``Y`` based on the comparison of their coordinate rings.
"""
function issubset(
    X::Spec{BRT, BRET, RT, RET, MST1}, 
    Y::Spec{BRT, BRET, RT, RET, MST2}
  ) where {BRT, BRET, RT, RET, MST1<:MPolyPowersOfElement{BRT, BRET, RT, RET}, MST2<:MPolyPowersOfElement{BRT, BRET, RT, RET}}
  R = base_ring(OO(X))
  R == base_ring(OO(Y)) || error("schemes can not be compared")
  UX = inverted_set(OO(X))
  UY = inverted_set(OO(Y))
  if !issubset(UY, UX) 
    # check whether the inverted elements in Y are units anyway
    for a in denominators(UY)
      isunit(OO(X)(a)) || return false
    end
  end
  J = localized_ring(OO(X))(modulus(OO(Y)))
  return issubset(J, localized_modulus(OO(X)))
end

function ==(X::T, Y::T) where {T<:Spec}
  return X === Y
end

function canonically_isomorphic(
    X::Spec{BRT, BRET, RT, RET, MST1}, 
    Y::Spec{BRT, BRET, RT, RET, MST2}
  ) where {BRT, BRET, RT, RET, MST1<:MPolyPowersOfElement{BRT, BRET, RT, RET}, MST2<:MPolyPowersOfElement{BRT, BRET, RT, RET}}
  isempty(X) && isempty(Y) && return true
  base_ring(OO(X)) == base_ring(OO(Y)) || return false
  return issubset(X, Y) && issubset(Y, X)
end

function canonically_isomorphic(X::Spec, Y::EmptyScheme)
  return issubset(X, Y)
end

canonically_isomorphic(X::EmptyScheme, Y::Spec) = canonically_isomorphic(Y, X)

Base.isempty(X::Spec) = iszero(one(OO(X)))

@Markdown.doc """
    is_open_embedding(
      X::Spec{BRT, BRET, RT, RET, MST1}, 
      Y::Spec{BRT, BRET, RT, RET, MST2}
    ) where {BRT, BRET, RT, RET, MST1<:MPolyPowersOfElement{BRT, BRET, RT, RET}, MST2<:MPolyPowersOfElement{BRT, BRET, RT, RET}}

Checks whether ``X`` is openly embedded in ``Y``.
"""
function is_open_embedding(
    X::Spec{BRT, BRET, RT, RET, MST1}, 
    Y::Spec{BRT, BRET, RT, RET, MST2}
  ) where {BRT, BRET, RT, RET, MST1<:MPolyPowersOfElement{BRT, BRET, RT, RET}, MST2<:MPolyPowersOfElement{BRT, BRET, RT, RET}}
  R = base_ring(OO(X))
  R == base_ring(OO(Y)) || return false
  UX = inverted_set(OO(X))
  UY = inverted_set(OO(Y))
  issubset(UY, UX) || return false
  J = localized_ring(OO(X))(modulus(OO(Y)))
  return localized_modulus(OO(X)) == J 
end

@Markdown.doc """
    is_closed_embedding(
      X::Spec{BRT, BRET, RT, RET, MST1}, 
      Y::Spec{BRT, BRET, RT, RET, MST2}
    ) where {BRT, BRET, RT, RET, MST1<:MPolyPowersOfElement{BRT, BRET, RT, RET}, MST2<:MPolyPowersOfElement{BRT, BRET, RT, RET}}

Checks whether ``X`` is closed embedded in ``Y``.
"""
function is_closed_embedding(
    X::Spec{BRT, BRET, RT, RET, MST1}, 
    Y::Spec{BRT, BRET, RT, RET, MST2}
  ) where {BRT, BRET, RT, RET, MST1<:MPolyPowersOfElement{BRT, BRET, RT, RET}, MST2<:MPolyPowersOfElement{BRT, BRET, RT, RET}}
  R = base_ring(OO(X))
  R == base_ring(OO(Y)) || return false
  inverted_set(OO(X)) == inverted_set(OO(Y)) || return false
  J = localized_ring(OO(X))(modulus(OO(Y)))
  return issubset(J, localized_modulus(OO(X)))
end

### set operations

intersect(E::EmptyScheme{BRT, BRET}, X::Scheme{BRT, BRET}) where {BRT, BRET} = E
intersect(X::Scheme{BRT, BRET}, E::EmptyScheme{BRT, BRET}) where {BRT, BRET} = E
intersect(X::EmptyScheme{BRT, BRET}, E::EmptyScheme{BRT, BRET}) where {BRT, BRET} = E

function Base.intersect(
    X::Spec{BRT, BRET, RT, RET, MST1}, 
    Y::Spec{BRT, BRET, RT, RET, MST2}
  ) where {BRT, BRET, RT, RET, MST1<:MPolyPowersOfElement{BRT, BRET, RT, RET}, MST2<:MPolyPowersOfElement{BRT, BRET, RT, RET}}
  base_ring(OO(X)) == base_ring(OO(Y)) || error("schemes can not be intersected")
  issubset(X, Y) && return X
  issubset(Y, X) && return Y
  UX = inverted_set(OO(X))
  UY = inverted_set(OO(Y))
  U = UX*UY
  IX = modulus(OO(X))
  IY = modulus(OO(Y))
  I = IX + IY
  R = base_ring(OO(X))
  L = MPolyQuoLocalizedRing(R, I, U)
  one(localized_ring(L)) in localized_modulus(L) && return EmptyScheme(coefficient_ring(R))
  return Spec(L)
end

### compute the closure of X in Y
@Markdown.doc """
    closure(
      X::Spec{BRT, BRET, RT, RET, MST1}, 
      Y::Spec{BRT, BRET, RT, RET, MST2}
    ) where {BRT, BRET, RT, RET, MST1<:MPolyPowersOfElement{BRT, BRET, RT, RET}, MST2<:MPolyPowersOfElement{BRT, BRET, RT, RET}}

Returns the closure of ``X`` in ``Y``.
"""
function closure(
    X::Spec{BRT, BRET, RT, RET, MST1}, 
    Y::Spec{BRT, BRET, RT, RET, MST2}
  ) where {BRT, BRET, RT, RET, MST1<:MPolyPowersOfElement{BRT, BRET, RT, RET}, MST2<:MPolyPowersOfElement{BRT, BRET, RT, RET}}
  issubset(X, Y) || error("the first argument is not a subset of the second")
  is_closed_embedding(X, Y) && return X
  W = Localization(inverted_set(OO(X))*inverted_set(OO(Y)))
  I = ideal(W, W.(gens(modulus(OO(X)))))
  lbpa = groebner_basis(I) # takes care of the saturation
  R = base_ring(OO(Y))
  return Spec(MPolyQuoLocalizedRing(R, ideal(R, numerator.(oscar_gens(lbpa))), inverted_set(OO(Y))))
end



########################################################################
# Morphisms of affine schemes                                      #
########################################################################

@attributes mutable struct SpecMor{DomainType<:Spec, CodomainType<:Spec, PullbackType<:MPolyQuoLocalizedRingHom}
  domain::DomainType
  codomain::CodomainType
  pullback::PullbackType

  function SpecMor(
      X::DomainType,
      Y::CodomainType,
      pullback::PullbackType
    ) where {DomainType<:Spec, CodomainType<:Spec, PullbackType<:MPolyQuoLocalizedRingHom}
    OO(X) == codomain(pullback) || error("the coordinate ring of the domain does not coincide with the codomain of the pullback")
    OO(Y) == domain(pullback) || error("the coordinate ring of the codomain does not coincide with the domain of the pullback")
    return new{DomainType, CodomainType, PullbackType}(X, Y, pullback)
  end
end

function morphism_type(::Type{SpecType1}, ::Type{SpecType2}) where {SpecType1<:Spec, SpecType2<:Spec}
  return SpecMor{SpecType1, SpecType1, morphism_type(ring_type(SpecType2), ring_type(SpecType1))}
end

morphism_type(X::Spec, Y::Spec) = morphism_type(typeof(X), typeof(Y))


### getter functions
pullback(phi::SpecMor) = phi.pullback
domain(phi::SpecMor) = phi.domain
codomain(phi::SpecMor) = phi.codomain

### additional constructors
function SpecMor(
      X::Spec{BRT, BRET, RT, RET, MST1},
      Y::Spec{BRT, BRET, RT, RET, MST2},
      f::Vector
  ) where {BRT, BRET, RT, RET, MST1, MST2}
  return SpecMor(X, Y, MPolyQuoLocalizedRingHom(OO(Y), OO(X), f))
end

identity_map(X::Spec) = SpecMor(X, X, gens(base_ring(OO(X))))
inclusion_map(X::T, Y::T) where {T<:Spec} = SpecMor(X, Y, gens(base_ring(OO(Y))))

function restrict(f::SpecMor, U::Spec, V::Spec; check::Bool=true)
  if check
    issubset(U, domain(f)) || error("second argument does not lay in the domain of the map")
    issubset(V, codomain(f)) || error("third argument does not lay in the codomain of the map")
    issubset(U, preimage(f, V)) || error("the image of the restriction is not contained in the restricted codomain")
  end
  return SpecMor(U, V, images(pullback(f)))
end

function compose(f::SpecMorType, g::SpecMorType) where {SpecMorType<:SpecMor}
  codomain(f) == domain(g) || error("Morphisms can not be composed")
  return SpecMor(domain(f), codomain(g), compose(pullback(g), pullback(f)))
end

function ==(f::SpecMorType, g::SpecMorType) where {SpecMorType<:SpecMor}
  X = domain(f)
  X == domain(g) || return false
  codomain(f) == codomain(g) || return false
  OO(X).(images(pullback(f))) == OO(X).(images(pullback(g))) || return false
  return true
end

### functionality
function preimage(
    phi::SpecMor{Spec{BRT, BRET, RT, RET, MST1}, Spec{BRT, BRET, RT, RET, MST2}},
    Z::Spec{BRT, BRET, RT, RET, MST3}
  ) where {BRT, BRET, RT, RET, MST1<:MPolyPowersOfElement{BRT, BRET, RT, RET}, MST2<:MPolyPowersOfElement{BRT, BRET, RT, RET}, MST3<:MPolyPowersOfElement{BRT, BRET, RT, RET}}
  X = domain(phi)
  Y = codomain(phi)
  issubset(Z, Y) || (Z = intersect(Y, Z))
  IZ = modulus(OO(Z))
  a = denominators(inverted_set(OO(Z)))
  R = base_ring(OO(X))
  f = pullback(phi)
  new_units = MPolyPowersOfElement(R, [lifted_numerator(f(d)) for d in a])
  new_gens = lifted_numerator.(f.(gens(IZ)))
  return Spec(MPolyQuoLocalizedRing(R, ideal(R, new_gens) + modulus(OO(X)), inverted_set(OO(X))*new_units))
end

function is_isomorphism(f::SpecMor)
  is_isomorphism(pullback(f)) || return false
  set_attribute!(f, :inverse, SpecMor(codomain(f), domain(f), inverse(pullback(f))))
  return true
end

function is_inverse_of(f::S, g::T) where {S<:SpecMor, T<:SpecMor}
  return is_isomorphism(f) && (inverse(f) == g)
end

function inverse(f::SpecMor) 
  if !has_attribute(f, :inverse) 
    is_isomorphism(f) || error("the given morphism is not an isomorphism")
  end
  return get_attribute(f, :inverse)::morphism_type(codomain(f), domain(f))
end

@Markdown.doc """
    product(X::Spec{BRT, BRET, RT, RET, MST}, Y::Spec{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST<:MPolyPowersOfElement}
    
Returns a triple ``(X×Y, p₁, p₂)`` consisting of the product ``X×Y`` and the two projections 
``p₁ : X×Y → X`` and ``p₂ : X×Y → Y``.
"""
function product(X::Spec{BRT, BRET, RT, RET, MST}, Y::Spec{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST<:MPolyPowersOfElement}
  K = OO(X)
  L = OO(Y) 
  V = localized_ring(K)
  W = localized_ring(L)
  R = base_ring(K)
  S = base_ring(L)
  k = base_ring(R)
  k == base_ring(S) || error("varieties are not defined over the same field")

  m = length(gens(R))
  n = length(gens(S))
  RS, z = PolynomialRing(k, vcat(symbols(R), symbols(S)))
  inc1 = AlgebraHomomorphism(R, RS, gens(RS)[1:m])
  inc2 = AlgebraHomomorphism(S, RS, gens(RS)[m+1:m+n])
  #pr1 = AlgebraHomomorphism(RS, R, vcat(gens(R), [zero(R) for i in 1:n]))
  #pr2 = AlgebraHomomorphism(RS, S, vcat([zero(S) for i in 1:m], gens(S)))
  IX = ideal(RS, inc1.(gens(modulus(OO(X)))))
  IY = ideal(RS, inc2.(gens(modulus(OO(Y)))))
  UX = MPolyPowersOfElement(RS, inc1.(denominators(inverted_set(OO(X)))))
  UY = MPolyPowersOfElement(RS, inc2.(denominators(inverted_set(OO(Y)))))
  XxY = Spec(RS, IX + IY, UX*UY)
  pr1 = SpecMor(XxY, X, gens(RS)[1:m])
  pr2 = SpecMor(XxY, Y, gens(RS)[m+1:m+n])
  return XxY, pr1, pr2
end

  
function graph(f::SpecMor{<:Spec{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement}, 
                          <:Spec{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement}}
  )
  X = domain(f)
  Y = codomain(f)
  XxY, prX, prY = product(X, Y)
  pb_X = pullback(prX)
  pb_Y = pullback(prY)
  pb_f = pullback(f)
  I = ideal(localized_ring(OO(XxY)), lift.(pb_X.(images(pb_f)) - pb_Y.(gens(OO(Y)))))
  G = subscheme(XxY, I)
  return G, restrict(prX, G, X), restrict(prY, G, Y)
end

function partition_of_unity(X::Spec{BRT, BRET, RT, RET, MST},
    f::Vector{RET}
  ) where {BRT, BRET, RT, RET, MST}
  error("not implemented")
end

#function hash(X::Spec, u::UInt)
#  r = 0x753087204385757820349592852
#  return xor(r, hash(OO(X), u))
#end
  
function affine_space(kk::BRT, n::Int; variable_name="x") where {BRT<:Ring}
  R, _ = PolynomialRing(kk, [ variable_name * "$i" for i in 1:n])
  return Spec(R)
end

function affine_space(kk::BRT, var_symbols::Vector{Symbol}) where {BRT<:Ring}
  R, _ = PolynomialRing(kk, var_symbols)
  return Spec(R)
end
