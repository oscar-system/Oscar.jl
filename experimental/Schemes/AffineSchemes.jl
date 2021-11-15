import AbstractAlgebra.Ring
import Base: intersect

export Spec, OO
export ⊂

export is_open_embedding, is_closed_embedding, hypersurface_complement, subscheme, name_of, set_name!
export closure 

export SpecHom
export pullback, domain, codomain, preimage

@Markdown.doc """
Scheme{ BaseRingType <: Ring }

A scheme over a ring ``k`` of type `BaseRingType`.
"""
abstract type Scheme{BaseRingType<:Ring, BaseRingElemType<:RingElement} end

struct EmptyScheme{BaseRingType, BaseRingElemType}<:Scheme{BaseRingType, BaseRingElemType} 
  k::BaseRingType
  function EmptyScheme(k::BaseRingType) where {BaseRingType<:Ring}
    return new{BaseRingType, elem_type(k)}(k)
  end
end

@doc Markdown.doc"""
AffineScheme{BaseRingType<:Ring, BaseRingElemType<:RingElement, RingType<:MPolyRing, RingElemType<:MPolyElem} <: Scheme{BaseRingType} 

An affine scheme over a base ring ``k`` of type `BaseRingType` and
elements of type `BaseRingElemType`, given by a ring ``R/I`` with ``R``
a polynomial ring of type `RingType` and elements of type `RingElemType`.
"""
abstract type AffineScheme{BaseRingType<:Ring, BaseRingElemType<:RingElement, RingType<:MPolyRing, RingElemType<:MPolyElem, MultSetType<:AbsMPolyMultSet} <: Scheme{BaseRingType, BaseRingElemType} end

@doc Markdown.doc"""
Spec{BRT, BRET, RT, RET} <: AffineScheme{BRT, BRET, RT, RET}

An affine scheme ``X = Spec (R/I)[S⁻¹]`` with ``R = k[x₁,…,xₙ]`` a free 
polynomial algebra of type `RT` over a base ring ``k`` of type 
`BRT` and ``I ⊂ R`` a finitely generated ideal 
with elements of type `RET`.
"""
mutable struct Spec{BRT, BRET, RT, RET, MST} <: AffineScheme{BRT, BRET, RT, RET, MST}
  # the basic fields 
  OO::MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MST}
  # fields for caching
  name::String # the name of this scheme for printing

  function Spec(OO::MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST}
    return new{BRT, BRET, RT, RET, MST}(OO)
  end
end


### Getter functions

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

### Copy constructor
Spec(X::Spec) = Spec(OO(X))

### internal copy routines
Base.deepcopy_internal(X::Spec, dict::IdDict) = Spec(deepcopy_internal(OO(X), dict))

### additional constructors 
Spec(R::MPolyRing) = Spec(MPolyQuoLocalizedRing(R))

Spec(Q::MPolyQuo) = Spec(MPolyQuoLocalizedRing(Q))

Spec(W::MPolyLocalizedRing) = Spec(MPolyQuoLocalizedRing(W))

### closed subschemes defined by ideals
function subscheme(X::Spec{BRT, BRET, RT, RET, MST}, I::MPolyLocalizedIdeal{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST}
  OO(X) == base_ring(I) || error("ideal does not live in the correct ring")
  return Spec(quo(OO(X), I))
end

function subscheme(X::Spec{BRT, BRET, RT, RET, MST}, I::MPolyIdeal{RET}) where {BRT, BRET, RT, RET, MST}
  base_ring(OO(X)) == base_ring(I) || error("ideal does not live in the correct ring")
  return Spec(quo(OO(X), I))
end
  
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
function hypersurface_complement(X::Spec{BRT, BRET, RT, RET, MST}, f::RET) where {BRT, BRET, RT, RET, MST<:MPolyPowersOfElement{BRT, BRET, RT, RET}}
  R = base_ring(OO(X))
  parent(f) == R || error("the element does not belong to the correct ring")
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
⊂(
  X::Scheme{BRT, BRET}, 
  Y::Scheme{BRT, BRET}
 ) where {BRT, BRET} = issubset(X, Y)

issubset(X::EmptyScheme{BRT, BRET}, Y::Scheme{BRT, BRET}) where {BRT, BRET} = true

function issubset(
    X::Spec{BRT, BRET, RT, RET, MST1}, 
    Y::Spec{BRT, BRET, RT, RET, MST2}
  ) where {BRT, BRET, RT, RET, MST1<:MPolyPowersOfElement{BRT, BRET, RT, RET}, MST2<:MPolyPowersOfElement{BRT, BRET, RT, RET}}
  R = base_ring(OO(X))
  R == base_ring(OO(Y)) || return false
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
# Homomorphisms of affine schemes                                      #
########################################################################

mutable struct SpecHom{BRT, BRET, RT, RET, MST1, MST2}
  domain::Spec{BRT, BRET, RT, RET, MST2}
  codomain::Spec{BRT, BRET, RT, RET, MST1}
  pullback::MPolyQuoLocalizedRingHom{BRT, BRET, RT, RET, MST1, MST2}

  function SpecHom(
      X::Spec{BRT, BRET, RT, RET, MST2},
      Y::Spec{BRT, BRET, RT, RET, MST1},
      pullback::MPolyQuoLocalizedRingHom{BRT, BRET, RT, RET, MST1, MST2}
    ) where {BRT, BRET, RT, RET, MST1, MST2}
    OO(X) == codomain(pullback) || error("the coordinate ring of the domain does not coincide with the codomain of the pullback")
    OO(Y) == domain(pullback) || error("the coordinate ring of the codomain does not coincide with the domain of the pullback")
    return new{BRT, BRET, RT, RET, MST1, MST2}(X, Y, pullback)
  end
end

### getter functions
pullback(phi::SpecHom) = phi.pullback
domain(phi::SpecHom) = phi.domain
codomain(phi::SpecHom) = phi.codomain

### functionality
function preimage(
    phi::SpecHom{BRT, BRET, RT, RET, MST1, MST2}, 
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
  return Spec(MPolyQuoLocalizedRing(R, ideal(R, new_gens), inverted_set(OO(X))*new_units))
end


