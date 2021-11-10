import AbstractAlgebra.Ring, Oscar.AlgHom, Oscar.compose, AbstractAlgebra.Generic.Frac
import Base: ∘
import Oscar: base_ring
import AbstractAlgebra: FPModule, FPModuleElem
import Base.copy

using Oscar.Multiindices
using Oscar.Misc

export AffineScheme, Spec, SpecPrincipalOpen
export base_ring, ambient_ring, defining_ideal, imgs_frac, pullback, pullback_from_parent, pullback_from_root, inclusion_in_parent, inclusion_in_root, set_name!, inverted_element, identity_morphism, denoms, inverses, divide_by_units
export affine_space, localize, subscheme

export AffSchMorphism
export domain, codomain, pullback

@Markdown.doc """
Scheme{ BaseRingType <: Ring }

A scheme over a ring ``k`` of type `BaseRingType`.
"""
abstract type Scheme{ S <: Ring } end

@doc Markdown.doc"""
    AffineScheme{BaseRingType, RingType <: MPolyRing, RingElemType <: MPolyElem} <: Scheme{BaseRingType}

An affine scheme over a base ring ``k`` of type `BaseRingType`, given by a ring ``R/I`` with 
``R`` a polynomial ring of type `RingType` and elements of type `RingElemType`.
"""
abstract type AffineScheme{BaseRingType<:Ring, BaseRingElemType<:RingElement, RingType<:MPolyRing, RingElemType<:MPolyElem} <: Scheme{BaseRingType} end

@Markdown.doc """
SchemeMorphism{BaseRingType<:Ring}

Morphism of Schemes over a ring ``k`` of type `BaseRingType`.
"""
abstract type SchemeMorphism{BaseRingType<:Ring} end

@doc Markdown.doc"""
Spec{BaseRingType, RingType, RingElemType} <: AffineScheme{BaseRingType, RingType, RingElemType}

An affine scheme ``X = Spec R/I`` with ``R = k[x₁,…,xₙ]`` a free 
polynomial algebra of type `RingType` over a base ring ``k`` of type 
`BaseRingType` and ``I ⊂ R`` a finitely generated ideal 
with elements of type ``RingElemType``.
"""
mutable struct Spec{BRT, BRET, RT, RET} <: AffineScheme{BRT, BRET, RT, RET}
  # the basic fields 
  OO::MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}}
  # fields for caching
  name::String # the name of this scheme for printing

  function Spec(OO::MPolyQuoLocalizedRing{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST<:AbsMultSet{RT, RET}}
    return new{BRT, BRET, RT, RET}(OO)
  end
end


### Getter functions

OO(X::Spec) = X.OO

function name(X::Spec) 
  if isdefined(X, :name)
    return X.name
  end
  return "unnamed affine variety"
end

function set_name!(X::Spec, name::String) 
  X.name = name
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
function subscheme(X::Spec{BRT, BRET, RT, RET}, I::MPolyLocalizedIdeal{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST}
  OO(X) == base_ring(I) || error("ideal does not live in the correct ring")
  return Spec(quo(OO(X), I))
end

function subscheme(X::Spec{BRT, BRET, RT, RET}, I::MPolyIdeal{RET}) where {BRT, BRET, RT, RET}
  base_ring(OO(X)) == base_ring(I) || error("ideal does not live in the correct ring")
  return Spec(quo(OO(X), I))
end
  
function subscheme(X::Spec{BRT, BRET, RT, RET}, f::RET) where {BRT, BRET, RT, RET}
  R = base_ring(OO(X))
  I = ideal(R, f)
  return subscheme(X, I)
end

function subscheme(X::Spec{BRT, BRET, RT, RET}, f::Vector{RET}) where {BRT, BRET, RT, RET}
  R = base_ring(OO(X))
  I = ideal(R, f)
  return subscheme(X, I)
end

function subscheme(X::Spec{BRT, BRET, RT, RET}, f::BRET) where {BRT, BRET, RT, RET}
  R = base_ring(OO(X))
  I = ideal(R, R(f))
  return subscheme(X, I)
end

function subscheme(X::Spec{BRT, BRET, RT, RET}, f::Vector{BRET}) where {BRT, BRET, RT, RET}
  R = base_ring(OO(X))
  I = ideal(R, R.(f))
  return subscheme(X, I)
end

### open subschemes defined by hypersurfaces
function hypersurface_complement(X::Spec{BRT, BRET, RT, RET}, f::RET) where {BRT, BRET, RT, RET}
  R = base_ring(OO(X))
  parent(f) == R || error("the element does not belong to the correct ring")
  return Spec(Localization(OO(X), MPolyPowersOfElement(f)))
end
