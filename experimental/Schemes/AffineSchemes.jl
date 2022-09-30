import AbstractAlgebra.Ring
import Base: intersect

export Scheme, AbsSpec
export Spec, OO, defining_ideal, ambient_ring
export spec_type, ring_type
export base_ring_type, base_ring_elem_type, poly_type, poly_ring_type, mult_set_type, ring_type
export affine_space, empty_spec
export EmptyScheme

export is_open_embedding, is_closed_embedding, is_canonically_isomorphic, hypersurface_complement, subscheme, name_of, set_name!
export closure, product

export SpecMor, morphism_type
export pullback, domain, codomain, preimage, restrict, graph, identity_map, inclusion_map, is_isomorphism, is_inverse_of, is_identity_map, lift_map

export strict_modulus

export simplify

@Markdown.doc """
    Scheme{BaseRingType<:Ring} 

A scheme over a ring ``𝕜`` of type `BaseRingType`.
"""
abstract type Scheme{BaseRingType} end

@Markdown.doc """
    SchemeMor{DomainType, CodomainType, MorphismType, BaseMorType}

A morphism of schemes ``f : X → Y`` of type `MorphismType` with 
``X`` of type `DomainType` and ``Y`` of type `CodomainType`. 

When ``X`` and ``Y`` are defined over schemes ``BX`` and ``BY`` other 
than ``Spec(𝕜)``, `BaseMorType` is the type of the underlying 
morphism ``BX → BY``; otherwise, it can be set to `Nothing`.
"""
abstract type SchemeMor{
                        DomainType, 
                        CodomainType, 
                        MorphismType,
                        BaseMorType
                       } <: Hecke.Map{
                                      DomainType, 
                                      CodomainType, 
                                      SetMap, 
                                      MorphismType
                                     } 
end

struct EmptyScheme{BaseRingType}<:Scheme{BaseRingType} 
  k::BaseRingType
  function EmptyScheme(k::BaseRingType) where {BaseRingType<:Ring}
    return new{BaseRingType}(k)
  end
end

is_empty(X::EmptyScheme) = true

########################################################################
#
# Interface for affine schemes
#
########################################################################

@Markdown.doc """
    AbsSpec{BaseRingType, RingType<:Ring}

An affine scheme ``X = Spec(R)`` with ``R`` of type `RingType` over 
a ring ``𝕜`` of type `BaseRingType`.
"""
abstract type AbsSpec{BaseRingType, RingType<:Ring} <: Scheme{BaseRingType} end

### essential getter methods

@Markdown.doc """
    OO(X::AbsSpec) 

On an affine scheme ``X = Spec(R)`` this returns the ring `R`.
"""
function OO(X::AbsSpec{BRT, RT}) where {BRT, RT} 
  OO(underlying_scheme(X))::RT
end

@Markdown.doc """
    ambient_ring(X::AbsSpec)

On an affine scheme ``X = Spec(R)`` over ``𝕜`` this returns a 
polynomial ring ``P = 𝕜[x₁,…,xₙ]`` with natural coercion 
``P → R`` and the property that for every other (commutative) 
ring ``S`` and any homomorphism ``φ : R → S`` there is a morphism 
``ψ : P → S`` factoring through ``φ`` and such that ``φ`` 
is uniquely determined by ``ψ``.
"""
function ambient_ring(X::AbsSpec)
  return ambient_ring(underlying_scheme(X))::MPolyRing
end

### type getters
ring_type(::Type{SpecType}) where {BRT, RT, SpecType<:AbsSpec{BRT, RT}} = RT
ring_type(X::AbsSpec) = ring_type(typeof(X))

base_ring_type(::Type{SpecType}) where {BRT, RT, SpecType<:AbsSpec{BRT, RT}} = BRT
base_ring_type(X::AbsSpec) = base_ring_type(typeof(X))
base_ring_elem_type(::Type{SpecType}) where {BRT, RT, SpecType<:AbsSpec{BRT, RT}} = elem_type(BRT)
base_ring_elem_type(X::AbsSpec) = base_ring_elem_type(typeof(X))

### generically derived getters
@Markdown.doc """
    base_ring(X::AbsSpec) 

On an affine scheme ``X/𝕜`` over ``𝕜`` this returns the ring `𝕜`.
"""
function base_ring(X::AbsSpec{BRT, RT}) where {BRT, RT}
 return base_ring(underlying_scheme(X))::BRT
end

### constructors
#
# Note that these default to the plain type Spec and must be overwritten 
# if something more sophisticated should be returned!

@Markdown.doc """
    subscheme(X::AbsSpec, I::Ideal)

For a scheme ``X = Spec(R)`` and an ideal ``I ⊂ 𝒪(X)`` 
this returns the closed subscheme defined by ``I``.
"""
function subscheme(X::AbsSpec, I::Ideal)
  base_ring(I) == OO(X) || return subscheme(X, ideal(OO(X), OO(X).(gens(I)))) # this will throw if coercion is not possible
  return Spec(quo(OO(X), I)[1])
end
subscheme(X::AbsSpec, f::RingElem) = subscheme(X, ideal(OO(X), [f]))
subscheme(X::AbsSpec, f::Vector{<:RingElem}) = subscheme(X, ideal(OO(X), f))

### open subschemes defined by complements of hypersurfaces
@Markdown.doc """
    hypersurface_complement(X::AbsSpec, f::RingElem)

For a scheme ``X = Spec(R)`` and an element ``f ∈ R`` 
this returns the open subscheme ``U = Spec(R[f⁻¹]) = X ∖ V(f)`` 
defined by the complement of the vanishing 
locus of ``f``.
"""
function hypersurface_complement(X::AbsSpec, f::RingElem)
  return hypersurface_complement(underlying_scheme(X), f)::AbsSpec
end

function hypersurface_complement(X::SpecType, f::RingElem) where {SpecType<:AbsSpec{<:Any, <:MPolyQuoLocalizedRing}}
  parent(f) == OO(X) || return hypersurface_complement(X, OO(X)(f))
  h = lifted_numerator(f)
  U = MPolyPowersOfElement(h)
  W, _ = Localization(OO(X), U)
  return Spec(W)
end

function hypersurface_complement(X::SpecType, f::RingElem) where {SpecType<:AbsSpec{<:Any, <:MPolyLocalizedRing}}
  parent(f) == OO(X) || return hypersurface_complement(X, OO(X)(f))
  h = numerator(f)
  U = MPolyPowersOfElement(h)
  W, _ = Localization(OO(X), U)
  return Spec(W)
end

function hypersurface_complement(X::SpecType, f::RingElem) where {SpecType<:AbsSpec{<:Any, <:MPolyRing}}
  parent(f) == OO(X) || return hypersurface_complement(X, OO(X)(f))
  U = MPolyPowersOfElement(f)
  W, _ = Localization(OO(X), U)
  return Spec(W)
end

function hypersurface_complement(X::SpecType, f::RingElem) where {SpecType<:AbsSpec{<:Any, <:MPolyQuo}}
  parent(f) == OO(X) || return hypersurface_complement(X, OO(X)(f))
  U = MPolyPowersOfElement(lift(f))
  W, _ = Localization(OO(X), U)
  return Spec(W)
end

function hypersurface_complement(X::SpecType, f::Vector{<:RingElem}) where {SpecType<:AbsSpec{<:Any, <:MPolyQuoLocalizedRing}}
  all(x->(parent(x) == OO(X)), f) || return hypersurface_complement(X, OO(X).(f))
  h = lifted_numerator.(f)
  U = MPolyPowersOfElement(ambient_ring(X), h)
  W, _ = Localization(OO(X), U)
  return Spec(W)
end

function hypersurface_complement(X::SpecType, f::Vector{<:RingElem}) where {SpecType<:AbsSpec{<:Any, <:MPolyLocalizedRing}}
  all(x->(parent(x) == OO(X)), f) || return hypersurface_complement(X, OO(X).(f))
  h = numerator.(f)
  U = MPolyPowersOfElement(ambient_ring(X), h)
  W, _ = Localization(OO(X), U)
  return Spec(W)
end

function hypersurface_complement(X::SpecType, f::Vector{<:RingElem}) where {SpecType<:AbsSpec{<:Any, <:MPolyRing}}
  all(x->(parent(x) == OO(X)), f) || return hypersurface_complement(X, OO(X).(f))
  U = MPolyPowersOfElement(ambient_ring(X), f)
  W, _ = Localization(OO(X), U)
  return Spec(W)
end

function hypersurface_complement(X::SpecType, f::Vector{<:RingElem}) where {SpecType<:AbsSpec{<:Any, <:MPolyQuo}}
  all(x->(parent(x) == OO(X)), f) || return hypersurface_complement(X, OO(X).(f))
  U = MPolyPowersOfElement(ambient_ring(X), lift.(f))
  W, _ = Localization(OO(X), U)
  return Spec(W)
end

# Shortcut notation for the rest of the file
StdSpec = AbsSpec{<:Ring, <:MPolyQuoLocalizedRing{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement}}

@Markdown.doc """
    Spec{BaseRingType, RingType}

An affine scheme ``X = Spec(R)`` with ``R`` a Noetherian ring of type `RingType`
over a base ring ``𝕜`` of type `BaseRingType`.
"""
@attributes mutable struct Spec{BaseRingType, RingType} <: AbsSpec{BaseRingType, RingType}
  # the basic fields 
  OO::RingType
  kk::BaseRingType

  function Spec(OO::MPolyQuoLocalizedRing) 
    kk = coefficient_ring(base_ring(OO))
    return new{typeof(kk), typeof(OO)}(OO, kk)
  end
  function Spec(OO::MPolyLocalizedRing) 
    kk = coefficient_ring(base_ring(OO))
    return new{typeof(kk), typeof(OO)}(OO, kk)
  end
  function Spec(OO::MPolyRing) 
    kk = coefficient_ring(OO)
    return new{typeof(kk), typeof(OO)}(OO, kk)
  end
  function Spec(OO::MPolyQuo) 
    kk = coefficient_ring(base_ring(OO))
    return new{typeof(kk), typeof(OO)}(OO, kk)
  end

  function Spec(R::Ring)
    return new{typeof(ZZ), typeof(R)}(R, ZZ)
  end

  function Spec(kk::Ring, R::Ring)
    return new{typeof(kk), typeof(R)}(R, kk)
  end

  function Spec(kk::Field)
    return new{typeof(kk), typeof(kk)}(kk, kk)
  end
end

### Type getters

ring_type(::Type{Spec{BRT, RT}}) where {BRT, RT} = RT
ring_type(X::Spec) = ring_type(typeof(X))
base_ring_type(::Type{Spec{BRT, RT}}) where {BRT, RT} = BRT
base_ring_type(X::Spec) = base_ring_type(typeof(X))

### type constructors

### Getter functions

OO(X::Spec) = X.OO
base_ring(X::Spec) = X.kk

### constructors
standard_spec(X::AbsSpec{<:Any, <:MPolyRing}) = Spec(MPolyQuoLocalizedRing(OO(X), ideal(OO(X), [zero(OO(X))]), units_of(OO(X))))

function standard_spec(X::AbsSpec{<:Any, <:MPolyQuo}) 
  A = OO(X)
  R = base_ring(A)
  return Spec(MPolyQuoLocalizedRing(R, modulus(A), units_of(R)))
end

standard_spec(X::AbsSpec{<:Any, <:MPolyLocalizedRing}) = Spec(MPolyQuoLocalizedRing(base_ring(OO(X)), ideal(base_ring(OO(X)), [zero(base_ring(OO(X)))]), inverted_set(OO(X))))
standard_spec(X::AbsSpec{<:Any, <:MPolyQuoLocalizedRing}) = Spec(OO(X))

ambient_ring(X::Spec{<:Any, <:MPolyRing}) = OO(X)
ambient_ring(X::Spec{<:Any, <:MPolyQuo}) = base_ring(OO(X))
ambient_ring(X::Spec{<:Any, <:MPolyLocalizedRing}) = base_ring(OO(X))
ambient_ring(X::Spec{<:Any, <:MPolyQuoLocalizedRing}) = base_ring(OO(X))
ambient_ring(X::Spec{T, T}) where {T<:Field} = base_ring(X)



@attr String function name(X::Spec)
  return "unnamed affine variety"
end

function set_name!(X::Spec, name::String) 
  return set_attribute!(X, :name, name)
end

function Base.show(io::IO, X::Spec) 
  if isdefined(X, :name)
    print(io, name_of(X))
    return
  end
  print(io, "Spec of $(OO(X))")
end

@attr defining_ideal(X::Spec{<:Any, <:MPolyRing}) = ideal(OO(X), [zero(OO(X))])
defining_ideal(X::Spec{<:Any, <:MPolyQuo}) = modulus(OO(X))
@attr defining_ideal(X::Spec{<:Any, <:MPolyLocalizedRing}) = ideal(OO(X), [zero(OO(X))])
defining_ideal(X::Spec{<:Any, <:MPolyQuoLocalizedRing}) = modulus(OO(X))

### Copy constructor
Spec(X::Spec) = Spec(OO(X))

### internal copy routines
Base.deepcopy_internal(X::Spec, dict::IdDict) = Spec(deepcopy_internal(OO(X), dict))

### additional constructors 
Spec(R::MPolyRing, I::MPolyIdeal, U::AbsMPolyMultSet) = Spec(MPolyQuoLocalizedRing(R, I, U))
Spec(R::MPolyRing, U::AbsMPolyMultSet) = Spec(Localization(R, U)[1])
Spec(R::MPolyRing, I::MPolyIdeal) = Spec(quo(R, I)[1])

### closed subschemes defined by elements
function subscheme(X::Spec, f::Vector{<:RingElem})
  all(x->(parent(x) == OO(X)), f) || return subscheme(X, OO(X).(f))
  return subscheme(X, ideal(OO(X), f))
end

### testing containment
issubset(X::EmptyScheme{BRT}, Y::Scheme{BRT}) where {BRT} = true

function issubset(Y::AbsSpec{BRT, <:Any}, X::EmptyScheme{BRT}) where {BRT}
  return iszero(one(OO(Y)))
end

@Markdown.doc """
    issubset(X::AbsSpec, Y::AbsSpec)

Checks whether ``X`` is a subset of ``Y`` based on the comparison of their coordinate rings.
"""
function issubset(X::AbsSpec, Y::AbsSpec)
  error("method `issubset(X, Y)` not implemented for `X` of type $(typeof(X)) and `Y` of type $(typeof(Y))")
end

function issubset(
    X::AbsSpec{BRT, RT}, 
    Y::AbsSpec{BRT, RT}
  ) where {BRT, RT<:MPolyRing}
  return OO(X) === OO(Y)
end

function issubset(
    X::AbsSpec{BRT, <:MPolyRing}, 
    Y::AbsSpec{BRT, <:MPolyQuo}
  ) where {BRT}
  R = OO(X)
  R == ambient_ring(Y) || error("schemes can not be compared")
  return iszero(modulus(OO(Y)))
end

function issubset(
    X::AbsSpec{BRT, <:MPolyRing},
    Y::AbsSpec{BRT, <:MPolyLocalizedRing}
  ) where {BRT}
  R = OO(X)
  R == base_ring(OO(Y)) || error("schemes can not be compared")
  return issubset(inverted_set(OO(Y)), units_of(R))
end

function issubset(
    X::AbsSpec{BRT, <:MPolyRing},
    Y::AbsSpec{BRT, <:MPolyQuoLocalizedRing}
  ) where {BRT}
  R = OO(X)
  R == base_ring(OO(Y)) || error("schemes can not be compared")
  return issubset(inverted_set(OO(Y)), units_of(R)) && iszero(defining_ideal(Y))
end

########################################################################
# MPolyQuo in first argument
########################################################################

function issubset(
    X::AbsSpec{BRT, <:MPolyQuo},
    Y::AbsSpec{BRT, <:MPolyRing}
  ) where {BRT}
  R = ambient_ring(Y)
  R == ambient_ring(X) || error("schemes can not be compared")
  return true
end

function issubset(
    X::AbsSpec{BRT, RT},
    Y::AbsSpec{BRT, RT}
  ) where {BRT, RT<:MPolyQuo}
  R = base_ring(OO(X))
  R == base_ring(OO(Y)) || error("schemes can not be compared")
  return issubset(defining_ideal(Y), defining_ideal(X))
end

function issubset(
    X::AbsSpec{BRT, <:MPolyQuo},
    Y::AbsSpec{BRT, <:MPolyLocalizedRing{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement}}
  ) where {BRT}
  R = base_ring(OO(X))
  R == base_ring(OO(Y)) || error("schemes can not be compared")
  all(x->isunit(OO(X)(x)), denominators(inverted_set(OO(Y)))) || return false
  return iszero(localized_ring(OO(Y))(modulus(OO(X))))
end

function issubset(
    X::AbsSpec{BRT, <:MPolyQuo},
    Y::AbsSpec{BRT, <:MPolyQuoLocalizedRing{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement}}
  ) where {BRT}
  R = base_ring(OO(X))
  R == base_ring(OO(Y)) || error("schemes can not be compared")
  all(x->isunit(OO(X)(x)), denominators(inverted_set(OO(Y)))) || return false
  return issubset(localized_modulus(OO(Y)), localized_ring(OO(Y))(modulus(OO(X))))
end

########################################################################
# MPolyLocalizedRing in first argument
########################################################################

function issubset(
    X::AbsSpec{BRT, <:MPolyLocalizedRing},
    Y::AbsSpec{BRT, <:MPolyRing}
  ) where {BRT}
  R = OO(Y)
  R == base_ring(OO(X)) || error("schemes can not be compared")
  return true
end

function issubset(
    X::AbsSpec{BRT, <:MPolyLocalizedRing},
    Y::AbsSpec{BRT, <:MPolyQuo}
  ) where {BRT}
  R = ambient_ring(Y)
  R == base_ring(OO(X)) || error("schemes can not be compared")
  return all(x->(iszero(OO(X)(x))), gens(modulus(OO(Y))))
end

function issubset(
    X::AbsSpec{BRT, RT},
    Y::AbsSpec{BRT, RT}
  ) where {BRT, RT<:MPolyLocalizedRing}
  R = base_ring(OO(X))
  R == base_ring(OO(Y)) || error("schemes can not be compared")
  UX = inverted_set(OO(X))
  UY = inverted_set(OO(Y))
  return issubset(UY, UX)
end

function issubset(
    X::AbsSpec{BRT, <:MPolyLocalizedRing},
    Y::AbsSpec{BRT, <:MPolyQuoLocalizedRing}
  ) where {BRT}
  R = base_ring(OO(X))
  R == base_ring(OO(Y)) || error("schemes can not be compared")
  UX = inverted_set(OO(X))
  UY = inverted_set(OO(Y))
  return issubset(UY, UX) && iszero(localized_modulus(OO(Y)))
end

########################################################################
# MPolyQuoLocalizedRing in first argument
########################################################################
function issubset(
    X::AbsSpec{BRT, <:MPolyQuoLocalizedRing},
    Y::AbsSpec{BRT, <:MPolyRing}
  ) where {BRT}
  R = OO(Y)
  R == base_ring(OO(X)) || error("schemes can not be compared")
  return true
end

function issubset(
    X::AbsSpec{BRT, <:MPolyQuoLocalizedRing},
    Y::AbsSpec{BRT, <:MPolyQuo}
  ) where {BRT}
  R = base_ring(OO(Y))
  R == base_ring(OO(X)) || error("schemes can not be compared")
  L = localized_ring(OO(X))
  return issubset(L(modulus(OO(Y))), localized_modulus(OO(X)))
end

function issubset(
    X::AbsSpec{BRT, <:MPolyQuoLocalizedRing},
    Y::AbsSpec{BRT, <:MPolyLocalizedRing{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement}}
  ) where {BRT}
  R = base_ring(OO(X))
  R == base_ring(OO(Y)) || error("schemes can not be compared")
  UX = inverted_set(OO(X))
  UY = inverted_set(OO(Y))
  if !issubset(UY, UX)
    # check whether the inverted elements in Y are units anyway
    for a in denominators(UY)
      is_unit(OO(X)(a)) || return false
    end
  end
  return true
end

function issubset(
    X::AbsSpec{BRT, RT},
    Y::AbsSpec{BRT, RT}
  ) where {BRT, RT<:MPolyQuoLocalizedRing{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement}}
  R = base_ring(OO(X))
  R == base_ring(OO(Y)) || error("schemes can not be compared")
  UX = inverted_set(OO(X))
  UY = inverted_set(OO(Y))
  if !issubset(UY, UX)
    # check whether the inverted elements in Y are units anyway
    for a in denominators(UY)
      is_unit(OO(X)(a)) || return false
    end
  end
  J = localized_ring(OO(X))(modulus(OO(Y)))
  return issubset(J, localized_modulus(OO(X)))
end

# TODO: Add further cross-type comparison methods as needed.

function ==(X::T, Y::T) where {T<:Spec}
  return X === Y
end

function is_canonically_isomorphic(X::AbsSpec, Y::AbsSpec)
  X === Y && return true
  return issubset(X, Y) && issubset(Y, X)
end

function is_canonically_isomorphic(X::AbsSpec, Y::EmptyScheme)
  return issubset(X, Y)
end

is_canonically_isomorphic(X::EmptyScheme, Y::Spec) = is_canonically_isomorphic(Y, X)

Base.isempty(X::AbsSpec) = iszero(one(OO(X)))

@Markdown.doc """
    is_open_embedding(X::AbsSpec, Y::AbsSpec)

Checks whether ``X`` is openly embedded in ``Y``.
"""
function is_open_embedding(X::AbsSpec, Y::AbsSpec)
  return is_open_embedding(standard_spec(X), standard_spec(Y))
end

function is_open_embedding(
    X::Spec{BRT, RT},
    Y::Spec{BRT, RT}
  ) where {BRT, RT<:MPolyQuoLocalizedRing{<:Any, <:Any, <:Any, <:Any,
                                          <:MPolyPowersOfElement}}
  R = base_ring(OO(X))
  R == base_ring(OO(Y)) || return false
  UX = inverted_set(OO(X))
  UY = inverted_set(OO(Y))
  issubset(UY, UX) || return false
  J = localized_ring(OO(X))(modulus(OO(Y)))
  return localized_modulus(OO(X)) == J
end

function is_open_embedding(
    X::Spec{BRT, <:MPolyQuoLocalizedRing},
    Y::Spec{BRT, <:MPolyRing}
  ) where {BRT}
  return OO(Y) == base_ring(OO(X)) && all(iszero, gens(modulus(OO(X))))
end

#TODO: Add more cross-type methods as needed.

@Markdown.doc """
    is_closed_embedding(X::AbsSpec, Y::AbsSpec)

Checks whether ``X`` is closed embedded in ``Y``.
"""
function is_closed_embedding(X::AbsSpec, Y::AbsSpec)
  error("`is_closed_embedding(X, Y)` not implemented for X of type $(typeof(X)) and Y of type $(typeof(Y))")
end

function is_closed_embedding(
    X::Spec{BRT, RT},
    Y::Spec{BRT, RT}
  ) where {BRT, RT<:MPolyQuoLocalizedRing{<:Any, <:Any, <:Any, <:Any,
                                          <:MPolyPowersOfElement}}
  R = base_ring(OO(X))
  R == base_ring(OO(Y)) || return false
  inverted_set(OO(X)) == inverted_set(OO(Y)) || return false
  J = localized_ring(OO(X))(modulus(OO(Y)))
  return issubset(J, localized_modulus(OO(X)))
end

function is_closed_embedding(
    X::Spec{BRT, <:MPolyQuo},
    Y::Spec{BRT, <:MPolyRing}
  ) where {BRT}
  return OO(Y) == base_ring(OO(X))
end

function is_closed_embedding(
    X::Spec{BRT, <:MPolyQuo},
    Y::Spec{BRT, <:RT}
  ) where {BRT, RT<:MPolyQuoLocalizedRing{<:Any, <:Any, <:Any, <:Any,
                                          <:MPolyPowersOfElement}}
  R = base_ring(OO(X))
  R == base_ring(OO(Y)) || error("schemes can not be compared")
  all(x->(isunit(OO(X)(x))), denominators(inverted_set(OO(Y)))) || return false
  return issubset(localized_modulus(OO(Y)), localized_ring(OO(Y))(modulus(OO(X))))
end

#TODO: Add more cross-type methods as needed.

### set operations

Base.intersect(E::EmptyScheme{BRT}, X::Scheme{BRT}) where {BRT} = E
Base.intersect(X::Scheme{BRT}, E::EmptyScheme{BRT}) where {BRT} = E
Base.intersect(X::EmptyScheme{BRT}, E::EmptyScheme{BRT}) where {BRT} = E

### For Specs of MPolyRings
# TODO  intersect X,Y for X<Y should return a copy of X with === ambient_rings
# Spec(X) does not apply for instance to principal open subsets hence a change
# is necessary
function Base.intersect(
    X::AbsSpec{BRT, <:MPolyRing},
    Y::AbsSpec{BRT, <:MPolyRing}
  ) where {BRT}
  R = OO(X)
  R == OO(Y) || error("schemes can not be compared")
  return Spec(X)
end

function Base.intersect(
    X::AbsSpec{BRT, <:MPolyRing},
    Y::AbsSpec{BRT, <:MPolyQuo}
  ) where {BRT}
  R = OO(X)
  R == base_ring(OO(Y)) || error("schemes can not be compared")
  return Spec(Y)
end

function Base.intersect(
    X::AbsSpec{BRT, <:MPolyRing},
    Y::AbsSpec{BRT, <:MPolyLocalizedRing}
  ) where {BRT}
  R = OO(X)
  R == base_ring(OO(Y)) || error("schemes can not be compared")
  return Spec(Y)
end

function Base.intersect(
    X::AbsSpec{BRT, <:MPolyRing},
    Y::AbsSpec{BRT, <:MPolyQuoLocalizedRing}
  ) where {BRT}
  @show "hi"
  R = OO(X)
  R == base_ring(OO(Y)) || error("schemes can not be compared")
  return Spec(Y)
end

function Base.intersect(
    Y::AbsSpec{BRT, <:Any},
    X::AbsSpec{BRT, <:MPolyRing}
  ) where {BRT}
  return intersect(X, Y)
end

### For Specs of MPolyQuos
function Base.intersect(
    X::AbsSpec{BRT, <:MPolyQuo},
    Y::AbsSpec{BRT, <:MPolyQuo}
  ) where {BRT}
  R = base_ring(OO(X))
  R == base_ring(OO(Y)) || error("schemes can not be compared")
  return Spec(quo(R, modulus(OO(X)) + modulus(OO(Y)))[1])
end

function Base.intersect(
    X::AbsSpec{BRT, <:MPolyQuo},
    Y::AbsSpec{BRT, <:MPolyLocalizedRing}
  ) where {BRT}
  R = base_ring(OO(X))
  R == base_ring(OO(Y)) || error("schemes can not be compared")
  return Spec(quo(OO(Y), OO(Y)(modulus(OO(X))))[1])
end

function Base.intersect(
    X::AbsSpec{BRT, <:MPolyQuo},
    Y::AbsSpec{BRT, <:MPolyQuoLocalizedRing}
  ) where {BRT}
  R = base_ring(OO(X))
  R == base_ring(OO(Y)) || error("schemes can not be compared")
  return Spec(quo(OO(Y), OO(Y)(modulus(OO(X))))[1])
end

function Base.intersect(
    Y::AbsSpec{BRT, <:Any},
    X::AbsSpec{BRT, <:MPolyQuo}
  ) where {BRT}
  return intersect(X, Y)
end

### For Specs of MPolyLocalizedRings
function Base.intersect(
    X::AbsSpec{BRT, <:MPolyLocalizedRing},
    Y::AbsSpec{BRT, <:MPolyLocalizedRing}
  ) where {BRT}
  R = base_ring(OO(X))
  R == base_ring(OO(Y)) || error("schemes can not be compared")
  return Spec(Localization(R, inverted_set(OO(X)) * inverted_set(OO(Y)))[1])
end

function Base.intersect(
    X::AbsSpec{BRT, <:MPolyLocalizedRing},
    Y::AbsSpec{BRT, <:MPolyQuoLocalizedRing}
  ) where {BRT}
  R = base_ring(OO(X))
  R == base_ring(OO(Y)) || error("schemes can not be compared")
  return Spec(R, modulus(OO(Y)), inverted_set(OO(X))*inverted_set(OO(Y)))
end

function Base.intersect(
    Y::AbsSpec{BRT, <:Any},
    X::AbsSpec{BRT, <:MPolyLocalizedRing}
  ) where {BRT}
  return intersect(X, Y)
end

### For Specs of MPolyQuoLocalizedRings
function Base.intersect(
    X::AbsSpec{BRT, <:MPolyQuoLocalizedRing},
    Y::AbsSpec{BRT, <:MPolyQuoLocalizedRing}
  ) where {BRT}
  R = base_ring(OO(X))
  R == base_ring(OO(Y)) || error("schemes can not be compared")
  Q, _ = quo(R, modulus(OO(X)) + modulus(OO(Y)))
  return Spec(R, modulus(OO(X)) + modulus(OO(Y)), 
              inverted_set(OO(X)) * inverted_set(OO(Y)))
end

#TODO: Add more cross-type methods as needed.

### compute the closure of X in Y
@Markdown.doc """
    closure(X::AbsSpec, Y::AbsSpec) 

Returns the closure of ``X`` in ``Y``.
"""
function closure(X::AbsSpec, Y::AbsSpec) 
  return closure(standard_spec(X), standard_spec(Y))
end

function closure(
    X::Spec{BRT, RT}, 
    Y::Spec{BRT, RT}
  ) where {BRT, RT<:MPolyQuoLocalizedRing{<:Any, <:Any, <:Any, <:Any, 
                                          <:MPolyPowersOfElement}}
  issubset(X, Y) || error("the first argument is not a subset of the second")
  is_closed_embedding(X, Y) && return X
  W, _ = Localization(inverted_set(OO(X))*inverted_set(OO(Y)))
  I = ideal(W, W.(gens(modulus(OO(X)))))
  Isat = saturated_ideal(I)
  R = base_ring(OO(Y))
  return Spec(MPolyQuoLocalizedRing(R, Isat, inverted_set(OO(Y))))
end

#TODO: Add more cross-type methods as needed.

########################################################################
# Morphisms of affine schemes                                      #
########################################################################

@Markdown.doc """
    AbsSpecMor{DomainType<:AbsSpec, 
               CodomainType<:AbsSpec, 
               PullbackType<:Hecke.Map,
               MorphismType, 
               BaseMorType
               }

Abstract type for morphisms ``f : X → Y`` of affine schemes where

  * ``X = Spec(S)`` is of type `DomainType`, 
  * ``Y = Spec(R)`` is of type `CodomainType`, 
  * ``f^* : R → S`` is a ring homomorphism of type `PullbackType`, 
  * ``f`` itself is of type `MorphismType` (required for the Map interface),
  * if ``f`` is defined over a morphism of base schemes ``BX → BY`` 
    (e.g. a field extension), then this base scheme morphism is of 
    type `BaseMorType`; otherwise, this can be set to `Nothing`.
"""
abstract type AbsSpecMor{
                         DomainType<:AbsSpec, 
                         CodomainType<:AbsSpec, 
                         PullbackType<:Hecke.Map,
                         MorphismType, 
                         BaseMorType
                        }<:SchemeMor{DomainType, CodomainType, MorphismType, BaseMorType}
end

underlying_morphism(f::AbsSpecMor) = error("`underlying_morphism(f)` not implemented for `f` of type $(typeof(f))")
@Markdown.doc """
    domain(f::AbsSpecMor)

On a morphism ``f : X → Y`` of affine schemes, this returns ``X``.
"""
domain(f::AbsSpecMor) = domain(underlying_morphism(f))

@Markdown.doc """
    codomain(f::AbsSpecMor)

On a morphism ``f : X → Y`` of affine schemes, this returns ``Y``.
"""
codomain(f::AbsSpecMor) = domain(underlying_morphism(f))

@Markdown.doc """
    pullback(f::AbsSpecMor)

On a morphism ``f : X → Y`` of affine schemes ``X = Spec(S)`` and 
``Y = Spec(R)``, this returns the ring homomorphism ``f^* : R → S``.
"""
pullback(f::AbsSpecMor) = pullback(underlying_morphism(f))

### Type getters
pullback_type(::Type{T}) where {DomType, CodType, PbType, T<:AbsSpecMor{DomType, CodType, PbType}} = PbType
pullback_type(f::AbsSpecMor) = pullback_type(typeof(f))
domain_type(::Type{T}) where {DomType, CodType, PbType, T<:AbsSpecMor{DomType, CodType, PbType}} = DomType
domain_type(f::AbsSpecMor) = domain_type(typeof(f))
codomain_type(::Type{T}) where {DomType, CodType, PbType, T<:AbsSpecMor{DomType, CodType, PbType}} = CodType
codomain_type(f::AbsSpecMor) = codomain_type(typeof(f))


@Markdown.doc """
    SpecMor{DomainType<:AbsSpec, 
            CodomainType<:AbsSpec, 
            PullbackType<:Hecke.Map
           }

A morphism ``f : X → Y`` of affine schemes ``X = Spec(S)`` of type 
`DomainType` and ``Y = Spec(R)`` of type `CodomainType`, both defined 
over the same `base_ring`, with underlying ring homomorphism 
``f^* : R → S`` of type `PullbackType`.
"""
@attributes mutable struct SpecMor{
                                   DomainType<:AbsSpec, 
                                   CodomainType<:AbsSpec, 
                                   PullbackType<:Hecke.Map
                                  } <: AbsSpecMor{DomainType, 
                                                  CodomainType, 
                                                  PullbackType, 
                                                  SpecMor, 
                                                  Nothing
                                                 }
  domain::DomainType
  codomain::CodomainType
  pullback::PullbackType

  function SpecMor(
      X::DomainType,
      Y::CodomainType,
      pullback::PullbackType;
      check::Bool=true
    ) where {DomainType<:AbsSpec, CodomainType<:AbsSpec, PullbackType<:Hecke.Map}
    OO(X) == codomain(pullback) || error("the coordinate ring of the domain does not coincide with the codomain of the pullback")
    OO(Y) == domain(pullback) || error("the coordinate ring of the codomain does not coincide with the domain of the pullback")
    if check
      # do some more expensive tests
    end
    return new{DomainType, CodomainType, PullbackType}(X, Y, pullback)
  end
end

function morphism_type(::Type{SpecType1}, ::Type{SpecType2}) where {SpecType1<:AbsSpec, SpecType2<:AbsSpec}
  return SpecMor{SpecType1, SpecType2, morphism_type(ring_type(SpecType2), ring_type(SpecType1))}
end

morphism_type(X::Spec, Y::Spec) = morphism_type(typeof(X), typeof(Y))
morphism_type(X::AbsSpec, Y::AbsSpec) = morphism_type(underlying_spec_type(typeof(X)), underlying_spec_type(typeof(Y)))


### getter functions
pullback(phi::SpecMor) = phi.pullback
domain(phi::SpecMor) = phi.domain
codomain(phi::SpecMor) = phi.codomain

### additional constructors
function SpecMor(
      X::AbsSpec,
      Y::AbsSpec{<:Ring, <:MPolyRing},
      f::Vector{<:RingElem};
      check::Bool=true
  )
  return SpecMor(X, Y, hom(OO(Y), OO(X), OO(X).(f)), check=check)
end

function SpecMor(
      X::AbsSpec,
      Y::AbsSpec,
      f::Vector{<:RingElem};
      check::Bool=true
  )
  return SpecMor(X, Y, hom(OO(Y), OO(X), OO(X).(f), check=check), check=check)
end

function SpecMor(
      X::AbsSpec,
      Y::AbsSpec,
      f::Vector;
      check::Bool=true
  )
  return SpecMor(X, Y, OO(X).(f), check=check)
end

identity_map(X::AbsSpec{<:Any, <:MPolyQuoLocalizedRing}) = SpecMor(X, X, hom(OO(X), OO(X), gens(base_ring(OO(X))), check=false))
identity_map(X::AbsSpec{<:Any, <:MPolyLocalizedRing}) = SpecMor(X, X, hom(OO(X), OO(X), gens(base_ring(OO(X))), check=false))
identity_map(X::AbsSpec{<:Any, <:MPolyRing}) = SpecMor(X, X, hom(OO(X), OO(X), gens(OO(X))))
identity_map(X::AbsSpec{<:Any, <:MPolyQuo}) = SpecMor(X, X, hom(OO(X), OO(X), gens(ambient_ring(X))))
inclusion_map(X::AbsSpec, Y::AbsSpec) = SpecMor(X, Y, gens(base_ring(OO(Y))))  # TODO: Remove
inclusion_morphism(X::AbsSpec, Y::AbsSpec; check::Bool=true) = SpecMor(X, Y, gens(ambient_ring(Y)), check=check)

function restrict(f::SpecMor, U::AbsSpec, V::AbsSpec; check::Bool=true)
  if check
    issubset(U, domain(f)) || error("second argument does not lay in the domain of the map")
    issubset(V, codomain(f)) || error("third argument does not lay in the codomain of the map")
    issubset(U, preimage(f, V)) || error("the image of the restriction is not contained in the restricted codomain")
  end
  return SpecMor(U, V, OO(U).(pullback(f).(gens(domain(pullback(f))))), check=check)
end

# TODO: Alias for compatibility. Needs to be cleaned up and removed.
restriction(f::SpecMor, U::AbsSpec, V::AbsSpec; check::Bool=true) = restrict(f,U,V,check=check)


function compose(f::AbsSpecMor, g::AbsSpecMor; check::Bool=true)
  codomain(f) == domain(g) || error("Morphisms can not be composed")
  return SpecMor(domain(f), codomain(g), compose(pullback(g), pullback(f)), check=check)
end

function ==(f::SpecMorType, g::SpecMorType) where {SpecMorType<:AbsSpecMor}
  X = domain(f)
  X == domain(g) || return false
  codomain(f) == codomain(g) || return false
  OO(X).(pullback(f).(gens(ambient_ring(codomain(f))))) == OO(X).(pullback(f).(gens(ambient_ring(codomain(g))))) || return false
  return true
end

### functionality
function preimage(
    phi::AbsSpecMor,
    Z::AbsSpec{<:Ring, <:MPolyQuoLocalizedRing{<:Any, <:Any, <:Any, <:Any, 
                                               <:MPolyPowersOfElement}};
    check::Bool=true
  )
  X = domain(phi)
  Y = codomain(phi)
  check && (issubset(Z, Y) || (Z = intersect(Y, Z)))
  IZ = modulus(OO(Z))
  a = denominators(inverted_set(OO(Z)))
  R = ambient_ring(X)
  f = pullback(phi)
  new_units = [lifted_numerator(f(d)) for d in a]
  new_gens = lifted_numerator.(f.(gens(IZ)))
  return hypersurface_complement(subscheme(X, new_gens), new_units)
end

function preimage(f::AbsSpecMor, Z::AbsSpec{<:Ring, <:MPolyRing}; check::Bool=true)
  OO(Z) == ambient_ring(codomain(f)) || error("schemes can not be compared")
  return subscheme(domain(f), ideal(OO(domain(f)), [zero(OO(domain(f)))]))
end

function preimage(f::AbsSpecMor, 
    Z::AbsSpec{<:Ring, <:MPolyLocalizedRing{<:Any, <:Any, <:Any, <:Any, 
                                            <:MPolyPowersOfElement}}; 
    check::Bool=true)
  return hypersurface_complement(domain(f), pullback(f).(denominators(inverted_set(OO(Z)))))
end

function preimage(f::AbsSpecMor, Z::AbsSpec; check::Bool=true)
  pbf = pullback(f)
  R = OO(codomain(f))
  S = OO(domain(f))
  I = ideal(R, gens(modulus(OO(Z))))
  J = ideal(S, pbf.(gens(I)))
  return subscheme(domain(f), J)
end

@attr Bool function is_isomorphism(f::AbsSpecMor)
  has_attribute(f, :inverse) && return true
  is_isomorphism(pullback(f)) || return false
  set_attribute!(f, :inverse, SpecMor(codomain(f), domain(f), inverse(pullback(f))))
  return true
end

function is_inverse_of(f::S, g::T) where {S<:AbsSpecMor, T<:AbsSpecMor}
  return is_isomorphism(f) && (inverse(f) == g)
end

is_identity_map(f::AbsSpecMor) = (domain(f) == codomain(f)) && all(x->(pullback(f)(x) == x), gens(OO(domain(f))))

@attr AbsSpecMor function inverse(f::AbsSpecMor) 
  is_isomorphism(f) || error("the given morphism is not an isomorphism")
  return get_attribute(f, :inverse)::morphism_type(codomain(f), domain(f))
end

@Markdown.doc """
    product(X::AbsSpec, Y::AbsSpec)
    
Returns a triple ``(X×Y, p₁, p₂)`` consisting of the product ``X×Y`` over 
the common base ring ``𝕜`` and the two projections ``p₁ : X×Y → X`` and 
``p₂ : X×Y → Y``.
"""
function product(X::AbsSpec, Y::AbsSpec;
    change_var_names_to::Vector{String}=["", ""]
  )
  base_ring(X) == base_ring(Y) || error("schemes are not defined over the same base ring")
  Xstd = standard_spec(X)
  Ystd = standard_spec(Y)
  XxY, prX, prY = product(Xstd, Ystd, change_var_names_to=change_var_names_to)
  return XxY, compose(prX, SpecMor(Xstd, X, gens(OO(Xstd)))), compose(prY, SpecMor(Ystd, Y, gens(OO(Ystd))))
end

function product(X::StdSpec, Y::StdSpec;
    change_var_names_to::Vector{String}=["", ""]
  )
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
  new_symb = Symbol[]
  if length(change_var_names_to[1]) == 0
    new_symb = symbols(R)
  else 
    new_symb = Symbol.([change_var_names_to[1]*"$i" for i in 1:ngens(R)])
  end
  if length(change_var_names_to[2]) == 0
    new_symb = vcat(new_symb, symbols(S))
  else 
    new_symb = vcat(new_symb, Symbol.([change_var_names_to[2]*"$i" for i in 1:ngens(S)]))
  end
  RS, z = PolynomialRing(k, new_symb)
  inc1 = hom(R, RS, gens(RS)[1:m])
  inc2 = hom(S, RS, gens(RS)[m+1:m+n])
  IX = ideal(RS, inc1.(gens(modulus(OO(X)))))
  IY = ideal(RS, inc2.(gens(modulus(OO(Y)))))
  UX = MPolyPowersOfElement(RS, inc1.(denominators(inverted_set(OO(X)))))
  UY = MPolyPowersOfElement(RS, inc2.(denominators(inverted_set(OO(Y)))))
  XxY = Spec(RS, IX + IY, UX*UY)
  pr1 = SpecMor(XxY, X, gens(RS)[1:m], check=false)
  pr2 = SpecMor(XxY, Y, gens(RS)[m+1:m+n], check=false)
  return XxY, pr1, pr2
end

function graph(f::AbsSpecMor) 
  X = standard_spec(domain(f))
  Y = standard_spec(codomain(f))
  fres = restrict(f, X, Y)
  G, prX, prY = graph(fres)
  return G, compose(prX, SpecMor(X, domain(f), gens(OO(X)))), compose(prY, SpecMor(Y, codomain(f), gens(OO(Y))))
end

function graph(f::AbsSpecMor{SpecType, SpecType}) where {SpecType<:StdSpec}
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

@Markdown.doc """
    fiber_product(f::SpecMor, g::SpecMor)

For morphisms ``f : Y → X`` and ``g : Z → X`` return the fiber 
product ``Y ×ₓ Z`` over ``X`` together with its two canonical 
projections.
"""
function fiber_product(
    f::SpecMor{SpecType, SpecType, <:Any}, 
    g::SpecMor{SpecType, SpecType, <:Any}
  ) where {SpecType<:StdSpec}
  Y = domain(f)
  X = codomain(f)
  X == codomain(g) || error("maps need to have the same codomain")
  Z = domain(g)
  YxZ, pY, pZ = product(Y, Z)
  RX = base_ring(OO(X))
  #I = ideal(OO(YxZ), [pullback(pY)(pullback(f)(x)) - pullback(pZ)(pullback(g)(x)) for x in gens(RX)])
  W = subscheme(YxZ, [pullback(pY)(pullback(f)(x)) - pullback(pZ)(pullback(g)(x)) for x in gens(RX)])
  return W, restrict(pY, W, Y, check=false), restrict(pZ, W, Z, check=false)
end

function affine_space(kk::BRT, n::Int; variable_name="x") where {BRT<:Ring}
  R, _ = PolynomialRing(kk, [ variable_name * "$i" for i in 1:n])
  return Spec(R)
end

function affine_space(kk::BRT, var_symbols::Vector{Symbol}) where {BRT<:Ring}
  R, _ = PolynomialRing(kk, var_symbols)
  return Spec(R)
end

@attr function dim(X::AbsSpec{<:Ring, <:MPolyQuoLocalizedRing})
  return dim(saturated_ideal(localized_modulus(OO(X))))
end

@attr function dim(X::AbsSpec{<:Ring, <:MPolyLocalizedRing})
  # the following line is supposed to refer the problem to the 
  # algebra side of the problem 
  return dim(ideal(ambient_ring(X), [zero(ambient_ring(X))]))
end

@attr function dim(X::AbsSpec{<:Ring, <:MPolyRing})
  return dim(ideal(ambient_ring(X), [zero(ambient_ring(X))]))
end

@attr function dim(X::AbsSpec{<:Ring, <:MPolyQuo})
  return dim(modulus(OO(X)))
end

@attr function codim(X::Spec)
  return dim(ideal(ambient_ring(X), [zero(ambient_ring(X))])) - dim(X)
end

strict_modulus(X::Spec) = saturated_ideal(localized_modulus(OO(X)))

function simplify(X::Spec)
  L, f, g = simplify(OO(X))
  Y = Spec(L)
  YtoX = SpecMor(Y, X, f)
  XtoY = SpecMor(X, Y, g)
  set_attribute!(YtoX, :inverse, XtoY)
  set_attribute!(XtoY, :inverse, YtoX)
  return Y, YtoX, XtoY
end

function is_non_zero_divisor(f::RingElem, X::AbsSpec{<:Ring, <:MPolyRing})
  return !iszero(OO(X)(f))
end

function is_non_zero_divisor(f::RingElem, X::AbsSpec{<:Ring, <:MPolyQuo})
  I = modulus(OO(X))
  J = ideal(OO(X), f)
  return I == quotient(I, J)
end

function is_non_zero_divisor(f::RingElem, X::AbsSpec{<:Ring, <:MPolyLocalizedRing})
  return !iszero(OO(X)(f))
end

function is_non_zero_divisor(f::RingElem, X::AbsSpec{<:Ring, <:MPolyQuoLocalizedRing})
  I = ideal(OO(X), [zero(OO(X))])
  zero_ideal = Oscar.pre_image_ideal(I)
  J = Oscar.pre_image_ideal(ideal(OO(X), [f]))
  Q = quotient(zero_ideal, J)
  return zero_ideal == Q 
end
