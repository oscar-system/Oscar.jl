export Scheme, SchemeMor
export EmptyScheme
export AbsSpec
export Spec
export AbsSpecMor
export SpecMor
export SpecOpen
export SpecOpenRing, SpecOpenRingElem

### Abstract type for arbitrary schemes ###############################
@Markdown.doc """
    Scheme{BaseRingType<:Ring} 

A scheme over a ring ``𝕜`` of type `BaseRingType`.
"""
abstract type Scheme{BaseRingType} end

### Abstract type for morphisms of arbitrary schemes ##################
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

### The empty scheme over a base ring #################################
struct EmptyScheme{BaseRingType}<:Scheme{BaseRingType} 
  k::BaseRingType
  function EmptyScheme(k::BaseRingType) where {BaseRingType<:Ring}
    return new{BaseRingType}(k)
  end
end

### Abstract affine schemes ###########################################
@Markdown.doc """
    AbsSpec{BaseRingType, RingType<:Ring}

An affine scheme ``X = Spec(R)`` with ``R`` of type `RingType` over 
a ring ``𝕜`` of type `BaseRingType`.
"""
abstract type AbsSpec{BaseRingType, RingType<:Ring} <: Scheme{BaseRingType} end

### Basic concrete type for affine schemes ############################
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

########################################################################
# Abstract morphisms of affine schemes                                 #
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

########################################################################
# Minimal concrete type for morphisms of affine schemes                #
########################################################################
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

########################################################################
# Type for Zariski-open subsets of affine schemes                      #
########################################################################
@Markdown.doc """
    SpecOpen{SpecType, BRT} <: Scheme{BRT}

Zariski open subset ``U`` of an affine scheme ``X = Spec(R)``. 
This stores a list of generators ``f₁,…,fᵣ`` of an ideal 
``I`` defining the complement ``Z = X ∖ U``. 
The scheme ``X`` is referred to as the *ambient scheme* and 
the list ``f₁,…,fᵣ`` as the *generators* for ``U``.
"""
@attributes mutable struct SpecOpen{SpecType, BRT} <: Scheme{BRT}
  X::SpecType # the ambient scheme
  gens::Vector # a list of functions defining the complement of the open subset

  # fields used for caching
  name::String
  patches::Vector{AbsSpec}
  intersections::Dict{Tuple{Int, Int}, AbsSpec}
  complement::AbsSpec
  complement_ideal::Ideal
  ring_of_functions::Ring

  function SpecOpen(
      X::SpecType, 
      f::Vector{RET}; 
      name::String="", 
      check::Bool=true
    ) where {SpecType<:AbsSpec, RET<:RingElem}
    for a in f
      parent(a) == ambient_ring(X) || error("element does not belong to the correct ring")
      if check
        !isempty(X) && iszero(OO(X)(a)) && error("generators must not be zero")
      end
    end
    U = new{SpecType, typeof(base_ring(X))}(X, f)
    U.intersections = Dict{Tuple{Int, Int}, SpecType}()
    length(name) > 0 && set_name!(U, name)
    return U
  end
end

########################################################################
# Rings of regular functions on Zariski open sets of affine schemes    #
########################################################################
@Markdown.doc """
    SpecOpenRing{SpecType, OpenType}

The ring of regular functions ``𝒪(X, U)`` on an open subset ``U`` of an 
affine scheme ``X``.

 * `SpecType` is the type of the affine scheme ``X`` on which 
this sheaf is defined;
 * `OpenType` is the type of the (Zariski) open subsets of ``U``.
"""
mutable struct SpecOpenRing{SpecType, OpenType} <: Ring
  scheme::SpecType
  domain::OpenType

  function SpecOpenRing(
      X::SpecType, 
      U::OpenType
    ) where {SpecType<:AbsSpec, OpenType<:SpecOpen}
    issubset(U, X) || error("open set does not lay in the scheme")
    return new{SpecType, OpenType}(X, U)
  end
end

########################################################################
# Elements of SpecOpenRings                                            #
########################################################################
@Markdown.doc """
    SpecOpenRingElem{SpecOpenType}

An element ``f ∈ 𝒪(X, U)`` of the ring of regular functions on 
an open set ``U`` of an affine scheme ``X``.

The type parameter `SpecOpenType` is the type of the open set
``U`` of ``X``.
"""
mutable struct SpecOpenRingElem{
      SpecOpenRingType<:SpecOpenRing
    } <: RingElem
  parent::SpecOpenRingType
  restrictions::Vector{<:RingElem}

  function SpecOpenRingElem(
      R::SpecOpenRingType,
      f::Vector{<:RingElem};
      check::Bool=true
    ) where {
        SpecOpenRingType<:SpecOpenRing
    }
    n = length(f)
    U = domain(R)
    n == length(affine_patches(U)) || error("the number of restrictions does not coincide with the number of affine patches")
    g = [OO(U[i])(f[i]) for i in 1:n] # will throw if conversion is not possible
    if check
      for i in 1:n-1
        for j in i+1:n
          W = U[i,j]
          OO(W)(f[i], check=false) == OO(W)(f[j], check=false) || error("elements are not compatible on overlap")
        end
      end
    end
    return new{SpecOpenRingType}(R, g)
  end
end

