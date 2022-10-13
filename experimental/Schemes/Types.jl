export Scheme, SchemeMor
export EmptyScheme
export AbsSpec, Spec
export AbsSpecMor, SpecMor
export OpenInclusion, ClosedEmbedding
export SpecOpen, SpecOpenMor
export SpecOpenRing, SpecOpenRingElem
export AbsGlueing, Glueing
export SimpleGlueing
export AbsProjectiveScheme, ProjectiveScheme
export ProjectiveSchemeMor
export Covering, CoveringMorphism
export AbsCoveredScheme, CoveredScheme
export AbsCoveredSchemeMorphism, CoveredSchemeMorphism
export VarietyFunctionField, VarietyFunctionFieldElem

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
# Special type for open inclusions of affine schemes                   #
########################################################################
@attributes mutable struct OpenInclusion{DomainType, CodomainType, PullbackType}<:AbsSpecMor{DomainType, CodomainType, PullbackType, OpenInclusion, Nothing}
  inc::SpecMor{DomainType, CodomainType, PullbackType}
  I::Ideal
  Z::Spec

  function OpenInclusion(f::AbsSpecMor, I::Ideal; check::Bool=true)
    U = domain(f)
    X = codomain(f)
    Z = subscheme(X, I)
    if check
      isempty(preimage(f, Z)) || error("image of the map is not contained in the complement of the vanishing locus of the ideal")
      #TODO: Do checks
    end
    return new{typeof(U), typeof(X), pullback_type(f)}(f, I, Z)
  end
end

########################################################################
# Special type for principal open subsets of affine schemes            #
########################################################################
@attributes mutable struct PrincipalOpenSubset{BRT, RT, AmbientType} <: AbsSpec{BRT, RT}
  X::AmbientType
  U::Spec{BRT, RT}
  f::RingElem
  inc::SpecMor

  function PrincipalOpenSubset(X::AbsSpec, f::RingElem)
    parent(f) == OO(X) || error("element does not belong to the correct ring")
    U = hypersurface_complement(X, [f])
    return new{base_ring_type(X), ring_type(U), typeof(X)}(X, U, f)
  end
  
  function PrincipalOpenSubset(X::AbsSpec, f::Vector{<:RingElem})
    all(x->(parent(x) == OO(X)), f) || return PrincipalOpenSubset(X, OO(X).(f))
    U = hypersurface_complement(X, f)
    return new{base_ring_type(X), ring_type(U), typeof(X)}(X, U, prod(f))
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
    U.intersections = Dict{Tuple{Int, Int}, AbsSpec}()
    length(name) > 0 && set_name!(U, name)
    return U
  end
  ### Conversion from PrincipalOpenSubsets
  function SpecOpen(U::PrincipalOpenSubset)
    X = ambient_scheme(U)
    h = complement_equation(U)
    V = new{typeof(X), typeof(base_ring(X))}(X, [lifted_numerator(h)])
    V.intersections = Dict{Tuple{Int, Int}, AbsSpec}()
    V.patches = [U]
    return V
  end
end

########################################################################
# Common type fo subsets of affine space                               #
########################################################################

SpecSubset = Union{SpecOpen,AbsSpec,PrincipalOpenSubset}


########################################################################
# Morphisms of Zariski-open subsets of affine schemes                  #
########################################################################
@Markdown.doc """
    SpecOpenMor{DomainType<:SpecOpen, CodomainType<:SpecOpen}

Morphisms ``f : U → V`` of open sets ``U ⊂ X`` and ``V ⊂ Y`` of affine schemes.
These are stored as morphisms ``fᵢ: Uᵢ→ Y`` on the affine patches 
``Uᵢ`` of ``U``.

The type parameters stand for the following: When ``X = Spec(R)`` and 
``Y = Spec(S)`` with ``R = (𝕜[x₁,…,xₘ]/I)[A⁻¹]`` and ``S = (𝕜[y₁,…,yₙ]/J)[B⁻¹]``
then 

 * `DomainType` is the type of the domain;
 * `CodomainType` is the type of the codomain;
affine patches of the domain to the affine ambient scheme of the codomain. 
"""
mutable struct SpecOpenMor{DomainType<:SpecOpen, 
                           CodomainType<:SpecOpen
                          }<:SchemeMor{DomainType, CodomainType, SpecOpenMor, Nothing}
  domain::DomainType
  codomain::CodomainType
  maps_on_patches::Vector{AbsSpecMor}

  # fields used for caching
  inverse::SpecOpenMor
  pullback::Hecke.Map

  function SpecOpenMor(
      U::DomainType,
      V::CodomainType,
      f::Vector{<:AbsSpecMor};
      check::Bool=true
    ) where {DomainType<:SpecOpen, CodomainType<:SpecOpen}
    Y = ambient(V)
    n = length(f)
    n == length(affine_patches(U)) || error("number of patches does not coincide with the number of maps")
    if check
      for i in 1:n
        domain(f[i]) === affine_patches(U)[i] || error("domain of definition of the map does not coincide with the patch")
        codomain(f[i]) === Y || error("codomain is not compatible")
      end
      for i in 1:n-1
	for j in i+1:n
	  A = intersect(domain(f[i]), domain(f[j]))
	  restrict(f[i], A, Y) == restrict(f[j], A, Y) || error("maps don't glue")
	end
      end
      for g in f
        is_empty(subscheme(domain(g), pullback(g).(gens(V)))) || error("image is not contained in the codomain")
      end
    end
    return new{DomainType, CodomainType}(U, V, f)
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

########################################################################
# Special Type for closed embeddings of affine schemes                 #
########################################################################
@Markdown.doc """
    ClosedEmbedding{DomainType, CodomainType, PullbackType}

A closed embedding ``f : X → Y`` of affine schemes ``X = Spec(S)`` 
into ``Y = Spec(R)`` such that ``S ≅ R/I`` via ``f`` for some 
ideal ``I ⊂ R``.
"""
@attributes mutable struct ClosedEmbedding{DomainType, 
                                           CodomainType, 
                                           PullbackType
                                          }<:AbsSpecMor{DomainType, 
                                                        CodomainType, 
                                                        PullbackType, 
                                                        ClosedEmbedding,
                                                        Nothing
                                                       }
  inc::SpecMor{DomainType, CodomainType, PullbackType}
  I::Ideal
  U::SpecOpen

  function ClosedEmbedding(X::AbsSpec, I::Ideal)
    base_ring(I) == OO(X) || error("ideal does not belong to the correct ring")
    Y = subscheme(X, I)
    inc = SpecMor(Y, X, hom(OO(X), OO(Y), gens(OO(Y))))
    return new{typeof(Y), typeof(X), pullback_type(inc)}(inc, I)
  end
end

########################################################################
# Abstract glueings for affine schemes                                 #
########################################################################
abstract type AbsGlueing{LeftSpecType<:AbsSpec, 
                         RightSpecType<:AbsSpec,
                         LeftOpenType<:Scheme, 
                         RightOpenType<:Scheme,
                         LeftMorType<:Hecke.Map,
                         RightMorType<:Hecke.Map
                        } end

########################################################################
# Concrete type for general glueings                                   #
########################################################################
@Markdown.doc """
    Glueing{SpecType<:Spec, OpenType<:SpecOpen, MorType<:SpecOpenMor}

Glueing of two affine schemes ``X ↩ U ≅ V ↪ Y`` along open subsets 
``U ⊂ X`` and ``V ⊂ Y via some isomorphism ``φ : U → V``.
"""
@attributes mutable struct Glueing{
                                   LeftSpecType<:AbsSpec, 
                                   RightSpecType<:AbsSpec,
                                   LeftOpenType<:SpecOpen, 
                                   RightOpenType<:SpecOpen,
                                   LeftMorType<:SpecOpenMor,
                                   RightMorType<:SpecOpenMor
                                  } <: AbsGlueing{
                                   LeftSpecType,
                                   RightSpecType,
                                   LeftOpenType,
                                   RightOpenType,
                                   LeftMorType,
                                   RightMorType
                                  }
  X::LeftSpecType	
  Y::RightSpecType
  U::LeftOpenType
  V::RightOpenType
  f::LeftMorType # f : U → V 
  g::RightMorType

  function Glueing(
      X::AbsSpec, Y::AbsSpec, f::SpecOpenMor, g::SpecOpenMor; check::Bool=true
    )
    ambient(domain(f)) === X || error("the domain of the glueing morphism is not an open subset of the first argument")
    ambient(codomain(f)) === Y || error("the codomain of the glueing morphism is not an open subset of the second argument")
    if check
      (domain(f) === codomain(g) &&
      domain(g) ===  codomain(f)) || error("maps can not be isomorphisms")
      compose(f, g) == identity_map(domain(f)) || error("glueing maps are not inverse of each other")
      compose(g, f) == identity_map(domain(g)) || error("glueing maps are not inverse of each other")
    end
    return new{typeof(X), typeof(Y), 
               typeof(domain(f)), typeof(domain(g)), 
               typeof(f), typeof(g)
              }(X, Y, domain(f), domain(g), f, g)
  end
end

########################################################################
# Special type for simple glueings of affine schemes along principal 
# open subsets   
#
# SimpleGlueing is for glueings X ↩ U ≅ V ↪ Y along principal 
# open subsets U ⊂ X and V ⊂ Y along identifications f : U ↔ V : g. 
# For general glueings it can not be guaranteed to have this setup, 
# but it is a situation often encountered and with significant 
# simplification of underlying algorithms in the background. 
# Hence, the special type.
########################################################################
@attributes mutable struct SimpleGlueing{LST<:AbsSpec, 
                                         RST<:AbsSpec, 
                                         LOT<:PrincipalOpenSubset, 
                                         ROT<:PrincipalOpenSubset, 
                                         LMT<:AbsSpecMor, 
                                         RMT<:AbsSpecMor
                                        } <: AbsGlueing{LST, RST, LOT, ROT, LMT, RMT} 
  X::LST
  Y::RST
  U::LOT
  V::ROT
  f::LMT
  g::RMT

  function SimpleGlueing(
      X::AbsSpec, Y::AbsSpec, 
      f::AbsSpecMor{<:PrincipalOpenSubset}, 
      g::AbsSpecMor{<:PrincipalOpenSubset};
      check::Bool=true
    )
    U = domain(f)
    V = domain(g)
    X === ambient_scheme(U) && Y === ambient_scheme(V) || error("schemes are not compatible")
    domain(f) === codomain(g) && domain(g) === codomain(f) || error("maps are not compatible")
    if check
      is_identity_map(compose(f, g)) || error("maps are not inverse to each other")
      is_identity_map(compose(g, f)) || error("maps are not inverse to each other")
    end
    set_attribute!(f, :inverse, g)
    set_attribute!(g, :inverse, f)
    return new{typeof(X), typeof(Y), 
               typeof(U), typeof(V), 
               typeof(f), typeof(g)
              }(X, Y, U, V, f, g)
  end
end

########################################################################
# Abstract projective schemes                                          #
########################################################################
abstract type AbsProjectiveScheme{BaseRingType, RingType} <: Scheme{BaseRingType} end

########################################################################
# Concrete type for projective schemes                                 #
########################################################################
@Markdown.doc """
    ProjectiveScheme{CoeffRingType, CoeffRingElemType, RingType, RingElemType}

Closed subschemes ``X ⊂ ℙʳ(A)`` of projective space of `fiber_dimension` ``r`` 
over a ring of coefficients ``A`` of type `CoeffRingType` with elements of 
type `CoeffRingElemType`. The subscheme ``X`` is given by means of a homogeneous 
ideal ``I`` in the graded ring ``A[s₀,…,sᵣ]`` and the latter is of type 
`RingType` with elements of type `RingElemType`.
"""
@attributes mutable struct ProjectiveScheme{CoeffRingType, CoeffRingElemType, RingType, RingElemType} <: AbsProjectiveScheme{CoeffRingType, RingType}
  A::CoeffRingType	# the base ring
  r::Int	# the relative dimension
  S::RingType   # A[s₀,…,sᵣ]
  I::MPolyIdeal{RingElemType} # generators for the defining ideal

  # fields used for caching
  C::Scheme # The affine cone of this scheme.
  Y::Scheme # the base scheme 
  projection_to_base::AbsSpecMor
  homog_coord::Vector # the homogeneous coordinates as functions on the affine cone

  function ProjectiveScheme(S::MPolyRing_dec)
    all(x->(total_degree(x) == 1), gens(S)) || error("ring is not standard graded")
    n = ngens(S)-1
    A = coefficient_ring(S)
    I = ideal(S, [zero(S)])
    return new{typeof(A), elem_type(A), typeof(S), elem_type(S)}(A, n, S, I)
  end

  function ProjectiveScheme(S::MPolyRing_dec, I::MPolyIdeal{T}) where {T<:RingElem}
    for f in gens(I)
      parent(f) == S || error("elements do not belong to the correct ring")
    end
    all(x->(total_degree(x) == 1), gens(S)) || error("ring is not standard graded")
    n = ngens(S)-1
    A = coefficient_ring(S)
    return new{typeof(A), elem_type(A), typeof(S), elem_type(S)}(A, n, S, I)
  end

  function ProjectiveScheme(Q::MPolyQuo{MPolyElem_dec{T, AbstractAlgebra.Generic.MPoly{T}}}) where {T}
    all(x->(total_degree(x) == 1), gens(S)) || error("ring is not standard graded")
    S = base_ring(Q)
    A = coefficient_ring(S)
    I = gens(modulus(Q))
    r = ngens(S)-1
    return new{typeof(A), elem_type(A), typeof(S), elem_type(S)}(A, r, S, I)
  end
end

########################################################################
# Morphisms of projective schemes                                      #
########################################################################
@Markdown.doc """
    ProjectiveSchemeMor

A morphism of projective schemes 

    ℙˢ(B)     ℙʳ(A)
      ∪         ∪
      P    →    Q
      ↓         ↓
   Spec(B) → Spec(A)
    
given by means of a commutative diagram of homomorphisms of 
graded rings 

  A[v₀,…,vᵣ] → B[u₀,…,uₛ]
      ↑            ↑
      A      →     B

If no morphism `A → B` of the base rings is specified, then 
both ``P`` and ``Q`` are assumed to be defined in relative projective 
space over the same ring with the identity on the base. 
"""
mutable struct ProjectiveSchemeMor{
    DomainType<:ProjectiveScheme, 
    CodomainType<:ProjectiveScheme, 
    PullbackType<:Hecke.Map, 
    BaseMorType
  } <: SchemeMor{DomainType, CodomainType,
                 ProjectiveSchemeMor, 
                 BaseMorType
                }
  domain::DomainType
  codomain::CodomainType
  pullback::PullbackType

  #fields for caching
  map_on_base_schemes::SchemeMor
  map_on_affine_cones::SchemeMor

  ### Simple morphism of projective schemes over the same base scheme
  function ProjectiveSchemeMor(
      P::DomainType,
      Q::CodomainType,
      f::PullbackType;
      check::Bool=true
    ) where {DomainType<:ProjectiveScheme, CodomainType<:ProjectiveScheme, PullbackType<:Map}
    T = ambient_ring(P)
    S = ambient_ring(Q)
    (S === domain(f) && T === codomain(f)) || error("pullback map incompatible")
    if check
      #TODO: Check map on ideals (not available yet)
    end
    return new{DomainType, CodomainType, PullbackType, Nothing}(P, Q, f)
  end

  ### complicated morphisms over a non-trivial morphism of base schemes
  function ProjectiveSchemeMor(
      P::DomainType,
      Q::CodomainType,
      f::PullbackType,
      h::BaseMorType;
      check::Bool=true
    ) where {DomainType<:ProjectiveScheme, 
             CodomainType<:ProjectiveScheme, 
             PullbackType<:Map,
             BaseMorType<:SchemeMor
            }
    T = ambient_ring(P)
    S = ambient_ring(Q)
    (S === domain(f) && T === codomain(f)) || error("pullback map incompatible")
    pbh = pullback(h)
    codomain(h) == coefficient_ring(T) || error("base scheme map not compatible")
    domain(h) == coefficient_ring(S) || error("base scheme map not compatible")
    if check
      T(pbh(one(domain(h)))) == f(S(one(domain(h)))) == one(T) || error("maps not compatible")
      coefficient_map(f) == pbh || error("maps not compatible")
    end
    return new{DomainType, CodomainType, PullbackType, BaseMorType}(P, Q, f, h)
  end
end

########################################################################
# Coverings for covered schemes                                        #
########################################################################
@Markdown.doc """
    Covering

A covering of a scheme ``X`` by affine patches ``Uᵢ`` which are glued 
along isomorphisms ``gᵢⱼ : Uᵢ⊃ Vᵢⱼ →  Vⱼᵢ ⊂ Uⱼ``.

**Note:** The distinction between the different affine patches of the scheme 
is made from their hashes. Thus, an affine scheme must not appear more than once 
in any covering!
"""
mutable struct Covering{BaseRingType}
  patches::Vector{<:AbsSpec} # the basic affine patches of X
  glueings::IdDict{Tuple{<:AbsSpec, <:AbsSpec}, <:AbsGlueing} # the glueings of the basic affine patches
  affine_refinements::IdDict{<:AbsSpec, <:Vector{<:Tuple{<:SpecOpen, Vector{<:RingElem}}}} # optional lists of refinements
      # of the basic affine patches.
      # These are stored as pairs (U, a) where U is a 'trivial' SpecOpen, 
      # meaning that its list of hypersurface equation (f₁,…,fᵣ) has empty 
      # intersection in the basic affine patch X and hence satisfies 
      # some equality 1 ≡ a₁⋅f₁ + a₂⋅f₂ + … + aᵣ⋅fᵣ on X. 
      # Since the coefficients aᵢ of this equality are crucial for computations, 
      # we store them in an extra tuple. 

  # fields for caching
  glueing_graph::Graph{Undirected}
  transition_graph::Graph{Undirected}
  edge_dict::Dict{Tuple{Int, Int}, Int}

  function Covering(
      patches::Vector{<:AbsSpec},
      glueings::IdDict{Tuple{<:AbsSpec, <:AbsSpec}, <:AbsGlueing};
      check::Bool=true,
      affine_refinements::IdDict{
          <:AbsSpec, 
          <:Vector{<:Tuple{<:SpecOpen, <:Vector{<:RingElem}}}
         }=IdDict{AbsSpec, Vector{Tuple{SpecOpen, Vector{RingElem}}}}()
    )
    n = length(patches)
    n > 0 || error("can not glue the empty scheme")
    kk = coefficient_ring(ambient_ring(patches[1]))
    for i in 2:n
      kk == coefficient_ring(base_ring(OO(patches[i]))) || error("schemes are not defined over the same base ring")
    end
    # Check that no patch appears twice
    for i in 1:n-1
      for j in i+1:n
        patches[i] === patches[j] && error("affine schemes must not appear twice among the patches")
      end
    end
    for (X, Y) in keys(glueings)
      X in patches || error("glueings are not compatible with the patches")
      Y in patches || error("glueings are not compatible with the patches")
      if haskey(glueings, (Y, X))
        if check
          inverse(glueings[(X, Y)]) == glueings[(Y, X)] || error("glueings are not inverse of each other")
        end
      else
        glueings[(Y, X)] = inverse(glueings[(X, Y)])
      end
    end

    # check the affine refinements
    for U in keys(affine_refinements)
      for (V, a) in affine_refinements[U]
        ambient(V) == U && error("the ambient scheme of the refinement of X must be X")
        U in patches && error("the ambient scheme of the refinement can not be found in the affine patches")
        if check
          isone(OO(U)(sum([c*g for (c, g) in zip(a, gens(U))]))) || error("the patch $V does not cover $U")
        end
      end
    end
    return new{base_ring_type(patches[1])}(patches, glueings, affine_refinements)
  end
end

########################################################################
# Morphisms of coverings                                               #
########################################################################
@Markdown.doc """
    CoveringMorphism{SpecType<:Spec, CoveringType<:Covering, SpecMorType<:SpecMor}

A morphism ``f : C → D`` of two coverings. For every patch ``U`` of ``C`` this 
provides a map `f[U']` of type `SpecMorType` from ``U' ⊂ U`` to 
some patch `codomain(f[U])` in `D` for some affine patches ``U'`` covering ``U``.

**Note:** For two affine patches ``U₁, U₂ ⊂ U`` the codomains of `f[U₁]` and `f[U₂]`
do not need to coincide! However, given the glueings in `C` and `D`, all affine maps 
have to coincide on their overlaps.
"""
mutable struct CoveringMorphism{DomainType<:Covering, CodomainType<:Covering, BaseMorType}
  domain::DomainType
  codomain::CodomainType
  morphisms::IdDict{<:AbsSpec, <:AbsSpecMor} # on a patch X of the domain covering, this
                                         # returns the morphism φ : X → Y to the corresponding 
                                         # patch Y of the codomain covering. 

  function CoveringMorphism(
      dom::DomainType, 
      cod::CodomainType, 
      mor::IdDict{<:AbsSpec, <:AbsSpecMor};
      check::Bool=true
    ) where {
             DomainType<:Covering,
             CodomainType<:Covering
            }
    # TODO: check domain/codomain compatibility
    # TODO: if check is true, check that all morphisms glue and that the domain patches 
    # cover the basic patches of `dom`.
    for U in keys(mor)
      U in dom || error("patch $U of the map not found in domain")
      codomain(mor[U]) in cod || error("codomain patch not found")
    end
    # check that the whole domain is covered
    for U in basic_patches(dom)
      if !haskey(mor, U)
        !haskey(affine_refinements(dom), U) || error("patch $U of the domain not covered")
        found = false
        for (V, a) in affine_refinements(dom)[U] 
          all(x->(haskey(mor, x)), affine_patches(V)) && (found = true)
        end
        !found && error("patch $U of the domain not covered")
      end
    end
    return new{DomainType, CodomainType, Nothing}(dom, cod, mor)
  end
end

########################################################################
# Abstract type for covered schemes                                    #
########################################################################
abstract type AbsCoveredScheme{BaseRingType} <: Scheme{BaseRingType} end

########################################################################
# A minimal implementation of AbsCoveredScheme                         #
########################################################################
@Markdown.doc """
    mutable struct CoveredScheme{
      CoveringType<:Covering, 
      CoveringMorphismType<:CoveringMorphism
    }

A covered scheme ``X`` given by means of at least one covering 
of type `CoveringType`. 

A scheme may posess several coverings which are partially ordered 
by refinement. Such refinements are special instances of `CoveringMorphism`

    ρ : C1 → C2

where for each patch ``U`` in `C1` the inclusion map ``ρ[U] : U → V`` 
into the corresponding patch ``V`` of `C2` is an open embedding for which 
both ``𝒪(U)`` and ``𝒪(V)`` have the same `base_ring` (so that they can be 
canonically compared). 
"""
@attributes mutable struct CoveredScheme{BaseRingType} <: AbsCoveredScheme{BaseRingType}
  coverings::Vector{<:Covering}
  refinements::Dict{<:Tuple{<:Covering, <:Covering}, <:CoveringMorphism}
  refinement_graph::Graph{Directed}
  kk::BaseRingType

  default_covering::Covering

  function CoveredScheme(coverings::Vector{<:Covering}, 
      refinements::Dict{Tuple{<:Covering, <:Covering}, <:CoveringMorphism}
    )
    # TODO: Check whether the refinements form a connected graph.
    BaseRingType = base_ring_type(coverings[1])
    all(x->(base_ring_type(x) == BaseRingType), coverings) || error("coverings are not compatible")
    X = new{BaseRingType}(coverings, refinements)
    X.default_covering = X.coverings[1]
    X.kk = base_ring(patches(coverings[1])[1])
    return X
  end
end

########################################################################
# Morphisms of covered schemes                                         #
########################################################################
abstract type AbsCoveredSchemeMorphism{
    DomainType<:CoveredScheme, 
    CodomainType<:CoveredScheme, 
    BaseMorphismType,
    CoveredSchemeMorphismType
   } <: SchemeMor{DomainType, CodomainType, CoveredSchemeMorphismType, BaseMorphismType}
end

########################################################################
# Concrete minimal type for morphisms of covered schemes               #
########################################################################
@attributes mutable struct CoveredSchemeMorphism{
    DomainType<:AbsCoveredScheme, 
    CodomainType<:AbsCoveredScheme, 
    BaseMorphismType
   } <: AbsCoveredSchemeMorphism{
                                 DomainType, 
                                 CodomainType, 
                                 CoveredSchemeMorphism, 
                                 BaseMorphismType
                                }
  X::DomainType
  Y::CodomainType
  f::CoveringMorphism

  function CoveredSchemeMorphism(
      X::DomainType, 
      Y::CodomainType, 
      f::CoveringMorphism{<:Any, <:Any, BaseMorType};
      check::Bool=true
    ) where {
             DomainType<:CoveredScheme, 
             CodomainType<:CoveredScheme,
             BaseMorType
            }
    domain(f) in coverings(X) || error("covering not found in domain")
    codomain(f) in coverings(Y) || error("covering not found in codomain")
    return new{DomainType, CodomainType, BaseMorType}(X, Y, f)
  end
end

########################################################################
# Rational functions on irreducible varieties                          #
########################################################################

mutable struct VarietyFunctionField{BaseRingType<:Field, 
                                    FracFieldType<:AbstractAlgebra.Generic.FracField,
                                    CoveredSchemeType<:AbsCoveredScheme,
                                    SpecType<:AbsSpec
                                   } <: Field
  kk::BaseRingType
  X::CoveredSchemeType
  U::SpecType  # representative patch to represent rational functions
  KK::FracFieldType

  function VarietyFunctionField(
      X::AbsCoveredScheme; 
      check::Bool=true,
      representative_patch::AbsSpec=default_covering(X)[1]
    )
    check && (is_irreducible(X) || error("variety is not irreducible"))
    representative_patch in default_covering(X) || error("representative patch not found")
    KK = FractionField(ambient_ring(representative_patch))
    kk = base_ring(X)
    return new{typeof(kk), typeof(KK), typeof(X), typeof(representative_patch)}(kk, X, representative_patch, KK)
  end
end

########################################################################
# Elements of VarietyFunctionFields                                    #
########################################################################
mutable struct VarietyFunctionFieldElem{FracType<:AbstractAlgebra.Generic.Frac, 
                                        ParentType<:VarietyFunctionField
                                       }
  KK::ParentType
  f::FracType

  function VarietyFunctionFieldElem(
      KK::VarietyFunctionField,
      f::AbstractAlgebra.Generic.Frac;
      check::Bool=true
    )
    representative_field(KK) == parent(f) || error("element does not have the correct parent")
    return new{typeof(f), typeof(KK)}(KK, f)
  end

  function VarietyFunctionFieldElem(
      KK::VarietyFunctionField,
      a::RingElem, b::RingElem;
      check::Bool=true
    )
    R = parent(a) 
    R == parent(b) || error("parent rings not compatible")
    R == base_ring(representative_field(KK))
    f = representative_field(KK)(a, b)
    return new{typeof(f), typeof(KK)}(KK, f)
  end
end

