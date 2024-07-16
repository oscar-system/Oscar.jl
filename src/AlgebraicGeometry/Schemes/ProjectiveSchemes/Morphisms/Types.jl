########################################################################
# Abstract type for morphisms of projective schemes for which the
# generic interface is defined.
########################################################################
abstract type AbsProjectiveSchemeMorphism{
    DomainType,
    CodomainType,
    SelfType, # The concrete type itself as required by the generic `Map` implementation
    BaseMorType
  } <: SchemeMor{DomainType, CodomainType,
                 SelfType,
                 BaseMorType
                }
end

########################################################################
# Morphisms of projective schemes                                      #
########################################################################
@doc raw"""
    ProjectiveSchemeMor

A morphism of projective schemes
```
     ℙˢ(B)     ℙʳ(A)
       ∪         ∪
       P    →    Q
       ↓         ↓
    Spec(B) → Spec(A)
```
given by means of a commutative diagram of homomorphisms of
graded rings
```
    A[v₀,…,vᵣ] → B[u₀,…,uₛ]
        ↑            ↑
        A      →     B
```
If no morphism `A → B` of the base rings is specified, then
both ``P`` and ``Q`` are assumed to be defined in relative projective
space over the same ring with the identity on the base.
"""
@attributes mutable struct ProjectiveSchemeMor{
    DomainType<:AbsProjectiveScheme,
    CodomainType<:AbsProjectiveScheme,
    PullbackType<:Map,
    BaseMorType
  } <: AbsProjectiveSchemeMorphism{DomainType, CodomainType,
                 ProjectiveSchemeMor,
                 BaseMorType
                }
  domain::DomainType
  codomain::CodomainType
  pullback::PullbackType
  base_ring_morphism::BaseMorType

  #fields for caching
  map_on_base_schemes::SchemeMor
  map_on_affine_cones::SchemeMor

  ### Simple morphism of projective schemes over the same base scheme
  function ProjectiveSchemeMor(
      P::DomainType,
      Q::CodomainType,
      f::PullbackType;
      check::Bool=true
    ) where {DomainType<:AbsProjectiveScheme,
             CodomainType<:AbsProjectiveScheme,
             PullbackType<:Map
            }
    T = homogeneous_coordinate_ring(P)
    S = homogeneous_coordinate_ring(Q)
    (S === domain(f) && T === codomain(f)) || error("pullback map incompatible")
    @check begin
      #TODO: Check map on ideals (not available yet)
      true
    end
    # TODO: Can we make this type stable? Or is it already?
    if _has_coefficient_map(f)
      return new{DomainType, CodomainType, PullbackType, typeof(coefficient_map(f))}(P, Q, f, coefficient_map(f))
    else
      return new{DomainType, CodomainType, PullbackType, Nothing}(P, Q, f)
    end
  end

  ### complicated morphisms over a non-trivial morphism of base schemes
  function ProjectiveSchemeMor(
      P::DomainType,
      Q::CodomainType,
      f::PullbackType,
      h::BaseMorType;
      check::Bool=true
    ) where {DomainType<:AbsProjectiveScheme,
             CodomainType<:AbsProjectiveScheme,
             PullbackType<:Map,
             BaseMorType<:SchemeMor
            }
    T = homogeneous_coordinate_ring(P)
    S = homogeneous_coordinate_ring(Q)
    (S === domain(f) && T === codomain(f)) || error("pullback map incompatible")
    pbh = pullback(h)
    OO(domain(h)) == coefficient_ring(T) || error("base scheme map not compatible")
    OO(codomain(h)) == coefficient_ring(S) || error("base scheme map not compatible")
    @check T(pbh(one(OO(codomain(h))))) == f(S(one(OO(codomain(h))))) == one(T) "maps not compatible"
    @check coefficient_map(f) == pbh "maps not compatible"
    return new{
               DomainType, CodomainType, PullbackType, typeof(coefficient_map(f))
              }(P, Q, f, coefficient_map(f), h)
  end
end

@attributes mutable struct ProjectiveClosedEmbedding{
    DomainType<:AbsProjectiveScheme,
    CodomainType<:AbsProjectiveScheme,
    PullbackType<:Map,
    BaseMorType, 
    IdealType<:Ideal
  } <: AbsProjectiveSchemeMorphism{DomainType, CodomainType,
                 ProjectiveClosedEmbedding,
                 BaseMorType
                }
  underlying_morphism::ProjectiveSchemeMor{DomainType, CodomainType, PullbackType, Nothing}
  ideal_of_image::IdealType

  function ProjectiveClosedEmbedding(
      P::DomainType,
      I::IdealType,
      check::Bool=true
    ) where {DomainType<:AbsProjectiveScheme, IdealType<:Ideal}
    S = homogeneous_coordinate_ring(P)
    @req base_ring(I) === S "ideal must be defined in the homogeneous coordinate ring of the scheme"
    T, pr = quo(S, I)
    Q = proj(T)
    f = ProjectiveSchemeMor(Q, P, pr, check=false)
    return new{typeof(Q), DomainType, typeof(pr), Nothing, IdealType}(f, I)
  end

  function ProjectiveClosedEmbedding(
      f::ProjectiveSchemeMor,
      I::Ideal;
      check::Bool=true
    )
    Y = codomain(f)
    SY = homogeneous_coordinate_ring(Y)
    ambient_coordinate_ring(Y) === ambient_coordinate_ring(domain(f)) || error("ambient coordinate rings are not compatible")
    base_ring(I) === SY || error("ideal does not belong to the correct ring")
    @check begin
      pbf = pullback(f)
      kernel(pbf) == I || error("ideal does not coincide with the kernel of the pullback")
    end
    return new{typeof(domain(f)), typeof(Y), typeof(pullback(f)), Nothing, typeof(I)}(f, I)
  end
end


########################################################################
# Abstract type for rational maps of projective varieties
########################################################################
abstract type AbsRationalMap{
    DomainType<:AbsProjectiveScheme,
    CodomainType<:AbsProjectiveScheme,
    SelfType, # The concrete type itself as required by the generic `Map` implementation
  } <: SchemeMor{DomainType, CodomainType,
                 SelfType,
                 Nothing
                }
end

########################################################################
# Concrete rational maps of projective varieties                       #
########################################################################
@doc raw"""
    RationalMap

A rational map of projective varieties over a field 𝕜
```
     ℙˢ     ℙʳ
     ∪      ∪
     P  →   Q
```
given by means of a commutative diagram of homomorphisms of
their `homogeneous_coordinate_rings` 
```
    𝕜[u₀,…,uₛ]/I ← 𝕜[v₀,…,vᵣ]/J
```
"""
@attributes mutable struct RationalMap{
    DomainType<:AbsProjectiveScheme,
    CodomainType<:AbsProjectiveScheme,
    PullbackType<:Map,
  } <: AbsRationalMap{
                      DomainType, CodomainType,
                      RationalMap
                     }
  domain::DomainType
  codomain::CodomainType
  pullback::PullbackType

  # Fields for caching
  graph_ring::Tuple{<:MPolyQuoRing, <:Map, <:Map}

  ### Simple morphism of projective schemes over the same base scheme
  function RationalMap(
      P::DomainType,
      Q::CodomainType,
      f::PullbackType;
      check::Bool=true
    ) where {DomainType<:AbsProjectiveScheme,
             CodomainType<:AbsProjectiveScheme,
             PullbackType<:Map
            }
    T = homogeneous_coordinate_ring(P)
    S = homogeneous_coordinate_ring(Q)
    (S === domain(f) && T === codomain(f)) || error("pullback map incompatible")
    @check begin
      is_irreducible(P) || error("domain must be irreducible")
      is_irreducible(Q) || error("codomain must be irreducible")
      is_reduced(P) || error("domain must be reduced")
      is_reduced(Q) || error("codomain must be reduced")
      true
    end
    return new{DomainType, CodomainType, PullbackType}(P, Q, f)
  end
end


