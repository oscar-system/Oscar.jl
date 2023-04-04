########################################################################
# Morphisms of projective schemes                                      #
########################################################################
@doc raw"""
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
@attributes mutable struct ProjectiveSchemeMor{
    DomainType<:AbsProjectiveScheme,
    CodomainType<:AbsProjectiveScheme,
    PullbackType<:Hecke.Map,
    BaseMorType
  } <: SchemeMor{DomainType, CodomainType,
                 ProjectiveSchemeMor,
                 BaseMorType
                }
  domain::DomainType
  codomain::CodomainType
  pullback::PullbackType
  base_ring_morphism::Hecke.Map

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
    if check
      #TODO: Check map on ideals (not available yet)
    end
    return new{DomainType, CodomainType, PullbackType, Nothing}(P, Q, f)
  end
  
  ### Morphisms with an underlying base change
  function ProjectiveSchemeMor(
      P::DomainType,
      Q::CodomainType,
      f::PullbackType;
      check::Bool=true
    ) where {DomainType<:AbsProjectiveScheme,
             CodomainType<:AbsProjectiveScheme,
             PullbackType<:MPolyAnyMap{<:Any, <:Any, <:Map}
            }
    T = homogeneous_coordinate_ring(P)
    S = homogeneous_coordinate_ring(Q)
    (S === domain(f) && T === codomain(f)) || error("pullback map incompatible")
    if check
      #TODO: Check map on ideals (not available yet)
    end
    return new{DomainType, CodomainType, PullbackType, Nothing}(P, Q, f, coefficient_map(f))
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
    if check
      T(pbh(one(OO(codomain(h))))) == f(S(one(OO(codomain(h))))) == one(T) || error("maps not compatible")
      coefficient_map(f) == pbh || error("maps not compatible")
    end
    return new{DomainType, CodomainType, PullbackType, BaseMorType}(P, Q, f, coefficient_map(f), h)
  end
end

