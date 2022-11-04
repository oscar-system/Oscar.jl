export AbsCoveredScheme, CoveredScheme

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
  function CoveredScheme(kk::Ring)
    res = new{typeof(kk)}()
    res.kk = kk
    return res
  end
end

