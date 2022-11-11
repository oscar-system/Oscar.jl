export AbsCoveredScheme, CoveredScheme

########################################################################
# Abstract type for covered schemes                                    #
########################################################################
@Markdown.doc """
    AbsCoveredScheme{BaseRingType}

An abstract scheme ``X`` over some `base_ring` ``𝕜`` of type 
`BaseRingType`, given by means of affine charts and their glueings.
"""
abstract type AbsCoveredScheme{BaseRingType} <: Scheme{BaseRingType} end

########################################################################
# A minimal implementation of AbsCoveredScheme                         #
########################################################################
@Markdown.doc """
    CoveredScheme{BaseRingType}

A covered scheme ``X`` given by means of at least one `Covering`.

A scheme may possess several coverings which are partially ordered
by refinement. Use `default_covering(X)` to obtain one covering of ``X``.
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

