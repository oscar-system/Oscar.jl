########################################################################
# Concrete type for projective schemes                                 #
########################################################################
@doc raw"""
    ProjectiveScheme{CoeffRingType, RingType}

Closed subschemes ``X ⊂ ℙʳ(A)`` of projective space of `fiber_dimension` ``r``
over a ring of coefficients ``A`` of type `CoeffRingType`.
The subscheme ``X`` is given by means of a homogeneous
ideal ``I`` in the graded ring ``A[s₀,…,sᵣ]`` and the latter is of type
`RingType`.
"""
@attributes mutable struct ProjectiveScheme{CoeffRingType, RingType} <: AbsProjectiveScheme{CoeffRingType, RingType}
  A::CoeffRingType	# the base ring
  r::Int	# the relative dimension
  S::RingType   # A[s₀,…,sᵣ]/I

  # fields used for caching
  C::Scheme # The affine cone of this scheme.
  Y::Scheme # the base scheme
  projection_to_base::SchemeMor
  homog_coord::Vector # the homogeneous coordinates as functions on the affine cone
  homogenization_cache::IdDict # TODO make it a WeakIdDict once available
  dehomogenization_cache::IdDict

  function ProjectiveScheme(S::MPolyDecRing)
    all(x->(total_degree(x) == 1), gens(S)) || error("ring is not standard graded")
    n = ngens(S)-1
    A = coefficient_ring(S)
    return new{typeof(A), typeof(S)}(A, n, S)
  end

  function ProjectiveScheme(S::MPolyDecRing, I::MPolyIdeal{T}) where {T<:RingElem}
    for f in gens(I)
      parent(f) == S || error("elements do not belong to the correct ring")
    end
    all(x->(total_degree(x) == 1), gens(S)) || error("ring is not standard graded")
    n = ngens(S)-1
    A = coefficient_ring(S)
    Q, _ = quo(S, I)
    return new{typeof(A), typeof(Q)}(A, n, Q)
  end

  function ProjectiveScheme(Q::MPolyQuoRing{MPolyDecRingElem{T, PT}}) where {T, PT<:MPolyRingElem{T}}
    # Test disabled because `total_degree` does not work at the moment.
    #all(x->(total_degree(x) == 1), gens(Q)) || error("ring is not standard graded") 
    S = base_ring(Q)
    all(x->(total_degree(x) == 1), gens(S)) || error("ring is not standard graded")
    A = coefficient_ring(S)
    r = ngens(S)-1
    result = new{typeof(A), typeof(Q)}(A, r, Q)
    return result
  end
end



