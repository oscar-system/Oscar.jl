########################################################################
# Concrete type for projective schemes                                 #
########################################################################
@doc raw"""
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
    return new{typeof(A), elem_type(A), typeof(S), elem_type(S)}(A, n, S)
  end

  function ProjectiveScheme(S::MPolyDecRing, I::MPolyIdeal{T}) where {T<:RingElem}
    for f in gens(I)
      parent(f) == S || error("elements do not belong to the correct ring")
    end
    all(x->(total_degree(x) == 1), gens(S)) || error("ring is not standard graded")
    n = ngens(S)-1
    A = coefficient_ring(S)
    Q, _ = quo(S, I)
    return new{typeof(A), elem_type(A), typeof(Q), elem_type(Q)}(A, n, Q)
  end

  function ProjectiveScheme(Q::MPolyQuoRing{MPolyDecRingElem{T, PT}}) where {T, PT<:MPolyRingElem{T}}
    # Test disabled because `total_degree` does not work at the moment.
    #all(x->(total_degree(x) == 1), gens(Q)) || error("ring is not standard graded") 
    S = base_ring(Q)
    all(x->(total_degree(x) == 1), gens(S)) || error("ring is not standard graded")
    A = coefficient_ring(S)
    r = ngens(S)-1
    result = new{typeof(A), elem_type(A), typeof(Q), elem_type(Q)}(A, r, Q)
    return result
  end
end

function Base.show(io::IO, P::AbsProjectiveScheme{<:Any, <:MPolyQuoRing})
  print(io, "subscheme of ℙ^$(relative_ambient_dimension(P))_{$(base_ring(P))} defined as the zero locus of  $(defining_ideal(P))")
end

function Base.show(io::IO, P::AbsProjectiveScheme{<:Any, <:MPolyDecRing}) 
  print(io, "ℙ^$(relative_ambient_dimension(P))_{$(base_ring(P))}")
end

