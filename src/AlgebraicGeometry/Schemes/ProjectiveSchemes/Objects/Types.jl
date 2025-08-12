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

# Examples
```jldoctest
julia> S, _ = QQ[:x, :y, :z];

julia> Sgr, _ = grade(S);

julia> P = proj(Sgr)
Projective space of dimension 2
  over rational field
with homogeneous coordinates [x, y, z]

julia> (x, y, z) = gens(Sgr);

julia> I = ideal(Sgr, x^3 + y^3 + z^3); # a hyperplane section

julia> Q, _ = quo(Sgr, I);

julia> C = proj(Q)
Projective scheme
  over rational field
defined by ideal (x^3 + y^3 + z^3)
```
"""
@attributes mutable struct ProjectiveScheme{CoeffRingType, RingType} <: AbsProjectiveScheme{CoeffRingType, RingType}
  A::CoeffRingType  # the base ring
  r::Int        # the relative dimension
  S::RingType   # A[s₀,…,sᵣ]/I

  # fields used for caching
  C::Scheme # The affine cone of this scheme.
  Y::Scheme # the base scheme
  projection_to_base::SchemeMor
  homog_coord::Vector # the homogeneous coordinates as functions on the affine cone
  homogenization_cache::IdDict # TODO make it a WeakIdDict once available
  dehomogenization_cache::IdDict

  function ProjectiveScheme(S::MPolyDecRing)
    is_standard_graded(S) || error("ring is not standard graded")
    n = ngens(S)-1
    A = coefficient_ring(S)
    return new{typeof(A), typeof(S)}(A, n, S)
  end

  function ProjectiveScheme(S::MPolyDecRing, I::MPolyIdeal{T}) where {T<:RingElem}
    for f in gens(I)
      parent(f) == S || error("elements do not belong to the correct ring")
    end
    is_standard_graded(S) || error("ring is not standard graded")
    n = ngens(S)-1
    A = coefficient_ring(S)
    Q, _ = quo(S, I)
    return new{typeof(A), typeof(Q)}(A, n, Q)
  end

  function ProjectiveScheme(Q::MPolyQuoRing{MPolyDecRingElem{T, PT}}) where {T, PT<:MPolyRingElem{T}}
    S = base_ring(Q)
    is_standard_graded(S) || error("ring is not standard graded")
    A = coefficient_ring(S)
    r = ngens(S)-1
    result = new{typeof(A), typeof(Q)}(A, r, Q)
    return result
  end
end



