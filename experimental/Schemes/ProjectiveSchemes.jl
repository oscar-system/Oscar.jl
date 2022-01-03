import Oscar.AlgHom

export ProjectiveSpace, base_ring, fiber_dimension, homogeneous_coordinate_ring, gens, getindex
export MorphismOfProjectiveSpaces, domain, codomain, graded_pullback

@Markdown.doc """
ProjectiveSpace{CoeffRingType, CoeffRingElemType}

Projective space ``ℙ^r_A`` of `fiber_dimension` ``r`` over a ring of 
coefficients ``A`` of type `CoeffRingType` with elements of type `CoeffRingElemType`.
"""
mutable struct ProjectiveSpace{CoeffRingType, CoeffRingElemType}
  A::CoeffRingType	# the base ring
  n::Int	# the relative dimension
  S::MPolyRing_dec{CoeffRingElemType, AbstractAlgebra.Generic.MPolyRing{CoeffRingElemType}}

  function ProjectiveSpace(A::CoeffRingType, n::Int; var_name::String="s") where {CoeffRingType<:Ring}
    S, _ = PolynomialRing(A, [var_name*"$i" for i in 0:n])
    return new{CoeffRingType, elem_type{A}}(A, n, grade(S, [1 for i in 0:n ]))
  end

  function ProjectiveSpace(S::MPolyRing_dec)
    n = ngens(S)-1
    A = coefficient_ring(S)
    return new{typeof(A), elem_type(A)}(A, n, S)
  end
end

@Markdown.doc """
base_ring(P::ProjectiveSpace)

On ``ℙ^r_A`` this returns ``A``.
"""
base_ring(P::ProjectiveSpace) = P.A

@Markdown.doc """
fiber_dimension(P::ProjectiveSpace)

On ``ℙ^r_A`` this returns ``r``.
"""
fiber_dimension(P::ProjectiveSpace) = P.n

@Markdown.doc """
homogeneous_coordinate_ring(P::ProjectiveSpace)

On ``ℙ^r_A`` this returns ``A[s₀,…,sᵣ]``.
"""
homogeneous_coordinate_ring(P::ProjectiveSpace) = P.S

@Markdown.doc """
gens(P::ProjectiveSpace)

On ``ℙ^r_A`` this returns a vector with the homogeneous 
coordinates ``[s₀,…,sᵣ]`` as entries.
"""
gens(P::ProjectiveSpace) = gens(homogeneous_coordinate_ring(P))

@Markdown.doc """
getindex(P::ProjectiveSpace, i::Int)

On ``ℙ^r_A`` this returns the ``i``-th homogeneous coordinate
``sᵢ``, starting from ``i=0``.
"""
getindex(P::ProjectiveSpace, i::Int) = homogeneous_coordinate_ring(P)[i-1]

function Base.show(io::IO, P::ProjectiveSpace) 
  print(io, "ℙ^$(fiber_dimension(P))_{$(base_ring(P))}")
end

@Markdown.doc """
change_of_rings(P::ProjectiveSpace, phi)

Given ``ℙ^r_A`` and a ring homomorphism ``φ : A → B``, 
return ``ℙ^r_B`` together with its morphism ``ℙ^r_B → ℙ^r_A``.
"""
function change_of_rings(P::ProjectiveSpace, phi)
  error("not implemented")
end

@Markdown.doc """
MorphismOfProjectiveSpaces{CoeffRingType, CoeffRingElemType}

A morphism ``ℙ^s_A → ℙ^r_A`` with ``A`` of type `CoeffRingType`, 
given by means of a homomorphism of graded rings 
``A[v₀,…,vᵣ] → A[u₀,…,uₛ]``.
"""
mutable struct MorphismOfProjectiveSpaces{CoeffRingType, CoeffRingElemType}
  domain::ProjectiveSpace{CoeffRingType, CoeffRingElemType}
  codomain::ProjectiveSpace{CoeffRingType, CoeffRingElemType}
  graded_pullback::AlgHom{CoeffRingElemType}

  function MorphismOfProjectiveSpaces(
      P::ProjectiveSpace{CoeffRingType, CoeffRingElemType},
      Q::ProjectiveSpace{CoeffRingType, CoeffRingElemType},
      imgs::Vector{MPolyElem_dec{CoeffRingElemType, AbstractAlgebra.Generic.MPoly{CoeffRingElemType}}}
    ) where {CoeffRingType, CoeffRingElemType}
    A = base_ring(P)
    A == base_ring(Q) || error("projective spaces do not have the same base ring")
    m = fiber_dimension(P)
    n = fiber_dimension(Q)
    length(imgs) == n+1 || error("the number of images does not agree with the number of variables")
    for i in i:n+1
      parent(imgs[i]) == homogeneous_coordinate_ring(P) || error("elements do not belong to the correct ring")
    end
    phi = AlgebraHomomorphism(homogeneous_coordinate_ring(Q), 
			      homogeneous_coordinate_ring(P),
			      imgs
			      )
    return new{CoeffRingType, CoeffRingElemType}(P, Q, phi)
  end
end

domain(phi::MorphismOfProjectiveSpaces) = phi.domain
codomain(phi::MorphismOfProjectiveSpaces) = phi.codomain
graded_pullback(phi::MorphismOfProjectiveSpaces) = phi.graded_pullback
