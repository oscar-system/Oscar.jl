import Oscar.AlgHom

export ProjectiveSpace, base_ring, fiber_dimension, homogeneous_coordinate_ring, gens, getindex
export MorphismOfProjectiveSpaces, domain, codomain, graded_pullback

@Markdown.doc """
    ProjectiveSpace{CoeffRingType, CoeffRingElemType, RingType, RingElemType}

Projective space ``ℙ^r_A`` of `fiber_dimension` ``r`` over a ring of 
coefficients ``A`` of type `CoeffRingType` with elements of type `CoeffRingElemType`.
"""
mutable struct ProjectiveSpace{CoeffRingType, CoeffRingElemType, RingType, RingElemType}
  A::CoeffRingType	# the base ring
  r::Int	# the relative dimension
  S::RingType
  I::MPolyIdeal{RingElemType}

  function ProjectiveSpace(A::CoeffRingType, r::Int; var_name::String="s") where {CoeffRingType<:Ring}
    R, _ = PolynomialRing(A, [var_name*"$i" for i in 0:n])
    S = grade(S, [1 for i in 0:n ])
    I = ideal(S, [zero(S)])
    return new{CoeffRingType, elem_type{A}, typeof(S), elem_type(S)}(A, r, S, I)
  end

  function ProjectiveSpace(S::MPolyRing_dec)
    #TODO: Check that all weights are equal to 1
    n = ngens(S)-1
    A = coefficient_ring(S)
    I = ideal(S, [zero(S)])
    return new{typeof(A), elem_type(A), typeof(S), elem_type(S)}(A, n, S, I)
  end

  function ProjectiveSpace(Q::MPolyQuo{MPolyElem_dec{T, AbstractAlgebra.Generic.MPoly{T}}}) where {T}
    #TODO: Check that all weights are equal to 1
    S = base_ring(Q)
    A = coefficient_ring(S)
    I = modulus(Q)
    r = ngens(S)-1
    return new{typeof(A), elem_type(A), typeof(S), elem_type(S)}(A, r, S, I)
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
fiber_dimension(P::ProjectiveSpace) = P.r

@Markdown.doc """
    homogeneous_coordinate_ring(P::ProjectiveSpace)

On ``ℙ^r_A`` this returns ``A[s₀,…,sᵣ]``.
"""
homogeneous_coordinate_ring(P::ProjectiveSpace) = P.S

@Markdown.doc """
    homogeneous_coordinates(P::ProjectiveSpace)

On ``ℙ^r_A`` this returns a vector with the homogeneous 
coordinates ``[s₀,…,sᵣ]`` as entries.
"""
homogeneous_coordinates(P::ProjectiveSpace) = gens(homogeneous_coordinate_ring(P))

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
mutable struct MorphismOfProjectiveSpaces{CRT, CRET, RT, RET}
  domain::ProjectiveSpace{CRT, CRET, RT, RET}
  codomain::ProjectiveSpace{CRT, CRET, RT, RET}
  graded_pullback::AlgHom{CRET}

  function MorphismOfProjectiveSpaces(
      P::ProjectiveSpace{CRT, CRET, RT, RET},
      Q::ProjectiveSpace{CRT, CRET, RT, RET},
      imgs::Vector{RET}
    ) where {CRT, CRET, RT, RET}
    A = base_ring(P)
    A == base_ring(Q) || error("projective spaces do not have the same base ring")
    m = fiber_dimension(P)
    n = fiber_dimension(Q)
    length(imgs) == n+1 || error("the number of images does not agree with the number of variables")
    for i in 1:n+1
      parent(imgs[i]) == homogeneous_coordinate_ring(P) || error("elements do not belong to the correct ring")
    end
    phi = AlgebraHomomorphism(homogeneous_coordinate_ring(Q), 
			      homogeneous_coordinate_ring(P),
			      imgs
			      )
    return new{CRT, CRET, RT, RET}(P, Q, phi)
  end
end

domain(phi::MorphismOfProjectiveSpaces) = phi.domain
codomain(phi::MorphismOfProjectiveSpaces) = phi.codomain
graded_pullback(phi::MorphismOfProjectiveSpaces) = phi.graded_pullback


mutable struct CoherentSheafOnProjectiveSpace{CRT, CRET, RT, RET} 
  variety::ProjectiveSpace{CRT, CRET, RT, RET}
  M

  function CoherentSheafOnProjectiveSpace(X::ProjectiveSpace{CRT, CRET, RT, RET}, M) where {CRT, CRET, RT, RET}
    return new{CRT, CRET, RT, RET}(X, M)
  end
end

function CoherentSheafOnProjectiveSpace(M) 
  S = base_ring(M)
  X = ProjectiveSpace(S)
  return CoherentSheafOnProjectiveSpace(X, M)
end

