import Oscar.AlgHom

export ProjectiveScheme, base_ring, fiber_dimension, homogeneous_coordinate_ring, gens, getindex
export projective_space, subscheme
export projection_to_base, affine_cone, base_scheme, homogeneous_coordinates, convert_to_fraction, convert_to_homog_polys
export MorphismOfProjectiveSchemes, domain, codomain, images_of_variables, map_on_affine_cones, is_well_defined

@Markdown.doc """
    ProjectiveScheme{CoeffRingType, CoeffRingElemType, RingType, RingElemType}

Closed subschemes ``X ⊂ ℙʳ(A)`` of projective space of ``fiber_dimension`` ``r`` 
over a ring of coefficients ``A`` of type `CoeffRingType` with elements of 
type `CoeffRingElemType`. The subscheme ``X`` is given by means of a homogeneous 
ideal ``I`` in the graded ring ``A[s₀,…,sᵣ]``.
"""
mutable struct ProjectiveScheme{CoeffRingType, CoeffRingElemType, RingType, RingElemType}
  A::CoeffRingType	# the base ring
  r::Int	# the relative dimension
  S::RingType
  I::Vector{RingElemType} # generators for the defining ideal
  #TODO: Once MPolyIdeal is finally generic, use that instead of storing the generators.

  # fields used for caching
  C::Spec # The affine cone of this scheme.
  Y::Spec # the base scheme 
  projection_to_base::SpecMor
  homog_coord::Vector # the homogeneous coordinates as functions on the affine cone

  function ProjectiveScheme(S::MPolyRing_dec)
    #TODO: Check that all weights are equal to 1
    n = ngens(S)-1
    A = coefficient_ring(S)
    I = [zero(S)]
    return new{typeof(A), elem_type(A), typeof(S), elem_type(S)}(A, n, S, I)
  end

  function ProjectiveScheme(S::MPolyRing_dec, I::Vector{T}) where {T<:RingElem}
    for f in I
      parent(f) == S || error("elements do not belong to the correct ring")
    end
    #TODO: Check that all weights are equal to 1
    n = ngens(S)-1
    A = coefficient_ring(S)
    return new{typeof(A), elem_type(A), typeof(S), elem_type(S)}(A, n, S, I)
  end

  function ProjectiveScheme(Q::MPolyQuo{MPolyElem_dec{T, AbstractAlgebra.Generic.MPoly{T}}}) where {T}
    #TODO: Check that all weights are equal to 1
    S = base_ring(Q)
    A = coefficient_ring(S)
    I = gens(modulus(Q))
    r = ngens(S)-1
    return new{typeof(A), elem_type(A), typeof(S), elem_type(S)}(A, r, S, I)
  end
end

@Markdown.doc """
    base_ring(X::ProjectiveScheme)

On ``X ⊂ ℙʳ(A)`` this returns ``A``.
"""
base_ring(P::ProjectiveScheme) = P.A

@Markdown.doc """
    base_scheme(X::ProjectiveScheme{CRT, CRET, RT, RET}) where {CRT<:MPolyQuoLocalizedRing, CRET, RT, RET}

Return the base scheme ``Y`` for ``X ⊂ ℙʳ×ₖ Y → Y`` with ``Y`` defined over a field ``𝕜``.
"""
function base_scheme(X::ProjectiveScheme{CRT, CRET, RT, RET}) where {CRT<:Union{MPolyRing, MPolyQuoLocalizedRing}, CRET, RT, RET}
  if !isdefined(X, :Y)
    X.Y = Spec(base_ring(X))
  end
  return X.Y
end

function projection_to_base(X::ProjectiveScheme{CRT, CRET, RT, RET}) where {CRT<:MPolyRing, CRET, RT, RET}
  if !isdefined(X, :projection_to_base)
    affine_cone(X)
  end
  return X.projection_to_base
end

function projection_to_base(X::ProjectiveScheme{CRT, CRET, RT, RET}) where {CRT<:MPolyQuoLocalizedRing, CRET, RT, RET}
  if !isdefined(X, :projection_to_base)
    affine_cone(X)
  end
  return X.projection_to_base
end

@Markdown.doc """
    fiber_dimension(X::ProjectiveScheme)

On ``X ⊂ ℙʳ(A)`` this returns ``r``.
"""
fiber_dimension(P::ProjectiveScheme) = P.r

@Markdown.doc """
    homogeneous_coordinate_ring(X::ProjectiveScheme)

On ``X ⊂ ℙʳ(A)`` this returns ``A[s₀,…,sᵣ]``.
"""
homogeneous_coordinate_ring(P::ProjectiveScheme) = P.S

@Markdown.doc """
    homogeneous_coordinates(X::ProjectiveScheme)

On ``X ⊂ ℙʳ(A)`` this returns a vector with the homogeneous 
coordinates ``[s₀,…,sᵣ]`` as entries where each one of the 
``sᵢ`` is a function on the `affine cone` of ``X``.
"""
function homogeneous_coordinates(P::ProjectiveScheme)
  if !isdefined(P, :homog_coord)
    affine_cone(P)
  end
  return P.homog_coord
end

homogeneous_coordinate(P::ProjectiveScheme, i::Int) = homogeneous_coordinates(P)[i]

@Markdown.doc """
    generators_of_defining_ideal(X::ProjectiveScheme) = X.I

On ``X ⊂ ℙʳ(A)`` this returns a set of homogeneous generators of the 
ideal ``I ⊂ A[s₀,…,sᵣ]`` defining ``X``.
"""
generators_of_defining_ideal(X::ProjectiveScheme) = X.I

function Base.show(io::IO, P::ProjectiveScheme) 
  print(io, "subscheme of ℙ^$(fiber_dimension(P))_{$(base_ring(P))} defined as the zero locus of  $(generators_of_defining_ideal(P))")
end

function projective_space(A::CoeffRingType, r::Int; var_name::String="s") where {CoeffRingType<:Ring}
  R, _ = PolynomialRing(A, [var_name*"$i" for i in 0:n])
  S = grade(S, [1 for i in 0:n ])
  I = ideal(S, [zero(S)])
  return ProjectiveScheme(S, I)
end

function subscheme(P::ProjectiveScheme, I::Vector{T}) where {T<:RingElem}
  S = homogeneous_coordinate_ring(P)
  for f in I
    parent(f) == S || error("elements do not belong to the correct ring")
  end
  return ProjectiveScheme(S, vcat(I, generators_of_defining_ideal(P)))
end

@Markdown.doc """
    convert_to_fraction(X::ProjectiveScheme, f)

Convert a homogeneous polynomial ``f`` to an element of the ring of 
regular functions on the `affine_cone` of ``X``.
"""
function convert_to_fraction(X::ProjectiveScheme{CRT, CRET, RT, RET}, f::RET) where {CRT<:MPolyQuoLocalizedRing, CRET, RT, RET}
  return evaluate(map_coefficients(pullback(projection_to_base(X)), f), homogeneous_coordinates(X))
end

function convert_to_fraction(X::ProjectiveScheme{CRT, CRET, RT, RET}, f::RET) where {CRT<:MPolyRing, CRET, RT, RET}
  return evaluate(map_coefficients(pullback(projection_to_base(X)), f), homogeneous_coordinates(X))
end

function convert_to_fraction(X::ProjectiveScheme{CRT, CRET, RT, RET}, f::RET) where {CRT<:AbstractAlgebra.Field, CRET, RT, RET}
  return evaluate(f, homogeneous_coordinates(X))
end
@Markdown.doc """
    convert_to_homog_polys(X::ProjectiveScheme, f) 

Convert a regular function ``f = a/b`` on some open subset of the affine 
cone of ``X`` to a pair of homogeneous polynomials ``(p, q)``, lifting 
``a`` and ``b``, respectively.
"""
function convert_to_homog_polys(X::ProjectiveScheme{CRT, CRET, RT, RET}, f::MPolyLocalizedRingElem) where {CRT<:MPolyRing, CRET, RT, RET}
  S = homogeneous_coordinate_ring(X)
  A = base_ring(S)
  a = numerator(f)
  b = denominator(f)
  a_nested = renest(S, a)
  b_nested = renest(S, b)
  return (a_nested, b_nested)
end

function convert_to_homog_polys(X::ProjectiveScheme{CRT, CRET, RT, RET}, f::MPolyQuoLocalizedRingElem) where {CRT<:MPolyRing, CRET, RT, RET}
  return convert_to_homog_polys(X, lift(f))
end

function convert_to_homog_polys(X::ProjectiveScheme{CRT, CRET, RT, RET}, f::MPolyLocalizedRingElem) where {CRT<:MPolyQuoLocalizedRing, CRET, RT, RET}
  S = homogeneous_coordinate_ring(X)
  A = base_ring(S)
  R = base_ring(A)
  T = PolynomialRing(R, symbols(S))
  a = numerator(f)
  b = denominator(f)
  a_nested = renest(T, a)
  b_nested = renest(T, b)
  function coeff_map(a::T) where {T<:MPolyElem} 
    return A(evaluate(a, gens(base_ring(A))))
  end
  a_evaluated = evaluate(map_coefficients(coeff_map, a_nested), gens(S))
  b_evaluated = evaluate(map_coefficients(coeff_map, b_nested), gens(S))
  return (a_evaluated, b_evaluated)
end

function convert_to_homog_polys(X::ProjectiveScheme{CRT, CRET, RT, RET}, f::MPolyQuoLocalizedRingElem) where {CRT<:MPolyQuoLocalizedRing, CRET, RT, RET}
  return convert_to_homog_polys(X, lift(f))
end
    
function convert_to_homog_polys(X::ProjectiveScheme{CRT, CRET, RT, RET}, f::MPolyLocalizedRingElem) where {CRT<:AbstractAlgebra.Field, CRET, RT, RET}
  S = homogeneous_coordinate_ring(X)
  return (evaluate(numerator(f), gens(S)), evaluate(denominator(f), gens(S)))
end

function convert_to_homog_polys(X::ProjectiveScheme{CRT, CRET, RT, RET}, f::MPolyQuoLocalizedRingElem) where {CRT<:AbstractAlgebra.Field, CRET, RT, RET}
  return convert_to_homog_polys(X, lift(f))
end
    
### This is a temporary fix that needs to be addressed in AbstractAlgebra, issue #1105
Generic.ordering(S::MPolyRing_dec) = :lex

@Markdown.doc """
    affine_cone(X::ProjectiveScheme) 

On ``X ⊂ ℙʳ(𝕜)`` this returns the affine cone ``C(X)⊂ 𝕜ʳ⁺¹`` and similar 
in the relative situation.
"""
function affine_cone(X::ProjectiveScheme{CRT, CRET, RT, RET}) where {CRT<:MPolyRing, CRET, RT, RET}
  if !isdefined(X, :C)
    A = base_ring(X)
    Y = Spec(A)
    X.Y = Y
    kk = base_ring(A)
    F = affine_space(kk, symbols(homogeneous_coordinate_ring(X)))
    C, pr_fiber, pr_base = product(F, Y)
    X.homog_coord = lift.([pullback(pr_fiber)(u) for u in gens(OO(F))])

    function complicated_evaluation(g::MPolyElem_dec)
      R = OO(Y)
      g1 = map_coefficients(R, g)
      return evaluate(map_coefficients(pullback(pr_base), g1), X.homog_coord)
    end

    I = ideal(OO(C), [complicated_evaluation(g) for g in generators_of_defining_ideal(X)])
    X.C = subscheme(C, I)
    X.projection_to_base = restrict(pr_base, X.C, Y)
  end
  return X.C
end

function affine_cone(X::ProjectiveScheme{CRT, CRET, RT, RET}) where {CRT<:MPolyQuoLocalizedRing, CRET, RT, RET}
  if !isdefined(X, :C)
    A = base_ring(X)
    Y = Spec(A)
    X.Y = Y
    R = base_ring(A)
    kk = coefficient_ring(R)
    F = affine_space(kk, symbols(homogeneous_coordinate_ring(X)))
    C, pr_fiber, pr_base = product(F, Y)
    X.homog_coord = lift.([pullback(pr_fiber)(u) for u in gens(OO(F))])

    function complicated_evaluation(g::MPolyElem_dec)
      R = OO(Y)
      g1 = map_coefficients(R, g)
      return evaluate(map_coefficients(pullback(pr_base), g1), X.homog_coord)
    end

    I = ideal(OO(C), [complicated_evaluation(g) for g in generators_of_defining_ideal(X)])
    X.C = subscheme(C, I)
    X.projection_to_base = restrict(pr_base, X.C, Y)
  end
  return X.C
end
    
function affine_cone(X::ProjectiveScheme{CRT, CRET, RT, RET}) where {CRT<:AbstractAlgebra.Field, CRET, RT, RET}
  if !isdefined(X, :C)
    kk = base_ring(X)
    C = affine_space(kk, symbols(homogeneous_coordinate_ring(X)))
    X.homog_coord = gens(OO(C))
    I = ideal(OO(C), [convert_to_fraction(X, g) for g in generators_of_defining_ideal(X)])
    X.C = subscheme(C, I)
  end
  return X.C
end
    
@Markdown.doc """
    change_of_rings(P::ProjectiveScheme, phi)

Given ``ℙʳ(A)`` and a ring homomorphism ``φ : A → B``, 
return ``ℙʳ(B)`` together with its morphism ``ℙʳ(B)→ ℙʳ(A)``.
"""
function change_of_rings(P::ProjectiveScheme, phi)
  error("not implemented")
end

@Markdown.doc """
    MorphismOfProjectiveSchemes{CoeffRingType, CoeffRingElemType}

A morphism ``ℙˢ(B) → ℙʳ(A)`` with ``A`` of type `CoeffRingType`, 
given by means of a homomorphism of graded rings 
``A[v₀,…,vᵣ] → A[u₀,…,uₛ]``.
"""
mutable struct MorphismOfProjectiveSchemes{CRT, CRET, RT, RET}
  domain::ProjectiveScheme{CRT, CRET, RT, RET}
  codomain::ProjectiveScheme{CRT, CRET, RT, RET}
  images_of_variables::Vector{RET}

  #fields for caching
  map_on_affine_cones::SpecMor

  function MorphismOfProjectiveSchemes(
      P::ProjectiveScheme{CRT, CRET, RT, RET},
      Q::ProjectiveScheme{CRT, CRET, RT, RET},
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
    # TODO: Once graded rings are half way functional, implement an 
    # optional check for being well defined. 
    return new{CRT, CRET, RT, RET}(P, Q, imgs)
  end
end

domain(phi::MorphismOfProjectiveSchemes) = phi.domain
codomain(phi::MorphismOfProjectiveSchemes) = phi.codomain
images_of_variables(phi::MorphismOfProjectiveSchemes) = phi.images_of_variables

function map_on_affine_cones(phi::MorphismOfProjectiveSchemes{CRT, CRET, RT, RET}) where {CRT<:MPolyQuoLocalizedRing, CRET, RT, RET}
  if !isdefined(phi, :map_on_affine_cones)
    A = base_ring(domain(phi))
    S = homogeneous_coordinate_ring(codomain(phi))
    T = homogeneous_coordinate_ring(domain(phi))
    P = domain(phi)
    Q = codomain(phi)
    pb_P = pullback(projection_to_base(P))
    pb_Q = pullback(projection_to_base(Q))
    imgs_base = pb_P.(gens(A))
    imgs_fiber = [convert_to_fraction(P, g) for g in images_of_variables(phi)]
    phi.map_on_affine_cones = SpecMor(affine_cone(P), affine_cone(Q), vcat(imgs_fiber, imgs_base))
  end
  return phi.map_on_affine_cones
end

function map_on_affine_cones(phi::MorphismOfProjectiveSchemes{CRT, CRET, RT, RET}) where {CRT<:MPolyRing, CRET, RT, RET}
  if !isdefined(phi, :map_on_affine_cones)
    Y = base_scheme(domain(phi))
    A = OO(Y)
    S = homogeneous_coordinate_ring(codomain(phi))
    T = homogeneous_coordinate_ring(domain(phi))
    P = domain(phi)
    Q = codomain(phi)
    pb_P = pullback(projection_to_base(P))
    pb_Q = pullback(projection_to_base(Q))
    imgs_base = pb_P.(gens(A))
    imgs_fiber = [convert_to_fraction(P, g) for g in images_of_variables(phi)]
    phi.map_on_affine_cones = SpecMor(affine_cone(P), affine_cone(Q), vcat(imgs_fiber, imgs_base))
  end
  return phi.map_on_affine_cones
end
    
function map_on_affine_cones(phi::MorphismOfProjectiveSchemes{CRT, CRET, RT, RET}) where {CRT<:AbstractAlgebra.Field, CRET, RT, RET}
  if !isdefined(phi, :map_on_affine_cones)
    S = homogeneous_coordinate_ring(codomain(phi))
    T = homogeneous_coordinate_ring(domain(phi))
    P = domain(phi)
    Q = codomain(phi)
    imgs_fiber = [convert_to_fraction(P, g) for g in images_of_variables(phi)]
    phi.map_on_affine_cones = SpecMor(affine_cone(P), affine_cone(Q), imgs_fiber)
  end
  return phi.map_on_affine_cones
end

function is_well_defined(phi::MorphismOfProjectiveSchemes) 
  CP = affine_cone(domain(phi))
  CQ = affine_cone(codomain(phi))
  return issubset(CP, preimage(map_on_affine_cones(phi), CQ))
end

mutable struct CoherentSheafOnProjectiveSpace{CRT, CRET, RT, RET} 
  variety::ProjectiveScheme{CRT, CRET, RT, RET}
  M # the graded module representing the sheaf 

  function CoherentSheafOnProjectiveSpace(X::ProjectiveScheme{CRT, CRET, RT, RET}, M) where {CRT, CRET, RT, RET}
    return new{CRT, CRET, RT, RET}(X, M)
  end
end

function CoherentSheafOnProjectiveSpace(M) 
  S = base_ring(M)
  X = ProjectiveScheme(S)
  return CoherentSheafOnProjectiveScheme(X, M)
end

