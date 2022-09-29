export ProjectiveScheme, base_ring, fiber_dimension, homogeneous_coordinate_ring, gens, getindex, affine_patch_type
export projective_scheme_type, affine_patch_type, base_ring_type, base_scheme_type, morphism_type
export projective_space, subscheme
export projection_to_base, affine_cone, set_base_scheme!, base_scheme, homogeneous_coordinates, homog_to_frac, as_covered_scheme, covered_projection_to_base, dehomogenize
export ProjectiveSchemeMor, domain, codomain, images_of_variables, map_on_affine_cones, is_well_defined, poly_to_homog, frac_to_homog_pair
export fiber_product, inclusion_map, identity_map

export ==

abstract type AbsProjectiveScheme{BaseRingType, RingType} <: Scheme{BaseRingType} end

function base_ring(P::AbsProjectiveScheme) 
  return base_ring(underlying_scheme(P))
end

@Markdown.doc """
    ambient_ring(P::AbsProjectiveScheme)

On a projective scheme ``P = Proj(S)`` with ``S = P/I`` 
for a standard graded polynomial ring ``P`` and a 
homogeneous ideal ``I`` this returns ``P``.

**Note:** This is preferred over the homogeneous coordinate 
ring ``S`` since quotient rings ``P/I`` can not be expected 
to be fully functional over arbitrary coefficient rings.
"""
function ambient_ring(P::AbsProjectiveScheme)
  return ambient_ring(underlying_scheme(P))
end

@attr AbsSpec function base_scheme(P::AbsProjectiveScheme)
  return base_scheme(underlying_scheme(P))
end

@Markdown.doc """
    affine_cone(X::ProjectiveScheme) 

On ``X ⊂ ℙʳ(𝕜)`` this returns the affine cone ``C(X)⊂ 𝕜ʳ⁺¹`` and similar 
in the relative situation.
"""
function affine_cone(P::AbsProjectiveScheme)
  return affine_cone(underlying_scheme(P))
end

@Markdown.doc """
    homog_to_frac(X::ProjectiveScheme) 

Returns a map that converts a polynomial in the 
`ambient_ring` of `X` into a function on the 
`affine_cone` of `X`.
"""
function homog_to_frac(P::AbsProjectiveScheme)
  return homog_to_frac(underlying_scheme(P))
end

@Markdown.doc """
    poly_to_homog(X::ProjectiveScheme)

Return a map that converts an element of the `base_ring` of the
ring of functions `OO` of the `affine_cone` of `X` into 
an element of the `ambient_ring` of `X`.
"""
function poly_to_homog(P::AbsProjectiveScheme)
  return poly_to_homog(underlying_scheme(P))
end

@Markdown.doc """
    frac_to_homog_pair(X::ProjectiveScheme)

Return a map that converts an element ``f = p/q`` of the ring of 
functions `OO` of the `affine_cone` of `X` into a pair 
``(a, b)`` of elements of the `ambient_ring` of `X`
corresponding to ``p`` and ``q``, respectively.
"""
function frac_to_homog_pair(P::AbsProjectiveScheme)
  return frac_to_homog_pair(underlying_scheme(P))
end


@Markdown.doc """
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
  S::RingType   # A[s₀,…,sᵣ]
  I::MPolyIdeal{RingElemType} # generators for the defining ideal

  # fields used for caching
  C::Scheme # The affine cone of this scheme.
  Y::Scheme # the base scheme 
  projection_to_base::AbsSpecMor
  homog_coord::Vector # the homogeneous coordinates as functions on the affine cone

  function ProjectiveScheme(S::MPolyRing_dec)
    all(x->(total_degree(x) == 1), gens(S)) || error("ring is not standard graded")
    n = ngens(S)-1
    A = coefficient_ring(S)
    I = ideal(S, [zero(S)])
    return new{typeof(A), elem_type(A), typeof(S), elem_type(S)}(A, n, S, I)
  end

  function ProjectiveScheme(S::MPolyRing_dec, I::MPolyIdeal{T}) where {T<:RingElem}
    for f in gens(I)
      parent(f) == S || error("elements do not belong to the correct ring")
    end
    all(x->(total_degree(x) == 1), gens(S)) || error("ring is not standard graded")
    n = ngens(S)-1
    A = coefficient_ring(S)
    return new{typeof(A), elem_type(A), typeof(S), elem_type(S)}(A, n, S, I)
  end

  function ProjectiveScheme(Q::MPolyQuo{MPolyElem_dec{T, AbstractAlgebra.Generic.MPoly{T}}}) where {T}
    all(x->(total_degree(x) == 1), gens(S)) || error("ring is not standard graded")
    S = base_ring(Q)
    A = coefficient_ring(S)
    I = gens(modulus(Q))
    r = ngens(S)-1
    return new{typeof(A), elem_type(A), typeof(S), elem_type(S)}(A, r, S, I)
  end
end

### type getters & constructors
projective_scheme_type(A::T) where {T<:AbstractAlgebra.Ring} = projective_scheme_type(typeof(A))
projective_scheme_type(::Type{T}) where {T<:AbstractAlgebra.Ring} = 
ProjectiveScheme{T, elem_type(T), mpoly_dec_ring_type(mpoly_ring_type(T)), mpoly_dec_type(mpoly_ring_type(T))}

base_ring_type(P::ProjectiveScheme{S, T, U, V}) where {S, T, U, V} = S
base_ring_type(::Type{ProjectiveScheme{S, T, U, V}}) where {S, T, U, V} = S

ring_type(P::ProjectiveScheme{S, T, U, V}) where {S, T, U, V} = U
ring_type(::Type{ProjectiveScheme{S, T, U, V}}) where {S, T, U, V} = U

base_scheme_type(P::ProjectiveScheme{S, T, U, V}) where {S, T, U, V} = spec_type(S)
base_scheme_type(::Type{ProjectiveScheme{S, T, U, V}}) where {S, T, U, V} = spec_type(S)

### type constructors 

# the type of a relative projective scheme over a given base scheme
projective_scheme_type(X::T) where {T<:AbsSpec} = projective_scheme_type(ring_type(T))
projective_scheme_type(::Type{T}) where {T<:AbsSpec} = projective_scheme_type(ring_type(T))


@Markdown.doc """
    base_ring(X::ProjectiveScheme)

On ``X ⊂ ℙʳ(A)`` this returns ``A``.
"""
base_ring(P::ProjectiveScheme) = P.A

@Markdown.doc """
    base_scheme(X::ProjectiveScheme{CRT, CRET, RT, RET}) where {CRT<:MPolyQuoLocalizedRing, CRET, RT, RET}

Return the base scheme ``Y`` for ``X ⊂ ℙʳ×ₖ Y → Y`` with ``Y`` defined over a field ``𝕜``.
"""
function base_scheme(X::ProjectiveScheme{CRT, CRET, RT, RET}) where {CRT<:Ring, CRET, RT, RET}
  if !isdefined(X, :Y)
    X.Y = Spec(base_ring(X))
  end
  return X.Y
end

function base_scheme(X::ProjectiveScheme{<:SpecOpenRing}) 
  return domain(base_ring(X))
end

function set_base_scheme!(
    P::ProjectiveScheme{CRT, CRET, RT, RET}, 
    X::AbsSpec
  ) where {CRT<:Ring, CRET, RT, RET}
  OO(X) == base_ring(P) || error("schemes are not compatible")
  P.Y = X
  return P
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

function projection_to_base(X::ProjectiveScheme{CRT, CRET, RT, RET}) where {CRT<:SpecOpenRing, CRET, RT, RET}
  if !isdefined(X, :projection_to_base)
    affine_cone(X)
  end
  return get_attribute(X, :projection_to_base)
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
ambient_ring(P::ProjectiveScheme) = homogeneous_coordinate_ring(P)

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
    defining_ideal(X::ProjectiveScheme)

On ``X ⊂ ℙʳ(A)`` this returns the homogeneous 
ideal ``I ⊂ A[s₀,…,sᵣ]`` defining ``X``.
"""
defining_ideal(X::ProjectiveScheme) = X.I

function Base.show(io::IO, P::ProjectiveScheme) 
  print(io, "subscheme of ℙ^$(fiber_dimension(P))_{$(base_ring(P))} defined as the zero locus of  $(defining_ideal(P))")
end

original_ring(S::MPolyRing_dec) = S.R

function subscheme(P::ProjectiveScheme, f::RingElemType) where {RingElemType<:MPolyElem_dec}
  S = homogeneous_coordinate_ring(P)
  parent(f) == S || error("ring element does not belong to the correct ring")
  Q = ProjectiveScheme(S, ideal(S, vcat(gens(defining_ideal(P)), [f])))
  if isdefined(P, :Y) 
    set_base_scheme!(Q, base_scheme(P))
  end
  return Q
end

function subscheme(P::ProjectiveScheme, f::Vector{RingElemType}) where {RingElemType<:MPolyElem_dec}
  length(f) == 0 && return P #TODO: Replace P by an honest copy!
  S = homogeneous_coordinate_ring(P)
  for i in 1:length(f)
    parent(f[i]) == S || error("ring element does not belong to the correct ring")
  end
  Q = ProjectiveScheme(S, ideal(S, vcat(gens(defining_ideal(P)),f)))
  if isdefined(P, :Y) 
    set_base_scheme!(Q, base_scheme(P))
  end
  return Q
end

function subscheme(P::ProjectiveScheme, I::MPolyIdeal{T}) where {T<:RingElem}
  S = homogeneous_coordinate_ring(P)
  base_ring(I) == S || error("ideal does not belong to the correct ring")
  Q = ProjectiveScheme(S, ideal(S, vcat(gens(I), gens(defining_ideal(P)))))
  if isdefined(P, :Y) 
    set_base_scheme!(Q, base_scheme(P))
  end
  return Q
end

function projective_space(A::CoeffRingType, var_symb::Vector{Symbol}) where {CoeffRingType<:Ring}
  n = length(var_symb)
  R, _ = PolynomialRing(A, var_symb)
  S, _ = grade(R, [1 for i in 1:n ])
  I = ideal(S, [zero(S)])
  return ProjectiveScheme(S, I)
end

projective_space(A::CoeffRingType, var_names::Vector{String}) where {CoeffRingType<:Ring} = projective_space(A, Symbol.(var_names))


function projective_space(A::CoeffRingType, r::Int; var_name::String="s") where {CoeffRingType<:Ring}
  R, _ = PolynomialRing(A, [var_name*"$i" for i in 0:r])
  S, _ = grade(R, [1 for i in 0:r ])
  I = ideal(S, [zero(S)])
  return ProjectiveScheme(S, I)
end

function projective_space(W::Spec, r::Int; var_name::String="s") 
  P = projective_space(OO(W), r, var_name=var_name)
  set_base_scheme!(P, W)
  return P
end

function projective_space(W::Spec, var_names::Vector{Symbol}) 
  P = projective_space(OO(W), var_names)
  set_base_scheme!(P, W)
  return P
end

function projective_space(W::Spec, var_names::Vector{String}) 
  P = projective_space(OO(W), var_names)
  set_base_scheme!(P, W)
  return P
end

function homog_to_frac(X::ProjectiveScheme) 
  if !has_attribute(X, :homog_to_frac)
    affine_cone(X)
  end
  return get_attribute(X, :homog_to_frac)
end

function poly_to_homog(X::ProjectiveScheme)
  if !has_attribute(X, :poly_to_homog)
    affine_cone(X)
  end
  return get_attribute(X, :poly_to_homog)
end

function frac_to_homog_pair(X::ProjectiveScheme)
  if !has_attribute(X, :frac_to_homog_pair)
    affine_cone(X)
  end
  return get_attribute(X, :frac_to_homog_pair)
end

    
### This is a temporary fix that needs to be addressed in AbstractAlgebra, issue #1105
Generic.ordering(S::MPolyRing_dec) = :degrevlex

function affine_cone(X::ProjectiveScheme{CRT, CRET, RT, RET}) where {CRT<:MPolyRing, CRET, RT, RET}
  if !isdefined(X, :C)
    A = base_ring(X)
    Y = Spec(A)
    X.Y = Y
    kk = base_ring(A)
    F = affine_space(kk, symbols(homogeneous_coordinate_ring(X)))
    C, pr_fiber, pr_base = product(F, Y)
    X.homog_coord = lift.([pullback(pr_fiber)(u) for u in gens(OO(F))])

    S = homogeneous_coordinate_ring(X)
    # use the new mapping types for polynomial rings.
    inner_help_map = hom(A, OO(C), [pullback(pr_base)(x) for x in gens(OO(Y))])
    help_map = hom(S, OO(C), inner_help_map, [pullback(pr_fiber)(y) for y in gens(OO(F))])

    # use the map to convert ideals:
    #I = ideal(OO(C), [help_map(g) for g in gens(defining_ideal(X))])
    I = help_map(defining_ideal(X))
    CX = subscheme(C, pre_image_ideal(I))
    set_attribute!(X, :affine_cone, CX)
    X.C = get_attribute(X, :affine_cone)
    pr_base_res = restrict(pr_base, CX, Y, check=false)
    pr_fiber_res = restrict(pr_fiber, CX, F, check=false)

    # store the various conversion maps
    set_attribute!(X, :homog_to_frac, 
                    hom(S, OO(CX), 
                          hom(A, OO(CX), [pullback(pr_base_res)(x) for x in gens(OO(Y))]),
                          [pullback(pr_fiber_res)(y) for y in gens(OO(F))]
                       )
                  )
    pth = hom(base_ring(OO(CX)), S, vcat(gens(S), S.(gens(A))))
    set_attribute!(X, :poly_to_homog, pth)
    set_attribute!(X, :frac_to_homog_pair, (f -> (pth(lifted_numerator(OO(CX)(f))), pth(lifted_denominator(OO(CX)(f))))))
    X.projection_to_base = pr_base_res
  end
  return X.C
end

function affine_cone(X::ProjectiveScheme{CRT, CRET, RT, RET}) where {CRT<:MPolyQuo, CRET, RT, RET}
  if !isdefined(X, :C)
    A = base_ring(X)
    R = base_ring(A)
    Y = base_scheme(X)
    kk = base_ring(R)
    F = affine_space(kk, symbols(homogeneous_coordinate_ring(X)))
    C, pr_fiber, pr_base = product(F, Y)
    X.homog_coord = lift.([pullback(pr_fiber)(u) for u in gens(OO(F))])

    S = homogeneous_coordinate_ring(X)
    # use the new mapping types for polynomial rings.
    inner_help_map = hom(A, OO(C), [pullback(pr_base)(x) for x in gens(OO(Y))])
    help_map = hom(S, OO(C), inner_help_map, [pullback(pr_fiber)(y) for y in gens(OO(F))])

    # use the map to convert ideals:
    #I = ideal(OO(C), [help_map(g) for g in gens(defining_ideal(X))])
    I = help_map(defining_ideal(X))
    CX = subscheme(C, I)
    set_attribute!(X, :affine_cone, CX)
    X.C = get_attribute(X, :affine_cone)
    pr_base_res = restrict(pr_base, CX, Y, check=false)
    pr_fiber_res = restrict(pr_fiber, CX, F, check=false)

    # store the various conversion maps
    set_attribute!(X, :homog_to_frac, 
                    hom(S, OO(CX), 
                          hom(A, OO(CX), [pullback(pr_base_res)(x) for x in gens(OO(Y))]),
                          [pullback(pr_fiber_res)(y) for y in gens(OO(F))]
                       )
                  )
    pth = hom(base_ring(OO(CX)), S, vcat(gens(S), S.(gens(A))))
    set_attribute!(X, :poly_to_homog, pth)
    set_attribute!(X, :frac_to_homog_pair, (f -> (pth(lifted_numerator(OO(CX)(f))), pth(lifted_denominator(OO(CX)(f))))))
    X.projection_to_base = pr_base_res
  end
  return X.C
end

function (f::MPolyAnyMap{<:MPolyRing, <:AbstractAlgebra.NCRing})(I::MPolyIdeal)
  return ideal(codomain(f), [f(g) for g in gens(I)])
end

function affine_cone(X::ProjectiveScheme{CRT, CRET, RT, RET}) where {CRT<:MPolyLocalizedRing, CRET, RT, RET}
  if !isdefined(X, :C)
    A = base_ring(X)
    Y = base_scheme(X)
    R = base_ring(A)
    kk = coefficient_ring(R)
    F = affine_space(kk, symbols(homogeneous_coordinate_ring(X)))
    C, pr_fiber, pr_base = product(F, Y)
    X.homog_coord = lift.([pullback(pr_fiber)(u) for u in gens(OO(F))])
    S = homogeneous_coordinate_ring(X)

    # store the various conversion maps
    help_map = hom(S, OO(C), 
                   (x -> pullback(pr_base)(x)),
                   [pullback(pr_fiber)(y) for y in gens(OO(F))]
                  )

    I = help_map(defining_ideal(X))
    CX = subscheme(C, pre_image_ideal(I))
    pr_base_res = restrict(pr_base, CX, Y, check=false)
    pr_fiber_res = restrict(pr_fiber, CX, F, check=false)

    set_attribute!(X, :homog_to_frac, 
                    hom(S, OO(CX), 
                        pullback(pr_base_res),
                        [pullback(pr_fiber_res)(y) for y in gens(OO(F))]
                       )
                  )
    pth = hom(base_ring(OO(CX)), S, vcat(gens(S), S.(gens(A))))
    set_attribute!(X, :poly_to_homog, pth)
    set_attribute!(X, :frac_to_homog_pair, (f -> (pth(lifted_numerator(OO(CX)(f))), pth(lifted_numerator(OO(CX)(f))))))
    X.C = CX
    X.projection_to_base = pr_base_res
  end
  return X.C
end
    
function affine_cone(X::ProjectiveScheme{CRT, CRET, RT, RET}) where {CRT<:MPolyQuoLocalizedRing, CRET, RT, RET}
  if !isdefined(X, :C)
    A = base_ring(X)
    Y = base_scheme(X)
    R = base_ring(A)
    kk = coefficient_ring(R)
    F = affine_space(kk, symbols(homogeneous_coordinate_ring(X)))
    C, pr_fiber, pr_base = product(F, Y)
    X.homog_coord = lift.([pullback(pr_fiber)(u) for u in gens(OO(F))])
    S = homogeneous_coordinate_ring(X)

    # store the various conversion maps
    help_map = hom(S, OO(C), 
                   (x -> pullback(pr_base)(x)),
                   [pullback(pr_fiber)(y) for y in gens(OO(F))]
                  )

    I = help_map(defining_ideal(X))
    CX = subscheme(C, pre_image_ideal(I))
    pr_base_res = restrict(pr_base, CX, Y, check=false)
    pr_fiber_res = restrict(pr_fiber, CX, F, check=false)

    set_attribute!(X, :homog_to_frac, 
                    hom(S, OO(CX), 
                        pullback(pr_base_res),
                        [pullback(pr_fiber_res)(y) for y in gens(OO(F))]
                       )
                  )
    pth = hom(base_ring(OO(CX)), S, vcat(gens(S), S.(gens(A))))
    set_attribute!(X, :poly_to_homog, pth)
    set_attribute!(X, :frac_to_homog_pair, (f -> (pth(lifted_numerator(OO(CX)(f))), pth(lifted_numerator(OO(CX)(f))))))
    X.C = CX
    X.projection_to_base = pr_base_res
  end
  return X.C
end
    
# assure compatibility with generic code for MPolyQuos:
lift(f::MPolyElem) = f

function affine_cone(X::ProjectiveScheme{CRT, CRET, RT, RET}) where {CRT<:AbstractAlgebra.Ring, CRET, RT, RET}
  if !isdefined(X, :C)
    kk = base_ring(X)
    C = affine_space(kk, symbols(homogeneous_coordinate_ring(X)))
    X.homog_coord = gens(OO(C))
    S = homogeneous_coordinate_ring(X)
    help_map = hom(S, OO(C), gens(OO(C)))
    I = help_map(defining_ideal(X))
    CX = subscheme(C, I)

    # store the various conversion maps
    set_attribute!(X, :homog_to_frac, hom(S, OO(CX), gens(OO(CX))))
    pth = hom(base_ring(OO(CX)), S, gens(S))
    set_attribute!(X, :poly_to_homog, pth)
    set_attribute!(X, :frac_to_homog_pair, (f -> (pth(lift(f)), one(S))))
    X.C = CX
  end
  return X.C
end

function affine_cone(X::ProjectiveScheme{CRT, CRET, RT, RET}) where {CRT<:SpecOpenRing, CRET, RT, RET}
  if !isdefined(X, :C)
    S = homogeneous_coordinate_ring(X)
    B = coefficient_ring(S)
    Y = scheme(B)
    U = domain(B)
    R = base_ring(OO(Y))
    kk = base_ring(R)
    F = affine_space(kk, symbols(homogeneous_coordinate_ring(X)))
    C, pr_base, pr_fiber = product(U, F)
    X.homog_coord = [pullback(pr_fiber)(u) 
                           for u in OO(codomain(pr_fiber)).(gens(OO(F)))]
    phi = hom(S, OO(C), pullback(pr_base), X.homog_coord)
    g = phi.(gens(defining_ideal(X)))
    CX = subscheme(C, g)
    X.C = CX

    set_attribute!(X, :homog_to_frac, compose(phi, restriction_map(C, CX)))
    set_attribute!(X, :base_scheme, U)
    set_attribute!(X, :projection_to_base, restrict(pr_base, CX, U, check=false))
  end
  return X.C
end
    
@Markdown.doc """
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
mutable struct ProjectiveSchemeMor{
    DomainType<:ProjectiveScheme, 
    CodomainType<:ProjectiveScheme, 
    PullbackType<:Hecke.Map, 
    BaseMorType
  } <: SchemeMor{DomainType, CodomainType,
                 ProjectiveSchemeMor, 
                 BaseMorType
                }
  domain::DomainType
  codomain::CodomainType
  pullback::PullbackType

  #fields for caching
  map_on_base_schemes::SchemeMor
  map_on_affine_cones::SchemeMor

  ### Simple morphism of projective schemes over the same base scheme
  function ProjectiveSchemeMor(
      P::DomainType,
      Q::CodomainType,
      f::PullbackType;
      check::Bool=true
    ) where {DomainType<:ProjectiveScheme, CodomainType<:ProjectiveScheme, PullbackType<:Map}
    T = homogeneous_coordinate_ring(P)
    S = homogeneous_coordinate_ring(Q)
    (S === domain(f) && T === codomain(f)) || error("pullback map incompatible")
    if check
      #TODO: Check map on ideals (not available yet)
    end
    return new{DomainType, CodomainType, PullbackType, Nothing}(P, Q, f)
  end

  ### complicated morphisms over a non-trivial morphism of base schemes
  function ProjectiveSchemeMor(
      P::DomainType,
      Q::CodomainType,
      f::PullbackType,
      h::BaseMorType;
      check::Bool=true
    ) where {DomainType<:ProjectiveScheme, 
             CodomainType<:ProjectiveScheme, 
             PullbackType<:Map,
             BaseMorType<:SchemeMor
            }
    T = homogeneous_coordinate_ring(P)
    S = homogeneous_coordinate_ring(Q)
    (S === domain(f) && T === codomain(f)) || error("pullback map incompatible")
    pbh = pullback(h)
    codomain(h) == coefficient_ring(T) || error("base scheme map not compatible")
    domain(h) == coefficient_ring(S) || error("base scheme map not compatible")
    if check
      T(pbh(one(domain(h)))) == f(S(one(domain(h)))) == one(T) || error("maps not compatible")
      coefficient_map(f) == pbh || error("maps not compatible")
    end
    return new{DomainType, CodomainType, PullbackType, BaseMorType}(P, Q, f, h)
  end
end

### type getters
morphism_type(P::S, Q::T) where {S<:ProjectiveScheme, T<:ProjectiveScheme} = morphism_type(S, T)
morphism_type(::Type{S}, ::Type{T}) where {S<:ProjectiveScheme, T<:ProjectiveScheme} = ProjectiveSchemeMor{S, T, MPolyAnyMap{ring_type(T), ring_type(S), morphism_type(base_ring_type(T), base_ring_type(S)), elem_type(ring_type(T))}}

morphism_type(P::S) where {S<:ProjectiveScheme} = morphism_type(S, S)
morphism_type(::Type{S}) where {S<:ProjectiveScheme} = morphism_type(S, S)

### getters 
domain(phi::ProjectiveSchemeMor) = phi.domain
codomain(phi::ProjectiveSchemeMor) = phi.codomain
pullback(phi::ProjectiveSchemeMor) = phi.pullback
base_ring_morphism(phi::ProjectiveSchemeMor) = coefficient_map(pullback(phi))

### additional constructors
function ProjectiveSchemeMor(X::T, Y::T, a::Vector{RET}) where {T<:ProjectiveScheme, RET<:MPolyElem_dec}
  base_ring(X) === base_ring(Y) || error("projective schemes must be defined over the same base ring")
  Q = homogeneous_coordinate_ring(X)
  P = homogeneous_coordinate_ring(Y)
  return ProjectiveSchemeMor(X, Y, hom(P, Q, a))
end

# in case we have honest base schemes, also make the map of schemes available
function base_map(phi::ProjectiveSchemeMor{<:ProjectiveScheme{<:MPolyQuoLocalizedRing}})
  if !isdefined(phi, :map_on_base_schemes)
    phi.map_on_base_schemes = SpecMor(base_scheme(domain(phi)), base_scheme(codomain(phi)), coefficient_map(pullback(phi)))
  end
  return phi.map_on_base_schemes::SchemeMor
end

function map_on_affine_cones(phi::ProjectiveSchemeMor{<:ProjectiveScheme{<:MPolyQuoLocalizedRing}}) 
  if !isdefined(phi, :map_on_affine_cones)
    A = base_ring(domain(phi))
    S = homogeneous_coordinate_ring(codomain(phi))
    T = homogeneous_coordinate_ring(domain(phi))
    P = domain(phi)
    Q = codomain(phi)
    pb_P = pullback(projection_to_base(P))
    pb_Q = pullback(projection_to_base(Q))
    imgs_base = pb_P.(gens(A))
    imgs_fiber = [homog_to_frac(P)(g) for g in pullback(phi).(gens(S))]
    phi.map_on_affine_cones = SpecMor(affine_cone(P), affine_cone(Q), vcat(imgs_fiber, imgs_base))
  end
  return phi.map_on_affine_cones::AbsSpecMor
end

function map_on_affine_cones(phi::ProjectiveSchemeMor{<:ProjectiveScheme{<:MPolyRing}})
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
    imgs_fiber = [homog_to_frac(P)(g) for g in pullback(phi).(gens(S))]
    phi.map_on_affine_cones = SpecMor(affine_cone(P), affine_cone(Q), vcat(imgs_fiber, imgs_base))
  end
  return phi.map_on_affine_cones::AbsSpecMor
end
    
function map_on_affine_cones(phi::ProjectiveSchemeMor{<:ProjectiveScheme{<:AbstractAlgebra.Ring}})
  if !isdefined(phi, :map_on_affine_cones)
    S = homogeneous_coordinate_ring(codomain(phi))
    T = homogeneous_coordinate_ring(domain(phi))
    P = domain(phi)
    Q = codomain(phi)
    imgs_fiber = [homog_to_frac(P)(g) for g in pullback(phi).(gens(S))]
    phi.map_on_affine_cones = SpecMor(affine_cone(P), affine_cone(Q), imgs_fiber)
  end
  return phi.map_on_affine_cones::AbsSpecMor
end

function map_on_affine_cones(phi::ProjectiveSchemeMor{<:ProjectiveScheme{<:SpecOpenRing}, <:ProjectiveScheme{<:SpecOpenRing}})
  if !isdefined(phi, :map_on_affine_cones)
    X = domain(phi)
    CX = affine_cone(X)
    P = ambient(CX)
    BX = base_scheme(X)
    BP = ambient(BX)
    Y = codomain(phi)
    CY = affine_cone(Y)
    Q = ambient(CY)
    BY = base_scheme(Y)
    BQ = ambient(BY)
    fiber_coord_imgs = homog_to_frac(X).(pullback(phi).(gens(homogeneous_coordinate_ring(Y)))) # elements in OO(CX)
    base_coord_imgs = pullback(phi).(pullback(projection_to_base(Y)).(gens(OO(BY))))
    coord_imgs = vcat(fiber_coord_imgs, base_coord_imgs)
    phi.map_on_affine_cones = SpecOpenMor(CX, CY, 
                                          [SpecMor(U, Q, (f->restriction_map(U, f)).(coord_imgs)) for U in CX], check=false)
  end
  return phi.map_on_affine_cones::SpecOpenMor
end

function is_well_defined(phi::ProjectiveSchemeMor) 
  CP = affine_cone(domain(phi))
  CQ = affine_cone(codomain(phi))
  return issubset(CP, preimage(map_on_affine_cones(phi), CQ))
end

function compose(f::T, g::T) where {T<:ProjectiveSchemeMor}
  return ProjectiveSchemeMor(domain(f), codomain(g), compose(pullback(g), pullback(f)))
end

function ==(f::ProjectiveSchemeMor, g::ProjectiveSchemeMor) 
  domain(f) === domain(g) || return false
  codomain(f) === codomain(g) || return false
  for s in gens(homogeneous_coordinate_ring(codomain(f)))
    pullback(f)(s) - pullback(g)(s) in defining_ideal(domain(f)) || return false
  end
  return true
end

function ==(f::ProjectiveSchemeMor{<:ProjectiveScheme{<:MPolyQuoLocalizedRing}}, 
            g::ProjectiveSchemeMor{<:ProjectiveScheme{<:MPolyQuoLocalizedRing}}) 
  domain(f) === domain(g) || return false
  codomain(f) === codomain(g) || return false
  return map_on_affine_cones(f) == map_on_affine_cones(g)
end

### additional constructors

function fiber_product(f::Hecke.Map{DomType, CodType}, P::ProjectiveScheme{DomType}) where {DomType<:Ring, CodType<:Ring}
  R = base_ring(P) 
  R == domain(f) || error("rings not compatible")
  Rnew = codomain(f)
  S = homogeneous_coordinate_ring(P)
  Qambient = projective_space(Rnew, symbols(S))
  Snew = homogeneous_coordinate_ring(Qambient)
  phi = hom(S, Snew, f, gens(Snew))
  Q = subscheme(Qambient, phi(defining_ideal(P)))
  return Q, ProjectiveSchemeMor(Q, P, hom(S, homogeneous_coordinate_ring(Q), f, gens(homogeneous_coordinate_ring(Q))))
end

function fiber_product(f::AbsSpecMor, P::ProjectiveScheme{<:MPolyQuoLocalizedRing})
  codomain(f) == base_scheme(P) || error("codomain and base_scheme are incompatible")
  X = domain(f)
  Y = codomain(f)
  Q_ambient = projective_space(X, symbols(homogeneous_coordinate_ring(P)))
  help_map = hom(homogeneous_coordinate_ring(P),
                 homogeneous_coordinate_ring(Q_ambient),
                 pullback(f),
                 gens(homogeneous_coordinate_ring(Q_ambient))
                )
  I = help_map(defining_ideal(P))
  Q = subscheme(Q_ambient, pre_image_ideal(I))
  return Q, ProjectiveSchemeMor(Q, P, 
                                hom(homogeneous_coordinate_ring(P),
                                    homogeneous_coordinate_ring(Q),
                                    pullback(f),
                                    gens(homogeneous_coordinate_ring(Q))
                                   )
                               )
end

fiber_product(X::AbsSpec, P::ProjectiveScheme{<:MPolyQuoLocalizedRing}) = fiber_product(inclusion_map(X, base_scheme(P)), P)

### canonical map constructors

@Markdown.doc """
    inclusion_map(P::T, Q::T)

Assuming that ``P ⊂ Q`` is a subscheme, both proper over an inclusion of 
their base schemes, this returns the associated `ProjectiveSchemeMor`.
"""
function inclusion_map(P::T, Q::T) where {T<:ProjectiveScheme{<:MPolyQuoLocalizedRing}}
  X = base_scheme(P)
  Y = base_scheme(Q)
  f = inclusion_map(X, Y) # will throw if X and Y are not compatible
  return ProjectiveSchemeMor(P, Q, 
                             hom(homogeneous_coordinate_ring(Q),
                                 homogeneous_coordinate_ring(P),
                                 pullback(f), 
                                 gens(homogeneous_coordinate_ring(P))
                                )
                            )
end

function inclusion_map(P::T, Q::T) where {T<:ProjectiveScheme{<:AbstractAlgebra.Ring}}
  A = base_ring(Q)
  B = base_ring(P)
  A === B || error("can not compare schemes for non-equal base rings") # TODO: Extend by check for canonical maps, once they are available
  return ProjectiveSchemeMor(P, Q, 
                             hom(homogeneous_coordinate_ring(Q),
                                 homogeneous_coordinate_ring(P),
                                 gens(homogeneous_coordinate_ring(P))
                                )
                            )
end

identity_map(P::ProjectiveScheme) = ProjectiveSchemeMor(P, P, 
                                                        hom(homogeneous_coordinate_ring(P),
                                                            homogeneous_coordinate_ring(P),
                                                            gens(homogeneous_coordinate_ring(P))
                                                           )
                                                       )

@attr function as_covered_scheme(P::ProjectiveScheme)
    C = standard_covering(P) 
    X = CoveredScheme(C)
    return X
end

function covered_projection_to_base(X::ProjectiveScheme{<:MPolyQuoLocalizedRing})
  if !has_attribute(X, :covered_projection_to_base) 
    C = standard_covering(X)
  end
  return get_attribute(X, :covered_projection_to_base) # TODO: establish type assertion here!
end

function covered_projection_to_base(X::ProjectiveScheme{<:MPolyLocalizedRing})
  if !has_attribute(X, :covered_projection_to_base) 
    C = standard_covering(X)
  end
  return get_attribute(X, :covered_projection_to_base) # TODO: establish type assertion here!
end

function covered_projection_to_base(X::ProjectiveScheme{<:MPolyQuo})
  if !has_attribute(X, :covered_projection_to_base) 
    C = standard_covering(X)
  end
  return get_attribute(X, :covered_projection_to_base) # TODO: establish type assertion here!
end

function covered_projection_to_base(X::ProjectiveScheme{<:MPolyRing})
  if !has_attribute(X, :covered_projection_to_base) 
    C = standard_covering(X)
  end
  return get_attribute(X, :covered_projection_to_base) # TODO: establish type assertion here!
end


function dehomogenize(
    X::ProjectiveScheme{CRT}, 
    i::Int
  ) where {
    CRT<:MPolyQuoLocalizedRing
  }
  i in 0:fiber_dimension(X) || error("the given integer is not in the admissible range")
  S = homogeneous_coordinate_ring(X)
  C = standard_covering(X)
  U = C[i+1]
  p = covered_projection_to_base(X)
  s = vcat(gens(OO(U))[1:i], [one(OO(U))], gens(OO(U))[i+1:fiber_dimension(X)])
  return hom(S, OO(U), pullback(p[U]), s)
end

function dehomogenize(
    X::ProjectiveScheme{CRT}, 
    i::Int
  ) where {
    CRT<:MPolyLocalizedRing
  }
  i in 0:fiber_dimension(X) || error("the given integer is not in the admissible range")
  S = homogeneous_coordinate_ring(X)
  C = standard_covering(X)
  U = C[i+1]
  p = covered_projection_to_base(X)
  s = vcat(gens(OO(U))[1:i], [one(OO(U))], gens(OO(U))[i+1:fiber_dimension(X)])
  return hom(S, OO(U), pullback(p[U]), s)
end

function dehomogenize(
    X::ProjectiveScheme{CRT}, 
    i::Int
  ) where {
    CRT<:MPolyRing
  }
  i in 0:fiber_dimension(X) || error("the given integer is not in the admissible range")
  S = homogeneous_coordinate_ring(X)
  C = standard_covering(X)
  U = C[i+1]
  p = covered_projection_to_base(X)
  s = vcat(gens(OO(U))[1:i], [one(OO(U))], gens(OO(U))[i+1:fiber_dimension(X)])
  return hom(S, OO(U), pullback(p[U]), s)
end

function dehomogenize(
    X::ProjectiveScheme{CRT}, 
    i::Int
  ) where {
    CRT<:MPolyQuo
  }
  i in 0:fiber_dimension(X) || error("the given integer is not in the admissible range")
  S = homogeneous_coordinate_ring(X)
  C = standard_covering(X)
  U = C[i+1]
  p = covered_projection_to_base(X)
  s = vcat(gens(OO(U))[1:i], [one(OO(U))], gens(OO(U))[i+1:fiber_dimension(X)])
  return hom(S, OO(U), pullback(p[U]), s)
end

function getindex(X::ProjectiveScheme, U::Spec)
  Xcov = as_covered_scheme(X)
  for C in coverings(Xcov)
    for j in 1:npatches(C)
      if U === C[j] 
        return C, j
      end
    end
  end
  return nothing, 0
end

function dehomogenize(
    X::ProjectiveScheme{CRT}, 
    U::AbsSpec
  ) where {
    CRT<:MPolyQuoLocalizedRing
  }
  # look up U in the coverings of X
  cover_of_U, index_of_U = X[U]
  Xcov = as_covered_scheme(X)
  S = homogeneous_coordinate_ring(X)

  s = Vector{elem_type(OO(U))}()
  if cover_of_U === standard_covering(X)
    S = homogeneous_coordinate_ring(X)
    C = standard_covering(X)
    p = covered_projection_to_base(X)
    s = vcat(gens(OO(U))[1:index_of_U-1], [one(OO(U))], gens(OO(U))[index_of_U:fiber_dimension(X)])
    return hom(S, OO(U), pullback(p[U]), s)
  else
    ref = Xcov[cover_of_U, standard_covering(X)]
    V = codomain(ref[U])
    index_of_V = standard_covering(X)[V]
    t = vcat(gens(OO(V))[1:index_of_V-1], [one(OO(V))], gens(OO(V))[index_of_V:fiber_dimension(X)])
    s = pullback(ref[U]).(t)
    pb = compose(ref[U], covered_projection_to_base(X)[V])
    return hom(S, OO(U), pullback(pb), s)
  end
end

function dehomogenize(
    X::ProjectiveScheme{CRT}, 
    i::Int
  ) where {
    CRT<:AbstractAlgebra.Ring
  }
  i in 0:fiber_dimension(X) || error("the given integer is not in the admissible range")
  S = homogeneous_coordinate_ring(X)
  C = standard_covering(X)
  U = C[i+1]
  s = vcat(gens(OO(U))[1:i], [one(OO(U))], gens(OO(U))[i+1:fiber_dimension(X)])
  return hom(S, OO(U), s)
end

### Hack for a detour to speed up mapping of elements 
# This is terribly slow in all kinds of quotient rings 
# because of massive checks for `iszero` due to memory 
# management.
function (f::Oscar.MPolyAnyMap{<:MPolyRing, <:MPolyQuoLocalizedRing, <:Nothing})(a::MPolyElem)
  if !has_attribute(f, :lifted_map)
    S = domain(f)
    W = codomain(f)
    L = localized_ring(W)
    g = hom(S, L, lift.(f.img_gens))
    set_attribute!(f, :lifted_map, g)
  end
  g = get_attribute(f, :lifted_map)
  return codomain(f)(g(a), check=false)
end

function (f::Oscar.MPolyAnyMap{<:MPolyRing, <:MPolyQuoLocalizedRing, <:MPolyQuoLocalizedRingHom})(a::MPolyElem)
  if !has_attribute(f, :lifted_map)
    S = domain(f)
    W = codomain(f)
    L = localized_ring(W)
    g = hom(S, L, x -> lift(f.coeff_map(x)), lift.(f.img_gens))
    set_attribute!(f, :lifted_map, g)
  end
  g = get_attribute(f, :lifted_map)
  return codomain(f)(g(a), check=false)
end
