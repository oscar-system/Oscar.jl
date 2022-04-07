export ProjectiveScheme, base_ring, fiber_dimension, homog_poly_ring, gens, getindex, affine_patch_type
export projective_scheme_type, affine_patch_type, base_ring_type, base_scheme_type, morphism_type
export projective_space, subscheme
export projection_to_base, affine_cone, set_base_scheme!, base_scheme, homogeneous_coordinates, homog_to_frac, as_covered_scheme, covered_projection_to_base, dehomogenize
export ProjectiveSchemeMor, domain, codomain, images_of_variables, map_on_affine_cones, is_well_defined, poly_to_homog, frac_to_homog_pair
export fiber_product, inclusion_map

export ==

@Markdown.doc """
    ProjectiveScheme{CoeffRingType, CoeffRingElemType, RingType, RingElemType}

Closed subschemes ``X ⊂ ℙʳ(A)`` of projective space of `fiber_dimension` ``r`` 
over a ring of coefficients ``A`` of type `CoeffRingType` with elements of 
type `CoeffRingElemType`. The subscheme ``X`` is given by means of a homogeneous 
ideal ``I`` in the graded ring ``A[s₀,…,sᵣ]`` and the latter is of type 
`RingType` with elements of type `RingElemType`.
"""
@attributes mutable struct ProjectiveScheme{CoeffRingType, CoeffRingElemType, RingType, RingElemType}
  A::CoeffRingType	# the base ring
  r::Int	# the relative dimension
  S::RingType   # A[s₀,…,sᵣ]
  I::MPolyIdeal{RingElemType} # generators for the defining ideal
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
    I = ideal(S, [zero(S)])
    return new{typeof(A), elem_type(A), typeof(S), elem_type(S)}(A, n, S, I)
  end

  function ProjectiveScheme(S::MPolyRing_dec, I::MPolyIdeal{T}) where {T<:RingElem}
    for f in gens(I)
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
projective_scheme_type(X::T) where {T<:Spec} = projective_scheme_type(ring_type(T))
projective_scheme_type(::Type{T}) where {T<:Spec} = projective_scheme_type(ring_type(T))

# the type of the affine cone for a projective scheme
affine_cone_type(P::ProjectiveScheme{CRT}) where {CRT<:AbstractAlgebra.Ring} = spec_type(CRT)
affine_cone_type(::Type{ProjectiveScheme{CRT}}) where {CRT<:AbstractAlgebra.Ring} = spec_type(CRT)

# again, this is the default assuming localizations at hypersurfaces
affine_cone_type(P::ProjectiveScheme{CRT}) where {CRT<:MPolyQuoLocalizedRing} = spec_type(CRT)
affine_cone_type(::Type{ProjectiveScheme{CRT}}) where {CRT<:MPolyQuoLocalizedRing} = spec_type(CRT)

# Other localization types in the base will lead to mixed localizations
# **Warning:** the methods for taking products of multiplicative sets are not type-stable themselves!
# So the following is a heuristic for what should happen.
affine_cone_type(
    P::ProjectiveScheme{MPolyQuoLocalizedRing{S, T, U, V, W}}
  ) where {
    S, T, U, V, W<:MPolyComplementOfPrimeIdeal
  } = spec_type(
    MPolyQuoLocalizedRing{
      S, T, U, V, 
      MPolyProductOfMultSets{
        S, T, U, V
      }
    }
  )
affine_cone_type(
    ::Type{ProjectiveScheme{MPolyQuoLocalizedRing{S, T, U, V, W}}}
  ) where {
    S, T, U, V, W<:MPolyComplementOfPrimeIdeal
  } = spec_type(
    MPolyQuoLocalizedRing{
      S, T, U, V, 
      MPolyProductOfMultSets{
        S, T, U, V
      }
    }
  )
affine_cone_type(
    P::ProjectiveScheme{MPolyQuoLocalizedRing{S, T, U, V, W}}
  ) where {
    S, T, U, V, W<:MPolyComplementOfKPointIdeal
  } = spec_type(
    MPolyQuoLocalizedRing{
      S, T, U, V, 
      MPolyProductOfMultSets{
        S, T, U, V
      }
    }
  )
affine_cone_type(
    ::Type{ProjectiveScheme{MPolyQuoLocalizedRing{S, T, U, V, W}}}
  ) where {
    S, T, U, V, W<:MPolyComplementOfKPointIdeal
  } = spec_type(
    MPolyQuoLocalizedRing{
      S, T, U, V, 
      MPolyProductOfMultSets{
        S, T, U, V
      }
    }
  )


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

function set_base_scheme!(P::ProjectiveScheme{CRT, CRET, RT, RET}, X::Spec) where {CRT<:MPolyQuoLocalizedRing, CRET, RT, RET}
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

@Markdown.doc """
    fiber_dimension(X::ProjectiveScheme)

On ``X ⊂ ℙʳ(A)`` this returns ``r``.
"""
fiber_dimension(P::ProjectiveScheme) = P.r

@Markdown.doc """
    homog_poly_ring(X::ProjectiveScheme)

On ``X ⊂ ℙʳ(A)`` this returns ``A[s₀,…,sᵣ]``.
"""
homog_poly_ring(P::ProjectiveScheme) = P.S

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

function affine_patch_type(X::ProjectiveScheme{CRT, CRET, RT, RET}) where {CRT<:AbstractAlgebra.Ring, CRET, RT, RET}
  return Spec{typeof(base_ring(X)), 
              elem_type(base_ring(X)), 
              typeof(original_ring(homog_poly_ring(X))), 
              elem_type(original_ring(homog_poly_ring(X))), 
              MPolyPowersOfElement{
                                   typeof(base_ring(X)), 
                                   elem_type(base_ring(X)), 
                                   typeof(original_ring(homog_poly_ring(X))), 
                                   elem_type(original_ring(homog_poly_ring(X)))
                                  }
             }
end

function affine_patch_type(X::ProjectiveScheme{CRT, CRET, RT, RET}) where {CRT<:MPolyRing, CRET, RT, RET}
  return Spec{typeof(coefficient_ring(base_ring(X))), 
              elem_type(coefficient_ring(base_ring(X))), 
              typeof(base_ring(X)), 
              elem_type(base_ring(X)),
              MPolyPowersOfElement{
                                   typeof(coefficient_ring(base_ring(X))), 
                                   elem_type(coefficient_ring(base_ring(X))), 
                                   typeof(original_ring(homog_poly_ring(X))), 
                                   elem_type(original_ring(homog_poly_ring(X)))
                                  }
             }
end

# TODO: This only supports the localizations at hypersurfaces for now. 
# In the future this will have to be modified as in the commented line; 
# but then caching on the type should be used to avoid the computations of 
# the product of multiplicative sets. 
function affine_patch_type(X::ProjectiveScheme{CRT, CRET, RT, RET}) where {CRT<:MPolyQuoLocalizedRing, CRET, RT, RET}
  Y = base_scheme(X)
  return Spec{typeof(coefficient_ring(base_ring(OO(Y)))), 
              elem_type(coefficient_ring(base_ring(OO(Y)))), 
              typeof(base_ring(OO(Y))), 
              elem_type(base_ring(OO(Y))),
              MPolyPowersOfElement{
                                   typeof(coefficient_ring(base_ring(OO(Y)))), 
                                   elem_type(coefficient_ring(base_ring(OO(Y)))), 
                                   typeof(base_ring(OO(Y))), 
                                   elem_type(base_ring(OO(Y)))
                                  }
              #typeof(inverted_set(OO(affine_cone(X)))*Localization(OO(affine_cone(X)), prod(lifted_numerator.(homogeneous_coordinates(X)))))
             }
end

function subscheme(P::ProjectiveScheme, f::RingElemType) where {RingElemType<:MPolyElem_dec}
  S = homog_poly_ring(P)
  parent(f) == S || error("ring element does not belong to the correct ring")
  Q = ProjectiveScheme(S, ideal(S, vcat(gens(defining_ideal(P)), [f])))
  if isdefined(P, :Y) 
    set_base_scheme!(Q, base_scheme(P))
  end
  return Q
end

function subscheme(P::ProjectiveScheme, f::Vector{RingElemType}) where {RingElemType<:MPolyElem_dec}
  S = homog_poly_ring(P)
  length(f) == 0 && return P #TODO: Replace P by an honest copy!
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
  S = homog_poly_ring(P)
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

@Markdown.doc """
    homog_to_frac(X::ProjectiveScheme) 

Returns a map that converts a polynomial in the 
`homog_poly_ring` of `X` into a function on the 
`affine_cone` of `X`.
"""
function homog_to_frac(X::ProjectiveScheme) 
  if !has_attribute(X, :homog_to_frac)
    affine_cone(X)
  end
  return get_attribute(X, :homog_to_frac)
end

@Markdown.doc """
    poly_to_homog(X::ProjectiveScheme)

Returns a map that converts an element of the `base_ring` of 
ring of functions `OO` of the `affine_cone` of `X` into 
an element of the `homog_poly_ring` of `X`.
"""
function poly_to_homog(X::ProjectiveScheme)
  if !has_attribute(X, :poly_to_homog)
    affine_cone(X)
  end
  return get_attribute(X, :poly_to_homog)
end

@Markdown.doc """
    function frac_to_homog_pair(X::ProjectiveScheme)

Returns a map that converts an element ``f = p/q`` of the ring of 
functions `OO` of the `affine_cone` of `X` into a pair 
``(a, b)`` of elements of the `homog_poly_ring` of `X`
corresponding to ``p`` and ``q``, respectively.
"""
function frac_to_homog_pair(X::ProjectiveScheme)
  if !has_attribute(X, :frac_to_homog_pair)
    affine_cone(X)
  end
  return get_attribute(X, :frac_to_homog_pair)
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
    F = affine_space(kk, symbols(homog_poly_ring(X)))
    C, pr_fiber, pr_base = product(F, Y)
    X.homog_coord = lift.([pullback(pr_fiber)(u) for u in gens(OO(F))])

    S = homog_poly_ring(X)
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

function affine_cone(X::ProjectiveScheme{CRT, CRET, RT, RET}) where {CRT<:MPolyQuoLocalizedRing, CRET, RT, RET}
  if !isdefined(X, :C)
    A = base_ring(X)
    Y = base_scheme(X)
    R = base_ring(A)
    kk = coefficient_ring(R)
    F = affine_space(kk, symbols(homog_poly_ring(X)))
    C, pr_fiber, pr_base = product(F, Y)
    X.homog_coord = lift.([pullback(pr_fiber)(u) for u in gens(OO(F))])
    S = homog_poly_ring(X)

    # store the various conversion maps
    help_map = hom(S, OO(C), 
                   (x -> pullback(pr_base)(x)),
                   [pullback(pr_fiber)(y) for y in gens(OO(F))]
                  )

    I = help_map(defining_ideal(X))
    CX = subscheme(C, I)
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
    
function affine_cone(X::ProjectiveScheme{CRT, CRET, RT, RET}) where {CRT<:AbstractAlgebra.Ring, CRET, RT, RET}
  if !isdefined(X, :C)
    kk = base_ring(X)
    C = affine_space(kk, symbols(homog_poly_ring(X)))
    X.homog_coord = gens(OO(C))
    S = homog_poly_ring(X)
    help_map = hom(S, OO(C), gens(OO(C)))
    I = help_map(defining_ideal(X))
    CX = subscheme(C, I)

    # store the various conversion maps
    set_attribute!(X, :homog_to_frac, hom(S, OO(CX), gens(OO(CX))))
    pth = hom(base_ring(OO(CX)), S, gens(S))
    set_attribute!(X, :poly_to_homog, pth)
    set_attribute!(X, :frac_to_homog_pair, (f -> (pth(lifted_numerator(OO(CX)(f))), pth(lifted_numerator(OO(CX)(f))))))
    X.C = CX
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
    PullbackType
  }
  domain::DomainType
  codomain::CodomainType
  pullback::PullbackType

  #fields for caching
  map_on_base_schemes::SpecMor
  map_on_affine_cones::SpecMor

  function ProjectiveSchemeMor(
      P::DomainType,
      Q::CodomainType,
      f::PullbackType;
      check::Bool=true
    ) where {DomainType<:ProjectiveScheme, CodomainType<:ProjectiveScheme, PullbackType<:Map}
    T = homog_poly_ring(P)
    S = homog_poly_ring(Q)
    (S === domain(f) && T === codomain(f)) || error("pullback map incompatible")
    if check
      #TODO: Check map on ideals (not available yet)
    end
    return new{DomainType, CodomainType, PullbackType}(P, Q, f)
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
  Q = homog_poly_ring(X)
  P = homog_poly_ring(Y)
  return ProjectiveSchemeMor(X, Y, hom(P, Q, a))
end

# in case we have honest base schemes, also make the map of schemes available
function base_map(phi::ProjectiveSchemeMor{<:ProjectiveScheme{<:MPolyQuoLocalizedRing}})
  if !isdefined(phi, :map_on_base_schemes)
    phi.map_on_base_schemes = SpecMor(base_scheme(domain(phi)), base_scheme(codomain(phi)), coefficient_map(pullback(phi)))
  end
  return phi.map_on_base_schemes::morphism_type(affine_patch_type(domain(phi)), affine_patch_type(codomain(phi)))
end

function map_on_affine_cones(phi::ProjectiveSchemeMor{<:ProjectiveScheme{<:MPolyQuoLocalizedRing}}) 
  if !isdefined(phi, :map_on_affine_cones)
    A = base_ring(domain(phi))
    S = homog_poly_ring(codomain(phi))
    T = homog_poly_ring(domain(phi))
    P = domain(phi)
    Q = codomain(phi)
    pb_P = pullback(projection_to_base(P))
    pb_Q = pullback(projection_to_base(Q))
    imgs_base = pb_P.(gens(A))
    imgs_fiber = [homog_to_frac(P)(g) for g in pullback(phi).(gens(S))]
    phi.map_on_affine_cones = SpecMor(affine_cone(P), affine_cone(Q), vcat(imgs_fiber, imgs_base))
  end
  return phi.map_on_affine_cones
end

function map_on_affine_cones(phi::ProjectiveSchemeMor{<:ProjectiveScheme{<:MPolyRing}})
  if !isdefined(phi, :map_on_affine_cones)
    Y = base_scheme(domain(phi))
    A = OO(Y)
    S = homog_poly_ring(codomain(phi))
    T = homog_poly_ring(domain(phi))
    P = domain(phi)
    Q = codomain(phi)
    pb_P = pullback(projection_to_base(P))
    pb_Q = pullback(projection_to_base(Q))
    imgs_base = pb_P.(gens(A))
    imgs_fiber = [homog_to_frac(P)(g) for g in pullback(phi).(gens(S))]
    phi.map_on_affine_cones = SpecMor(affine_cone(P), affine_cone(Q), vcat(imgs_fiber, imgs_base))
  end
  return phi.map_on_affine_cones
end
    
function map_on_affine_cones(phi::ProjectiveSchemeMor{<:ProjectiveScheme{<:AbstractAlgebra.Ring}})
  if !isdefined(phi, :map_on_affine_cones)
    S = homog_poly_ring(codomain(phi))
    T = homog_poly_ring(domain(phi))
    P = domain(phi)
    Q = codomain(phi)
    imgs_fiber = [homog_to_frac(P)(g) for g in pullback(phi).(gens(S))]
    phi.map_on_affine_cones = SpecMor(affine_cone(P), affine_cone(Q), imgs_fiber)
  end
  return phi.map_on_affine_cones
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
  for s in gens(homog_poly_ring(codomain(f)))
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

function fiber_product(f::SpecMor, P::ProjectiveScheme{<:MPolyQuoLocalizedRing})
  codomain(f) == base_scheme(P) || error("codomain and base_scheme are incompatible")
  X = domain(f)
  Y = codomain(f)
  Q_ambient = projective_space(X, symbols(homog_poly_ring(P)))
  help_map = hom(homog_poly_ring(P), 
                 homog_poly_ring(Q_ambient),
                 pullback(f),
                 gens(homog_poly_ring(Q_ambient))
                )
  I = help_map(defining_ideal(P))
  Q = subscheme(Q_ambient, I)
  return Q, ProjectiveSchemeMor(Q, P, 
                                hom(homog_poly_ring(P), 
                                    homog_poly_ring(Q),
                                    pullback(f),
                                    gens(homog_poly_ring(Q))
                                   )
                               )
end

fiber_product(X::Spec, P::ProjectiveScheme{<:MPolyQuoLocalizedRing}) = fiber_product(inclusion_map(X, base_scheme(P)), P)

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
                             hom(homog_poly_ring(Q),
                                 homog_poly_ring(P),
                                 pullback(f), 
                                 gens(homog_poly_ring(P))
                                )
                            )
end

function inclusion_map(P::T, Q::T) where {T<:ProjectiveScheme{<:AbstractAlgebra.Ring}}
  A = base_ring(Q)
  B = base_ring(P)
  A === B || error("can not compare schemes for non-equal base rings") # TODO: Extend by check for canonical maps, once they are available
  return ProjectiveSchemeMor(P, Q, 
                             hom(homog_poly_ring(Q),
                                 homog_poly_ring(P),
                                 gens(homog_poly_ring(P))
                                )
                            )
end

function as_covered_scheme(P::ProjectiveScheme)
  if !has_attribute(P, :as_covered_scheme) 
    C = standard_covering(P) 
    X = CoveredScheme(C)
    set_attribute!(P, :as_covered_scheme, X)
  end
  return get_attribute(P, :as_covered_scheme)
end

function covered_projection_to_base(X::ProjectiveScheme{<:MPolyQuoLocalizedRing})
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
  S = homog_poly_ring(X)
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
    U::SpecType
  ) where {
    CRT<:MPolyQuoLocalizedRing,
    SpecType<:Spec
  }
  # look up U in the coverings of X
  cover_of_U, index_of_U = X[U]
  Xcov = as_covered_scheme(X)
  S = homog_poly_ring(X)

  s = Vector{elem_type(OO(U))}()
  if cover_of_U === standard_covering(X)
    S = homog_poly_ring(X)
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
  S = homog_poly_ring(X)
  C = standard_covering(X)
  U = C[i+1]
  s = vcat(gens(OO(U))[1:i], [one(OO(U))], gens(OO(U))[i+1:fiber_dimension(X)])
  return hom(S, OO(U), s)
end

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
