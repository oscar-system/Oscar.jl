import Oscar.AlgHom

export ProjectiveScheme, base_ring, fiber_dimension, homogeneous_coordinate_ring, gens, getindex, affine_patch_type
export projective_scheme_type, affine_patch_type
export projective_space, subscheme
export projection_to_base, affine_cone, set_base_scheme!, base_scheme, homogeneous_coordinates, frac_to_homog, homog_to_frac, as_covered_scheme, covered_projection_to_base, dehomogenize
export ProjectiveSchemeMor, domain, codomain, images_of_variables, map_on_affine_cones, is_well_defined

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
    @show I
    @show typeof(I)
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

### type getters
#affine_patch_type(X::ProjectiveScheme{R, S, T, U}) where {R<:AbstractAlgebra.Field, S, T, U} = 

### type constructors 

# the type of a relative projective scheme over a given base scheme
projective_scheme_type(X::Spec{RT, RET, BRT, BRET, MST}) where {RT, RET, BRT, BRET, MST} = ProjectiveScheme{MPolyQuoLocalizedRing{RT, RET, BRT, BRET, MST}, MPolyQuoLocalizedRingElem{RT, RET, BRT, BRET, MST}, BRT, BRET}
projective_scheme_type(::Type{Spec{RT, RET, BRT, BRET, MST}}) where {RT, RET, BRT, BRET, MST} = ProjectiveScheme{MPolyQuoLocalizedRing{RT, RET, BRT, BRET, MST}, MPolyQuoLocalizedRingElem{RT, RET, BRT, BRET, MST}, BRT, BRET}

# the type of the affine cone for a projective scheme
affine_cone_type(P::ProjectiveScheme{CRT}) where {CRT<:AbstractAlgebra.Field} = spec_type(CRT)
affine_cone_type(::Type{ProjectiveScheme{CRT}}) where {CRT<:AbstractAlgebra.Field} = spec_type(CRT)

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
    defining_ideal(X::ProjectiveScheme)

On ``X ⊂ ℙʳ(A)`` this returns the homogeneous 
ideal ``I ⊂ A[s₀,…,sᵣ]`` defining ``X``.
"""
defining_ideal(X::ProjectiveScheme) = X.I

function Base.show(io::IO, P::ProjectiveScheme) 
  print(io, "subscheme of ℙ^$(fiber_dimension(P))_{$(base_ring(P))} defined as the zero locus of  $(defining_ideal(P))")
end

original_ring(S::MPolyRing_dec) = S.R

function affine_patch_type(X::ProjectiveScheme{CRT, CRET, RT, RET}) where {CRT<:AbstractAlgebra.Field, CRET, RT, RET}
  return Spec{typeof(base_ring(X)), 
              elem_type(base_ring(X)), 
              typeof(original_ring(homogeneous_coordinate_ring(X))), 
              elem_type(original_ring(homogeneous_coordinate_ring(X))), 
              MPolyPowersOfElement{
                                   typeof(base_ring(X)), 
                                   elem_type(base_ring(X)), 
                                   typeof(original_ring(homogeneous_coordinate_ring(X))), 
                                   elem_type(original_ring(homogeneous_coordinate_ring(X)))
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
                                   typeof(original_ring(homogeneous_coordinate_ring(X))), 
                                   elem_type(original_ring(homogeneous_coordinate_ring(X)))
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
  S = homogeneous_coordinate_ring(P)
  parent(f) == S || error("ring element does not belong to the correct ring")
  return ProjectiveScheme(S, defining_ideal(P) + ideal(S, [f]))
end

function subscheme(P::ProjectiveScheme, f::Vector{RingElemType}) where {RingElemType<:MPolyElem_dec}
  S = homogeneous_coordinate_ring(P)
  length(f) == 1 && return P #TODO: Replace P by an honest copy!
  for i in 1:length(f)
    parent(f) == S || error("ring element does not belong to the correct ring")
  end
  return ProjectiveScheme(S, defining_ideal(P) + ideal(S, f))
end

function subscheme(P::ProjectiveScheme, I::MPolyIdeal{T}) where {T<:RingElem}
  S = homogeneous_coordinate_ring(P)
  base_ring(I) == S || error("ideal does not belong to the correct ring")
  return ProjectiveScheme(S, I + defining_ideal(P))
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
  P = projective_space(OO(W), var_name)
  set_base_scheme!(P, W)
  return P
end

function projective_space(W::Spec, var_names::Vector{String}) 
  P = projective_space(OO(W), var_name)
  set_base_scheme!(P, W)
  return P
end

@Markdown.doc """
    homog_to_frac(X::ProjectiveScheme{CRT, CRET, RT, RET}, f::RET) where {CRT<:MPolyQuoLocalizedRing, CRET, RT, RET}

Convert a homogeneous polynomial ``f`` to an element of the ring of 
regular functions on the `affine_cone` of ``X``.
"""
function homog_to_frac(X::ProjectiveScheme{CRT, CRET, RT, RET}, f::RET) where {CRT<:MPolyQuoLocalizedRing, CRET, RT, RET}
  return evaluate(map_coefficients(pullback(projection_to_base(X)), f), homogeneous_coordinates(X))
end

function homog_to_frac(X::ProjectiveScheme{CRT, CRET, RT, RET}, f::RET) where {CRT<:MPolyRing, CRET, RT, RET}
  return evaluate(map_coefficients(pullback(projection_to_base(X)), f), homogeneous_coordinates(X))
end

function homog_to_frac(X::ProjectiveScheme{CRT, CRET, RT, RET}, f::RET) where {CRT<:AbstractAlgebra.Field, CRET, RT, RET}
  return evaluate(f, homogeneous_coordinates(X))
end

function homog_to_frac(X::ProjectiveScheme) 
  if !has_attribute(X, :homog_to_frac)
    affine_cone(X)
  end
  #TODO: insert type assertion here!
  return get_attribute(X, :homog_to_frac)
end

@Markdown.doc """
    frac_to_homog(X::ProjectiveScheme, f::T) where {T<:RingElem}

Convert a regular function ``f = a/b`` on some open subset of the affine 
cone of ``X`` to a pair of homogeneous polynomials ``(p, q)``, lifting 
``a`` and ``b``, respectively.
"""
function frac_to_homog(X::ProjectiveScheme{CRT, CRET, RT, RET}, f::CRET) where {CRT<:MPolyRing, CRET, RT, RET}
  S = homogeneous_coordinate_ring(X)
  A = base_ring(S)
  f_nested = zero(S)
  if parent(f) == A
    f_nested = S(f)
  else
    f_nested = renest(S, f)
  end
  return f_nested
end

function frac_to_homog(X::ProjectiveScheme{CRT, CRET, RT, RET}, f::MPolyLocalizedRingElem) where {CRT<:MPolyRing, CRET, RT, RET}
  S = homogeneous_coordinate_ring(X)
  A = base_ring(S)
  a = numerator(f)
  b = denominator(f)
  a_nested = renest(S, a)
  b_nested = renest(S, b)
  return (a_nested, b_nested)
end

function frac_to_homog(X::ProjectiveScheme{CRT, CRET, RT, RET}, f::MPolyQuoLocalizedRingElem) where {CRT<:MPolyRing, CRET, RT, RET}
  return frac_to_homog(X, lift(f))
end

function frac_to_homog(X::ProjectiveScheme{CRT, CRET, RT, RET}, f::Vector{CRET}) where {CRT<:MPolyRing, CRET, RT, RET}
  length(f) == 0 && return elem_type(S)[]
  S = homogeneous_coordinate_ring(X)
  A = base_ring(S)
  f_nested = elem_type(S)[]
  if parent(f[1]) == A
    f_nested = S.(f)
  else
    f_nested = [renest(S, a) for a in f]
  end
  return f_nested
end

function frac_to_homog(X::ProjectiveScheme{CRT, CRET, RT, RET}, f::Vector{MPolyLocalizedRingElem}) where {CRT<:MPolyRing, CRET, RT, RET}
  S = homogeneous_coordinate_ring(X)
  A = base_ring(S)
  a = numerator.(f)
  b = denominator.(f)
  a_nested = [renest(S, p) for p in a]
  b_nested = [renest(S, q) for q in b]
  return [(a_nested[i], b_nested[i]) for i in 1:length(a_nested)]
end

function frac_to_homog(X::ProjectiveScheme{CRT, CRET, RT, RET}, f::Vector{MPolyQuoLocalizedRingElem}) where {CRT<:MPolyRing, CRET, RT, RET}
  return frac_to_homog(X, lift.(f))
end

### projective space over an honest base scheme
function frac_to_homog(X::ProjectiveScheme{CRT, CRET, RT, RET}, f::PolyType) where {CRT<:MPolyQuoLocalizedRing, PolyType<:MPolyElem, CRET, RT, RET}
  S = homogeneous_coordinate_ring(X)
  A = base_ring(S)
  R = base_ring(A)
  T, _ = PolynomialRing(R, symbols(S))
  f_nested = zero(T)
  if parent(f) == R
    f_nested = T(f)
  else
    f_nested = renest(T, f)
  end
  function coeff_map(a::T) where {T<:MPolyElem} 
    return A(evaluate(a, gens(base_ring(A))))
  end
  f_evaluated = evaluate(map_coefficients(coeff_map, f_nested), gens(S))
  return f_evaluated
end

function frac_to_homog(X::ProjectiveScheme{CRT, CRET, RT, RET}, f::MPolyLocalizedRingElem) where {CRT<:MPolyQuoLocalizedRing, CRET, RT, RET}
  S = homogeneous_coordinate_ring(X)
  A = base_ring(S)
  R = base_ring(A)
  T, _ = PolynomialRing(R, symbols(S))
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

function frac_to_homog(X::ProjectiveScheme{CRT, CRET, RT, RET}, f::MPolyQuoLocalizedRingElem) where {CRT<:MPolyQuoLocalizedRing, CRET, RT, RET}
  return frac_to_homog(X, lift(f))
end
    
function frac_to_homog(X::ProjectiveScheme{CRT, CRET, RT, RET}, f::Vector{PolyType}) where {CRT<:MPolyQuoLocalizedRing, PolyType<:MPolyElem, CRET, RT, RET}
  length(f) == 0 && return elem_type(S)[]
  S = homogeneous_coordinate_ring(X)
  A = base_ring(S)
  R = base_ring(A)
  T, _ = PolynomialRing(R, symbols(S))
  if parent(f[1]) == R
    f_nested = [T(a) for a in f]
  else
    f_nested = [renest(T, a) for a in f]
  end
  function coeff_map(a::T) where {T<:MPolyElem} 
    return A(evaluate(a, gens(base_ring(A))))
  end
  f_evaluated = [evaluate(map_coefficients(coeff_map, a), gens(S)) for a in f_nested]
  return f_evaluated
end

function frac_to_homog(X::ProjectiveScheme{CRT, CRET, RT, RET}, f::Vector{FractionType}) where {CRT<:MPolyQuoLocalizedRing, FractionType<:MPolyLocalizedRingElem, CRET, RT, RET}
  S = homogeneous_coordinate_ring(X)
  A = base_ring(S)
  R = base_ring(A)
  T, _ = PolynomialRing(R, symbols(S))
  a = numerator.(f)
  b = denominator.(f)
  a_nested = [renest(S, p) for p in a]
  b_nested = [renest(S, q) for q in b]
  function coeff_map(a::T) where {T<:MPolyElem} 
    return A(evaluate(a, gens(base_ring(A))))
  end
  a_evaluated = [evaluate(map_coefficients(coeff_map, p), gens(S)) for p in a_nested]
  b_evaluated = [evaluate(map_coefficients(coeff_map, q), gens(S)) for q in b_nested]
  return [(a_evaluated[i], b_evaluated[i]) for i in 1:length(a_evaluated)]
end

function frac_to_homog(X::ProjectiveScheme{CRT, CRET, RT, RET}, f::Vector{T}) where {CRT<:MPolyQuoLocalizedRing, T<:MPolyQuoLocalizedRingElem, CRET, RT, RET}
  return frac_to_homog(X, lift.(f))
end
    
### projective space over a field
function frac_to_homog(X::ProjectiveScheme{CRT, CRET, RT, RET}, f::MPolyLocalizedRingElem) where {CRT<:AbstractAlgebra.Field, CRET, RT, RET}
  S = homogeneous_coordinate_ring(X)
  return (evaluate(numerator(f), gens(S)), evaluate(denominator(f), gens(S)))
end

function frac_to_homog(X::ProjectiveScheme{CRT, CRET, RT, RET}, f::MPolyQuoLocalizedRingElem) where {CRT<:AbstractAlgebra.Field, CRET, RT, RET}
  return frac_to_homog(X, lift(f))
end
    
function frac_to_homog(X::ProjectiveScheme{CRT, CRET, RT, RET}, f::Vector{T}) where {CRT<:AbstractAlgebra.Field, T<:MPolyLocalizedRingElem, CRET, RT, RET}
  S = homogeneous_coordinate_ring(X)
  return [(evaluate(numerator(g), gens(S)), evaluate(denominator(g), gens(S))) for g in f]
end

function frac_to_homog(X::ProjectiveScheme{CRT, CRET, RT, RET}, f::Vector{T}) where {CRT<:AbstractAlgebra.Field, T<:MPolyQuoLocalizedRingElem, CRET, RT, RET}
  return frac_to_homog(X, lift.(f))
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

    S = homogeneous_coordinate_ring(X)
    inner_help_map = hom(A, OO(C), [pr_base(x) for x in gens(OO(Y))])
    help_map = hom(S, OO(C), inner_help_map, [pr_fiber(y) for y in gens(OO(F))])
    I = ideal(OO(C), [help_map(g) for g in gens(defining_ideal(X))])
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

    I = ideal(OO(C), [complicated_evaluation(g) for g in gens(defining_ideal(X))])
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
    I = ideal(OO(C), [homog_to_frac(X, g) for g in gens(defining_ideal(X))])
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
    BasePullbackType, 
    RingElemType<:MPolyElem
  }
  domain::DomainType
  codomain::CodomainType
  images_of_variables::Vector{RingElemType}
  pullback_on_base::BasePullbackType

  #fields for caching
  map_on_base_schemes::SpecMor
  map_on_affine_cones::SpecMor

  # constructor for morphisms over the same base
  function ProjectiveSchemeMor(
      P::DomainType,
      Q::CodomainType,
      imgs::Vector{RET};
      check::Bool=true
    ) where {DomainType, CodomainType, BasePullbackType, RET}
    A = base_ring(P)
    A == base_ring(Q) || error("projective schemes are not defined over the same ring")
    m = fiber_dimension(P)
    n = fiber_dimension(Q)
    length(imgs) == n+1 || error("the number of images does not agree with the number of variables")
    for i in 1:n+1
      parent(imgs[i]) == homogeneous_coordinate_ring(P) || error("elements do not belong to the correct ring")
    end
    return new{DomainType, CodomainType, Any, RET}(P, Q, imgs)
  end

  # constructor for morphisms with base change
  function ProjectiveSchemeMor(
      P::DomainType,
      Q::CodomainType,
      f::BasePullbackType,
      imgs::Vector{RET};
      check::Bool=true
    ) where {DomainType, CodomainType, BasePullbackType, RET}
    B = base_ring(P)
    A = base_ring(Q)
    domain(f) == B || error("pullback is not compatible with the given base rings")
    codomain(f) == A || error("pullback is not compatible with the given base rings")
    m = fiber_dimension(P)
    n = fiber_dimension(Q)
    length(imgs) == n+1 || error("the number of images does not agree with the number of variables")
    for i in 1:n+1
      parent(imgs[i]) == homogeneous_coordinate_ring(P) || error("elements do not belong to the correct ring")
    end
    return new{DomainType, CodomainType, BasePullbackType, RET}(P, Q, imgs, f)
  end
end

domain(phi::ProjectiveSchemeMor) = phi.domain
codomain(phi::ProjectiveSchemeMor) = phi.codomain
images_of_variables(phi::ProjectiveSchemeMor) = phi.images_of_variables
pullback_on_base(phi::ProjectiveSchemeMor) = phi.pullback_on_base

# in case we have honest base schemes, also make the map of schemes available
function base_map(phi::ProjectiveSchemeMor{ProjectiveScheme{CRT}}) where {CRT<:MPolyQuoLocalizedRing}
  if !isdefined(phi, :map_on_base_schemes)
    phi.map_on_base_schemes = SpecMor(base_scheme(domain(phi)), base_scheme(codomain(phi)), pullback_on_base)
  end
  return phi.map_on_base_schemes::morphism_type(affine_patch_type(domain(phi)), affine_patch_type(codomain(phi)))
end

### dispatch for morphisms over the same base scheme/ring
function map_on_affine_cones(phi::ProjectiveSchemeMor{ProjectiveScheme{CRT}}) where {CRT<:MPolyQuoLocalizedRing}
  if !isdefined(phi, :map_on_affine_cones)
    A = base_ring(domain(phi))
    S = homogeneous_coordinate_ring(codomain(phi))
    T = homogeneous_coordinate_ring(domain(phi))
    P = domain(phi)
    Q = codomain(phi)
    pb_P = pullback(projection_to_base(P))
    pb_Q = pullback(projection_to_base(Q))
    imgs_base = pb_P.(gens(A))
    imgs_fiber = [homog_to_frac(P, g) for g in images_of_variables(phi)]
    phi.map_on_affine_cones = SpecMor(affine_cone(P), affine_cone(Q), vcat(imgs_fiber, imgs_base))
  end
  return phi.map_on_affine_cones
end

function map_on_affine_cones(phi::ProjectiveSchemeMor{ProjectiveScheme{CRT}}) where {CRT<:MPolyRing}
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
    imgs_fiber = [homog_to_frac(P, g) for g in images_of_variables(phi)]
    phi.map_on_affine_cones = SpecMor(affine_cone(P), affine_cone(Q), vcat(imgs_fiber, imgs_base))
  end
  return phi.map_on_affine_cones
end
    
function map_on_affine_cones(phi::ProjectiveSchemeMor{ProjectiveScheme{CRT}}) where {CRT<:AbstractAlgebra.Field}
  if !isdefined(phi, :map_on_affine_cones)
    S = homogeneous_coordinate_ring(codomain(phi))
    T = homogeneous_coordinate_ring(domain(phi))
    P = domain(phi)
    Q = codomain(phi)
    imgs_fiber = [homog_to_frac(P, g) for g in images_of_variables(phi)]
    phi.map_on_affine_cones = SpecMor(affine_cone(P), affine_cone(Q), imgs_fiber)
  end
  return phi.map_on_affine_cones
end


### dispatch for morphisms over a given morphism of bases
function map_on_affine_cones(
    phi::ProjectiveSchemeMor{ProjectiveScheme{CRT}, CDT, BMT}
  ) where {
    CRT<:MPolyQuoLocalizedRing,
    CDT,
    BMT<:AbstractAlgebra.Map
  }
  if !isdefined(phi, :map_on_affine_cones)
    A = base_ring(domain(phi))
    S = homogeneous_coordinate_ring(codomain(phi))
    T = homogeneous_coordinate_ring(domain(phi))
    P = domain(phi)
    Q = codomain(phi)
    pb_Q = pullback(projection_to_base(Q))
    # also employ the morphism of the base implicitly
    imgs_base = pb_Q.(pullback_on_base(phi).(gens(A)))
    imgs_fiber = [homog_to_frac(P, g) for g in images_of_variables(phi)]
    phi.map_on_affine_cones = SpecMor(affine_cone(P), affine_cone(Q), vcat(imgs_fiber, imgs_base))
  end
  return phi.map_on_affine_cones::morphism_type{affine_patch_type(domain(phi)), affine_patch_type(codomain(phi))}
end

function map_on_affine_cones(
    phi::ProjectiveSchemeMor{ProjectiveScheme{CRT}, CDT, BMT}
  ) where {
    CRT<:MPolyRing,
    CDT,
    BMT<:AbstractAlgebra.Map
  }
  if !isdefined(phi, :map_on_affine_cones)
    Y = base_scheme(domain(phi))
    A = OO(Y)
    S = homogeneous_coordinate_ring(codomain(phi))
    T = homogeneous_coordinate_ring(domain(phi))
    P = domain(phi)
    Q = codomain(phi)
    pb_Q = pullback(projection_to_base(Q))
    # also employ the morphism of the base implicitly
    imgs_base = pb_Q.(pullback_on_base(phi).(gens(A)))
    imgs_fiber = [homog_to_frac(P, g) for g in images_of_variables(phi)]
    phi.map_on_affine_cones = SpecMor(affine_cone(P), affine_cone(Q), vcat(imgs_fiber, imgs_base))
  end
  return phi.map_on_affine_cones
end
    
function map_on_affine_cones(
    phi::ProjectiveSchemeMor{ProjectiveScheme{CRT}, CDT, BMT}
  ) where {
    CRT<:Field,
    CDT,
    BMT<:AbstractAlgebra.Map
  }
  error("not implemented")
end

function is_well_defined(phi::ProjectiveSchemeMor) 
  CP = affine_cone(domain(phi))
  CQ = affine_cone(codomain(phi))
  return issubset(CP, preimage(map_on_affine_cones(phi), CQ))
end

function as_covered_scheme(P::ProjectiveScheme)
  if has_attribute(P, :as_covered_scheme) 
    return get_attribute(P, :as_covered_scheme)
  end
  C = standard_covering(P) 
  X = CoveredScheme(C)
  set_attribute!(P, :as_covered_scheme, X)
  return X
end

function covered_projection_to_base(X::ProjectiveSchemeMor{ProjectiveScheme{CRT}}) where {CRT<:MPolyQuoLocalizedRing}
  if has_attribute(X, :covered_projection_to_base) 
    return get_attribute(X, :covered_projection_to_base)
  end
  C = standard_covering(X)
  return get_attribute(X, :covered_projection_to_base) # TODO: establish type assertion here!
end


function dehomogenize(
    X::ProjectiveScheme{CRT}, 
    f::RingElemType, 
    i::Int
  ) where {
    CRT<:MPolyQuoLocalizedRing, 
    RingElemType<:MPolyElem_dec
  }
  i in 0:fiber_dimension(X) || error("the given integer is not in the admissible range")
  S = homogeneous_coordinate_ring(X)
  parent(f) == S || error("the given polynomial does not have the correct parent")
  C = standard_covering(X)
  U = C[i+1]
  p = covered_projection_to_base(X)
  s = vcat(gens(OO(U))[1:i], [one(OO(U))], gens(OO(U))[i+1:fiber_dimension(X)])
  return evaluate(map_coefficients(pullback(p[U]), f), s)
end

function dehomogenize(
    X::ProjectiveScheme{CRT}, 
    f::Vector{RingElemType}, 
    i::Int
  ) where {
    CRT<:MPolyQuoLocalizedRing, 
    RingElemType<:MPolyElem_dec
  }
  i in 0:fiber_dimension(X) || error("the given integer is not in the admissible range")
  S = homogeneous_coordinate_ring(X)
  for a in f 
    parent(a) == S || error("the given polynomial does not have the correct parent")
  end
  C = standard_covering(X)
  U = C[i+1]
  p = covered_projection_to_base(X)
  s = vcat(gens(OO(U))[1:i], [one(OO(U))], gens(OO(U))[i+1:fiber_dimension(X)])
  return [evaluate(map_coefficients(pullback(p[U]), a), s) for a in f]
end

# TODO: Write the other dispatch for the dehomogenization routines

function fiber_product(f::SpecMor, P::ProjectiveScheme) 
  X = domain(f)
  Y = codomain(f)
  Y == base_scheme(P) || error("map and scheme are not compatible")
  R, _ = PolynomialRing(OO(X), symbols(homogeneous_coordinate_ring(P)))
  S, s = grade(R, [1 for i in 1:nvars(R)])
  g = [map_coefficients(pullback(f), g, parent=S) for g in gens(defining_ideal(P))]
  Q = ProjectiveScheme(S, ideal(S, g))
  f_up = ProjectiveSchemeMor
  return
end

