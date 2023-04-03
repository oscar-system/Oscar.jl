export homogeneous_coordinate_ring

########################################################################
# Interface for abstract projective schemes                            #
########################################################################

function base_ring(P::AbsProjectiveScheme) 
  return base_ring(underlying_scheme(P))
end

@doc Markdown.doc"""
    ambient_coordinate_ring(P::AbsProjectiveScheme)

On a projective scheme ``P = Proj(S)`` with ``S = P/I`` 
for a standard graded polynomial ring ``P`` and a 
homogeneous ideal ``I`` this returns ``P``.
"""
function ambient_coordinate_ring(P::AbsProjectiveScheme)
  return ambient_coordinate_ring(underlying_scheme(P))
end

@doc Markdown.doc"""
    homogeneous_coordinate_ring(P::AbsProjectiveScheme)

On a projective scheme ``P = Proj(S)`` for a standard 
graded finitely generated algebra ``S`` this returns ``S``.
"""
function homogeneous_coordinate_ring(P::AbsProjectiveScheme)
  return homogeneous_coordinate_ring(underlying_scheme(P))
end

@attr AbsSpec function base_scheme(P::AbsProjectiveScheme)
  return base_scheme(underlying_scheme(P))
end

@doc Markdown.doc"""
    affine_cone(X::ProjectiveScheme) 

On ``X = Proj(S) ⊂ ℙʳ_𝕜`` this returns a pair `(C, f)` where ``C = C(X) ⊂ 𝕜ʳ⁺¹`` 
is the affine cone of ``X`` and ``f : S → 𝒪(C)`` is the morphism of rings 
from the `homogeneous_coordinate_ring` to the `coordinate_ring` of the affine cone.
"""
@attr function affine_cone(
    P::AbsProjectiveScheme{RT}
  ) where {RT<:Union{MPolyRing, MPolyQuoRing, MPolyQuoLocRing, MPolyLocRing}}
  S = homogeneous_coordinate_ring(P)
  phi = RingFlattening(S)
  A = codomain(phi)
  C = Spec(A)
  B = base_scheme(P)
  P.projection_to_base = SpecMor(C, B, hom(OO(B), OO(C), gens(OO(C))[ngens(S)+1:end], check=false), check=false)
  return C, phi
end

@attr function affine_cone(
    P::AbsProjectiveScheme{RT, <:MPolyQuoRing}
  ) where {RT<:Field}
  S = homogeneous_coordinate_ring(P)
  PS = base_ring(S)
  PP = forget_grading(PS) # the ungraded polynomial ring
  I = modulus(S)
  II = forget_grading(I)
  SS, _ = quo(PP, II)
  phi = hom(S, SS, gens(SS))
  C = Spec(SS)
  return C, phi
end

@attr function affine_cone(
    P::AbsProjectiveScheme{RT, <:MPolyDecRing}
  ) where {RT<:Field}
  S = homogeneous_coordinate_ring(P)
  PP = forget_grading(S) # the ungraded polynomial ring
  phi = hom(S, PP, gens(PP))
  C = Spec(PP)
  return C, phi
end

########################################################################
# Methods for the concrete minimal instance                            #
########################################################################

@doc Markdown.doc"""
    base_ring(X::ProjectiveScheme)

On ``X ⊂ ℙʳ_A`` this returns ``A``.
"""
base_ring(P::ProjectiveScheme) = P.A

@doc Markdown.doc"""
    base_scheme(X::ProjectiveScheme{CRT, CRET, RT, RET}) where {CRT<:MPolyQuoLocRing, CRET, RT, RET}

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
    X::Union{<:AbsSpec, <:SpecOpen}
  ) where {CRT<:Ring, CRET, RT, RET}
  OO(X) === base_ring(P) || error("schemes are not compatible")
  P.Y = X
  return P
end

function projection_to_base(X::ProjectiveScheme{CRT, CRET, RT, RET}) where {CRT<:Union{<:MPolyRing, <:MPolyQuoRing, <:MPolyLocRing, <:MPolyQuoLocRing, <:SpecOpenRing}, CRET, RT, RET}
  if !isdefined(X, :projection_to_base)
    affine_cone(X)
  end
  return X.projection_to_base
end


@doc Markdown.doc"""
    relative_ambient_dimension(X::ProjectiveScheme)

On ``X ⊂ ℙʳ_A`` this returns ``r``.
"""
relative_ambient_dimension(P::ProjectiveScheme) = P.r

@doc Markdown.doc"""
    homogeneous_coordinate_ring(X::ProjectiveScheme)

On ``X ⊂ ℙʳ_A`` this returns ``A[s₀,…,sᵣ]``.
"""
homogeneous_coordinate_ring(P::ProjectiveScheme) = P.S

ambient_coordinate_ring(P::ProjectiveScheme{<:Any, <:Any, <:MPolyQuoRing}) = base_ring(homogeneous_coordinate_ring(P))
ambient_coordinate_ring(P::ProjectiveScheme{<:Any, <:Any, <:MPolyDecRing}) = homogeneous_coordinate_ring(P)

@doc Markdown.doc"""
    homogeneous_coordinates(X::ProjectiveScheme)

Return the generators of the homogeneous coordinate ring of ``X``.
"""
function homogeneous_coordinates(X::ProjectiveScheme)
  return gens(homogeneous_coordinate_ring(X))
end

### TODO: Replace by the map of generators.
@doc Markdown.doc"""
    homogeneous_coordinates_on_affine_cone(X::ProjectiveScheme)

On ``X ⊂ ℙʳ_A`` this returns a vector with the homogeneous
coordinates ``[s₀,…,sᵣ]`` as entries where each one of the 
``sᵢ`` is a function on the `affine cone` of ``X``.
"""
function homogeneous_coordinates_on_affine_cone(P::ProjectiveScheme)
  if !isdefined(P, :homog_coord)
    C, f = affine_cone(P)
    P.homog_coord = f.(gens(homogeneous_coordinate_ring(P)))
  end
  return P.homog_coord
end

homogeneous_coordinate_on_affine_cone(P::ProjectiveScheme, i::Int) = homogeneous_coordinates_on_affine_cone(P)[i]

@doc Markdown.doc"""
    defining_ideal(X::AbsProjectiveScheme)

On ``X ⊂ ℙʳ_A`` this returns the homogeneous
ideal ``I ⊂ A[s₀,…,sᵣ]`` defining ``X``.
"""
defining_ideal(X::AbsProjectiveScheme{<:Any, <:MPolyDecRing}) = ideal(homogeneous_coordinate_ring(X), Vector{elem_type(homogeneous_coordinate_ring(X))}())
defining_ideal(X::AbsProjectiveScheme{<:Any, <:MPolyQuoRing}) = modulus(homogeneous_coordinate_ring(X))

### type getters 
projective_scheme_type(A::T) where {T<:AbstractAlgebra.Ring} = projective_scheme_type(typeof(A))
projective_scheme_type(::Type{T}) where {T<:AbstractAlgebra.Ring} = 
ProjectiveScheme{T, elem_type(T), mpoly_dec_ring_type(mpoly_ring_type(T)), mpoly_dec_type(mpoly_ring_type(T))}

base_ring_type(P::ProjectiveScheme) = base_ring_type(typeof(P))
base_ring_type(::Type{ProjectiveScheme{S, T, U, V}}) where {S, T, U, V} = S

ring_type(P::ProjectiveScheme) = ring_type(typeof(P))
ring_type(::Type{ProjectiveScheme{S, T, U, V}}) where {S, T, U, V} = U

### type constructors 

# the type of a relative projective scheme over a given base scheme
projective_scheme_type(X::AbsSpec) = projective_scheme_type(typeof(X))
projective_scheme_type(::Type{T}) where {T<:AbsSpec} = projective_scheme_type(ring_type(T))

@doc Markdown.doc"""
    affine_cone(X::AbsProjectiveScheme) -> AbsSpec

Return the affine cone of `X`.
"""
@attr function affine_cone(
    X::AbsProjectiveScheme{CRT, RT}
  ) where {
           CRT<:SpecOpenRing, 
           RT<:MPolyRing 
          }
  S = ambient_coordinate_ring(X)
  B = coefficient_ring(S)
  Y = scheme(B)
  U = domain(B)
  R = base_ring(OO(Y))
  kk = base_ring(R)
  F = affine_space(kk, symbols(ambient_coordinate_ring(X)))
  C, pr_base, pr_fiber = product(U, F)
  X.homog_coord = [pullback(pr_fiber)(u) 
                   for u in OO(codomain(pr_fiber)).(gens(OO(F)))]
  phi = hom(S, OO(C), pullback(pr_base), X.homog_coord)
  g = phi.(gens(defining_ideal(X)))
  CX = subscheme(C, g)
  X.C = CX

  psi = compose(phi, restriction_map(C, CX))
  set_attribute!(X, :base_scheme, U)
  X.projection_to_base = restrict(pr_base, CX, U, check=false)
  return X.C, psi 
end

@doc Markdown.doc"""
    affine_cone(X::AbsProjectiveScheme{<:SpecOpenRing}) -> SpecOpen

Return the affine_cone of `X`.

Note that if the base scheme is not affine, then the affine cone is not affine.
"""
@attr function affine_cone(
    X::AbsProjectiveScheme{CRT, RT}
  ) where {
           CRT<:SpecOpenRing,
           RT<:MPolyQuoRing
          }
  P = ambient_coordinate_ring(X)
  S = homogeneous_coordinate_ring(X)
  B = coefficient_ring(P)
  Y = scheme(B)
  U = domain(B)
  R = base_ring(OO(Y))
  kk = base_ring(R)
  F = affine_space(kk, symbols(ambient_coordinate_ring(X)))
  C, pr_base, pr_fiber = product(U, F)
  homog_coord = [pullback(pr_fiber)(u) 
                 for u in OO(codomain(pr_fiber)).(gens(OO(F)))]
  phi = hom(P, OO(C), pullback(pr_base), homog_coord)
  g = phi.(gens(modulus(S)))
  CX = subscheme(C, g)
  pr_base_res = restrict(pr_base, CX, codomain(pr_base), check=true)
  X.C = CX
  X.homog_coord = OO(CX).(homog_coord)

  #psi = hom(S, OO(CX), pullback(pr_base), OO(CX).(X.homog_coord), check=false)

  psi = compose(phi, restriction_map(C, CX))
  psi_res = hom(S, OO(CX), pullback(pr_base_res), X.homog_coord, check=false)
  set_attribute!(X, :base_scheme, U)
  X.projection_to_base = restrict(pr_base, CX, U, check=false)
  return X.C, psi_res 
end

# Basic functionality required for Warham
@attr Int function dim(P::AbsProjectiveScheme{<:Field})
  return dim(defining_ideal(P))-1
end

@attr QQPolyRingElem function hilbert_polynomial(P::AbsProjectiveScheme{<:Field})
  return hilbert_polynomial(homogeneous_coordinate_ring(P))
end

@attr ZZRingElem function degree(P::AbsProjectiveScheme{<:Field})
  return degree(homogeneous_coordinate_ring(P))
end

@attr QQFieldElem function arithmetic_genus(P::AbsProjectiveScheme{<:Field})
  h = hilbert_polynomial(P)
  return (-1)^dim(P) * (first(coefficients(h)) - 1)
end

@attr Bool function is_smooth(P::AbsProjectiveScheme)
  return is_smooth(covered_scheme(P))
end

@doc Markdown.doc"""
    covered_scheme(P::ProjectiveScheme)
    
Return a `CoveredScheme` ``X`` isomorphic to `P` with standard affine charts given by dehomogenization. 

Use `dehomogenization_map` with `U` one of the `affine_charts` of ``X`` to
obtain the dehomogenization map from the `homogeneous_coordinate_ring` of `P`
to the `coordinate_ring` of `U`.

# Examples
```jldoctest
julia> P = projective_space(QQ, 2);

julia> Pcov = covered_scheme(P)
covered scheme with 3 affine patches in its default covering
```
"""
@attr AbsCoveredScheme function covered_scheme(P::ProjectiveScheme)
    C = standard_covering(P) 
    X = CoveredScheme(C)
    return X
end

@attr function covered_projection_to_base(X::ProjectiveScheme{<:Union{<:MPolyQuoLocRing, <:MPolyLocRing, <:MPolyQuoRing, <:MPolyRing}})
  if !has_attribute(X, :covering_projection_to_base) 
    C = standard_covering(X)
  end
  covering_projection = get_attribute(X, :covering_projection_to_base)::CoveringMorphism
  projection = CoveredSchemeMorphism(covered_scheme(X), CoveredScheme(codomain(covering_projection)), covering_projection)
end

