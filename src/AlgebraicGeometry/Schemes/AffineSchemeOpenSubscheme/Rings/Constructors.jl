
########################################################################
# Constructors for AffineSchemeOpenSubschemeRing                                        #
########################################################################

@doc raw"""
    OO(U::AffineSchemeOpenSubscheme) -> AffineSchemeOpenSubschemeRing

Given a Zariski open subset `U` of an affine scheme `X`, return the ring
`𝒪(X, U)` of regular functions on `U`.

# Examples
```jldoctest
julia> P, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> I = ideal([x^3-y^2*z]);

julia> A = spec(P)
Spectrum
  of multivariate polynomial ring in 3 variables x, y, z
    over rational field

julia> Y = spec(P, I)
Spectrum
  of quotient
    of multivariate polynomial ring in 3 variables x, y, z
      over rational field
    by ideal (x^3 - y^2*z)

julia> U = complement(A, Y)
Open subset
  of affine 3-space
complement to V(x^3 - y^2*z)

julia> OO(U)
Ring of regular functions
  on complement to V(x^3 - y^2*z) in affine scheme with coordinates [x, y, z]

julia> one(OO(U))
Regular function
  on open subset
    of affine scheme with coordinates [x, y, z]
  complement to V(x^3 - y^2*z)
  covered by 1 affine patch
    1: [x, y, z]   AA^3 \ scheme(x^3 - y^2*z)
with restriction
  patch 1: 1
```
"""
function OO(U::AffineSchemeOpenSubscheme)
  if !isdefined(U, :ring_of_functions)
    U.ring_of_functions = AffineSchemeOpenSubschemeRing(ambient_scheme(U), U, check=false)
  end
  return U.ring_of_functions::AffineSchemeOpenSubschemeRing
  #return U.ring_of_functions::AffineSchemeOpenSubschemeRing{affine_patch_type(U), typeof(U)}
end

########################################################################
# Constructors for AffineSchemeOpenSubschemeRingElem                                    #
########################################################################

########################################################################
# Coercion                                                             #
########################################################################
(R::AffineSchemeOpenSubschemeRing)(f::RingElem; check::Bool=true) = AffineSchemeOpenSubschemeRingElem(R, [OO(U)(f, check=check) for U in affine_patches(domain(R))])
(R::AffineSchemeOpenSubschemeRing)(f::MPolyQuoLocRingElem; check::Bool=true) = AffineSchemeOpenSubschemeRingElem(R, [_cast_fraction(OO(U),lifted_numerator(f), lifted_denominator(f), check=check) for U in affine_patches(domain(R))], check=false)

(R::AffineSchemeOpenSubschemeRing)(f::Vector{T}; check::Bool=true) where {T<:RingElem} = AffineSchemeOpenSubschemeRingElem(R, [OO(domain(R)[i])(f[i], check=check) for i in 1:length(f)])

function (R::AffineSchemeOpenSubschemeRing)(f::AffineSchemeOpenSubschemeRingElem; check::Bool=true)
  parent(f) === R && return f
  return AffineSchemeOpenSubschemeRingElem(R, [restrict(f, U, check=check) for U in affine_patches(domain(R))])
end

########################################################################
# Additional constructors for the Ring interface                       #
########################################################################

one(R::AffineSchemeOpenSubschemeRing) = AffineSchemeOpenSubschemeRingElem(R, [one(OO(U)) for U in affine_patches(domain(R))], check=false)
zero(R::AffineSchemeOpenSubschemeRing) = AffineSchemeOpenSubschemeRingElem(R, [zero(OO(U)) for U in affine_patches(domain(R))], check=false)
(R::AffineSchemeOpenSubschemeRing)() = zero(R)
(R::AffineSchemeOpenSubschemeRing)(a::Integer) = AffineSchemeOpenSubschemeRingElem(R, [OO(U)(a) for U in affine_patches(domain(R))], check=false)
(R::AffineSchemeOpenSubschemeRing)(a::Int64) = AffineSchemeOpenSubschemeRingElem(R, [OO(U)(a) for U in affine_patches(domain(R))], check=false)
(R::AffineSchemeOpenSubschemeRing)(a::ZZRingElem) = AffineSchemeOpenSubschemeRingElem(R, [OO(U)(a) for U in affine_patches(domain(R))], check=false)

########################################################################
# Copying                                                              #
########################################################################
function Base.deepcopy_internal(f::AffineSchemeOpenSubschemeRingElem, dict::IdDict)
  return AffineSchemeOpenSubschemeRingElem(parent(f), copy(restrictions(f)), check=false)
end

########################################################################
# Maximal extensions of rational functions on affine schemes           #
########################################################################
@doc raw"""
    maximal_extension(X::AffineScheme, f::AbstractAlgebra.Generic.FracFieldElem)

Return the maximal extension of the restriction of ``f``
to a rational function on ``X`` on a maximal domain of
definition ``U ⊂ X``.

**Note:** When ``X = Spec(R)`` with ``R = (𝕜[x₁,…,xₙ]/I)[S⁻¹]``,
the numerator and denominator of ``f`` have to be elements of
the ring ``𝕜[x₁,…,xₙ]``.
"""
function maximal_extension(
    X::AbsAffineScheme{<:Ring, <:MPolyLocRing},
    f::AbstractAlgebra.Generic.FracFieldElem{RET}
  ) where {RET<:MPolyRingElem}

  a = numerator(f)
  b = denominator(f)
  g = gcd(a, b)
  if !isone(g)
    a = divexact(a, g)
    b = divexact(b, g)
    f = parent(f)(a,b)
  end
  W = OO(X)
  U = AffineSchemeOpenSubscheme(X, [b])
  g = [OO(V)(f) for V in affine_patches(U)]
  R = AffineSchemeOpenSubschemeRing(X, U)
  return R(g)
end

function maximal_extension(
    X::AbsAffineScheme{<:Ring, <:MPolyQuoLocRing},
    f::AbstractAlgebra.Generic.FracFieldElem{RET}
  ) where {RET<:RingElem}

  a = numerator(f)
  b = denominator(f)
  W = localized_ring(OO(X))
  I = quotient(ideal(W, b) + modulus(OO(X)), ideal(W, a))
  U = AffineSchemeOpenSubscheme(X, I)
  g = [OO(V)(f) for V in affine_patches(U)]
  R = AffineSchemeOpenSubschemeRing(X, U)
  return R(g)
end

function maximal_extension(
    X::AbsAffineScheme{<:Ring, <:MPolyQuoRing},
    f::AbstractAlgebra.Generic.FracFieldElem{RET}
  ) where {RET<:RingElem}
  a = numerator(f)
  b = denominator(f)
  W = ambient_coordinate_ring(X)
  I = quotient(ideal(W, b) + modulus(OO(X)), ideal(W, a))
  U = AffineSchemeOpenSubscheme(X, I)
  g = [OO(V)(f) for V in affine_patches(U)]
  R = AffineSchemeOpenSubschemeRing(X, U)
  return R(g)
end

function maximal_extension(
    X::AbsAffineScheme{<:Ring, <:MPolyRing},
    f::AbstractAlgebra.Generic.FracFieldElem{RET}
  ) where {RET<:RingElem}
  a = numerator(f)
  b = denominator(f)
  W = ambient_coordinate_ring(X)
  I = quotient(ideal(W, b), ideal(W, a))
  U = AffineSchemeOpenSubscheme(X, I)
  g = [OO(V)(f) for V in affine_patches(U)]
  R = AffineSchemeOpenSubschemeRing(X, U)
  return R(g)
end

@doc raw"""
    maximal_extension(X::AffineScheme, f::Vector{AbstractAlgebra.Generic.FracFieldElem})

Return the extension of the restriction of the ``fᵢ`` as a
set of rational functions on ``X`` as *regular* functions to a
common maximal domain of definition ``U ⊂ X``.

**Note:** When ``X = Spec(R)`` with ``R = (𝕜[x₁,…,xₙ]/I)[S⁻¹]``,
the numerators and denominators of the entries of ``f`` have to
be elements of the ring ``𝕜[x₁,…,xₙ]``.
"""
function maximal_extension(
    X::AbsAffineScheme{<:Ring, <:AbsLocalizedRing},
    f::Vector{AbstractAlgebra.Generic.FracFieldElem{RET}}
  ) where {RET<:RingElem}
  if length(f) == 0
    return AffineSchemeOpenSubscheme(X), Vector{structure_sheaf_elem_type(X)}()
  end
  R = base_ring(parent(f[1]))
  for a in f
    R == base_ring(parent(a)) || error("fractions do not belong to the same ring")
  end
  R == ambient_coordinate_ring(X) || error("fractions do not belong to the base ring of the scheme")
  a = numerator.(f)
  b = denominator.(f)
  W = OO(X)
  I = ideal(W, one(W))
  for p in f
    I = intersect(quotient(ideal(W, denominator(p)), ideal(W, numerator(p))), I)
  end
  U = AffineSchemeOpenSubscheme(X, I)
  S = AffineSchemeOpenSubschemeRing(X, U)
  # TODO: For some reason, the type of the inner vector is not inferred if it has no entries.
  # Investigate why? Type instability?
  return U, [AffineSchemeOpenSubschemeRingElem(S, (elem_type(OO(X))[OO(V)(a) for V in affine_patches(U)])) for a in f]
end

function maximal_extension(
    X::AbsAffineScheme{<:Ring, <:MPolyRing},
    f::Vector{AbstractAlgebra.Generic.FracFieldElem{RET}}
  ) where {RET<:RingElem}
  if length(f) == 0
    return AffineSchemeOpenSubscheme(X), Vector{structure_sheaf_elem_type(X)}()
  end
  R = base_ring(parent(f[1]))
  for a in f
    R == base_ring(parent(a)) || error("fractions do not belong to the same ring")
  end
  R == ambient_coordinate_ring(X) || error("fractions do not belong to the base ring of the scheme")
  W = ambient_coordinate_ring(X)
  I = ideal(W, one(W))
  for p in f
    I = intersect(quotient(ideal(W, denominator(p)), ideal(W, numerator(p))), I)
  end
  U = AffineSchemeOpenSubscheme(X, I)
  S = AffineSchemeOpenSubschemeRing(X, U)
  # TODO: For some reason, the type of the inner vector is not inferred if it has no entries.
  # Investigate why? Type instability?
  return U, [AffineSchemeOpenSubschemeRingElem(S, [OO(V)(a) for V in affine_patches(U)]) for a in f]
end

function maximal_extension(
    X::AbsAffineScheme{<:Ring, <:MPolyQuoRing},
    f::Vector{AbstractAlgebra.Generic.FracFieldElem{RET}}
  ) where {RET<:RingElem}
  if length(f) == 0
    return AffineSchemeOpenSubscheme(X), Vector{structure_sheaf_elem_type(X)}()
  end
  R = base_ring(parent(f[1]))
  for a in f
    R == base_ring(parent(a)) || error("fractions do not belong to the same ring")
  end
  R == ambient_coordinate_ring(X) || error("fractions do not belong to the base ring of the scheme")
  W = ambient_coordinate_ring(X)
  I = ideal(W, one(W))
  for p in f
    I = intersect(quotient(ideal(W, denominator(p)), ideal(W, numerator(p))), I)
  end
  U = AffineSchemeOpenSubscheme(X, I)
  S = AffineSchemeOpenSubschemeRing(X, U)
  # TODO: For some reason, the type of the inner vector is not inferred if it has no entries.
  # Investigate why? Type instability?
  return U, [AffineSchemeOpenSubschemeRingElem(S, [OO(V)(a) for V in affine_patches(U)]) for a in f]
end
#TODO: implement the catchall versions of the above functions.

########################################################################
# Subscheme constructors                                               #
########################################################################
function subscheme(U::AffineSchemeOpenSubscheme, g::Vector{T}) where {T<:AffineSchemeOpenSubschemeRingElem}
  all(x->(parent(x)==OO(U)), g) || error("elements do not belong to the correct ring")
  X = ambient_scheme(U)
  gen_list = Vector{elem_type(OO(X))}()
  for f in g
    gen_list = vcat(gen_list, OO(X).([lifted_numerator(f[i]) for i in 1:ngens(U)]))
  end
  Z = subscheme(X, gen_list)
  return AffineSchemeOpenSubscheme(Z, complement_equations(U))
end

