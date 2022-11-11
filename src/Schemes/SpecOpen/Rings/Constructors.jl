export maximal_extension

########################################################################
# Constructors for SpecOpenRing                                        #
########################################################################

function OO(U::SpecOpen) 
  if !isdefined(U, :ring_of_functions) 
    U.ring_of_functions = SpecOpenRing(ambient_scheme(U), U)
  end
  return U.ring_of_functions::SpecOpenRing
  #return U.ring_of_functions::SpecOpenRing{affine_patch_type(U), typeof(U)}
end

########################################################################
# Constructors for SpecOpenRingElem                                    #
########################################################################

########################################################################
# Coercion                                                             #
########################################################################
(R::SpecOpenRing)(f::RingElem) = SpecOpenRingElem(R, [OO(U)(f) for U in affine_patches(domain(R))])
(R::SpecOpenRing)(f::MPolyQuoLocalizedRingElem) = SpecOpenRingElem(R, [OO(U)(lifted_numerator(f), lifted_denominator(f)) for U in affine_patches(domain(R))], check=false)

(R::SpecOpenRing)(f::Vector{T}) where {T<:RingElem} = SpecOpenRingElem(R, [OO(domain(R)[i])(f[i]) for i in 1:length(f)])

function (R::SpecOpenRing)(f::SpecOpenRingElem)
  parent(f) === R && return f
  return SpecOpenRingElem(R, [restrict(f, U) for U in affine_patches(domain(R))])
end

########################################################################
# Additional constructors for the Ring interface                       #
########################################################################

one(R::SpecOpenRing) = SpecOpenRingElem(R, [one(OO(U)) for U in affine_patches(domain(R))], check=false)
zero(R::SpecOpenRing) = SpecOpenRingElem(R, [zero(OO(U)) for U in affine_patches(domain(R))], check=false)
(R::SpecOpenRing)() = zero(R)
(R::SpecOpenRing)(a::Integer) = SpecOpenRingElem(R, [OO(U)(a) for U in affine_patches(domain(R))], check=false)
(R::SpecOpenRing)(a::Int64) = SpecOpenRingElem(R, [OO(U)(a) for U in affine_patches(domain(R))], check=false)
(R::SpecOpenRing)(a::fmpz) = SpecOpenRingElem(R, [OO(U)(a) for U in affine_patches(domain(R))], check=false)

########################################################################
# Copying                                                              #
########################################################################
function Base.deepcopy_internal(f::SpecOpenRingElem, dict::IdDict)
  return SpecOpenRingElem(parent(f), copy(restrictions(f)), check=false)
end

########################################################################
# Maximal extensions of rational functions on affine schemes           #
########################################################################
@Markdown.doc """
    maximal_extension(X::Spec, f::AbstractAlgebra.Generic.Frac)

Return the maximal extension of the restriction of ``f``
to a rational function on ``X`` on a maximal domain of 
definition ``U ⊂ X``. 

**Note:** When ``X = Spec(R)`` with ``R = (𝕜[x₁,…,xₙ]/I)[S⁻¹]``, 
the numerator and denominator of ``f`` have to be elements of 
the ring ``𝕜[x₁,…,xₙ]``.
"""
function maximal_extension(
    X::AbsSpec{<:Ring, <:MPolyLocalizedRing}, 
    f::AbstractAlgebra.Generic.Frac{RET}
  ) where {RET<:MPolyElem}

  a = numerator(f)
  b = denominator(f)
  g = gcd(a, b)
  if !isone(g)
    a = divexact(a, g)
    b = divexact(b, g)
    f = parent(f)(a,b)
  end
  W = OO(X)
  U = SpecOpen(X, [b])
  g = [OO(V)(f) for V in affine_patches(U)]
  R = SpecOpenRing(X, U)
  return R(g)
end

function maximal_extension(
    X::AbsSpec{<:Ring, <:MPolyQuoLocalizedRing}, 
    f::AbstractAlgebra.Generic.Frac{RET}
  ) where {RET<:RingElem}

  a = numerator(f)
  b = denominator(f)
  W = localized_ring(OO(X))
  I = quotient(ideal(W, b) + modulus(OO(X)), ideal(W, a))
  U = SpecOpen(X, I)
  g = [OO(V)(f) for V in affine_patches(U)]
  R = SpecOpenRing(X, U)
  return R(g)
end

function maximal_extension(
    X::AbsSpec{<:Ring, <:MPolyQuo}, 
    f::AbstractAlgebra.Generic.Frac{RET}
  ) where {RET<:RingElem}
  a = numerator(f)
  b = denominator(f)
  W = ambient_coordinate_ring(X)
  I = quotient(ideal(W, b) + modulus(OO(X)), ideal(W, a))
  U = SpecOpen(X, I)
  g = [OO(V)(f) for V in affine_patches(U)]
  R = SpecOpenRing(X, U)
  return R(g)
end

function maximal_extension(
    X::AbsSpec{<:Ring, <:MPolyRing}, 
    f::AbstractAlgebra.Generic.Frac{RET}
  ) where {RET<:RingElem}
  a = numerator(f)
  b = denominator(f)
  W = ambient_coordinate_ring(X)
  I = quotient(ideal(W, b), ideal(W, a))
  U = SpecOpen(X, I)
  g = [OO(V)(f) for V in affine_patches(U)]
  R = SpecOpenRing(X, U)
  return R(g)
end

@Markdown.doc """
    maximal_extension(X::Spec, f::Vector{AbstractAlgebra.Generic.Frac})

Return the extension of the restriction of the ``fᵢ`` as a
set of rational functions on ``X`` as *regular* functions to a 
common maximal domain of definition ``U ⊂ X``.

**Note:** When ``X = Spec(R)`` with ``R = (𝕜[x₁,…,xₙ]/I)[S⁻¹]``, 
the numerators and denominators of the entries of ``f`` have to 
be elements of the ring ``𝕜[x₁,…,xₙ]``.
"""
function maximal_extension(
    X::AbsSpec{<:Ring, <:AbsLocalizedRing}, 
    f::Vector{AbstractAlgebra.Generic.Frac{RET}}
  ) where {RET<:RingElem}
  if length(f) == 0
    return SpecOpen(X), Vector{structure_sheaf_elem_type(X)}()
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
  U = SpecOpen(X, I)
  S = SpecOpenRing(X, U)
  # TODO: For some reason, the type of the inner vector is not inferred if it has no entries. 
  # Investigate why? Type instability?
  return U, [SpecOpenRingElem(S, (elem_type(OO(X))[OO(V)(a) for V in affine_patches(U)])) for a in f]
end

function maximal_extension(
    X::AbsSpec{<:Ring, <:MPolyRing}, 
    f::Vector{AbstractAlgebra.Generic.Frac{RET}}
  ) where {RET<:RingElem}
  if length(f) == 0
    return SpecOpen(X), Vector{structure_sheaf_elem_type(X)}()
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
  U = SpecOpen(X, I)
  S = SpecOpenRing(X, U)
  # TODO: For some reason, the type of the inner vector is not inferred if it has no entries. 
  # Investigate why? Type instability?
  return U, [SpecOpenRingElem(S, [OO(V)(a) for V in affine_patches(U)]) for a in f]
end

function maximal_extension(
    X::AbsSpec{<:Ring, <:MPolyQuo}, 
    f::Vector{AbstractAlgebra.Generic.Frac{RET}}
  ) where {RET<:RingElem}
  if length(f) == 0
    return SpecOpen(X), Vector{structure_sheaf_elem_type(X)}()
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
  U = SpecOpen(X, I)
  S = SpecOpenRing(X, U)
  # TODO: For some reason, the type of the inner vector is not inferred if it has no entries. 
  # Investigate why? Type instability?
  return U, [SpecOpenRingElem(S, [OO(V)(a) for V in affine_patches(U)]) for a in f]
end
#TODO: implement the catchall versions of the above functions.

########################################################################
# Subscheme constructors                                               #
########################################################################
function subscheme(U::SpecOpen, g::Vector{T}) where {T<:SpecOpenRingElem}
  all(x->(parent(x)==OO(U)), g) || error("elements do not belong to the correct ring")
  X = ambient_scheme(U)
  Z = subscheme(X, vcat([[lifted_numerator(f[i]) for i in 1:ngens(U)] for f in g]...))
  return SpecOpen(Z, gens(U))
end
