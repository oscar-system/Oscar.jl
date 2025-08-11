################################################################################
#
#  Morphisms between affine algebra homomorphisms
#
################################################################################

const AffAlgHom = MPolyAnyMap{DT, CT, Nothing} where {T <: FieldElem,
                                                       U1 <: MPolyRingElem{T},
                                                       U2 <: MPolyRingElem{T}, # types in domain an codomain might differ: one might be decorated, the other not.
                                                       DT <: Union{MPolyRing{T}, MPolyQuoRing{U1}},
                                                       CT <: Union{MPolyRing{T}, MPolyQuoRing{U2}}}

affine_algebra_morphism_type(::Type{T}) where {T <: Union{MPolyRing, MPolyQuoRing}} = morphism_type(T, T)

affine_algebra_morphism_type(::T) where {T} = affine_algebra_morphism_type(T)

affine_algebra_morphism_type(::Type{S}, ::Type{T}) where {S <: Union{MPolyRing, MPolyQuoRing},
                                                          T <: Union{MPolyRing, MPolyQuoRing}} = morphism_type(S, T)

affine_algebra_morphism_type(R::S, U::T) where {S <: Ring, T} = affine_algebra_morphism_type(S, T)

################################################################################
#
#  Singular data stuff
#
################################################################################

@attr Any _singular_ring_domain(f::MPolyAnyMap) = singular_poly_ring(domain(f))

@attr Any _singular_ring_codomain(f::MPolyAnyMap) = singular_poly_ring(codomain(f))

@attr Any function _singular_algebra_morphism(f::MPolyAnyMap{<:MPolyRing, <:Union{MPolyRing, MPolyQuoRing}, Nothing})
  @assert coefficient_ring(domain(f)) === coefficient_ring(codomain(f))  "singular does not handle coefficient maps"
  DS = _singular_ring_domain(f)
  CS = _singular_ring_codomain(f)
  CSimgs = CS.(_images(f))
  return Singular.AlgebraHomomorphism(DS, CS, CSimgs)
end

################################################################################
#
#  Kernel
#
################################################################################

@doc raw"""
    kernel(F::AffAlgHom)

Return the kernel of `F`.
"""
@attr Any function kernel(f::AffAlgHom) # TODO: need some ideal_type(domain(f)) here :)
  C = codomain(f)
  return preimage(f, ideal(C, [zero(C)]))
end

##############################################################################
#
#  Injectivity
#
##############################################################################

@doc raw"""
    is_injective(F::AffAlgHom)

Return `true` if `F` is injective, `false` otherwise.
"""
function is_injective(F::AffAlgHom)
  return iszero(kernel(F))
end

################################################################################
#
#  Surjectivity
#
################################################################################

@doc raw"""
    is_surjective(F::AffAlgHom)

Return `true` if `F` is surjective, `false` otherwise.
"""
function is_surjective(F::AffAlgHom)
  return all(x -> has_preimage_with_preimage(F, x)[1], gens(codomain(F)))
end

################################################################################
#
#  Bijectivity
#
################################################################################

@doc raw"""
    is_bijective(F::AffAlgHom)

Return `true` if `F` is bijective, `false` otherwise.
"""
function is_bijective(F::AffAlgHom)
  return is_injective(F) && is_surjective(F)
end
 
################################################################################
#
#  Finiteness
#
################################################################################

@doc raw"""
    is_finite(F::AffAlgHom)

Return `true` if `F` is finite, `false` otherwise.
"""
function is_finite(F::AffAlgHom)
  # Use [GP08, Proposition 3.1.5]
  T, _, _, J = _groebner_data(F)
  n = ngens(codomain(F))
  o = lex(gens(T)[1:n])*induce(gens(T)[n + 1:end], default_ordering(domain(F)))
  gb = groebner_basis(J, ordering = o)

  # Check if for all i, powers of x_i occur as leading monomials
  b = falses(n)
  for f in gb
    exp = exponent_vector(leading_monomial(f, ordering = o), 1)
    inds = findall(!is_zero, exp)
    if length(inds) > 1 || inds[1] > n
      continue
    end
    b[inds[1]] = true
  end

  return all(b)
end

##############################################################################
#
#  Inverse of maps of affine algebras and preimages of elements
#
##############################################################################

@doc raw"""
    inverse(F::AffAlgHom)

If `F` is bijective, return its inverse.

# Examples
```jldoctest
julia> D1, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> D, _ = quo(D1, [y-x^2, z-x^3]);

julia> C, (t,) = polynomial_ring(QQ, [:t]);

julia> F = hom(D, C, [t, t^2, t^3]);

julia> is_bijective(F)
true

julia> G = inverse(F)
Ring homomorphism
  from multivariate polynomial ring in 1 variable over QQ
  to quotient of multivariate polynomial ring by ideal (-x^2 + y, -x^3 + z)
defined by
  t -> x

julia> G(t)
x
```
"""
function inverse(F::AffAlgHom)
  !is_injective(F) && error("Homomorphism is not injective")

  R = domain(F)
  S = codomain(F)

  preimgs = elem_type(R)[]
  for i in 1:ngens(S)
    fl, p = has_preimage_with_preimage(F, gen(S, i))
    fl || error("Homomorphism is not surjective")
    push!(preimgs, p)
  end

  return hom(S, R, preimgs)
end

function preimage_with_kernel(F::AffAlgHom, f::Union{MPolyRingElem, MPolyQuoRingElem})
  return preimage(F, f), kernel(F)
end

function preimage(F::AffAlgHom, f::Union{MPolyRingElem, MPolyQuoRingElem})
  fl, g = has_preimage_with_preimage(F, f)
  !fl && error("Element not contained in image")
  return g
end

function has_preimage_with_preimage(F::AffAlgHom, f::Union{MPolyRingElem, MPolyQuoRingElem})
  # Basically [GP09, p. 86, Solution 2]
  @req parent(f) === codomain(F) "Polynomial is not element of the codomain"

  R = domain(F)
  S = codomain(F)
  m = ngens(R)
  n = ngens(S)

  T, inc, pr, J = _groebner_data(F)
  o = induce(gens(T)[1:n], default_ordering(S))*induce(gens(T)[n + 1:end], default_ordering(R))
  nf = normal_form(inc(lift(f)), J, ordering = o)
  if isone(cmp(o, gen(T, n), leading_monomial(nf, ordering = o)))
    return true, pr(nf)
  end
  return false, zero(R)
end

@doc raw"""
    preimage(F::MPolyAnyMap, I::Ideal)

Return the preimage of the ideal `I` under `F`.
"""
function preimage(F::MPolyAnyMap, I::Ideal)
  # This generic routine does not work for maps where the domain is a quotient ring. 
  # error message: _singular_algebra_morphism(...) does not have a method for this.
  # Hence it has been split into two specialized methods below.
  error("not implemented")
end

function preimage(
    f::MPolyAnyMap{<:MPolyRing{T}, CT, Nothing}, 
    I::Union{MPolyIdeal, MPolyQuoIdeal}
  ) where {T <: RingElem,
           CT <: Union{MPolyRing{T}, MPolyQuoRing{<:MPolyRingElem{T}}}}
  return _preimage_via_singular(f, I)
end

function preimage(
    f::MPolyAnyMap{S, S, Nothing}, 
    I::MPolyIdeal
  ) where {S <: MPolyRing{<:RingElem}}
  # If the map is the inclusion of a subring in a subset of 
  # variables, use the `eliminate` command instead.
  _assert_has_maps_variables_to_variables!(f)
  success, ind = _maps_variables_to_variables(f)
  if success
    img_gens = elem_type(domain(f))[] # for the projection map
    elim_ind = Int[]
    for (i, x) in enumerate(gens(codomain(f)))
      k = findfirst(==(i), ind)
      if isnothing(k)
        push!(img_gens, zero(domain(f)))
        push!(elim_ind, i)
      else
        push!(img_gens, gen(domain(f), k::Int))
      end
    end
    J = eliminate(I, elim_ind)
    pr = hom(codomain(f), domain(f), img_gens)
    return pr(J)
  end
  return _preimage_via_singular(f, I)
end

function _preimage_via_singular(
    f::MPolyAnyMap{<:MPolyRing{T}, CT, Nothing}, 
    I::Union{MPolyIdeal, MPolyQuoIdeal}
  ) where {T <: RingElem,
           CT <: Union{MPolyRing{T}, MPolyQuoRing{<:MPolyRingElem{T}}}}
  @req base_ring(I) === codomain(f) "Parent mismatch"
  D = domain(f)
  salghom = _singular_algebra_morphism(f)
  CS = codomain(salghom)
  V = gens(I)
  Ix = Singular.Ideal(CS, CS.(V))
  prIx = Singular.preimage(salghom, Ix)
  return ideal(D, D.(gens(prIx)))
end

function preimage(
    f::MPolyAnyMap{<:MPolyQuoRing, CT}, 
    I::Union{MPolyIdeal, MPolyQuoIdeal}
  ) where {T <: RingElem,
           CT <: Union{MPolyRing{T}, MPolyQuoRing{<:MPolyRingElem{T}}}}
  @req base_ring(I) === codomain(f) "Parent mismatch"
  R = base_ring(domain(f))
  help_map = hom(R, domain(f), gens(domain(f)); check=false)
  g = compose(help_map, f)
  K = preimage(g, I)
  return ideal(domain(f), help_map.(gens(K)))
end

@attr Any function kernel(f::MPolyAnyMap{<:Union{MPolyRing, MPolyQuoRing}, <:Generic.LaurentMPolyWrapRing, <:Any, <:Any})
  C = codomain(f)
  return preimage(f, ideal(C, [zero(C)]))
end

# Let F: K[x]/I_1 -> K[y]/I_2, x_i \mapsto f_i .
# Construct the polynomial ring K[y, x], the natural maps K[x] -> K[y, x]
# and K[y, x] -> K[y], and the ideal I_2 + (y_i - f_i) in it.
# No actual Gröbner basis computation is done here, but computed bases are
# cached in the ideal.
function _groebner_data(F::AffAlgHom)
  R = domain(F)
  S = codomain(F)
  K = coefficient_ring(R)
  @req K === coefficient_ring(S) "Coefficient rings of domain and codomain must coincide"
  S2 = base_ring(modulus(S))
  m = ngens(R)
  n = ngens(S)
  J = get_attribute!(F, :groebner_data) do
    T, _ = polynomial_ring(K, m + n; cached = false)

    S2toT = hom(S2, T, [ gen(T, i) for i in 1:n ])

    fs = map(lift, _images(F))
    return S2toT(modulus(S)) + ideal(T, [ gen(T, n + i) - S2toT(fs[i]) for i in 1:m ])
  end::MPolyIdeal{mpoly_type(K)}
  T = base_ring(J)
  S2toT = hom(S2, T, [ gen(T, i) for i in 1:n ])
  TtoR = hom(T, R, append!(zeros(R, n), gens(R)))
  return T, S2toT, TtoR, J
end
