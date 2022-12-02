################################################################################
#
#  Morphisms between affine algebra homomorphisms
#
################################################################################

const AffAlgHom = MPolyAnyMap{DT, CT, Nothing} where {T <: FieldElem,
                                                       U <: MPolyElem{T},
                                                       DT <: Union{MPolyRing{T}, MPolyQuo{U}},
                                                       CT <: Union{MPolyRing{T}, MPolyQuo{U}}}

affine_algebra_morphism_type(::Type{T}) where {T <: Union{MPolyRing, MPolyQuo}} = morphism_type(T, T)

affine_algebra_morphism_type(::T) where {T} = affine_algebra_morphism_type(T)

affine_algebra_morphism_type(::Type{S}, ::Type{T}) where {S <: Union{MPolyRing, MPolyQuo},
                                                          T <: Union{MPolyRing, MPolyQuo}} = morphism_type(S, T)

affine_algebra_morphism_type(R::S, U::T) where {S <: Ring, T} = affine_algebra_morphism_type(S, T)

################################################################################
#
#  Singular data stuff
#
################################################################################

@attr Any _singular_ring_domain(f::MPolyAnyMap) = singular_poly_ring(domain(f))

@attr Any _singular_ring_codomain(f::MPolyAnyMap) = singular_poly_ring(codomain(f))

@attr Any function _singular_algebra_morphism(f::MPolyAnyMap)
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

@doc Markdown.doc"""
    kernel(F::AffAlgHom)

Return the kernel of `F`.
"""
function kernel(f::AffAlgHom)
  get_attribute!(f, :kernel) do
    C = codomain(f)
    return preimage(f, ideal(C, [zero(C)]))
  end # TODO: need some ideal_type(domain(f)) here :)
end

##############################################################################
#
#  Injectivity
#
##############################################################################

@doc Markdown.doc"""
    is_injective(F::AffAlgHom)

Return `true` if `F` is injective, `false` otherwise.
"""
function is_injective(F::AffAlgHom)
  iszero(kernel(F))
end

# Helper function related to the computation of surjectivity, preimage etc.
# Stores the necessary data in groebner_data, resp. groebner_data_lex
function _groebner_data(F::AffAlgHom, ord::Symbol)
  r = domain(F)
  s = codomain(F)
  n = ngens(r)
  m = ngens(s)
  return get_attribute!(F, ord) do
    (S, I, W, _) = _ring_helper(s, zero(s), _images(F))
    # Build auxiliary objects
    (T, inc, J) = _containment_helper(S, n, m, I, W, ord)
    # make sure to compute the NF with respect to the correct ordering  
    D = normal_form([gen(T, i) for i in 1:m], J,
                    ordering = first(collect(keys(J.gb))))
    A = [zero(r) for i in 1:m]
    B = [gen(r, i) for i in 1:n]
    pr = hom(T, r, vcat(A, B))
    groebner_data_lex = (T, inc, pr, J, D)
    return groebner_data_lex
  end
end

################################################################################
#
#  Surjectivity
#
################################################################################

@doc Markdown.doc"""
    is_surjective(F::AffAlgHom)

Return `true` if `F` is is_surjective, `false` otherwise.
"""
function is_surjective(F::AffAlgHom)
  # Compute data necessary for computation
  r = domain(F)
  s = codomain(F)
  n = ngens(r)
  m = ngens(s)
  (T, _, _, _, D) = _groebner_data(F, :degrevlex)

  # Check if map is surjective

  for i in 1:m
     if !(leading_monomial(D[i]) < gen(T, m))
        return false
     end
  end
  return true
end

################################################################################
#
#  Bijectivity
#
################################################################################

@doc Markdown.doc"""
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

@doc Markdown.doc"""
    isfinite(F::AffAlgHom)

Return `true` if `F` is finite, `false` otherwise.
"""
function isfinite(F::AffAlgHom)
  (T, _, _, J, _) = _groebner_data(F, :lex)
  G = collect(groebner_assure(J))
  # Find all elements with leading monomial which contains the 
  # variables x_i.
  s = codomain(F)
  m = ngens(s)
  L = map(AbstractAlgebra.leading_monomial, G)

  # Check if for all i, powers of x_i occur as leading monomials
  N = Vector{Int}()
  for i in 1:length(L)
     exp = exponent_vector(L[i], 1)
     f = findall(x->x!=0, exp) 
     length(f) == 1 && f[1] <= m && union!(N, f)
  end

  return length(N) == m
end

##############################################################################
#
#  Inverse of maps of affine algebras and preimages of elements
#
##############################################################################

@doc Markdown.doc"""
    inverse(F::AffAlgHom)

If `F` is bijective, return its inverse.

# Examples
```jldoctest
julia> D1, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"]);

julia> D, _ = quo(D1, [y-x^2, z-x^3])
(Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(-x^2 + y, -x^3 + z), Map from
Multivariate Polynomial Ring in x, y, z over Rational Field to Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(-x^2 + y, -x^3 + z) defined by a julia-function with inverse)

julia> C, (t,) = PolynomialRing(QQ, ["t"]);

julia> F = hom(D, C, [t, t^2, t^3]);

julia> is_bijective(F)
true

julia> G = inverse(F)
Map with following data
Domain:
=======
Multivariate Polynomial Ring in t over Rational Field
Codomain:
=========
Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(-x^2 + y, -x^3 + z)

julia> G(t)
x
```
"""
function inverse(F::AffAlgHom)
  !is_injective(F) && error("Homomorphism is not injective")
  !is_surjective(F) && error("Homomorphism is not surjective")

  # Compute inverse map via preimages of algebra generators
  r = domain(F)
  s = codomain(F)
  n = ngens(r)
  m = ngens(s)

  (T, _, pr, _, D) = _groebner_data(F, :degrevlex)
  psi = hom(s, r, [pr(D[i]) for i in 1:m])
  #psi.kernel = ideal(s, [zero(s)])
  return psi
end

function preimage(F::AffAlgHom, f::Union{MPolyElem, MPolyQuoElem})
  @assert parent(f) === codomain(F)
  return preimage_with_kernel(F, f)[1]
end

function preimage_with_kernel(F::AffAlgHom, f::Union{MPolyElem, MPolyQuoElem})
  @assert parent(f) === codomain(F)
  r = domain(F)
  s = codomain(F)
  n = ngens(r)
  m = ngens(s)

  (S, _, _, g) = _ring_helper(s, f, [zero(s)])
  (T, inc, pr, J, o) = _groebner_data(F, :degrevlex)
  D = normal_form([inc(g)], J, ordering = first(collect(keys(J.gb))))
  !(leading_monomial(D[1]) < gen(T, m)) && error("Element not contained in image")
  return (pr(D[1]), kernel(F))
end

@doc Markdown.doc"""
    preimage(F::AffAlgHom, I::U) where U <: Union{MPolyIdeal, MPolyQuoIdeal}

Return the preimage of the ideal `I` under `F`.
"""
function preimage(f::AffAlgHom, I::Union{MPolyIdeal, MPolyQuoIdeal})
  @req base_ring(I) === codomain(f) "Parent mismatch"
  D = domain(f)
  salghom = _singular_algebra_morphism(f)
  CS = codomain(salghom)
  V = gens(I)
  Ix = Singular.Ideal(CS, CS.(V))
  prIx = Singular.preimage(salghom, Ix)
  return ideal(D, D.(gens(prIx)))
end
