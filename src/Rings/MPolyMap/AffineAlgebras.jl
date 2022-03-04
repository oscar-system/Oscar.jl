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
# Legacy construct
#
################################################################################

#@doc Markdown.doc"""
#    AlgebraHomomorphism(D::U, C::W, V::Vector{X}) where 
#    {T, S <: MPolyElem{T},
#    U <: Union{MPolyRing{T}, MPolyQuo{S}},
#    W <: Union{MPolyRing{T}, MPolyQuo{S}},
#    X <: Union{S, MPolyQuoElem{S}}}
#   
#Create the algebra homomorphism $D \rightarrow C$ defined by sending the $i$th generator of $D$ to the $i$th element of $V$. 
#Allow types `MPolyRing` and `MPolyQuo` for $C$ and $D$ as well as entries of type `MPolyElem` and `MPolyQuoElem` for `X`.
#Alternatively, use `hom(D::U, C::W, V::Vector{X})`.
#
## Examples
#```jldoctest
#julia> D, (t,) = PolynomialRing(QQ, ["t"])
#(Multivariate Polynomial Ring in t over Rational Field, fmpq_mpoly[t])
#
#julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"])
#(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])
#
#julia> C, p = quo(R,  ideal(R, [x*y-1]))
#(Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x*y - 1), Map from
#Multivariate Polynomial Ring in x, y over Rational Field to Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x*y - 1) defined by a julia-function with inverse)
#
#julia> V = [p(y)]
#1-element Vector{MPolyQuoElem{fmpq_mpoly}}:
# y
#
#julia> P = hom(D, C, V)
#Algebra homomorphism with
#
#domain: Multivariate Polynomial Ring in t over Rational Field
#
#codomain: Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x*y - 1)
#
#defining images of generators: MPolyQuoElem{fmpq_mpoly}[y]
#```
#"""
#function AlgebraHomomorphism(D::U, C::W, V::Vector{X}) where 
#    {T, S <: MPolyElem{T},
#    U <: Union{MPolyRing{T}, MPolyQuo{S}},
#    W <: Union{MPolyRing{T}, MPolyQuo{S}},
#    X <: Union{S, MPolyQuoElem{S}}}
#   return hom(D, C, copy(V))
#end

################################################################################
#
#  Singular data stuff
#
################################################################################

@attr Any _singular_ring_domain(f::MPolyAnyMap) = singular_ring(domain(f))

@attr Any _singular_ring_codomain(f::MPolyAnyMap) = singular_ring(codomain(f))

@attr Any function _singular_algebra_morphism(f::MPolyAnyMap)
  DS = _singular_ring_domain(f)
  CS = _singular_ring_codomain(f)
  CSimgs = CS.(_images(f))
  return Singular.AlgebraHomomorphism(DS, CS, CSimgs)
end

################################################################################
#
#  Preimage of ideal
#
################################################################################

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

################################################################################
#
#  Kernel
#
################################################################################

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
    isinjective(F::AlgHom)

Return `true` if `F` is injective, `false` otherwise.
"""
function isinjective(F::AffAlgHom)
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
    # Build auxilliary objects
    (T, inc, J) = _containement_helper(S, n, m, I, W, ord)
    D = normal_form([gen(T, i) for i in 1:m], J)
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

function issurjective(F::AffAlgHom)
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

function isbijective(F::AffAlgHom)
  return isinjective(F) && issurjective(F)
end

################################################################################
#
#  Finiteness
#
################################################################################

@doc Markdown.doc"""
    isfinite(F::AlgHom)

Return `true` if `F` is finite, `false` otherwise.
"""
function isfinite(F::AffAlgHom)
  (T, _, _, J, _) = _groebner_data(F, :lex)
  G = collect(groebner_assure(J))
  # Find all elements with leading monomial which contains the 
  # variables x_i.
  s = codomain(F)
  m = ngens(s)
  L = leading_monomial.(G)

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

function inverse(F::AffAlgHom)
  !isinjective(F) && error("Homomorphism is not injective")
  !issurjective(F) && error("Homomorphism is not surjective")

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
  D = normal_form([inc(g)], J)
  !(leading_monomial(D[1]) < gen(T, m)) && error("Element not contained in image")
  return (pr(D[1]), kernel(F))
end
