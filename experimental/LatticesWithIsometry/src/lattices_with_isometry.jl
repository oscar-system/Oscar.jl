
###############################################################################
#
#  Accessors
#
###############################################################################

@doc raw"""
    lattice(Lf::ZZLatWithIsom) -> ZZLat

Given a lattice with isometry $(L, f)$, return the underlying lattice `L`.
"""
lattice(Lf::ZZLatWithIsom) = Lf.Lb

@doc raw"""
    isometry(Lf::ZZLatWithIsom) -> QQMatrix

Given a lattice with isometry $(L, f)$, return the underlying isometry `f`.
"""
isometry(Lf::ZZLatWithIsom) = Lf.f

@doc raw"""
   ambient_isometry(Lf::ZZLatWithIsom) -> QQMatrix

Given a lattice with isometry $(L, f)$, return an isometry of underlying isometry
of the ambient space of `L` inducing `f` on `L`
"""
ambient_isometry(Lf::ZZLatWithIsom) = Lf.f_ambient

@doc raw"""
    order_of_isometry(Lf::ZZLatWithIsom) -> Integer

Given a lattice with isometry $(L, f)$, return the order of the underlying
isometry `f`.
"""
order_of_isometry(Lf::ZZLatWithIsom) = Lf.n

###############################################################################
#
#  Attributes
#
###############################################################################

@doc raw"""
    rank(Lf::ZZLatWithIsom) -> Integer

Given a lattice with isometry $(L, f)$, return the rank of the underlying lattice
`L`.
"""
rank(Lf::ZZLatWithIsom) = rank(lattice(Lf))::Integer

@doc raw"""
    charpoly(Lf::ZZLatWithIsom) -> QQPolyRingElem

Given a lattice with isometry $(L, f)$, return the characteristic polynomial of the
underlying isometry `f`.
"""
charpoly(Lf::ZZLatWithIsom) = charpoly(isometry(Lf))::QQPolyRingElem

@doc raw"""
    minpoly(Lf::ZZLatWithIsom) -> QQPolyRingElem

Given a lattice with isometry $(L, f)$, return the minimal polynomial of the
underlying isometry `f`.
"""
minpoly(Lf::ZZLatWithIsom) = minpoly(isometry(Lf))::QQPolyRingElem

@doc raw"""
    genus(Lf::ZZLatWithIsom) -> ZZGenus

Given a lattice with isometry $(L, f)$, return the genus of the underlying 
lattice `L` (see [`genus(::ZZLat)`](@ref)).
"""
genus(Lf::ZZLatWithIsom) = genus(lattice(Lf))::ZZGenus

@doc raw"""
    ambient_space(Lf::ZZLatWithIsom) -> QuadSpace

Given a lattice with isometry $(L, f)$, return the ambient space of the underlying
lattice `L` (see [`ambient_space(::ZZLat)`](@ref)).
"""
ambient_space(Lf::ZZLatWithIsom) = ambient_space(lattice(Lf))::Hecke.QuadSpace{FlintRationalField, QQMatrix}

@doc raw"""
    basis_matrix(Lf::ZZLatWithIsom) -> QQMatrix

Given a lattice with isometry $(L, f)$, return the basis matrix of the underlying
lattice `L` (see [`basis_matrix(::ZZLat)`](@ref)).
"""
basis_matrix(Lf::ZZLatWithIsom) = basis_matrix(lattice(Lf))::QQMatrix

@doc raw"""
    gram_matrix(Lf::ZZLatWithIsom) -> QQMatrix

Given a lattice with isometry $(L, f)$ with basis matric `B` (see [`basis_matrix(Lf::ZZLatWithIsom)`](@ref))
inside the space $(V, \Phi)$ (see [`ambient_space(Lf::ZZLatWithIsom)`](@ref)), return the gram matrix
of lattice `L` associted to `B` with respect to $\Phi$.
"""
gram_matrix(Lf::ZZLatWithIsom) = gram_matrix(lattice(Lf))::QQMatrix

@doc raw"""
    rational_span(Lf::ZZLatWithIsom) -> QuadSpace

Given a lattice with isometry $(L, f)$, return the rational span $L \otimes \mathbb{Q}$
of the underlying lattice `L`.
"""
rational_span(Lf::ZZLatWithIsom) = rational_span(lattice(Lf))::Hecke.QuadSpace{FlintRationalField, QQMatrix}

@doc raw"""
    det(Lf::ZZLatWithIsom) -> QQFieldElem

Given a lattice with isometry $(L, f)$, return the determinant of the
underlying lattice `L` (see [`det(::ZZLat)`](@ref)).
"""
det(Lf::ZZLatWithIsom) = det(lattice(Lf))::QQFieldElem

@doc raw"""
    scale(Lf::ZZLatWithIsom) -> QQFieldElem

Given a lattice with isometry $(L, f)$, return the scale of the underlying
lattice `L` (see [`scale(::ZZLat)`](@ref)).
"""
scale(Lf::ZZLatWithIsom) = scale(lattice(Lf))::QQFieldElem

@doc raw"""
    norm(Lf::ZZLatWithIsom) -> QQFieldElem

Given a lattice with isometry $(L, f)$, return the norm of the underlying
lattice `L` (see [`norm(::ZZLat)`](@ref)).
"""
norm(Lf::ZZLatWithIsom) = norm(lattice(Lf))::QQFieldElem

@doc raw"""
    is_positive_definite(Lf::ZZLatWithIsom) -> Bool

Given a lattice with isometry $(L, f)$, return whether the underlying
lattice `L` is positive definite (see [`is_positive_definite(::ZZLat)`](@ref)).
"""
is_positive_definite(Lf::ZZLatWithIsom) = is_positive_definite(lattice(Lf))::Bool

@doc raw"""
    is_negative_definite(Lf::ZZLatWithIsom) -> Bool

Given a lattice with isometry $(L, f)$, return whether the underlying
lattice `L` is negative definite (see [`is_positive_definite(::ZZLat)`](@ref)).
"""
is_negative_definite(Lf::ZZLatWithIsom) = is_negative_definite(lattice(Lf))::Bool

@doc raw"""
    is_definite(Lf::ZZLatWithIsom) -> Bool

Given a lattice with isometry $(L, f)$, return whether the underlying
lattice `L` is definite (see [`is_definite(::ZZLat)`](@ref)).
"""
is_definite(Lf::ZZLatWithIsom) = is_definite(lattice(Lf))::Bool

@doc raw"""
    minimum(Lf::ZZLatWithIsom) -> QQFieldElem

Given a positive definite lattice with isometry $(L, f)$, return the minimum
of the underlying lattice `L` (see [`minimum(::ZZLat)`](@ref)).
"""
function minimum(Lf::ZZLatWithIsom)
  @req is_positive_definite(Lf) "Underlying lattice must be positive definite"
  return minimum(lattice(Lf))
end

@doc raw"""
    is_integral(Lf::ZZLatWithIsom) -> Bool

Given a lattice with isometry $(L, f)$, return whether the underlying lattice
is integral, i.e. whether its scale is an integer (see [`scale(::ZZLatWithIsom)`](@ref)).
"""
is_integral(Lf::ZZLatWithIsom) = is_integral(lattice(Lf))::Bool

@doc raw"""
    degree(Lf::ZZLatWithIsom) -> Int

Given a lattice with isometry $(L, f)$ inside the quadratic space $(V, \Phi)$,
return the dimension of `V` as a $\mathbb Q$ vector space.
"""
degree(Lf::ZZLatWithIsom) = degree(lattice(Lf))::Int

@doc raw"""
    is_even(Lf::ZZLatWithIsom) -> Bool

Given a lattice with isometry $(L, f)$, return whether the underlying lattice
`L` is even, i.e. whether its norm is an even integer ([`norn(::ZZLatWithIsom)`](@ref)).

Note that to be even, `L` must be integral (see [`is_integral(::ZZLat)`](@ref)).
"""
is_even(Lf::ZZLatWithIsom) = iseven(lattice(Lf))::Bool

@doc raw"""
    discriminant(Lf::ZZLatWithIsom) -> QQFieldElem

Given a lattice with isometry $(L, f)$, return the discriminant of the underlying
lattice `L` (see [`discriminant(::ZZLat)`](@ref)).
"""
discriminant(Lf::ZZLatWithIsom) = discriminant(lattice(Lf))::QQFieldElem

@doc raw"""
    signature_tuple(Lf::ZZLatWithIsom) -> Tuple{Int, Int, Int}

Given a lattice with isometry $(L, f)$, return the signature tuple of the
underlying lattice `L` (see [`signature_tuple(::ZZLat)`](@ref)).
"""
signature_tuple(Lf::ZZLatWithIsom) = signature_tuple(lattice(Lf))::Tuple{Int, Int, Int}

###############################################################################
#
#  Constructor
#
###############################################################################

@doc raw"""
    lattice_with_isometry(L::ZZLat, f::QQMatrix; check::Bool = true,
                                                ambient_representation = true)
				                                                             -> ZZLatWithIsom

Given a $\mathbb Z$-lattice `L` and a matrix `f`, if `f` defines an isometry
of `L` of order `n`, return the corresponding lattice with isometry pair $(L, f)$.

If `ambient_representation` is set to true, `f` is consider as an isometry of the
ambient space of `L` and the induced isometry on `L` is automatically computed.
Otherwise, an isometry of the ambient space of `L` is constructed, setting the identity
on the complement of the rational span of `L` if it is not of full rank.
"""
function lattice_with_isometry(L::ZZLat, f::QQMatrix; check::Bool = true,
                                                     ambient_representation::Bool = true)
  if rank(L) == 0
    return ZZLatWithIsom(L, matrix(QQ,0,0,[]), identity_matrix(QQ, degree(L)), -1)
  end

  if check
    @req det(f) != 0 "f is not invertible"
  end

  if ambient_representation
    f_ambient = f
    B = basis_matrix(L)
    ok, f = can_solve_with_solution(B, B*f_ambient, side = :left)
    @req ok "Isometry does not restrict to L"
  else
    V = ambient_space(L)
    B = basis_matrix(L)
    B2 = orthogonal_complement(V, B)
    C = vcat(B, B2)
    f_ambient = block_diagonal_matrix([f, identity_matrix(QQ, nrows(B2))])
    f_ambient = inv(C)*f_ambient*C
  end

  n = multiplicative_order(f)

  if check
    @req f*gram_matrix(L)*transpose(f) == gram_matrix(L) "f does not define an isometry of L"
    @req f_ambient*gram_matrix(ambient_space(L))*transpose(f_ambient) == gram_matrix(ambient_space(L)) "f_ambient is not an isometry of the ambient space of L"
    @hassert :ZZLatWithIsom 1 basis_matrix(L)*f_ambient == f*basis_matrix(L)
  end

  return ZZLatWithIsom(L, f, f_ambient, n)::ZZLatWithIsom
end


@doc raw"""
    lattice_with_isometry(L::ZZLat; neg::Bool = false) -> ZZLatWithIsom

Given a $\mathbb Z$-lattice `L` return the lattice with isometry pair $(L, f)$,
where `f` corresponds to the identity mapping of `L`.

If `neg` is set to `true`, then the isometry `f` is negative the identity of `L`.
"""
function lattice_with_isometry(L::ZZLat; neg::Bool = false)
  d = degree(L)
  f = identity_matrix(QQ, d)
  if neg
    f = -f
  end
  return lattice_with_isometry(L, f, check = false, ambient_representation = true)::ZZLatWithIsom
end

###############################################################################
#
#  Operations on lattices with isometry
#
###############################################################################

@doc raw"""
    rescale(Lf::ZZLatWithIsom, a::RationalUnion) -> ZZLatWithIsom

Given a lattice with isometry $(L, f)$ and a rational number `a`, return the lattice
with isometry $(L(a), f)$ (see [`rescale(::ZZLat, ::RationalUnion)`](@ref)).
"""
function rescale(Lf::ZZLatWithIsom, a::Hecke.RationalUnion)
  return lattice_with_isometry(rescale(lattice(Lf), a), ambient_isometry(Lf), check=false)
end

@doc raw"""
    dual(Lf::ZZLatWithIsom) -> ZZLatWithIsom

Given a lattice with isometry $(L, f)$ inside the space $(V, \Phi)$, such that `f` is
induced by an isometry `g` of $(V, \Phi)$, return the lattice with isometry $(L^{\vee}, h)$
where $L^{\vee}$ is the dual of `L` in $(V, \Phi)$ (see [`dual(::ZZLat)`](@ref)) and `h` is
induced by `g`.
"""
function dual(Lf::ZZLatWithIsom)
  @req is_integral(Lf) "Underlying lattice must be integral"
  return lattice_with_isometry(dual(lattice(Lf)), ambient_isometry(Lf), check = false)
end

@doc raw"""
    lll(Lf::ZZLatWithIsom) -> ZZLatWithIsom

Given a lattice with isometry $(L, f)$, return the same lattice with isometry with a different
basis matrix for `L` (see [`basis_matrix(::ZZLat)`](@ref)) obtained by performing an LLL-reduction
on the associated gram matrix of `L` (see [`gram_matrix(::ZZLat)`](@ref)).

Note that matrix representing the action of `f` on `L` changes but the global action
on the ambient space of `L` stays the same.
"""
function lll(Lf::ZZLatWithIsom)
  f = ambient_isometry(Lf)
  L2 = lll(lattice(Lf), same_ambient=true)
  return lattice_with_isometry(L2, f, ambient_representation = true)
end

@doc raw"""
    direct_sum(x::Vector{ZZLatWithIsom}) -> ZZLatWithIsom, Vector{AbstractSpaceMor}
    direct_sum(x::Vararg{ZZLatWithIsom}) -> ZZLatWithIsom, Vector{AbstractSpaceMor}

Given a collection of lattices with isometries $(L_1, f_1) \ldots, (L_n, f_n)$,
return the lattice with isometry $(L, f)$ together with the injections $L_i \to L$,
where `L` is the direct sum $L := L_1 \oplus \ldots \oplus L_n$ and `f` is the
isometry of `L` induced by the diagonal actions of the $f_i$'s.

For objects of type `ZZLatWithIsom`, finite direct sums and finite direct products
agree and they are therefore called biproducts.
If one wants to obtain $(L, f)$ as a direct product with the projections $L \to L_i$,
one should call `direct_product(x)`.
If one wants to obtain $(L, f)$ as a biproduct with the injections $L_i \to L$ and
the projections $L \to L_i$, one should call `biproduct(x)`.
"""
function direct_sum(x::Vector{ZZLatWithIsom})
  @req length(x) >= 2 "Input must consist of at least 2 lattices with isometries"
  W, inj = direct_sum(lattice.(x))
  f = block_diagonal_matrix(ambient_isometry.(x))
  return lattice_with_isometry(W, f, check=false), inj
end

direct_sum(x::Vararg{ZZLatWithIsom}) = direct_sum(collect(x))

@doc raw"""
    direct_product(x::Vector{ZZLatWithIsom}) -> ZZLatWithIsom, Vector{AbstractSpaceMor}
    direct_product(x::Vararg{ZZLatWithIsom}) -> ZZLatWithIsom, Vector{AbstractSpaceMor}

Given a collection of lattices with isometries $(L_1, f_1), \ldots, (L_n, f_n)$,
return the lattice with isometry $(L, f)$ together with the projections $L \to L_i$,
where `L` is the direct product $L := L_1 \times \ldots \times L_n$ and `f` is the
isometry of `L` induced by the diagonal actions of the $f_i$'s.

For objects of type `ZZLatWithIsom`, finite direct sums and finite direct products
agree and they are therefore called biproducts.
If one wants to obtain $(L, f)$ as a direct sum with the injections $L_i \to L$,
one should call `direct_sum(x)`.
If one wants to obtain $(L, f)$ as a biproduct with the injections $L_i \to L$ and
the projections $L \to L_i$, one should call `biproduct(x)`.
"""
function direct_product(x::Vector{ZZLatWithIsom})
  @req length(x) >= 2 "Input must consist of at least 2 lattices with isometries"
  W, proj = direct_product(lattice.(x))
  f = block_diagonal_matrix(ambient_isometry.(x))
  return lattice_with_isometry(W, f, check=false), proj
end

direct_product(x::Vararg{ZZLatWithIsom}) = direct_product(collect(x))

@doc raw"""
    biproduct(x::Vector{ZZLatWithIsom}) -> ZZLatWithIsom, Vector{AbstractSpaceMor}, Vector{AbstractSpaceMor}
    biproduct(x::Vararg{ZZLatWithIsom}) -> ZZLatWithIsom, Vector{AbstractSpaceMor}, Vector{AbstractSpaceMor}

Given a collection of lattices with isometries $(L_1, f_1), \ldots, (L_n, f_n)$,
return the lattice with isometry $(L, f)$ together with the injections
$L_i \to L$ and the projections $L \to L_i$, where `L` is the biproduct
$L := L_1 \oplus \ldots \oplus L_n$ and `f` is the isometry of `L` induced by the
diagonal actions of the $f_i$'s.

For objects of type `ZZLatWithIsom`, finite direct sums and finite direct products
agree and they are therefore called biproducts.
If one wants to obtain $(L, f)$ as a direct sum with the injections $L_i \to L$,
one should call `direct_sum(x)`.
If one wants to obtain $(L, f)$ as a direct product with the projections $L \to L_i$,
one should call `direct_product(x)`.
"""
function biproduct(x::Vector{ZZLatWithIsom})
  @req length(x) >= 2 "Input must consist of at least 2 lattices with isometries"
  W, inj, proj = biproduct(lattice.(x))
  f = block_diagonal_matrix(ambient_isometry.(x))
  return lattice_with_isometry(W, f, check=false), inj, proj
end

biproduct(x::Vararg{ZZLatWithIsom}) = biproduct(collect(x))

###############################################################################
#
#  Hermitian structure
#
###############################################################################

@doc raw"""
    is_of_hermitian_type(Lf::ZZLatWithIsom) -> Bool

Given a lattice with isometry $(L, f)$, return whether the minimal polynomial of
the underlying isometry `f` is (irreducible) cyclotomic.

Note that if $(L, f)$ is of hermitian type with `f` of order `n`, then `L` can
be seen as a hermitian lattice over the order $\mathbb{Z}[\zeta_n]$ where
$\zeta_n$ is a primitive $n$-th root of unity.
"""
function is_of_hermitian_type(Lf::ZZLatWithIsom)
  @req rank(Lf) > 0 "Underlying lattice must have positive rank"
  n = order_of_isometry(Lf)
  if n == -1 || !is_finite(n)
    return false
  end
  return is_cyclotomic_polynomial(minpoly(isometry(Lf)))
end

@doc raw"""
    hermitian_structure(Lf::ZZLatWithIsom) -> HermLat

Given a lattice with isometry $(L, f)$ such that the minimal polynomial of the
underlying isometry `f` is cyclotomic, return the hermitian structure of the
underlying lattice `L` over the $n$th cyclotomic field, where $n$ is the
order of `f`.

If it exists, the hermitian structure is stored.
"""
@attr HermLat function hermitian_structure(Lf::ZZLatWithIsom)
  @req is_of_hermitian_type(Lf) "Lf is not of hermitian type"
  f = isometry(Lf)

  H, res = Hecke.hermitian_structure_with_transfer_data(lattice(Lf), f, ambient_representation = false)

  set_attribute!(Lf, :transfer_data, res)
  return H
end

###############################################################################
#
#  Image centralizer in discriminant group
#
###############################################################################

@doc raw"""
    discriminant_group(Lf::ZZLatWithIsom) -> TorQuadModule, AutomorphismGroupElem

Given an integral lattice with isometry $(L, f)$, return the discriminant group `q`
of the underlying lattice `L` as well as this image of the underlying isometry
`f` inside $O(q)$.
"""
function discriminant_group(Lf::ZZLatWithIsom)
  @req is_integral(Lf) "Underlying lattice must be integral"
  L = lattice(Lf)
  f = ambient_isometry(Lf)
  q = discriminant_group(L)
  Oq = orthogonal_group(q)
  return (q, Oq(gens(matrix_group(f))[1], check = false))::Tuple{TorQuadModule, AutomorphismGroupElem{TorQuadModule}}
end

@doc raw"""
    image_centralizer_in_Oq(Lf::ZZLatWithIsom) -> AutomorphismGroup{TorQuadModule}

Given an integral lattice with isometry $(L, f)$, return the image $G_L$ in
$O(q_L, \bar{f})$ of the centralizer $O(L, f)$ of `f` in $O(L)$. Here $q_L$
denotes the discriminant group of `L` and $\bar{f}$ is the isometry of
$q_L$ induced by `f`.
"""
@attr AutomorphismGroup{TorQuadModule} function image_centralizer_in_Oq(Lf::ZZLatWithIsom)
  @req is_integral(Lf) "Underlying lattice must be integral"
  n = order_of_isometry(Lf)
  L = lattice(Lf)
  f = ambient_isometry(Lf)
  if (n in [1, -1]) || (isometry(Lf) == -identity_matrix(QQ, rank(L)))
    GL, _ = image_in_Oq(L)
  elseif is_definite(L) 
    OL = orthogonal_group(L)
    f = OL(f)
    UL = QQMatrix[matrix(OL(s)) for s in gens(centralizer(OL, f)[1])]
    qL = discriminant_group(L)
    UL = ZZMatrix[hom(qL, qL, elem_type(qL)[qL(lift(t)*g) for t in gens(qL)]).map_ab.map for g in UL]
    unique!(UL)
    GL = Oscar._orthogonal_group(qL, UL, check = false)
  elseif rank(L) == euler_phi(n)
    qL = discriminant_group(L)
    UL = ZZMatrix[hom(qL, qL, elem_type(qL)[qL(lift(t)*g) for t in gens(qL)]).map_ab.map for g in [-f^0, f]]
    unique!(UL)
    GL = Oscar._orthogonal_group(qL, UL, check = false)
  else
    @req is_of_hermitian_type(Lf) "Not yet implemented for indefinite lattices with isometry which are not of hermitian type"
    dets = Oscar._local_determinants_morphism(Lf)
    GL, _ = kernel(dets)
  end
  return GL::AutomorphismGroup{TorQuadModule}
end

###############################################################################
#
#  Signatures
#
###############################################################################

function _real_kernel_signatures(L::ZZLat, M) 
  C = base_ring(M)
  bL = basis_matrix(L)
  GL = gram_matrix(ambient_space(L))
  bLC = change_base_ring(C, bL)
  GLC = change_base_ring(C, GL)
  k, KC = left_kernel(M)
  newGC = KC*bLC*GLC*transpose(KC*bLC)

  newGC = Hecke._gram_schmidt(newGC, C)[1]
  diagC = diagonal(newGC)

  @hassert :ZZLatWithIsom 1 all(z -> isreal(z), diagC)
  @hassert :ZZLatWithIsom 1 all(z -> !iszero(z), diagC)

  k1 = count(z -> z > 0, diagC)
  k2 = length(diagC) - k1

  return (k1, k2)
end

@doc raw"""
    signatures(Lf::ZZLatWithIsom) -> Dict{Int, Tuple{Int, Int}}

Given a lattice with isometry $(L, f)$ where the minimal polynomial of `f`
is irreducible cyclotomic, return the signatures of $(L, f)$.

In this context, if we denote $z$ a primitive `n`-th root of unity, where `n`
is the order of `f`, then for each $1 \leq i \leq n/2$ such that $(i, n) = 1$,
the $i$-th signature of $(L, f)$ is given by the signatures of the real quadratic
form $\Ker(f + f^{-1} - z^i - z^{-i})$.
"""
function signatures(Lf::ZZLatWithIsom)
  @req is_of_hermitian_type(Lf) "Lf must be of hermitian type"
  L = lattice(Lf)
  f = isometry(Lf)
  n = order_of_isometry(Lf)
  C = CalciumField()
  eig = eigenvalues(f, QQBar)
  j = findfirst(z -> findfirst(k -> isone(z^k), 1:n) == n, eig)
  lambda = C(eig[j])
  Sq = [i for i in 1:div(n,2) if gcd(i,n) == 1]
  D = Dict{Integer, Tuple{Int64, Int64}}()
  fC = change_base_ring(C, f)
  for i in Sq
    M = fC + inv(fC) - lambda^i - lambda^(-i)
    D[i] = _real_kernel_signatures(L, M)
  end
  return D
end

###############################################################################
#
#  Kernel sublattices
#
###############################################################################

function _divides(k::IntExt, n::Int)
  is_finite(k) && return Hecke.divides(k, n)[1]
  return true
end

@doc raw"""
    kernel_lattice(Lf::ZZLatWithIsom, p::Union{fmpz_poly, QQPolyRingElem})
                                                         -> ZZLatWithIsom

Given a lattice with isometry $(L, f)$ and a polynomial `p` with rational
coefficients, return the sublattice $M := \ker(p(f))$ of the underlying lattice
`L` with isometry `f`, together with the restriction $f_{\mid M}$.
"""
kernel_lattice(Lf::ZZLatWithIsom, p::Union{ZZPolyRingElem, QQPolyRingElem})

function kernel_lattice(Lf::ZZLatWithIsom, p::QQPolyRingElem)
  n = order_of_isometry(Lf)
  L = lattice(Lf)
  f = isometry(Lf)
  M = p(f)
  d = denominator(M)
  k, K = left_kernel(change_base_ring(ZZ, d*M))
  L2 = lattice_in_same_ambient_space(L, K*basis_matrix(L))
  f2 = solve_left(change_base_ring(QQ, K), K*f)
  @hassert :ZZLatWithIsom 1 f2*gram_matrix(L2)*transpose(f2) == gram_matrix(L2)
  chi = parent(p)(collect(coefficients(minpoly(f2))))
  chif = parent(p)(collect(coefficients(minpoly(Lf))))
  _chi = gcd(p, chif)
  @hassert :ZZLatWithIsom 1 (rank(L2) == 0) || (chi == _chi)
  return lattice_with_isometry(L2, f2, ambient_representation = false)
end

kernel_lattice(Lf::ZZLatWithIsom, p::ZZPolyRingElem) = kernel_lattice(Lf, change_base_ring(QQ, p))

@doc raw"""
    kernel_lattice(Lf::ZZLatWithIsom, l::Integer) -> ZZLatWithIsom

Given a lattice with isometry $(L, f)$ and an integer `l`, return the kernel
lattice of $(L, f)$ associated to the `l`-th cyclotomic polynomial.
"""
function kernel_lattice(Lf::ZZLatWithIsom, l::Integer)
  @req _divides(order_of_isometry(Lf), l)[1] "l must divide the order of the underlying isometry"
  p = cyclotomic_polynomial(l)
  return kernel_lattice(Lf, p)
end

@doc raw"""
    invariant_lattice(Lf::ZZLatWithIsom) -> ZZLatWithIsom

Given a lattice with isometry $(L, f)$, return the invariant lattice $L^f$ of
$(L, f)$ together with the restriction of `f` to $L^f$ (which is the identity
in this case)
"""
invariant_lattice(Lf::ZZLatWithIsom) = kernel_lattice(Lf, 1)

@doc raw"""
    coinvariant_lattice(Lf::ZZLatWithIsom) -> ZZLatWithIsom

Given a lattice with isometry $(L, f)$, return the coinvariant lattice $L_f$ of
$(L, f)$ together with the restriction of `f` to $L_f$.

The coinvariant lattice $L_f$ of $(L, f)$ is the orthogonal complement in
`L` of the invariant lattice $L_f$.
"""
function coinvariant_lattice(Lf::ZZLatWithIsom)
  chi = minpoly(Lf)
  if chi(1) == 0
    R = parent(chi)
    x = gen(R)
    chi = divexact(chi, x-1)
  end
  return kernel_lattice(Lf, chi)
end

###############################################################################
#
#  Type
#
###############################################################################

@doc raw"""
    type(Lf::ZZLatWithIsom)
                      -> Dict{Int, Tuple{ <: Union{ZZGenus, HermGenus}, ZZGenus}}

Given a lattice with isometry $(L, f)$ with `f` of finite order `n`, return the
type of $(L, f)$.

In this context, the type is defined as follows: for each divisor `k` of `n`,
the `k`-type of $(L, f)$ is the tuple $(H_k, A_K)$ consisting of the genus
$H_k$ of the lattice $\Ker(\Phi_k(f))$ viewed as a hermitian $\mathbb{Z}[\zeta_k]$-
lattice (so a $\mathbb{Z}$-lattice for k= 1, 2) and of the genus $A_k$ of the
$\mathbb{Z}$-lattice $\Ker(f^k-1)$.
"""
@attr function type(Lf::ZZLatWithIsom)
  L = lattice(Lf)
  f = isometry(Lf)
  n = order_of_isometry(Lf)
  @req is_finite(n) "Isometry must be of finite order"
  divs = divisors(n)
  Qx = Hecke.Globals.Qx
  x = gen(Qx)
  t = Dict{Integer, Tuple}()
  for l in divs
    Hl = kernel_lattice(Lf, cyclotomic_polynomial(l))
    if !(order_of_isometry(Hl) in [-1,1,2])
      Hl = Hecke.hermitian_structure(lattice(Hl), isometry(Hl), check=false, ambient_representation=false)
    end
    Al = kernel_lattice(Lf, x^l-1)
    t[l] = (genus(Hl), genus(Al))
  end
  return t
end

@doc raw"""
    is_of_type(Lf::ZZLatWithIsom, t::Dict) -> Bool

Given a lattice with isometry $(L, f)$, return whether $(L, f)$ is of type `t`.
"""
function is_of_type(L::ZZLatWithIsom, t::Dict)
  @req is_finite(order_of_isometry(L)) "Type is defined only for finite order isometries"
  divs = sort(collect(keys(t)))
  x = gen(Hecke.Globals.Qx)
  for l in divs
    Hl = kernel_lattice(L, cyclotomic_polynomial(l))
    if !(order_of_isometry(Hl) in [-1, 1, 2])
      t[l][1] isa Hecke.HermGenus || return false
      Hl = Hecke.hermitian_structure(lattice(Hl), isometry(Hl), check=false, ambient_representation=false, E = base_field(t[l][1]))
    end
    genus(Hl) == t[l][1] || return false
    Al = kernel_lattice(L, x^l-1)
    genus(Al) == t[l][2] || return false
  end
  return true
end

@doc raw"""
    is_of_same_type(Lf::ZZLatWithIsom, Mg::ZZLatWithIsom) -> Bool

Given two lattices with isometry $(L, f)$ and $(M, g)$, return whether they are
of the same type.
"""
function is_of_same_type(L::ZZLatWithIsom, M::ZZLatWithIsom)
  @req is_finite(order_of_isometry(L)*order_of_isometry(M)) "Type is defined only for finite order isometries"
  order_of_isometry(L) != order_of_isometry(M) && return false
  genus(L) != genus(M) && return false
  return is_of_type(L, type(M))
end

@doc raw"""
    is_hermitian(t::Dict) -> Bool

Given a type `t` of lattices with isometry, return whether `t` is hermitian, i.e.
whether it defines the type of a hermitian lattice with isometry.
"""
function is_hermitian(t::Dict)
  ke = collect(keys(t))
  n = maximum(ke)
  return all(i -> rank(t[i][1]) == rank(t[i][2]) == 0, [i for i in ke if i != n])
end

###############################################################################
#
#  Useful
#
###############################################################################

function to_oscar(Lf::ZZLatWithIsom)
  L = lattice(Lf)
  f = ambient_isometry(Lf)
  println(stdout, "B = matrix(QQ, $(rank(L)), $(degree(L)), " , basis_matrix(L), " );")
  println(stdout, "G = matrix(QQ, $(degree(L)), $(degree(L)), ", gram_matrix(ambient_space(L)), " );")
  println(stdout, "L = integer_lattice(B, gram = G);")
  println(stdout, "f = matrix(QQ, $(degree(L)), $(degree(L)), ", f, " );")
  println(stdout, "Lf = lattice_with_isometry(L, f);")
end

