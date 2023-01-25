export ambient_isometry,
       coinvariant_lattice,
       hermitian_structure,
       image_centralizer_in_Oq,
       isometry,
       is_of_pure_type,
       is_of_same_type,
       is_of_type,
       is_pure,
       lattice_with_isometry,
       order_of_isometry,
       type
    
###############################################################################
#
#  String I/O
#
###############################################################################

function Base.show(io::IO,  Lf::LatticeWithIsometry)
  println(io, "Lattice of rank $(rank(Lf)) with isometry of order $(order_of_isometry(Lf))")
  println(io, "given by")
  print(IOContext(io, :compact => true), isometry(Lf))
end

###############################################################################
#
#  Accessors
#
###############################################################################

@doc Markdown.doc"""
    lattice(Lf::LatticeWithIsometry) -> ZLat

Given a lattice with isometry $(L, f)$, return the underlying lattice `L`.
"""
lattice(Lf::LatticeWithIsometry) = Lf.Lb

@doc Markdown.doc"""
    isometry(Lf::LatticeWithIsometry) -> fmpq_mat

Given a lattice with isometry $(L, f)$, return the underlying isometry `f`.
"""
isometry(Lf::LatticeWithIsometry) = Lf.f

@doc Markdown.doc"""
   ambient_isometry(Lf::LatticeWithIsometry) -> fmpq_mat

Given a lattice with isometry $(L, f)$, return an isometry of underlying isometry
of the ambient space of `L` inducing `f` on `L`
"""
ambient_isometry(Lf::LatticeWithIsometry) = Lf.f_ambient

@doc Markdown.doc"""
    order_of_isometry(Lf::LatticeWithIsometry) -> Integer

Given a lattice with isometry $(L, f)$, return the order of the underlying
isometry `f`.
"""
order_of_isometry(Lf::LatticeWithIsometry) = Lf.n

###############################################################################
#
#  Attributes
#
###############################################################################

@doc Markdown.doc"""
    rank(Lf::LatticeWithIsometry) -> Integer

Given a lattice with isometry $(L, f)$, return the rank of the underlying lattice
`L`.
"""
rank(Lf::LatticeWithIsometry) = rank(lattice(Lf))::Integer

@doc Markdown.doc"""
    charpoly(Lf::LatticeWithIsometry) -> fmpq_poly

Given a lattice with isometry $(L, f)$, return the characteristic polynomial of the
underlying isometry `f`.
"""
charpoly(Lf::LatticeWithIsometry) = charpoly(isometry(Lf))::fmpq_poly

@doc Markdown.doc"""
    minpoly(Lf::LatticeWithIsometry) -> fmpq_poly

Given a lattice with isometry $(L, f)$, return the minimal polynomial of the
underlying isometry `f`.
"""
minpoly(Lf::LatticeWithIsometry) = minpoly(isometry(Lf))::fmpq_poly

@doc Markdown.doc"""
    genus(Lf::LatticeWithIsometry) -> ZGenus

Given a lattice with isometry $(L, f)$, return the genus of the underlying 
lattice `L`.
"""
genus(Lf::LatticeWithIsometry) = begin; L = lattice(Lf); is_integral(L) ? genus(L)::ZGenus : error("Underlying lattice must be integral"); end

@doc Markdown.doc"""
    ambient_space(Lf::LatticeWithIsometry) -> QuadSpace

Given a lattice with isometry $(L, f)$, return the ambient space of the underlying
lattice `L`.
"""
ambient_space(Lf::LatticeWithIsometry) = ambient_space(lattice(Lf))::Hecke.QuadSpace{FlintRationalField, fmpq_mat}

###############################################################################
#
#  Constructor
#
###############################################################################

@doc Markdown.doc"""
    lattice_with_isometry(L::ZLat, f::fmpq_mat, n::IntExt; check::Bool = true,
                                                            ambient_representation = true)
				                                                                -> LatticeWithIsometry

Given a $\mathbb Z$-lattice `L` and a matrix `f`, if `f` defines an isometry
of `L` of order `n`, return the corresponding lattice with isometry pair $(L, f)$.

If `ambient_representation` is set to true, `f` is consider as an isometry of the
ambient space of `L` and the induced isometry on `L` is automatically computed.
Otherwise, an isometry of the ambient space of `L` is constructed, setting the identity
on the complement of the rational span of `L` if it is not of full rank.
"""
function lattice_with_isometry(L::ZLat, f::fmpq_mat, n::IntExt; check::Bool = true,
                                                                 ambient_representation::Bool = true)
  if rank(L) == 0
    return LatticewithIsometry(L, matrix(QQ,0,0,[]), -1)
  end

  if check
    @req det(f) != 0 "f is not invertible"
    m = multiplicative_order(f)
    @req n == m "The order of f is equal to $m, not $n"
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

  if check
    @req f*gram_matrix(L)*transpose(f) == gram_matrix(L) "f does not define an isometry of L"
    @req f_ambient*gram_matrix(ambient_space(L))*transpose(f_ambient) == gram_matrix(ambient_space(L)) "f_ambient is not an isometry of the ambient space of L"
    @assert basis_matrix(L)*f_ambient == f*basis_matrix(L)
  end

  return LatticeWithIsometry(L, f, f_ambient, n)::LatticeWithIsometry
end

@doc Markdown.doc"""
    lattice_with_isometry(L::ZLat, f::fmpq_mat; check::Bool = true,
                                                ambient_representation::Bool, = true)
                                                                  -> LatticeWithIsometry

Given a $\mathbb Z$-lattice `L` and a matrix `f`, if `f` defines an isometry
of `L`, return the corresponding lattice with isometry pair $(L, f)$.

If `ambient_representation` is set to true, `f` is consider as an isometry of the
ambient space of `L` and the induced isometry on `L` is automatically computed.
Otherwise, an isometry of the ambient space of `L` is constructed, setting the identity
on the complement of the rational span of `L` if it is not of full rank.
"""
function lattice_with_isometry(L::ZLat, f::fmpq_mat; check::Bool = true,
                                                     ambient_representation::Bool = true)
  if rank(L) == 0
    return LatticeWithIsometry(L, matrix(QQ,0,0,[]), identity_matrix(QQ, degree(L)), -1)
  end

  n = multiplicative_order(f)
  return lattice_with_isometry(L, f, n, check = check,
                               ambient_representation = ambient_representation)::LatticeWithIsometry
end

@doc Markdown.doc"""
    lattice_with_isometry(L::ZLat) -> LatticeWithIsometry

Given a $\mathbb Z$-lattice `L` return the lattice with isometry pair $(L, f)$,
where `f` corresponds to the identity mapping of `L`.
"""
function lattice_with_isometry(L::ZLat)
  d = degree(L)
  f = identity_matrix(QQ, d)
  return lattice_with_isometry(L, f, check = false, ambient_representation = true)::LatticeWithIsometry
end

###############################################################################
#
#  Hermitian structure
#
###############################################################################

@doc Markdown.doc"""
    hermitian_structure(Lf::LatticeWithIsometry; check::Bool = true) -> HermLat

Given a lattice with isometry $(L, f)$ such that the minimal polynomial of the
underlying isometry `f` is irreducible and cyclotomic, return the
hermitian structure of the underlying lattice `L` over the $n$th cyclotomic
field, where $n$ is the order of `f`.

If it exists, the hermitian structure is cached.
"""
@attr HermLat function hermitian_structure(Lf::LatticeWithIsometry)
  @req rank(Lf) > 0 "Lf must be of positive rank"
  f = isometry(Lf)
  n = order_of_isometry(Lf)

  @req n >= 3 "No hermitian structures for n smaller than 3"
  @req is_cyclotomic_polynomial(minpoly(f)) "The minimal polynomial must be irreducible and cyclotomic"

  return Oscar._hermitian_structure(lattice(Lf), f, n = n, check = false,
                                                     ambient_representation = false)
end

###############################################################################
#
#  Image centralizer in discriminant group
#
###############################################################################

@doc Markdown.doc"""
    discriminant_group(Lf::LatticeWithIsometry) -> TorQuadMod, AutomorphismGroupElem

Given an integral lattice with isometry $(L, f)$, return the discriminant group `q`
of the underlying lattice `L` as well as this image of the underlying isometry
`f` inside $O(q)$.
"""
function discriminant_group(Lf::LatticeWithIsometry)
  L = lattice(Lf)
  f = ambient_isometry(Lf)
  @req is_integral(L) "Underlying lattice must be integral"
  q = discriminant_group(L)
  Oq = orthogonal_group(q)
  return (q, Oq(gens(matrix_group(f))[1], check = false))::Tuple{TorQuadMod, AutomorphismGroupElem{TorQuadMod}}
end

@doc Markdown.doc"""
    image_centralizer_in_Oq(Lf::LatticeWithIsometry) -> AutomorphismGroup

Given an integral lattice with isometry $(L, f)$, return the image $G_L$ in
$O(q_L, \bar{f})$ of the centralizer $O(L, f)$ of `f` in $O(L)$. Here $q_L$
denotes the discriminant group of `L` and $\bar{f}$ is the isometry of
$q_L$ induced by `f`.
"""
@attr AutomorphismGroup{TorQuadMod} function image_centralizer_in_Oq(Lf::LatticeWithIsometry)
  n = order_of_isometry(Lf)
  L = lattice(Lf)
  f = ambient_isometry(Lf)
  @req is_integral(L) "Underlying lattice must be integral"
  if n in [1, -1]
    GL, _ = image_in_Oq(L)
  elseif is_definite(L)
    OL = orthogonal_group(L)
    f = OL(f)
    UL = fmpq_mat[matrix(OL(s)) for s in gens(centralizer(OL, f)[1])]
    qL = discriminant_group(L)
    UL = fmpz_mat[hom(qL, qL, elem_type(qL)[qL(lift(t)*g) for t in gens(qL)]).map_ab.map for g in UL]
    GL = Oscar._orthogonal_group(qL, UL)
  elseif rank(L) == euler_phi(n)
    qL = discriminant_group(L)
    UL = fmpz_mat[hom(qL, qL, elem_type(qL)[qL(lift(t)*g) for t in gens(qL)]).map_ab.map for g in [-f^0, f]]
    GL = Oscar._orthogonal_group(qL, UL)
  else
    qL, fqL = discriminant_group(Lf)
    OqL = orthogonal_group(qL)
    CdL, _ =  centralizer(OqL, fqL)
    GL, _ = sub(OqL, [OqL(s.X) for s in CdL])
  end
  return GL::AutomorphismGroup{TorQuadMod}
end

###############################################################################
#
#  Signatures
#
###############################################################################

function _real_kernel_signatures(L::ZLat, M) 
  C = base_ring(M)
  bL = basis_matrix(L)
  GL = gram_matrix(ambient_space(L))
  bLC = change_base_ring(C, bL)
  GLC = change_base_ring(C, GL)
  k, KC = left_kernel(M)
  newGC = KC*bLC*GLC*transpose(KC*bLC)

  newGC = Hecke._gram_schmidt(newGC, C)[1]
  diagC = diagonal(newGC)

  @assert all(z -> isreal(z), diagC)
  @assert all(z -> !iszero(z), diagC)

  k1 = count(z -> z > 0, diagC)
  k2 = length(diagC) - k1

  return (k1, k2)
end

@doc Markdown.doc"""
    signatures(Lf::LatticeWithIsometry) -> Dict{Int, Tuple{Int, Int}}

Given a lattice with isometry $(L, f)$ where the minimal polynomial of `f`
is irreducible cyclotomic, return the signatures of $(L, f)$.

In this context, if we denote $z$ a primitive `n`-th root of unity, where `n`
is the order of `f`, then for each $1 \leq i \leq n/2$ such that $(i, n) = 1$,
the $i$-th signature of $(L, f)$ is given by the signatures of the real quadratic
form $\Ker(f + f^{-1} - z^i - z^{-i})$.
"""
function signatures(Lf::LatticeWithIsometry)
  @req rank(Lf) != 0 "Signatures non available for the empty lattice"
  L = lattice(Lf)
  f = isometry(Lf)
  @req is_cyclotomic_polynomial(minpoly(f)) "Minimal polynomial must be irreducible and cyclotomic"
  n = order_of_isometry(Lf)
  @req divides(rank(L), euler_phi(n))[1] "The totient of the order of the underlying isometry must divide the rank of the underlying lattice"
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

divides(k::PosInf, n::Int) = true

@doc Markdown.doc"""
    kernel_lattice(Lf::LatticeWithIsometry, p::Union{fmpz_poly, fmpq_poly})
                                                         -> LatticeWithIsometry

Given a lattice with isometry $(L, f)$ and a polynomial `p` with rational
coefficients, return the sublattice $M := \ker(p(f))$ of the underlying lattice
`L` with isometry `f`, together with the restriction $f_{\mid M}$.
"""
kernel_lattice(Lf::LatticeWithIsometry, p::Union{fmpz_poly, fmpq_poly})

function kernel_lattice(Lf::LatticeWithIsometry, p::fmpq_poly)
  n = order_of_isometry(Lf)
  L = lattice(Lf)
  f = isometry(Lf)
  M = p(f)
  d = denominator(M)
  k, K = left_kernel(change_base_ring(ZZ, d*M))
  L2 = lattice_in_same_ambient_space(L, K*basis_matrix(L))
  f2 = solve_left(change_base_ring(QQ, K), K*f)
  @assert f2*gram_matrix(L2)*transpose(f2) == gram_matrix(L2)
  chi = parent(p)(collect(coefficients(minpoly(f2))))
  chif = parent(p)(collect(coefficients(minpoly(Lf))))
  _chi = gcd(p, chif)
  @assert (rank(L2) == 0) || (chi == _chi)
  return lattice_with_isometry(L2, f2, check = true, ambient_representation = false)
end

kernel_lattice(Lf::LatticeWithIsometry, p::fmpz_poly) = kernel_lattice(Lf, change_base_ring(QQ, p))

@doc Markdown.doc"""
    kernel_lattice(Lf::LatticeWithIsometry, l::Integer) -> LatticeWithIsometry

Given a lattice with isometry $(L, f)$ and an integer `l`, return the kernel
lattice of $(L, f)$ associated to the `l`-th cyclotomic polynomial.
"""
function kernel_lattice(Lf::LatticeWithIsometry, l::Integer)
  @req divides(order_of_isometry(Lf), l)[1] "l must divide the order of the underlying isometry"
  p = cyclotomic_polynomial(l)
  return kernel_lattice(Lf, p)
end

@doc Markdown.doc"""
    invariant_lattice(Lf::LatticeWithIsometry) -> LatticeWithIsometry

Given a lattice with isometry $(L, f)$, return the invariant lattice $L^f$ of
$(L, f)$ together with the restriction of `f` to $L^f$ (which is the identity
in this case)
"""
invariant_lattice(Lf::LatticeWithIsometry) = kernel_lattice(Lf, 1)

@doc Markdown.doc"""
    coinvariant_lattice(Lf::LatticeWithIsometry) -> LatticeWithIsometry

Given a lattice with isometry $(L, f)$, return the coinvariant lattice $L_f$ of
$(L, f)$ together with the restriction of `f` to $L_f$.

The coinvariant lattice $L_f$ of $(L, f)$ is the orthogonal complement in
`L` of the invariant lattice $L_f$.
"""
function coinvariant_lattice(Lf::LatticeWithIsometry)
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

@doc Markdown.doc"""
    type(Lf::LatticeWithIsometry)
                      -> Dict{Int, Tuple{ <: Union{ZGenus, GenusHerm}, ZGenus}}

Given a lattice with isometry $(L, f)$ with `f` of finite order `n`, return the
type of $(L, f)$.

In this context, the type is defined as follows: for each divisor `k` of `n`,
the `k`-type of $(L, f)$ is the tuple $(H_k, A_K)$ consisting of the genus
$H_k$ of the lattice $\Ker(\Phi_k(f))$ viewed as a hermitian $\mathbb{Z}[\zeta_k]$-
lattice (so a $\mathbb{Z}$-lattice for k= 1, 2) and of the genus $A_k$ of the
$\mathbb{Z}$-lattice $\Ker(f^k-1)$.
"""
@attr function type(Lf::LatticeWithIsometry)
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
      Hl = Oscar._hermitian_structure(lattice(Hl), isometry(Hl), n=order_of_isometry(Hl), check=false, ambient_representation=false)
    end
    Al = kernel_lattice(Lf, x^l-1)
    t[l] = (genus(Hl), genus(Al))
  end
  return t
end

@doc Markdown.doc"""
    is_of_type(Lf::LatticeWithIsometry, t::Dict) -> Bool

Given a lattice with isometry $(L, f)$, return whether $(L, f)$ is of type `t`.
"""
function is_of_type(L::LatticeWithIsometry, t::Dict)
  @req is_finite(order_of_isometry(L)) "Type is defined only for finite order isometries"
  divs = sort(collect(keys(t)))
  x = gen(Hecke.Globals.Qx)
  for l in divs
    Hl = kernel_lattice(L, cyclotomic_polynomial(l))
    if !(order_of_isometry(Hl) in [-1, 1, 2])
      t[l][1] isa Hecke.GenusHerm || return false
      Hl = Oscar._hermitian_structure(lattice(Hl), isometry(Hl), n=order_of_isometry(Hl), check=false, ambient_representation=false, E = base_field(t[l][1]))
    end
    genus(Hl) == t[l][1] || return false
    Al = kernel_lattice(L, x^l-1)
    genus(Al) == t[l][2] || return false
  end
  return true
end

@doc Markdown.doc"""
    is_of_same_type(Lf::LatticeWithIsometry, Mg::LatticeWithIsometry) -> Bool

Given two lattices with isometry $(L, f)$ and $(M, g)$, return whether they are
of the same type.
"""
function is_of_same_type(L::LatticeWithIsometry, M::LatticeWithIsometry)
  @req is_finite(order_of_isometry(L)*order_of_isometry(M)) "Type is defined only for finite order isometries"
  order_of_isometry(L) != order_of_isometry(M) && return false
  genus(L) != genus(M) && return false
  return is_of_type(L, type(M))
end

@doc Markdown.doc"""
    is_of_pure_type(Lf::LatticeWithIsometry) -> Bool

Given a lattice with isometry $(L, f)$, return whether the minimal polynomial
of `f` is irreducible cyclotomic.
"""
function is_of_pure_type(L::LatticeWithIsometry)
  @req is_finite(order_of_isometry(L)) "Type is defined only for finite order isometries"
  return is_cyclotomic_polynomial(minpoly(L))
end

@doc Markdown.doc"""
    is_pure(t::Dict) -> Bool

Given a type `t` of lattices with isometry, return whether `t` is pure, i.e.
whether it defines the type of lattice with isometry whose minimal polynomial
is irreducible cyclotomic.
"""
function is_pure(t::Dict)
  ke = collect(keys(t))
  n = maximum(ke)
  return all(i -> rank(t[i][1]) == rank(t[i][2]) == 0, [i for i in ke if i != n])
end

