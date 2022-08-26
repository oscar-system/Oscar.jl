export hermitian_structure, isometry, lattice_with_isometry, order_of_isometry,
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
#  Attributes
#
###############################################################################

@doc Markdown.doc"""
    lattice(Lf::LatticeWithIsometry) -> ZLat

Given a lattice with isometry `Lf`, return the underlying lattice `L`.
"""
lattice(Lf::LatticeWithIsometry) = Lf.Lb

@doc Markdown.doc"""
    isometry(Lf::LatticeWithIsometry) -> fmpq_mat

Given a lattice with isometry `Lf`, return the underlying isometry `f`.
"""
isometry(Lf::LatticeWithIsometry) = Lf.f

@doc Markdown.doc"""
    order_of_isometry(Lf::LatticeWithIsometry) -> Integer

Given a lattice with isometry `Lf`, return the order of the underlying
isometry `f`.
"""
order_of_isometry(Lf::LatticeWithIsometry) = Lf.n

###############################################################################
#
#  Predicates
#
###############################################################################

@doc Markdown.doc"""
    rank(Lf::LatticeWithIsometry) -> Integer

Given a lattice with isometry `Lf`, return the rank of the underlying lattice
`L`.
"""
rank(Lf::LatticeWithIsometry) = rank(lattice(Lf))

@doc Markdown.doc"""
    charpoly(Lf::LatticeWithIsometry) -> fmpq_poly

Given a lattice with isometry `Lf`, return the characteristic polynomial of the
underlying isometry `f`.
"""
charpoly(Lf::LatticeWithIsometry) = charpoly(isometry(Lf))

@doc Markdown.doc"""
    minpoly(Lf::LatticeWithIsometry) -> fmpq_poly

Given a lattice with isometry `Lf`, return the minimal polynomial of the
underlying isometry `f`.
"""
minpoly(Lf::LatticeWithIsometry) = minpoly(isometry(Lf))

@doc Markdown.doc"""
    genus(Lf::LatticeWithIsometry) -> ZGenus

Given a lattice with isometry `Lf`, return the genus of the underlying 
lattice `L`.

For now, in order for the genus to exist, the lattice must be integral.
"""
genus(Lf::LatticeWithIsometry) = begin; L = lattice(Lf); is_integral(L) ? genus(L) : error("Underlying lattice must be integral"); end


###############################################################################
#
#  Constructor
#
###############################################################################

@doc Markdown.doc"""
    lattice_with_isometry(L::ZLat, f::fmpq_mat, n::Integer; check::Bool = true)
				                      -> LatticeWithIsometry

Given a $\mathbb Z$-lattice `L` and a matrix `f`, if `f` defines an isometry
of `L` of finite order `n`, return the corresponding lattice with isometry
pair $(L, f)$.
"""
function lattice_with_isometry(L::ZLat, f::fmpq_mat, n::Integer; check::Bool = true)
  if rank(L) == 0
    return LatticewithIsometry(L, matrix(QQ,0,0,[]), -1)
  end

  if check
    @req det(f) != 0 "f is not invertible"
    G = gram_matrix(L)
    @req f*G*transpose(f) == G "f does not define an isometry of L"
    m = Oscar._exponent(f)
    @req m > 0 "f is not of finite exponent"
    @req n == m "The order of f is not equal to n"
  end

  return LatticeWithIsometry(L, f, n)
end

@doc Markdown.doc"""
    lattice_with_isometry(L::ZLat, f::fmpq_mat; check::Bool = true)
                                                      -> LatticeWithIsometry

Given a $\mathbb Z$-lattice `L` and a matrix `f`, if `f` defines an isometry
of `L` of finite order, return the corresponding lattice with isometry
pair $(L, f)$.
"""
function lattice_with_isometry(L::ZLat, f::fmpq_mat; check::Bool = true)
  if rank(L) == 0
    return LatticeWithIsometry(L, matrix(QQ,0,0,[]), -1)
  end

  if check
    @req det(f) != 0 "f is not invertibe"
    G = gram_matrix(L)
    @req f*G*transpose(f) == G "f does not define an isometry of L"
    m = Oscar._exponent(f)
    @req m > 0 "f is not of finite exponent"
  end

  n = Oscar._exponent(f)
  return lattice_with_isometry(L, f, Int(n), check = false)
end

@doc Markdown.doc"""
    lattice_with_isometry(L::ZLat) -> LatticeWithIsometry

Given a $\mathbb Z$-lattice `L` return the lattice with isometry pair $(L, f)$,
where `f` corresponds to the identity mapping of `L`.
"""
function lattice_with_isometry(L::ZLat)
  r = rank(L)
  f = identity_matrix(QQ, r)
  return lattice_with_isometry(L, f, check = false)
end

###############################################################################
#
#  Hermitian structure
#
###############################################################################

@doc Markdown.doc"""
    hermitian_structure(Lf::LatticeWithIsometry; check::Bool = true) -> HermLat

Given a lattice with isometry `Lf` such that the minimal polynomial of the
underlying isometry `f` is irreducible and cyclotomic, return the
hermitian structure of the underlying lattice `L` over the $n$th cyclotomic
field, where $n$ is the order of `f`.

If it exists, the hermitian structure is cached.
"""
@attr function hermitian_structure(Lf::LatticeWithIsometry)
  @req rank(Lf) > 0 "Lf must be of positive rank"
  f = isometry(Lf)
  n = order_of_isometry(Lf)
  
  @req n >= 3 "No hermitian structures for n smaller than 3"

  return inverse_trace_lattice(lattice(Lf), f, n = n, check = true)
end

###############################################################################
#
#  Image centralizer in discriminant group
#
###############################################################################

@doc Markdown.doc"""
    discriminant_group(Lf::LatticeWithIsometry) -> TorQuadMod, TorQuadModMor

Given an integral lattice with isometry `Lf`, return the discriminant group `q`
of the underlying lattice `L` as well as this image of the underlying isometry
`f` inside $O(q)$.
"""
function discriminant_group(Lf::LatticeWithIsometry)
  L = lattice(Lf)
  f = isometry(Lf)
  @req is_integral(L) "Underlying lattice must be integral"
  q = discriminant_group(L)
  Oq = orthogonal_group(q)
  return q, Oq(gens(matrix_group(f))[1])
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

function signatures(Lf::LatticeWithIsometry)
  @req rank(Lf) != 0 "Signatures non available for the empty lattice"
  L = lattice(Lf)
  f = isometry(Lf)
  @req Oscar._is_cyclotomic_polynomial(minpoly(f)) "Minimal polynomial must be irreducible and cyclotomic"
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

@doc Markdown.doc"""
    kernel_lattice(Lf::LatticeWithIsometry, p::fmpq_poly) -> LatticeWithIsometry

Given a lattice with isometry `Lf` and a polynomial `p` with rational
coefficients, return the sublattice $M := \ker(p(f))$ of the underlying lattice
`L` with isometry `f`, together with the restriction $f_{\mid M}$.
"""
function kernel_lattice(Lf::LatticeWithIsometry, p::fmpq_poly)
  n = order_of_isometry(Lf)
  L = lattice(Lf)
  f = isometry(Lf)
  M = p(f)
  k, K = left_kernel(change_base_ring(ZZ, M))
  L2 = Zlattice(gram = K*gram_matrix(L)*transpose(K))
  f2 = solve_left(change_base_ring(QQ, K), K*f)
  @assert f2*gram_matrix(L2)*transpose(f2) == gram_matrix(L2)
  chi = parent(p)(collect(coefficients(minpoly(f2))))
  chif = parent(p)(collect(coefficients(minpoly(Lf))))
  _chi = gcd(p, chif)
  @assert (rank(L2) == 0) || (chi == _chi)
  return lattice_with_isometry(L2, f2, check = false)
end

kernel_lattice(Lf::LatticeWithIsometry, p::fmpz_poly) = kernel_lattice(Lf, change_base_ring(QQ, p))

function kernel_lattice(Lf::LatticeWithIsometry, l::Integer)
  @req divides(order_of_isometry(Lf), l)[1] "l must divide the order of the underlying isometry"
  p = Oscar._cyclotomic_polynomial(l)
  return kernel_lattice(Lf, p)
end

###############################################################################
#
#  Type
#
###############################################################################

@attr function type(Lf::LatticeWithIsometry)
  L = lattice(Lf)
  f = isometry(Lf)
  n = order_of_isometry(Lf)
  divs = divisors(n)
  Qx = Hecke.Globals.Qx
  x = gen(Qx)
  t = Dict{Integer, Tuple{LatticeWithIsometry, LatticeWithIsometry}}()
  for l in divs
    Hl = kernel_lattice(Lf, Oscar._cyclotomic_polynomial(l))
    Al = kernel_lattice(Lf, x^l-1)
    t[l] = (Hl, Al)
  end
  return t
end

