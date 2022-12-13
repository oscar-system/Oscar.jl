export weyl_vector, separating_hyperplanes, walls, K3_surface_automorphism_group,
      adjacent_chamber, aut, has_zero_entropy, K3Chamber, chamber, borcherds_method


################################################################################
# Types
################################################################################

@doc Markdown.doc"""
    BorcherdsCtx

Contains all the data necessary to run Borcherds' method.

The assumptions are as follows:

- `L::ZLat`: even, hyperbolic, unimodular `Z`-lattice of rank n = 10, 18 or 26
  with basis given by the n x n standard basis

- `S::ZLat`: primitive sublattice

- `weyl_vector::fmpz_mat`: given in the basis of `L`

- `R::ZLat`: the orthogonal complement of `S` in `L`
  We assume that the basis of R consists of the last rank R standard basis vectors.
- `SS::ZLat`: lattice with standard basis and same gram matrix as `S`.

- `deltaR::Vector{fmpz_mat}`:

- `dualDeltaR::Vector{fmpz_mat}`:

- `prRdelta::Vector{Tuple{fmpq_mat,fmpq}}`:

- `membership_test`: takes the `s x s` matrix describing `g: S -> S` with respect
  to the basis of `S` and returns whether `g` lies in the group `G`.

- `prS::fmpq_mat`: `prS: L --> S^\vee` given with respect to the (standard)
  basis of `L` and the basis of `S`
"""
mutable struct BorcherdsCtx
  L::ZLat
  S::ZLat
  weyl_vector::fmpz_mat # given in the basis of L
  SS::ZLat
  R::ZLat
  deltaR::Vector{fmpz_mat}
  dualDeltaR::Vector{fmpz_mat}
  prRdelta::Vector{Tuple{fmpq_mat,fmpq}}
  membership_test
  gramL::fmpz_mat  # avoid a few conversions because gram_matrix(::ZLat) -> fmpq_mat
  gramS::fmpz_mat
  prS::fmpq_mat
  compute_OR::Bool
  # TODO: Store temporary variables for the computations
  # in order to make the core-functions adjacent_chamber and walls
  # as non-allocating as possible.
end

function Base.show(io::IOContext, d::BorcherdsCtx)
  print(io, "BorcherdsCtx: dim(L) = $(rank(d.L)), rank(S) = $(rank(d.S)), det(S) = $(det(d.S))")
end


@doc Markdown.doc"""
    BorcherdsCtx(L::ZLat, S::ZLat, compute_OR::Bool=true) -> BorcherdsCtx

Return the context for Borcherds' method.

# Arguments
- `L::ZLat`: an even, hyperbolic, unimodular `Z`-lattice of rank 10, 18, or 26
- `S::ZLat`: a primitive sublattice of `L` in the same ambient space
- if `compute_OR` is `false`, then `G` is the subgroup of the orthogonal group
  of `S` acting as $\pm 1$ on the discriminant group.
  If `compute_OR` is `true`, then `G` consists the subgroup consisting of
  isometries of `S` that can be extended to isometries of `L`.
"""
function BorcherdsCtx(L::ZLat, S::ZLat, weyl, compute_OR::Bool=true)
  r = rank(L)
  lw = (weyl*gram_matrix(L)*transpose(weyl))[1,1]
  if  r == 26
    @req lw == 0 "not a weyl vector"
  end
  if  r == 18
    @req lw == 620 "not a weyl vector"
  end
  if  r == 10
    @req lw == 1240 "not a weyl vector"
  end
  # transform L to have the standard basis
  # we assume that the basis of L is obtained by completing a basis of R
  # hence we can throw away the R coordinates of a weyl vector when projecting to S
  R = lll(Hecke.orthogonal_submodule(L, S))

  # the following completes the basis of R to a basis of L
  basisRL = solve_left(basis_matrix(L),basis_matrix(R))
  basisRL = change_base_ring(ZZ, basisRL)

  A, j = snf(abelian_group(basisRL))
  T = reduce(vcat, [j(i).coeff for i in gens(A)])
  basisL1 = vcat(T, basisRL)*basis_matrix(L)

  # carry the weyl vector along
  L1 = lattice(ambient_space(L), basisL1)
  weyl = change_base_ring(ZZ, solve_left(basisL1, weyl*basis_matrix(L)))
  basisSL1 = solve_left(basis_matrix(L1), basis_matrix(S))
  basisRL1 = solve_left(basis_matrix(L1), basis_matrix(R))

  # Assure that L has the standard basis.
  L = Zlattice(gram=gram_matrix(L1))
  V = ambient_space(L)
  # carry S along
  S = lattice(V, basisSL1)
  R = lattice(V, basisRL1)
  @req is_S_nondegenerate(L, S, change_base_ring(QQ,weyl)) "Weyl vector is S degenerate"

  SS = Zlattice(gram=gram_matrix(S))

  # precomputations

  @assert iszero(basis_matrix(R)[1:end,1:rank(S)])
  bSR = vcat(basis_matrix(S),basis_matrix(R))
  ibSR = inv(bSR)
  I = identity_matrix(QQ,degree(L))
  # prS: L --> S^\vee given with respect to the standard basis of L and the basis of S
  prS = ibSR*I[:,1:rank(S)]#*basis_matrix(S)
  @assert prS[rank(S)+1,:]==0


  if compute_OR
    dd = diagonal(gram_matrix(R))
    @vprint :K3Auto 2 "computing orthogonal group\n"
    OR = orthogonal_group(R)
    @vprint :K3Auto 2 "done\n"
    DR = discriminant_group(R)
    ODR = orthogonal_group(DR)
    imOR = [ODR(hom(DR,DR,[DR(lift(d)*f) for d in gens(DR)])) for f in gens(OR)]
    DS = discriminant_group(S)
    DSS = discriminant_group(SS)
    ODSS = orthogonal_group(DSS)
    orderimOR = order(sub(ODR,imOR)[1])
    @vprint :K3Auto 1 "[O(S):G] = $(order(ODSS)//orderimOR)\n"
    if order(ODR)== orderimOR
      membership_test = (g->true)
    else
      phiSS_S = hom(DSS,DS,[DS(lift(x)*basis_matrix(S)) for x in gens(DSS)])
      phi,i,j = glue_map(L,S,R)
      phi = phiSS_S*inv(i)*phi*j
      img,_ = sub(ODSS,[ODSS(phi*hom(g)*inv(phi)) for g in imOR])
      ds = degree(SS)
      membership_test = (g->ODSS(hom(DSS,DSS,[DSS(vec(matrix(QQ, 1, ds, lift(x))*g)) for x in gens(DSS)])) in img)
    end
  else
    membership_test(g) = is_pm1_on_discr(SS,g)
  end

  d = exponent(discriminant_group(S))
  Rdual = dual(R)
  sv = short_vectors(rescale(Rdual,-1), 2)
  # not storing the following for efficiency
  # append!(sv,[(-v[1],v[2]) for v in sv])
  # but for convenience we include zero
  T = typeof(sv).parameters[1].parameters[1].parameters[1]
  push!(sv,(zeros(T, rank(Rdual)), QQ(0)))
  rkR = rank(R)
  prRdelta = [(matrix(QQ, 1, rkR, v[1])*basis_matrix(Rdual),v[2]) for v in sv]
  gramL = change_base_ring(ZZ,gram_matrix(L))
  gramS = change_base_ring(ZZ,gram_matrix(S))
  deltaR = [change_base_ring(ZZ, matrix(QQ, 1, rkR, v[1])*basis_matrix(R)) for v in short_vectors(rescale(R,-1),2)]
  dualDeltaR = [gramL*transpose(r) for r in deltaR]
  return BorcherdsCtx(L, S, weyl, SS, R, deltaR, dualDeltaR, prRdelta, membership_test,
                      gramL, gramS, prS, compute_OR)
end

################################################################################
# Chambers
################################################################################

@doc Markdown.doc"""
    K3Chamber

The chamber in `S` induced from a Weyl vector in `L`.
"""
mutable struct K3Chamber
  weyl_vector::fmpz_mat
  # for v in walls, the corresponding half space is defined by the equation
  # x * gram_matrix(S)*v >= 0, further v is primitive in S (and, in contrast to Shimada, not S^\vee)
  walls::Vector{fmpz_mat}
  lengths::Vector{fmpq}  #
  B::fmpz_mat # QQ-basis consisting of rays #... why do we bother to save this?
  gramB::fmpz_mat # the basis matrix inferred from the QQ-basis
  parent_wall::fmpz_mat # for the spanning tree
  data::BorcherdsCtx
  fp::Matrix{Int} # fingerprint for the backtrack search
  #per::Vector{Int}  # permutation
  fp_diagonal::Vector{Int}

  # TODO: Be more memory efficient and store only the indices for the
  # basis matrix and the permutation.
  # I am not sure if storing gram matrix stuff in memory actually increases performance...

  function K3Chamber()
    return new()
  end
end

export chamber, BorcherdsCtx

@doc Markdown.doc"""
    chamber(data::BorcherdsCtx, weyl_vector::fmpz_mat, [parent_wall::fmpz_mat, walls::Vector{fmpz_mat}])

Return the chamber with the given Weyl vector.
"""
function chamber(data::BorcherdsCtx, weyl_vector::fmpz_mat, parent_wall::fmpz_mat=zero_matrix(ZZ, 0, 0))
  D = K3Chamber()
  D.weyl_vector = weyl_vector
  D.parent_wall = parent_wall
  D.data = data
  return D
end

function chamber(data::BorcherdsCtx, weyl_vector::fmpz_mat, parent_wall::fmpz_mat, walls::Vector{fmpz_mat})
  D = K3Chamber()
  D.weyl_vector = weyl_vector
  D.parent_wall = parent_wall
  D.data = data
  D.walls = walls
  return D
end

# needed to create sets of K3Chambers
function Base.hash(C::K3Chamber)
  return hash(C.weyl_vector[:,1:rank(C.data.S)])
end

# Two chambers are equal if and only if their Weyl vectors
# project to the same point in S
# By the choice of our coordinates this projection is determined
# by the first rank(S) coordinates.
function Base.:(==)(C::K3Chamber, D::K3Chamber)
  @req C.data===D.data "K3Chambers do not have the same context"
  return C.weyl_vector[:,1:rank(C.data.S)] == D.weyl_vector[:,1:rank(D.data.S)]
end

@doc Markdown.doc"""
    walls(D::K3Chamber) -> Vector{fmpz_mat}

Return the walls of the chamber `D`, i.e. its facets.

The corresponding half space of the wall defined by `v` in `walls(D)` is

```Math
\{x \in S \otimes \mathbb{R} |  \langle x,v \rangle  \geq 0\}.
```

`v` is given with respect to the basis of `S` and is primitive in `S`.
Note that Shimada follows a different convention and takes `v` primitive in `S^\vee`.
"""
function walls(D::K3Chamber)
  if !isdefined(D, :walls)
    D.walls = _walls_of_chamber(D.data, D.weyl_vector)
    @assert length(D.walls)>=rank(D.data.S) "$(D.weyl_vector)"
  end
  return D.walls
end

weyl_vector(D::K3Chamber) = D.weyl_vector

@doc Markdown.doc"""
    rays(D::K3Chamber)

Return the rays of the induced chamber `D`.

They are represented as primitive row vectors with respect to the basis of `S`.
"""
function rays(D::K3Chamber)
  r = reduce(vcat, walls(D), init=zero_matrix(ZZ,0,rank(D.data.SS)))
  rQ = change_base_ring(QQ, r) * gram_matrix(D.data.SS)
  C = positive_hull(rQ)
  Cd = polarize(C)
  L = rays(Cd)
  Lq = fmpq_mat[matrix(QQ,1,rank(D.data.SS),i) for i in L]
  # clear denominators
  Lz = fmpz_mat[change_base_ring(ZZ,i*denominator(i)) for i in Lq]
  # primitive in S
  Lz = fmpz_mat[divexact(i,gcd(vec(i))) for i in Lz]
  @hassert :K3Auto 2 all(all(x>=0 for x in vec(r*gram_matrix(D.data.SS)*transpose(i))) for i in Lz)
  return Lz
end

function Base.show(io::IO, c::K3Chamber)
  if isdefined(c,:walls)
    print(IOContext(io, :compact => true), "Chamber  in dimension $(length(walls(c)[1])) with $(length(walls(c))) walls")
  else
    print(IOContext(io, :compact => true), "Chamber: $(c.weyl_vector[1,1:rank(c.data.S)])")
  end
end

@doc Markdown.doc"""
    fingerprint(D::K3Chamber)

Return the fingerprint of this chamber.

The fingerprint is an invariant computed from the rays and their inner products.
"""
function fingerprint(D::K3Chamber)
  v = sum(walls(D))
  G = D.data.gramS
  m1 = (v*G*transpose(v))[1,1]
  m2 = [(a*G*transpose(a))[1,1] for a in walls(D)]
  sort!(m2)
  m3 = [(v*G*transpose(a))[1,1] for a in walls(D)]
  sort!(m3)
  m4 = fmpz[]
  for i in 1:length(walls(D))
    for j in 1:i-1
      push!(m4,(walls(D)[i]*G*transpose(walls(D)[j]))[1,1])
    end
  end
  sort!(m4)
  V = Dict{Tuple{fmpz,fmpz},Vector{fmpz_mat}}()
  for w in walls(D)
    i =  (v*G*transpose(w))[1,1]
    j =  (w*G*transpose(w))[1,1]
    if (i,j) in keys(V)
      push!(V[(i,j)],w)
    else
      V[(i,j)] = [w]
    end
  end
  #=
  # so far m5 was not needed to separate the O(S)-orbits
  m5 = []
  for i in keys(V)
    vi = sum(V[i])
    push!(m5, [i,sort!([(vi*G*transpose(j))[1,1] for j in walls(D)])])
  end
  sort!(m5)
  # So far we have only O(S)-invariants. There are also ways to produce G-invariants
  # by working the the images of the rays in the discriminant group and their
  # G-orbits. Perhaps one has to switch to S^\vee primitive vectors in this case.
  =#
  return (m1, m2, m3, m4)
end


"""
Compute the fingerprint defined by Plesken-Souvignier and change the basis
matrix and gram matrix accordingly.

It is computed from the walls and the gram matrix of `S`.

A permutation `per` for the
order of the basis-vectors is chosen
such that in every step the number of
possible continuations is minimal,
for j from per[i] to per[dim-1] the
value f[i][j] in the fingerprint f is
the number of vectors, which have the
same scalar product with the
basis-vectors per[0]...per[i-1] as the
basis-vector j and the same length as
this vector with respect to all
invariant forms
"""
function _fingerprint_backtrack!(D::K3Chamber)
  n = rank(D.data.S)
  V = walls(D)
  gramS = gram_matrix(D.data.S)
  B, indB = _find_basis(V, n)
  tmp = V[indB]
  deleteat!(V, indB)
  prepend!(V, tmp)
  lengths = fmpq[(v*gramS*transpose(v))[1,1] for v in V]
  D.lengths = lengths
  gramB = change_base_ring(ZZ, B*gramS*transpose(B))
  D.gramB = gramB
  D.B = B

  per = Vector{Int}(undef, n)
  for i in 1:n
    per[i] = i
  end

  fp = zeros(Int, n, n)

  # fp[1, i] = # vectors v such that v has same length as b_i for all forms
  for i in 1:n
    cvl = gramB[i,i]
    fp[1, i] = count(x->x==cvl, lengths)

  end

  for i in 1:(n - 1)
    # Find the minimal non-zero entry in the i-th row
    mini = i
    @inbounds for j in (i+1):n
      if fp[i, per[j]] < fp[i, per[mini]]
        mini = j
      end
    end

    per[mini], per[i] = per[i], per[mini]

    # Set entries below the minimal entry to zero
    @inbounds for j in (i + 1):n
      fp[j, per[i]] = 0
    end

    # Now compute row i + 1

    for j in (i + 1):n
      fp[i + 1, per[j]] = _possible(D, per, i, per[j])
    end
  end

  # Extract the diagonal

  res = Vector{Int}(undef, n)

  @inbounds for i in 1:n
    res[i] = fp[i, per[i]]
    @assert res[i]>0
  end




  #D.per = per
  D.B = B[per,:]
  D.gramB = gramB[per,per]
  D.fp = fp[:,per]
  D.fp_diagonal = res
end

@doc Markdown.doc"""
    _possible(D::K3Chamber, per, I, J) -> Int

Return the number of possible extensions of an `n`-partial isometry to
an `n+1`-partial one.
"""
function _possible(D::K3Chamber, per, I, J)
  vectors = walls(D)
  lengths = D.lengths
  gramB = D.gramB
  gramS = D.data.gramS
  n = length(vectors)
  @assert n == length(lengths)

  count = 0
  T = gramS*transpose(reduce(vcat,vectors[per[1:I]]))
  for j in 1:n
    lengthsj = lengths[j]
    vectorsj = vectors[j]
    good_scalar = true
    if lengthsj != gramB[J, J]
      continue
    end

    for i in 1:I
      if (vectorsj*T[:,i])[1,1] != gramB[J,per[i]]
        good_scalar = false
        break
      end
    end

    if !good_scalar
      continue
    end
    count = count + 1

    # length is correct
  end
  return count
end

################################################################################
# close vector functions
################################################################################

@doc Markdown.doc"""
    enumerate_quadratic_triple -> Vector{Tuple{Vector{Int}, fmpq}}

Return $\{x \in \mathbb Z^n : x Q x^T + 2xb^T + c <=0\}$.

#Input:
- `Q`: positive definite matrix
- `b`: row vector
- `c`: rational number
"""
function enumerate_quadratic_triple(Q, b, c; algorithm=:short_vectors, equal=false)
  if algorithm == :short_vectors
    L, p, dist = Hecke._convert_type(Q, b, QQ(c))
    #@vprint :K3Auto 1 ambient_space(L), basis_matrix(L), p, dist
    if equal
      cv = Hecke.close_vectors(L, vec(p), dist, dist, check=false)
    else
      cv = Hecke.close_vectors(L, vec(p), dist, check=false)
    end
  end
  return cv
end

@doc Markdown.doc"""
    short_vectors_affine(S::ZLat, v::MatrixElem, alpha, d)
    short_vectors_affine(gram::MatrixElem, v::MatrixElem, alpha, d)

Return the vectors of squared length `d` in the given affine hyperplane.

```Math
\{x \in S : x^2=d, x.v=\alpha \}.
```
The matrix version takes `S` with standard basis and the given gram matrix.

# Arguments
- `v`: row vector with $v^2 > 0$
- `S`: a hyperbolic `Z`-lattice

The output is given in the ambient representation.

The implementation is based on Algorithm 2.2 in [Shimada](@cite)
"""
function short_vectors_affine(S::ZLat, v::MatrixElem, alpha, d)
  alpha = QQ(alpha)
  gram = gram_matrix(S)
  tmp = v*gram_matrix(ambient_space(S))*transpose(basis_matrix(S))
  v_S = solve_left(gram_matrix(S),tmp)
  sol = short_vectors_affine(gram, v_S, alpha, d)
  B = basis_matrix(S)
  return [s*B for s in sol]
end

function short_vectors_affine(gram::MatrixElem, v::MatrixElem, alpha::fmpq, d)
  # find a solution <x,v> = alpha with x in L if it exists
  w = gram*transpose(v)
  tmp = FakeFmpqMat(w)
  wn = numerator(tmp)
  wd = denominator(tmp)
  b, x = can_solve_with_solution(transpose(wn), matrix(ZZ, 1, 1, [alpha*wd]))
  if !b
    return fmpq_mat[]
  end
  _, K = left_kernel(wn)
  # (x + y*K)*gram*(x + y*K) = x gram x + 2xGKy + y K G K y

  # now I want to formulate this as a cvp
  # (x +y K) gram (x+yK) ==d
  # (x
  GK = gram*transpose(K)
  Q = K * GK
  b = transpose(x) * GK
  c = (transpose(x)*gram*x)[1,1] - d
  # solve the quadratic triple
  Q = change_base_ring(QQ, Q)
  b = change_base_ring(QQ, transpose(b))
  cv = enumerate_quadratic_triple(-Q, -b,-QQ(c),equal=true)
  xt = transpose(x)
  cv = [xt+matrix(ZZ,1,nrows(Q),u[1])*K for u in cv]
  @hassert :K3Auto 1 all((v*gram*transpose(u))[1,1]==alpha for u in cv)
  @hassert :K3Auto 1 all((u*gram*transpose(u))[1,1]== d for u in cv)
  return cv #[u for u in cv if (u*gram*transpose(u))[1,1]==d]
end


@doc Markdown.doc"""
    separating_hyperplanes(S::ZLat, v::fmpq_mat, h::fmpq_mat, d)

Return $\{x in S : x^2=d, x.v>0, x.h<0\}$.

# Arguments
- `S`:  a hyperbolic lattice
- `d`: a negative integer
- `v`,`h`: vectors of positive square
"""
function separating_hyperplanes(S::ZLat, v::fmpq_mat, h::fmpq_mat, d)
  V = ambient_space(S)
  @hassert :K3Auto 1 inner_product(V,v,v)[1,1]>0
  @hassert :K3Auto 1 inner_product(V,h,h)[1,1]>0
  gram = gram_matrix(S)
  B = basis_matrix(S)
  vS = solve_left(B,v)
  hS = solve_left(B,h)
  return [a*B for a in separating_hyperplanes(gram,vS,hS,d)]
end

function separating_hyperplanes(gram::fmpq_mat, v::fmpq_mat, h::fmpq_mat, d)
  L = Zlattice(gram=gram)
  n = ncols(gram)
  ch = QQ((h*gram*transpose(h))[1,1])
  cv = QQ((h*gram*transpose(v))[1,1])
  b = basis_matrix(L)
  prW = reduce(vcat,[b[i,:] - (b[i,:]*gram*transpose(h))*ch^-1*h for i in 1:n])
  W = lattice(ambient_space(L), prW, isbasis=false)
  bW = basis_matrix(W)
  # set up the quadratic triple for SW
  gramW = gram_matrix(W)
  s = solve_left(bW, v*prW) * gramW
  Q = gramW + transpose(s)*s*ch*cv^-2

  @vprint :K3Auto 5 Q
  LQ = Zlattice(gram=-Q*denominator(Q))

  S = fmpq_mat[]
  h = change_base_ring(QQ, h)
  rho = abs(d)*ch^-1
  t,sqrtho = issquare_with_sqrt(rho)
  if t
    r = sqrtho*h
    if denominator(r)==1 && (r*gram*transpose(h))[1,1]>0 && (r*gram*transpose(v))[1,1] < 0
      push!(S,r)
    end
  end
  for (x,_) in short_vectors_iterator(LQ,  abs(d*denominator(Q)))
    rp = matrix(ZZ, 1, nrows(Q), x)*bW
    rho = abs(d - (rp*gram*transpose(rp))[1,1])*ch^-1
    t,rho = issquare_with_sqrt(rho)
    if !t
      continue
    end
    r = rho*h + rp
    if denominator(r)==1 && (r*gram*transpose(h))[1,1]>0 && (r*gram*transpose(v))[1,1] < 0
      push!(S,r)
    end
    r = rho*h - rp
    if denominator(r)==1 && (r*gram*transpose(h))[1,1]>0 && (r*gram*transpose(v))[1,1] < 0
      push!(S,r)
    end
  end
  return S
end

@doc Markdown.doc"""
    _find_basis(row_matrices::Vector, dim::Integer)

Return the first `dim` linearly independent vectors in row_matrices and their indices.

We assume that row_matrices consists of row vectors.
"""
function _find_basis(row_matrices::Vector, dim::Integer)
  @req length(row_matrices)>=dim > 0 "must contain at least a single vector"
  r = row_matrices[1]
  n = ncols(r)
  B = zero_matrix(base_ring(r), 0, n)
  rk = 0
  indices = Int[]
  for i in 1:length(row_matrices)
    r = row_matrices[i]
    Br = vcat(B, r)
    rk = rank(Br)
    if rk > nrows(B)
      B = Br
      push!(indices, i)
    end
    if rk == dim
      break
    end
  end
  @assert length(indices) == rk == dim
  return B, indices
end

_find_basis(row_matrices::Vector) = _find_basis(row_matrices, ncols(row_matrices[1]))

@doc Markdown.doc"""
    is_pm1_on_discr(S::ZLat, g::fmpz_mat) -> Bool

Return whether the isometry `g` of `S` acts as `+-1` on the discriminant group.
"""
function is_pm1_on_discr(S::ZLat, g::fmpz_mat)
  D = discriminant_group(S)
  imgs = [D(vec(matrix(QQ,1,rank(S),lift(d))*g)) for d in gens(D)]
  return all(imgs[i] == gens(D)[i] for i in 1:length(gens(D))) || all(imgs[i] == -gens(D)[i] for i in 1:length(gens(D)))
  # OD = orthogonal_group(D)
  # g1 = hom(D,D,[D(lift(d)*g) for d in gens(D)])
  # gg = OD(g1)
  # return isone(gg) || gg == OD(-matrix(one(OD)))
end

@doc Markdown.doc"""
    hom(D::K3Chamber, E::K3Chamber) -> Vector{fmpz_mat}

Return the set ``Hom_G(D, E)`` of elements of `G` mapping `D` to `E`.

The elements are represented with respect to the basis of `S`.
"""
Hecke.hom(D::K3Chamber, E::K3Chamber) = alg319(D, E)
#alg319(gram_matrix(D.data.SS), D.B,D.gramB, walls(D), walls(E), D.data.membership_test)

@doc Markdown.doc"""
    aut(D::K3Chamber, E::K3Chamber) -> Vector{fmpz_mat}

Return the stabilizer of ``E`` in ``G``.

The elements are represented with respect to the basis of `S`.
"""
aut(D::K3Chamber) = hom(D, D)

function alg319(D::K3Chamber, E::K3Chamber)
  if !isdefined(D,:B)
    _fingerprint_backtrack!(D)  # compute a favorable basis
  end
  gram_basis = D.gramB
  gram = D.data.gramS
  fp = D.fp_diagonal
  basis = D.B
  n = ncols(gram)
  raysD = walls(D)
  raysE = walls(E)
  partial_homs = [zero_matrix(ZZ, 0, n)]
  # breadth first search
  # Since we expect D and E to be isomorphic,
  # a depth first search with early abort could be more efficient.
  # for now this does not seem to be a bottleneck
  for i in 1:n
    @vprint :K3Auto 4 "level $(i-1), partial homs $(length(partial_homs)) \n"
    partial_homs_new = fmpz_mat[]
    for img in partial_homs
      extensions = fmpz_mat[]
      k = nrows(img)
      gi = gram*transpose(img)
      for r in raysE
        if (r*gram*transpose(r))[1,1] != gram_basis[k+1,k+1] || (k>0 && r*gi != gram_basis[k+1,1:k])
          continue
        end
        # now r has the correct inner products with what we need
        push!(extensions, vcat(img, r))
      end
      if fp[i] != length(extensions)
        continue
      end
      append!(partial_homs_new, extensions)
    end
    partial_homs = partial_homs_new
  end
  basisinv = inv(change_base_ring(QQ, basis))
  homs = fmpz_mat[]
  is_in_hom_D_E(fz) = all(r*fz in raysE for r in raysD)
  vE = sum(raysE) # center of mass of the dual cone
  vD = sum(raysD)
  for f in partial_homs
    f = basisinv*f
    if denominator(f)!=1
      continue
    end
    fz = change_base_ring(ZZ, f)
    if !D.data.membership_test(fz)
      continue
    end
    if !(vD*fz == vE)
      continue
    end
    # The center of mass is an interior point
    # Further it uniquely determines the chamber and is compatible with homomorphisms
    # This is basically Remark 3.20
    # -> but this is not true for the center of mass of the dual cone
    # hence this extra check
    if !is_in_hom_D_E(fz)
      continue
    end
    push!(homs, fz)
  end
  @hassert :K3Auto 1 all(f*gram*transpose(f)==gram for f in homs)
  return homs
end


# legacy worker for hom and aut without Plesken-Souvignier preprocessing
function alg319(gram::MatrixElem, raysD::Vector{fmpz_mat}, raysE::Vector{fmpz_mat}, membership_test)
  n = ncols(gram)
  partial_homs = [zero_matrix(ZZ, 0, n)]
  basis,_ = _find_basis(raysD, n)
  gram_basis = basis*gram*transpose(basis)
  return alg319(gram, basis, gram_basis, raysD, raysE, membership_test)
end

function alg319(gram::MatrixElem, basis::fmpz_mat, gram_basis::fmpq_mat, raysD::Vector{fmpz_mat}, raysE::Vector{fmpz_mat}, membership_test)
  n = ncols(gram)
  partial_homs = [zero_matrix(ZZ, 0, n)]
  # breadth first search
  # Since we expect D and E to be isomorphic,
  # a depth first search with early abort would be more efficient.
  # for now this does not seem to be a bottleneck
  for i in 1:n
    @vprint :K3Auto 4 "level $(i), partial homs $(length(partial_homs)) \n"
    partial_homs_new = fmpz_mat[]
    for img in partial_homs
      extensions = fmpz_mat[]
      k = nrows(img)
      gi = gram*transpose(img)
      for r in raysE
        if (r*gram*transpose(r))[1,1] != gram_basis[k+1,k+1] || (k>0 && r*gi != gram_basis[k+1,1:k])
          continue
        end
        # now r has the correct inner products with what we need
        push!(extensions, vcat(img,r))
      end
      append!(partial_homs_new, extensions)
    end
    partial_homs = partial_homs_new
  end
  basisinv = inv(change_base_ring(QQ, basis))
  homs = fmpz_mat[]
  is_in_hom_D_E(fz) = all(r*fz in raysE for r in raysD)
  vE = sum(raysE) # center of mass of the dual cone
  vD = sum(raysD)
  for f in partial_homs
    f = basisinv*f
    if denominator(f)!=1
      continue
    end
    fz = change_base_ring(ZZ,f)
    if !membership_test(fz)
      continue
    end
    if !(vD*fz == vE)
      continue
    end
    # The center of mass is an interior point
    # Further it uniquely determines the chamber and is compatible with homomorphisms
    # This is basically Remark 3.20
    # -> but this is not true for the center of mass of the dual cone
    if !is_in_hom_D_E(fz)
      continue
    end
    push!(homs, fz)
  end
  @hassert :K3Auto 1 all(f*gram*transpose(f)==gram for f in homs)
  return homs
end


@doc Markdown.doc"""
    _alg58(L::ZLat, S::ZLat, R::ZLat, prRdelta, w)

Compute Delta_w

Tuples (r_S, r) where r is an element of Delta_w and r_S is the
orthogonal projection of `r` to `S`.

Correponds to Algorithm 5.8 in [Shimada](@cite)
but this implementation is different.
"""
# legacy function needed for precomputations
function _alg58(L::ZLat, S::ZLat, R::ZLat, prRdelta, w)
  V = ambient_space(L)
  d = exponent(discriminant_group(S))
  @hassert :K3Auto 1 V == ambient_space(S)
  n_R = [QQ(i)//d for i in (-2*d+1):0 if mod(d*i,2)==0]
  Rdual = dual(R)
  Sdual = dual(S)
  rkR = rank(R)
  delta_w = fmpq_mat[]
  iB = inv(basis_matrix(L))
  for c in n_R
    cm = -c
    for (vr0,vsquare) in prRdelta
      if vsquare != cm
        continue
      end
      a0 = inner_product(V,w,vr0)[1,1]
      if c == 0
        VA = [(vr0,a0)]
      else
        VA = [(vr0,a0),(-vr0,-a0)]
      end
      for (vr,a) in VA
        Sdual_na = short_vectors_affine(Sdual, w, 1 - a, -2 - c)
        for vs in Sdual_na
          vv = vs +  vr
          if all(denominator(i)==1 for i in collect(vv*iB))
            push!(delta_w, vs)
          end
        end
      end
    end
  end
  return delta_w
end

function _alg58(L::ZLat, S::ZLat, R::ZLat, w::MatrixElem)
  Rdual = dual(R)
  sv = short_vectors(rescale(Rdual, -1), 2)
  # not storing the following for efficiency
  # append!(sv,[(-v[1],v[2]) for v in sv])
  # but for convenience we include zero
  T = typeof(sv).parameters[1].parameters[1].parameters[1]
  push!(sv,(zeros(T, rank(Rdual)), QQ(0)))
  rkR = rank(R)
  prRdelta = [(matrix(QQ, 1, rkR, v[1])*basis_matrix(Rdual),v[2]) for v in sv]
  return _alg58(L, S, R, prRdelta, w)
end

# the actual somewhat optimized implementation relying on short vector enumeration
function _alg58_short_vector(data::BorcherdsCtx, w::fmpz_mat)
  L = data.L
  V = ambient_space(L)
  S = data.S
  R = data.R
  wS = w*data.prS
  wSL = wS*basis_matrix(S)
  wL = gram_matrix(L)*transpose(w)
  wSsquare = (wS*data.gramS*transpose(wS))[1,1]
  W = lattice(V, wS*basis_matrix(S))
  N = orthogonal_submodule(S, W)
  # W + N + R < L of finite index
  svp_input = Tuple{fmpq,fmpq_mat,fmpq,Int}[]
  for (rR, rRsq) in data.prRdelta
    if rRsq==2
      continue
    end
    @inbounds rwS = (rR*wL)[1,1]
    alpha = 1 - rwS
    usq = alpha^2*wSsquare^-1 - rRsq
    sq = -2 - usq
    push!(svp_input, (alpha, rR,sq,1))
    alpha = 1 + rwS
    usq = alpha^2*wSsquare^-1 - rRsq
    sq = -2 - usq
    push!(svp_input, (alpha, rR, sq, -1))
  end
  @inbounds bounds = unique!([-i[3] for i in svp_input])
  Ndual = dual(N)
  G = -gram_matrix(Ndual)
  d = denominator(G)
  bounds = [i for i in bounds if divides(d,denominator(i))[1]]
  mi = minimum(bounds)
  ma = maximum(bounds)

  svN = Hecke._short_vectors_gram(Hecke.LatEnumCtx, G,mi,ma, fmpz)
  result = fmpq_mat[]
  # treat the special case of the zero vector by copy paste.
  if QQ(0) in bounds
    (rN,sqrN) = (zeros(Int64,rank(Ndual)),0)
    rN1 = zero_matrix(ZZ,1,degree(Ndual))
    found1 = false
    found2 = false
    sqrN = QQ(0)
    for (alpha, rR, sq, si) in svp_input
      if sqrN != sq
        continue
      end
      rr = alpha*wSsquare^-1*wSL + si*rR
      r = rr + rN1
      if !found1 && @inbounds all(denominator(r[1,i])==1 for i in 1:ncols(r))==1
        found1 = true
        push!(result, r*data.prS)
        break
      end
      r = rr - rN1
      if !found2 && @inbounds all(denominator(r[1,i])==1 for i in 1:ncols(r))==1
        found2 = true
        push!(result, r*data.prS)
        break
      end
      if found1 && found2
        break
      end
    end
  end
  for (rN, sqrN) in svN
    if !(sqrN in bounds)
      continue
    end
    rN1 = matrix(ZZ,1,rank(Ndual),rN)*basis_matrix(Ndual)
    found1 = false
    found2 = false
    sqrN = -sqrN
    for (alpha, rR, sq, si) in svp_input
      if sqrN != sq
        continue
      end
      rr = alpha*wSsquare^-1*wSL + si*rR
      r = rr + rN1
      if !found1 && @inbounds all(denominator(r[1,i])==1 for i in 1:ncols(r))==1
        found1 = true
        push!(result, r*data.prS)
        break
      end
      r = rr - rN1
      if !found2 && @inbounds all(denominator(r[1,i])==1 for i in 1:ncols(r))==1
        found2 = true
        push!(result, r*data.prS)
        break
      end
      if found1 && found2
        break
      end
    end
  end
  return result
end

@doc Markdown.doc"""
    _alg58_close_vector(data::BorcherdsCtx, w::fmpz_mat)

Return tuples (r_S, r) where r is an element of Delta_w and r_S is the
orthogonal projection of `r` to `S^\vee` given in the basis of S.

Corresponds to Algorithm 5.8 in [Shimada](@cite)
"""
function _alg58_close_vector(data::BorcherdsCtx, w::fmpz_mat)
  V = ambient_space(data.L)
  S = data.S
  d = exponent(discriminant_group(S))
  @hassert :K3Auto 1 V == ambient_space(S)
  @hassert :K3Auto 2 basis_matrix(data.L)==1
  n_R = [QQ(i)//d for i in (-2*d+1):0 if mod(d*i,2)==0]
  SSdual = dual(data.SS)
  delta_w = fmpq_mat[]
  wS = w*data.prS
  #wS = solve_left(gram_matrix(S),w*gram_matrix(V)*transpose(basis_matrix(S)))
  Vw = data.gramL*transpose(w)
  # since we do repeated cvp in the same lattice
  # we do the preprocessing here

  # collect the cvp inputs to avoid repeated calculation of the same cvp
  cvp_inputs = Dict{Tuple{fmpq,fmpq},Vector{fmpq_mat}}()
  for c in n_R
    kc = -2-c
    cm = -c
    for (vr,vsquare) in data.prRdelta
      if vsquare != cm
        continue
      end
      a = (vr*Vw)[1,1]  # TODO: could be improved by working in R
      key = (1-a,kc)
      if key in keys(cvp_inputs)
        push!(cvp_inputs[key], vr)
      else
        cvp_inputs[key] = [vr]
      end
      if c != 0
        key1 = (1+a, kc)
        if key1 in keys(cvp_inputs)
          push!(cvp_inputs[key1], -vr)
        else
          cvp_inputs[key1] = [-vr]
        end
      end
    end
  end

  # setup
  gram = gram_matrix(SSdual)
  #B = basis_matrix(S)
  #return [s*B for s in sol]

  # much of this code is copied from
  # short_vectors_affine(gram::MatrixElem, v::MatrixElem, alpha::fmpq, d)
  # to avoid repeated calculation of the same stuff e.g. K and Q
  # find a solution <x,v> = alpha with x in L if it exists
  ww = transpose(wS)
  tmp = FakeFmpqMat(ww)
  wn = numerator(tmp)
  wd = denominator(tmp)
  _, K = left_kernel(wn)
  K = lll!(K)  # perhaps doing this has no gain?
  # (x + y*K)*gram*(x + y*K) = x gram x + 2xGKy + y K G K y

  # now I want to formulate this as a cvp
  # (x +y K) gram (x+yK) ==d
  # (x
  KG = K*gram
  Q = -KG * transpose(K)
  Qi = inv(Q)
  N = Zlattice(gram=Q,check=false)
  V = ambient_space(N)

  #@show sum(length.(values(cvp_inputs)))
  tmp = zero_matrix(QQ,1,rank(SSdual))


  B = basis_matrix(SSdual)
  KB = K*B
  for (alpha, d) in keys(cvp_inputs)
    can_solve_i, x = can_solve_with_solution(transpose(wn), matrix(ZZ, 1, 1, [alpha*wd]))
    if !can_solve_i
      continue
    end
    # looks like premature optimization ...
    x = change_base_ring(QQ, x)
    b = KG*x
    transpose!(tmp, x)
    # c = (transpose(x)*gram*x)[1,1] - d
    c = (tmp*mul!(x,gram,x))[1,1] - d
    #cv = enumerate_quadratic_triple(Q,-b,-c,equal=true)
    mul!(b, Qi, b)
    #b = Qi*b
    v = vec(b)
    upperbound = inner_product(V,v,v) + c
    # solve the quadratic triple
    cv = close_vectors(N, v, upperbound, upperbound, check=false)
    mul!(tmp,tmp,B)
    #xtB = transpose(x)*B
    Sdual_na1 = [matrix(ZZ, 1, nrows(Q), u)*KB for (u,_) in cv]
    for v in Sdual_na1
      add!(v,v,tmp)
    end
    Sdual_na2 = [vs*basis_matrix(data.S) for vs in Sdual_na1]
    for i in 1:length(Sdual_na1)
      v = Sdual_na2[i]
      for vr in cvp_inputs[(alpha,d)]
        vv =  v +  vr
        if denominator(vv)==1
          push!(delta_w, Sdual_na1[i])
          break # delta_w is a set, hence we may break
        end
      end
    end
  end

  # the chamber should intersect the boundary only at the QQ-rational points
  @hassert :K3Auto 2 rank(S) == rank(reduce(vcat,[s for s in delta_w]))
  return delta_w
end


@doc Markdown.doc"""
    _walls_of_chamber(data::BorcherdsCtx, weyl_vector)

Return the walls of the L|S chamber induced by `weyl_vector`.

Corresponds Algorithm 5.11 in [Shimada](@cite) and calls Polymake.
"""
function _walls_of_chamber(data::BorcherdsCtx, weyl_vector, alg=:short)
  if alg==:short
    walls1 = _alg58_short_vector(data, weyl_vector)
  elseif alg==:close
    walls1 = _alg58_close_vector(data, weyl_vector)
  end
  if length(walls1)==rank(data.S)
    # shortcut which avoids calling Polymake
    d = rank(data.S)
    walls = Vector{fmpz_mat}(undef,d)
    for i in 1:d
      vs = numerator(FakeFmpqMat(walls1[i]))
      g = gcd(vec(vs))
      if g != 1
        vs = divexact(vs, g)
      end
      walls[i] = vs
    end
    return walls
  end
  i = zero_matrix(QQ, 0, degree(data.SS))
  D = reduce(vcat, (v for v in walls1), init=i)
  P = positive_hull(D)
  r = rays(P)
  d = length(r)
  walls = Vector{fmpz_mat}(undef,d)
  for i in 1:d
    v = matrix(QQ, 1, degree(data.SS), r[i])
    # rescale v to be primitive in S
    vs = numerator(FakeFmpqMat(v))
    g = gcd(vec(vs))
    if g!=1
      vs = divexact(vs, g)
    end
    walls[i] = vs
  end
  return walls
end

@doc Markdown.doc"""
    is_S_nondegenerate(L::ZLat, S::ZLat, w::fmpq_mat)

Return whether the ``L|S`` chamber defined by `w` is `S`-nondegenerate.

This is the case if and only if $C = C(w) \cap S \otimes \RR$ has the
expected dimension `dim(S)`.
"""
function is_S_nondegenerate(L::ZLat, S::ZLat, w::fmpq_mat)
  R = Hecke.orthogonal_submodule(L, S)
  Delta_w = _alg58(L, S, R, w)
  V = ambient_space(L)
  G = gram_matrix(V)
  BS = transpose(basis_matrix(S))
  prSDelta_w = [v * G * BS for v in Delta_w]
  i = zero_matrix(QQ, 0, rank(S))
  D = reduce(vcat, prSDelta_w, init=i)
  P = positive_hull(D)  # the dual cone of C
  # If P has a linear subspace, then its dual C is not of full dimension.
  return ispointed(P)
end

@doc Markdown.doc"""
    inner_point(L::ZLat, S::ZLat, w::fmpq_mat)
    inner_point(C::K3Chamber)

Return a reasonably small integer inner point of the given L|S chamber.
"""
function inner_point(L::ZLat, S::ZLat, w::fmpq_mat)
  R = Hecke.orthogonal_submodule(L, S)
  Delta_w = _alg58(L, S, R, w)
  V = ambient_space(L)
  G = gram_matrix(V)
  prSDelta_w = [v * G * transpose(basis_matrix(S)) for v in Delta_w]
  i = zero_matrix(QQ, 0, rank(S))
  D = reduce(vcat, prSDelta_w, init=i)
  P = positive_hull(D)  # dual to C
  @hassert :K3Auto 1 is_pointed(P) # we check S-nondegenerateness
  C = polarize(P)  # C
  facets = [matrix(QQ,rank(S),1, r) for r in rays(P)]
  H = [matrix(ZZ,1,rank(S),i) for i in hilbert_basis(C)]

  # add hilbert basis vectors until we are in the interior
  n = length(facets)
  h = H[1]
  j = 2
  for i in 1:n
    if (h*facets[i])[1,1] > 0
      continue
    end
    for k in 1:length(H)
      v = H[k]
      if (v*facets[i])[1,1]> 0
        h = h + v
        break
      end
    end
  end
  # confirm the computation
  hV = h*basis_matrix(S)
  @hassert :K3Auto 1 all(0 < inner_product(V, v, hV)[1,1] for v in Delta_w)
  return h
end

inner_point(C::K3Chamber) = inner_point(C.data.L, C.data.S, change_base_ring(QQ,C.weyl_vector))


@doc Markdown.doc"""
    unproject_wall(data::BorcherdsCtx, vS::fmpz_mat)

Return the (-2)-walls of L containing $v^{perp_S}$ but not all of ``S``.

Based on Algorithm 5.13 in [Shimada](@cite)

# Arguments
- `vS`: Given with respect to the basis of `S`.
"""
function unproject_wall(data::BorcherdsCtx, vS::fmpz_mat)
  d = gcd(vec(vS*data.gramS))
  v = QQ(1,d)*(vS*basis_matrix(data.S))  # primitive in Sdual
  vsq = QQ((vS*data.gramS*transpose(vS))[1,1],d^2)

  @hassert :K3Auto 1 vsq>=-2
  rkR = rank(data.R)
  Pv = fmpz_mat[]
  for alpha in 1:Int64(floor(sqrt(Float64(-2//vsq))))
    c = 2 + alpha^2*vsq
    alphav = alpha*v
    for (vr,cc) in data.prRdelta
      # probably we could speed up the for loop and compute
      # and compute the result without a membership test
      # by working modulo ZZ^n directly
      if cc != c
        continue
      end
      vv = alphav + vr
      if denominator(vv)==1
        push!(Pv,change_base_ring(ZZ,vv))
      end
      if c != 0
        vv = alphav - vr
        if denominator(vv)==1
          push!(Pv, change_base_ring(ZZ,vv))
        end
      end
    end
  end
  return Pv
end

@doc Markdown.doc"""
    adjacent_chamber(D::K3Chamber, v) -> K3Chamber

Return return the L|S chamber adjacent to `D` via the wall defined by `v`.
"""
function adjacent_chamber(D::K3Chamber, v)
  gramL = D.data.gramL
  dualDeltaR = D.data.dualDeltaR
  deltaR = D.data.deltaR
  dimL = ncols(gramL)
  Pv = unproject_wall(D.data, v)
  l = length(Pv)
  @hassert :K3Auto 1 length(Pv) == length(unique(Pv))
  a = 1000000
  @label getu
  a = 2*a
  rep = Array{Tuple{Int,fmpq,Bool}}(undef,l+length(dualDeltaR))
  u = matrix(ZZ, 1, dimL, rand(-a:a, dimL))
  Vw = gramL*transpose(D.weyl_vector)
  Vu = gramL*transpose(u)

  z = zero_matrix(ZZ,1,1)
  for i in 1:length(Pv)
    r = Pv[i]
    mul!(z,r, Vw)
    s = (r*Vu)[1,1]
    divexact!(s,z[1,1])
    if any(rep[j][2]==s for j in 1:i-1)
      @goto getu
    end
    rep[i] = (i, s, false)
  end

  for i in 1:length(deltaR)
    r = deltaR[i]
    mul!(z, r, Vw)
    s = (r*Vu)[1,1]
    divexact!(s,z[1,1])
    if any(rep[j][2]==s for j in 1:l+i-1)
      @goto getu
    end
    rep[l+i] = (i, s, true)
  end
  @hassert :K3Auto 2 length(unique([r[2] for r in rep]))==length(rep)
  sort!(rep, by=x->x[2])
  w = deepcopy(D.weyl_vector)
  tmp = zero_matrix(ZZ,ncols(w),1)
  for (i,s,indualDeltaR) in rep
    if indualDeltaR
      # saves a matrix multiplication
      rdual = dualDeltaR[i]
      r = deltaR[i]
      mul!(z, w, rdual)
      addmul!(w, r, z[1,1])
    else
      r = Pv[i]
      #g = (r*gramL*wt)[1,1]
      #w = w + g*r
      transpose!(tmp,w)
      mul!(tmp,gramL,tmp)
      mul!(z, r, tmp)
      addmul!(w,r, z[1,1])
    end
  end
  # both weyl vectors should lie in the positive cone.
  @assert ((D.weyl_vector)*D.data.gramL*transpose(w))[1,1]>0 "$(D.weyl_vector)    $v\n"
  return chamber(D.data, w, v)
end



@doc Markdown.doc"""
    K3_surface_automorphism_group(S::ZLat [, ample_class]) -> generators, rational curves, chambers

Compute the automorphism group of a very-general $S$-polarized K3 surface.

Further return representatives of the `Aut(X)`-orbits of (-2)-curves on `X` and
a fundamental domain for the action of Aut(X) on the set of nef L|S chambers.
This is almost a fundamental domain for Aut(X) on the nef cone.

Here very general means that `Num(X)` is isomorphic to `S` and the image of
$Aut(X) \to H^0(X,\Omega^2_X)$ is $ \pm 1$.

The function returns generators for the image of

\[f: Aut(X) \to O(Num(X)) \]

The output is represented with respect to  the basis of `S`.

Note that under our genericity assumptions the kernel of $f$ is of order at most $2$
and it is equal to $2$ if and only if $S$ is $2$-elementary.
If an ample class is given, then the generators returned preserve it.

This kind of computation can be very expensive. To print progress information
use `set_verbose_level(:K3Auto, 2)` or higher.

# Input
- `S`: a hyperbolic lattice
- `ample`: a row matrix or a vector given with respect to the ambient space of S.

"""
function K3_surface_automorphism_group(S::ZLat)
  return borcherds_method(S, 26, compute_OR=false)[2:end]
end

function K3_surface_automorphism_group(S::ZLat, ample_class::fmpq_mat)
  ample_classS = solve_left(basis_matrix(S), ample_class)
  L, S, weyl = borcherds_method_preprocessing(S, n, ample=ample_class)
  return borcherds_method(L, S, weyl, compute_OR=false)[2:end]
end

K3_surface_automorphism_group(S::ZLat, ample_class::Vector{fmpq}) = K3_surface_automorphism_group(S, matrix(QQ, 1, degree(S), ample_class))


function borcherds_method(S::ZLat, n::Integer; compute_OR=true, entropy_abort=false, max_nchambers=-1)
  @req n in [10,18,26] "n(=$(n)) must be one of 10,18 or 26"
  L, S, weyl = borcherds_method_preprocessing(S, n)
  return borcherds_method(L, S, weyl; compute_OR=compute_OR, entropy_abort=entropy_abort, max_nchambers=max_nchambers)
end

@doc Markdown.doc"""
    borcherds_method(S::ZLat, n::Integer; compute_OR=true, entropy_abort=false, max_nchambers=-1)
    borcherds_method(L::ZLat, S::ZLat, w::fmpq_mat; compute_OR=true, entropy_abort=false, max_nchambers=-1)

Compute the symmetry group of a weyl chamber up to finite index.

# Arguments
- `w`:  initial Weyl row vector represented with respect to the basis of `L`;
- `L`:  even, unimodular, hyperbolic lattice of rank n=10,18 or 26;
- `S`:  a primitive sublattice of `L`;
- `compute_OR=true`: if false take as `G` all isometries of `S` extending to `L`;
- `max_nchambers`: break the computation after `max_nchambers` are found;
- `entropy_abort` abort if an automorphism of positive entropy is found.
"""
function borcherds_method(L::ZLat, S::ZLat, w::fmpz_mat; compute_OR=true, entropy_abort=false, max_nchambers=-1)
  data = BorcherdsCtx(L, S, w, compute_OR)
  return borcherds_method(data, entropy_abort=entropy_abort, max_nchambers=max_nchambers)
end

function borcherds_method(data::BorcherdsCtx; entropy_abort::Bool, max_nchambers=-1)
  S = data.S
  # for G-sets
  F = FreeModule(ZZ,rank(S))
  # initialization
  chambers = Dict{UInt64,Vector{K3Chamber}}()
  explored = Set{K3Chamber}()
  D = chamber(data, data.weyl_vector, zero_matrix(ZZ, 1, rank(S)))
  waiting_list = [D]

  automorphisms = Set{fmpz_mat}()
  rational_curves = Set{fmpz_mat}()

  # gogo
  ncircles = 0
  ntry = 0
  nchambers = 0
  while length(waiting_list) > 0
    ntry = ntry + 1
    if mod(ntry, 100)==0
      @vprint :K3Auto 2 "largest bucket: $(maximum(length(i) for i in values(chambers))) "
      @vprint :K3Auto 1 "buckets: $(length(chambers)) explored: $(nchambers) unexplored: $(length(waiting_list)) generators: $(length(automorphisms))\n"
    end
    D = popfirst!(waiting_list)
    if D in explored
      continue
    end
    # check G-congruence
    fp = hash(fingerprint(D))  # this is the bottleneck - the computation of the walls.
    if !haskey(chambers, fp)
      chambers[fp] = K3Chamber[]
    end
    is_explored = false
    for E in chambers[fp]
      @vprint :K3Auto 4 "$(D.weyl_vector)    $(E.weyl_vector)\n"
      gg = hom(D, E)
      if length(gg) > 0
        # enough to add a single homomorphism
        if !(gg[1] in automorphisms)
          push!(automorphisms, gg[1])
          if entropy_abort
            C = lattice(rational_span(S),_common_invariant(automorphisms)[2])
            d = diagonal(rational_span(C))
            if 0 > maximum(push!([sign(i) for i in d],-1))
              @vprint :K3Auto 1 "entropy abort \n"
              return data, automorphisms, chambers, rational_curves, false
            end
          end
        end
        is_explored = true
        break
      end
    end
    if is_explored
      continue
    end
    push!(chambers[fp], D)
    push!(explored, D)
    nchambers = nchambers+1
    @vprint :K3Auto 3 "new weyl vector $(D.weyl_vector)\n"
    @vprint :K3Auto 3 "$D\n"

    autD = aut(D)
    autD = [a for a in autD if !isone(a)]
    # we need the orbits of the walls only
    if length(autD) > 0
      for f in autD
        push!(automorphisms, f)
      end
      @vprint :K3Auto 1 "Found a chamber with $(length(autD)) automorphisms\n"
      # compute the orbits
      @vprint :K3Auto 3 "computing orbits\n"
      Omega = [F(v) for v in walls(D)]
      W = gset(matrix_group(autD),Omega)
      vv = F(D.parent_wall)
      wallsDmodAutD = [representative(w).v for w in orbits(W) if !(vv in w)]
      @vprint :K3Auto 3 "done\n"
    else
      # the minus shouldn't be necessary ... but who knows?
      wallsDmodAutD = (v for v in walls(D) if !(v==D.parent_wall || -v==D.parent_wall))
    end
    # compute the adjacent chambers to be explored
    for v in wallsDmodAutD
      if -2 == (v*gram_matrix(S)*transpose(v))[1,1]
        # v comes from a rational curve
        push!(rational_curves, v)
        continue
      end
      Dv = adjacent_chamber(D, v)
      push!(waiting_list, Dv)
    end
    if max_nchambers != -1 && ntry > max_nchambers
      return data, collect(automorphisms), reduce(append!,values(chambers), init=K3Chamber[]), collect(rational_curves), true
    end
  end
  @vprint :K3Auto "$(length(automorphisms)) automorphism group generators\n"
  @vprint :K3Auto "$(nchambers) congruence classes of chambers \n"
  @vprint :K3Auto "$(length(rational_curves)) orbits of rational curves\n"
  return data, collect(automorphisms), reduce(append!,values(chambers), init=K3Chamber[]), collect(rational_curves), true
end

function _dist(V::Hecke.QuadSpace, r::fmpq_mat, h1::fmpq_mat, h2::fmpq_mat)
  if inner_product(V,h1-h2,r)!=0
    return inner_product(V, h1, r)[1,1]//inner_product(V, h1 - h2, r)[1,1]
  else
    return PosInf()
  end
end

################################################################################
# preprocessing
################################################################################



function chain_reflect(V::Hecke.QuadSpace, h1, h2, w, separating_walls::Vector{fmpq_mat})
  @hassert :K3Auto 1 inner_product(V,h1,h2)[1,1]>0
  @hassert :K3Auto 1 all(inner_product(V,h1,r)[1,1]>=0 for r in separating_walls)
  @hassert :K3Auto 1 all(inner_product(V,h2,r)[1,1]<=0 for r in separating_walls)
  di(r) = _dist(V, r, h1, h2)
  #sort!(separating_walls, by=di)
  separating_walls0 = deepcopy(separating_walls)
  while length(separating_walls0)>0
    _,i = findmax([di(r) for r in separating_walls0])
    r = separating_walls0[i]
    deleteat!(separating_walls0,i)
    if inner_product(V,h2,r)[1,1]>0
      continue
    end
    h2 = h2 + inner_product(V, h2, r)*r
    w = w + inner_product(V, w, r)*r
    separating_walls0 = [r for r in separating_walls0 if inner_product(V,h2,r)[1,1]<0]
    # should be decreasing
    # @vprint :K3Auto 1 length([s for s in separating_walls0 if 0>sign(inner_product(V,h2,s)[1,1])])
  end
  # confirm output .... since I did not yet prove this algorithm .. it looks a bit fishy
  @assert all(inner_product(V,h2,r)[1,1]>=0 for r in separating_walls)
  return h2, w
end

@doc Markdown.doc"""
    span_in_S(L, S, weyl) -> fmpq_mat

Return a basis matrix of the linear hull of $C(weyl) \cap S$.
with respect to the basis of `S`.
"""
function span_in_S(L, S, weyl)
  R = Hecke.orthogonal_submodule(L, S)
  V = ambient_space(L)
  Delta_w = _alg58(L, S, R, weyl)
  G = gram_matrix(V)
  BS = transpose(basis_matrix(S))
  prSDelta_w = [v*G*BS for v in Delta_w]

  i = zero_matrix(QQ, 0, rank(S))
  Cdual = positive_hull(reduce(vcat, prSDelta_w, init=i))

  spanC = linear_span(polarize(Cdual))
  N = length(spanC)
  if N==0
    M = zero_matrix(QQ, 0, rank(S))
  else
    M = reduce(vcat, (matrix(QQ, 1, rank(S), v.a) for v in spanC))
  end
  k, K = kernel(M)
  gensN = transpose(K)[1:k,:]
  return gensN
end

@doc Markdown.doc"""
    weyl_vector_non_degenerate(L::ZLat, S::ZLat, u0::fmpq_mat, weyl::fmpq_mat, ample0::fmpq_mat, perturbation_factor=1000)

Return an `S`-nondegenerate Weyl vector of `L`.

- u0: inner point of the chamber `C(weyl)` of `L`.
- ample0: an ample class... the weyl vector and u0 are moved to the same chamber first
- perturbation_factor: used to get a random point close to the ample vector. If it is (too) big, computations become harder.
"""
function weyl_vector_non_degenerate(L::ZLat, S::ZLat, u0::fmpq_mat, weyl::fmpq_mat,
                                    ample0::fmpq_mat, perturbation_factor=1000)
  V = ambient_space(L)
  ample = ample0
  u = u0

  @vprint :K3Auto 2 "calculating separating hyperplanes\n"
  separating_walls = separating_hyperplanes(L, u, ample, -2)
  @vprint :K3Auto 2 "moving weyl vector $(solve_left(basis_matrix(L),weyl)) towards the ample class\n"
  u, weyl = chain_reflect(V, ample, u, weyl, separating_walls)
  @vprint :K3Auto "new weyl: $(solve_left(basis_matrix(L),weyl)) \n"
  if is_S_nondegenerate(L,S,weyl)
    return weyl, u, ample
  end
  @vprint :K3Auto 2 "calculating QQDcapS\n"
  QQDcapS = lattice(V, span_in_S(L, S, weyl)*basis_matrix(S))

  N = Hecke.orthogonal_submodule(L, QQDcapS)
  N = lll(N)
  @vprint :K3Auto 2 "computing the relevant roots\n"
  @vprint :K3Auto 3 "$(gram_matrix(N))\n"

  sv = short_vectors(N^(-1//1), 2)
  @vprint :K3Auto 2 "done\n"
  relevant_roots = [matrix(QQ,1,rank(N),a[1])*basis_matrix(N) for a in sv]
  T = Hecke.orthogonal_submodule(S, QQDcapS)
  if rank(T)==0
    return weyl,u,u
  end
  @vprint :K3Auto 1 "degeneracy dimension of the chamber $(rank(T))\n"
  @label choose_h
  h = perturbation_factor*ample + matrix(QQ,1,rank(T),rand(-10:10,rank(T)))*basis_matrix(T)
  # roots orthogonal to S do not help. Therefore discard them.
  relevant_roots = [r for r in relevant_roots if inner_product(V,basis_matrix(S),r)!=0]
  if any(inner_product(V,h,r)==0 for r in relevant_roots)
    @goto choose_h
  end
  separating = fmpq_mat[r for r in relevant_roots if sign(inner_product(V, h, r)[1,1])*sign(inner_product(V, u, r)[1,1])<0]
  # fix signs
  for i in 1:length(separating)
    r = separating[i]
    if inner_product(V, u, r)[1,1]>0
      separating[i] = -r
    end
    r = separating[i]
  end
  @hassert :K3Auto 1 all(inner_product(V,h,r)[1,1] > 0 for r in separating)
  @hassert :K3Auto 1 all(inner_product(V,u,r)[1,1] < 0 for r in separating)

  u, weyl = chain_reflect(V, h, u, weyl, separating)
  @hassert :K3Auto 1 all(inner_product(V,u,r)[1,1] > 0 for r in separating)

  @assert is_S_nondegenerate(L, S, weyl)
  return weyl, u, h
end


@doc Markdown.doc"""
    weyl_vector(L::ZLat, U0::ZLat) -> weyl, u0

Return a Weyl vector of ``L`` as well as an interior point of the corresponding chamber.

For ``L`` of signature (1, 25) it uses the 23 holy constructions of the Leech lattice.

# Arguments
`L`: an even unimodular lattice of signature ``(1, 9)``, ``(1, 17)`` or ``(1, 25)``
`U0`: - a sublattice of `L` with gram matrix `[0 1; 1 -2]`
"""
# we can do the 24 constructions of the leech lattice
function weyl_vector(L::ZLat, U0::ZLat)
  @vprint :K3Auto 1 "computing an initial Weyl vector \n"
  @assert gram_matrix(U0) == QQ[0 1; 1 -2] "$(gram_matrix(U0))"
  V = ambient_space(L)
  U = U0
  R = Hecke.orthogonal_submodule(L,U)
  R0 = Hecke.orthogonal_submodule(L,U)
  if rank(L)==10
    E8 = R0
    # normalize the basis
    e8 = rescale(root_lattice(:E,8), -1)
    _, T = is_isometric_with_isometry(e8, E8, ambient_representation=false)
    E8 = lattice(V, T * basis_matrix(E8))
    B = vcat(basis_matrix(U), basis_matrix(E8))
    Bdual = inv(gram_matrix(V) * transpose(B))
    # this one does not have ample projection
    weyl = QQ[30 1 1 1 1 1 1 1 1 1] * Bdual
    @hassert :K3Auto 1 inner_product(V, weyl, weyl)[1,1] == 1240
    return weyl, weyl
  elseif rank(L) == 18
    # normalize the basis
    e8 = rescale(root_lattice(:E,8), -1)
    e8e8,_,_ = orthogonal_sum(e8, e8)
    while true
      R = Hecke.orthogonal_submodule(L,U)
      @vprint :K3Auto 1 "starting isometry test\n"
      isiso, T = is_isometric_with_isometry(e8e8, R, ambient_representation=false)
      @vprint :K3Auto 1 "done\n"
      if isiso
        E8E8 = R
        break
      end
      U = U0
      R = R0

      # compute a 2-neighbor
      v = zero_matrix(QQ,1,rank(R))
      while true
        v = matrix(QQ,1,rank(R),rand(0:1,rank(R)))
        if !iszero(v) && mod(numerator((v*gram_matrix(R)*transpose(v))[1,1]),4)==0
          break
        end
      end
      b = change_base_ring(ZZ, v*gram_matrix(R)*transpose(v)*1//4)
      A = change_base_ring(ZZ, gram_matrix(R)*transpose(v))
      b = change_base_ring(GF(2), b)
      A = change_base_ring(GF(2), A)
      x = lift(solve_left(A, b))
      v = (v + 2*x)*basis_matrix(R)
      @hassert :K3Auto 1 mod(inner_product(V,v,v)[1,1], 8)==0
      u = basis_matrix(U)
      f1 = u[1,:]
      e1 = u[2,:] + u[1,:]
      f2 = -inner_product(V, v, v)*1//4*f1 + 2*e1 + v
      @hassert :K3Auto 1 inner_product(V, f2, f2)==0

      e2 = find_section(L, f2)

      @hassert :K3Auto 2 inner_product(V, f2, e2)[1,1] == 1
      @hassert :K3Auto 2 inner_product(V,e2,e2)[1,1] == -2
      u = vcat(f2, e2)
      U = lattice(V, u)
      @hassert :K3Auto 1 gram_matrix(U) == QQ[0 1; 1 -2]
    end
    E8E8 = lattice(V, T * basis_matrix(E8E8))
    B = vcat(basis_matrix(U), basis_matrix(E8E8))
    Bdual = inv(gram_matrix(V) * transpose(B))
    # this one does not have ample projection
    weyl = QQ[30 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1] * Bdual
    @hassert :K3Auto 1 inner_product(V, weyl, weyl)[1,1] == 620
    return weyl, weyl
  elseif rank(L)==26
    while true
        R = lll(Hecke.orthogonal_submodule(L,U))
        m = minimum(rescale(R,-1))
        @vprint :K3Auto 1 "found a lattice of minimum $(m) \n"
        if m==4
          # R is isomorphic to the Leech lattice
          fu = basis_matrix(U)[1,:]
          zu = basis_matrix(U)[2,:]
          u0 = 3*fu+zu
          return fu,u0
        end
        @assert m == 2
        leech,v,h = leech_lattice(rescale(R,-1))
        # the leech lattice is the h-neighbor of R with respect to v
        # with the attached hyperbolic planes this can be engineered to give an isometry
        @hassert :K3Auto 1 mod(inner_product(V,v,v)[1,1],2*h^2)==0
        u = basis_matrix(U)
        f1 = u[1,:]
        e1 = u[2,:] + u[1,:]
        f2 = -inner_product(V, v, v)*1//(2*h)*f1 + h*e1 + v
        @hassert :K3Auto 1 inner_product(V, f2, f2)==0

        e2 = find_section(L,f2)
        u = vcat(f2,e2)
        U = lattice(V,u)
        @hassert :K3Auto 1 gram_matrix(U) == QQ[0 1; 1 -2]
      end
      E8E8E8 = lattice(V, T * basis_matrix(R))
      B = vcat(basis_matrix(U), basis_matrix(E8E8E8))
      Bdual = inv(gram_matrix(V) * transpose(B))
      # this one does not have ample projection
      weyl = QQ[30 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1] * Bdual
      @hassert :K3Auto 1 inner_product(V, weyl, weyl)[1,1] == 0
      s = find_section(L,weyl)
      u0 = 3*weyl + s
      @hassert :K3Auto 1 inner_product(V, u0,u0)[1,1] == 4
      @hassert :K3Auto 1 inner_product(V, u0,weyl)[1,1] == 1
      return weyl, u0
    end
  error("L must be even, hyperbolic unimodular of rank 10,18,26")
end


@doc Markdown.doc"""
    borcherds_method_preprocessing(S::ZLat, n::Integer; ample=nothing) -> ZLat, ZLat, fmpz_mat

Return an embedding of `S` into an even unimodular, hyperbolic lattice `L` of
rank `n=10, 18, 26` as well as an `S`-nondegenerate Weyl-vector.

The weyl vector is represented with respect to the basis of `L`.
"""
function borcherds_method_preprocessing(S::ZLat, n::Integer; ample=nothing)
  @req n in [10,18,26] "n must be one of 10, 18 or 26"
  # another example
  S = Zlattice(gram=gram_matrix(S))
  L,S,iS = embed_in_unimodular(S::ZLat, 1, n-1,primitive=true,even=true)
  V = ambient_space(L)
  if gram_matrix(S)[1:2,1:2] == QQ[0 1; 1 -2]
    U = lattice(V,basis_matrix(S)[1:2,:])
  else
    # find a hyperbolic plane
    G = gram_matrix(L)
    g1, u1 = lll_gram_indef_with_transform(change_base_ring(ZZ, G))
    # apparently we need to run lll two times to produce a zero
    g2, u2 = lll_gram_indef_with_transform(change_base_ring(ZZ, g1))
    u = u2*u1
    @assert g2 == u*G*transpose(u)
    B = u*basis_matrix(L)
    B = vcat(B[1,:],B[end,:]-B[1,:])
    if inner_product(V,B,B) != QQ[0 1; 1 -2]
      # find an isotropic vector
      if gram_matrix(S)[1,1]==0
        v = basis_matrix(S)[1,:]
      else
        b, v = Hecke._isisotropic_with_vector(gram_matrix(L))
        @assert b
        v = matrix(QQ, 1, degree(L), v)
        v = v*basis_matrix(L)
      end
      s = find_section(L, v)
      B = vcat(v, s)
    end
    U = lattice(V, B)
  end
  @assert gram_matrix(U) == QQ[0 1; 1 -2]
  weyl, u0 = weyl_vector(L, U)

  if ample isa Nothing
    @vprint :K3Auto 1 "searching a random ample vector in S\n"
    h = ample_class(S)
    hS = solve_left(basis_matrix(S),h)
    @vprint :K3Auto 1 "ample vector: $hS \n"
  else
    h = ample*basis_matrix(S)
  end
  # double check
  Q = orthogonal_submodule(S, lattice(V,h))
  @assert length(short_vectors(rescale(Q,-1),2))==0
  if (h*gram_matrix(ambient_space(L))*transpose(weyl))[1,1] < 0
    weyl = -weyl
    u0 = -u0
  end
  weyl1, u, hh = weyl_vector_non_degenerate(L, S, u0, weyl, h)


  weyl2 = change_base_ring(ZZ, solve_left(basis_matrix(L), weyl1))
  return L, S, weyl2
end


################################################################################
# ample
################################################################################

@doc Markdown.doc"""
    ample_class(S::ZLat)

Return an interior point of a Weyl chamber of the hyperbolic lattice `S`.

This means that its orthogonal complement $s^\perp$ does not contain any
``(-2)``-vectors. If `S` is the Picard lattice of a K3 surface, then the ample
cone is the (interior of) some Weyl chamber.
"""
function ample_class(S::ZLat)
  @req signature_tuple(S)[1] == 1 "lattice must be hyperbolic"
  V = ambient_space(S)
  u = basis_matrix(S)[1:2,:]
  if inner_product(V,u,u) == QQ[0 1; 1 -2]
    h = zero_matrix(QQ,1,rank(S))
    v = 3*u[1,:] + u[2,:]
    fudge = 1
    nt = 0
    while true
      a = floor(fudge//2)
      h = matrix(QQ,1,rank(S)-2,rand(-5-a:5+a,rank(S)-2))
      h = hcat(zero_matrix(QQ,1,2),h)*basis_matrix(S)
      b = inner_product(V,h,h)[1,1]
      bb = ZZ(ceil(sqrt(Float64(abs(b)))/2))+fudge
      h = h + bb*v
      @hassert :K3Auto 1 inner_product(V,h,h)[1,1]>0
      # confirm that h is in the interior of a weyl chamber,
      # i.e. check that Q does not contain any -2 vector and h^2>0
      Q = rescale(Hecke.orthogonal_submodule(S, lattice(V, h)),-1)
      @vprint :K3Auto 3 "testing ampleness $(inner_product(V,h,h)[1,1])\n"
      nt = nt+1
      if nt >10
        fudge = fudge+1
        nt = 0
      end
      sv = short_vectors(Q, 2)
      if length(sv)==0
        @vprint :K3Auto 1 "found ample class $(h)\n"
        return h
      end
      @vprint :K3Auto 3 "not ample\n"
    end
  end
  G = gram_matrix(S)
  D, B = Hecke._gram_schmidt(G,identity,true)
  i = findfirst(x->x, [d>0 for d in diagonal(D)])
  v = B[i,:]
  v = denominator(v)*v
  vsq = (v*gram_matrix(S)*transpose(v))[1,1]
  @assert vsq > 0
  # search ample
  ntry = 0
  R,x = PolynomialRing(QQ,"x")
  while true
    ntry = ntry+1
    range = 10 + floor(ntry//100)
    r = matrix(QQ, 1, rank(S), rand(-range:range,rank(S)))
    rsq = (r*gram_matrix(S)*transpose(r))[1,1]
    if rsq > 0
      h = r
    else
      rv = (r*gram_matrix(S)*transpose(v))[1,1]
      p = x^2*vsq + 2*x*rv + rsq
      rp = roots(p, CalciumQQBar)
      a = rp[1]
      b = rp[2]
      if a > b
        (a,b) = (b,a)
      end
      a = fmpz(floor(a))
      b = fmpz(ceil(b))
      if p(a) == 0  # catches the case of an integer root
        a = a -1
        @assert p(a) > 0
      end
      if p(b) == 0
        b = b + 1
        @assert p(b) > 0
      end
      if abs(a) > abs(b)
        alpha = b
      else
        alpha = a
      end
      h = alpha*v + r
      @assert (h*gram_matrix(S)*transpose(h))[1,1]>0
    end
    h = h*basis_matrix(S)
    Q = rescale(orthogonal_submodule(S, lattice(V,h)), -1)
    if length(short_vectors(Q, 2)) == 0
      return h
    end
  end
end


################################################################################
#  Elliptic Fibrations
################################################################################

@doc Markdown.doc"""
    fibration_type(NS::ZLat, f) -> rank, torsion, ADE_fibers

Return the Mordell-Weil rank, torsion subgroup and reducible singular fibers
of the Jacobian genus one fibration induced by the isotropic class `f`
on a K3 surface.

More precisely let $f \in NS$ with $f^2=0$.
Then $F = f^\perp / <f>$ is a definite lattice.
Let $R$ be its root sublattice. Its ADE type gives the reducible singular fibers of
the fibration. Let $\overline{R}$ be the primitive closure of $R$ in $NS$.
Then $\overline{R}/R$ is the torsion part of the Mordell-Weil group and finally
the Mordell-Weil rank is given by rank F - rank R.
"""
function fibration_type(NS::ZLat, f::fmpq_mat)
  p,z,n = signature_tuple(NS)
  @req p==1 && z==0 "lattice not hyperbolic"
  if rank(NS) == 2
    return 0
  end
  V = ambient_space(NS)
  # compute f^\perp / ZZf
  K = orthogonal_submodule(NS, lattice(V, f))
  fK = change_base_ring(ZZ, solve_left(basis_matrix(K), f))
  g = gcd(vec(collect(fK)))
  fK = divexact(fK, g)
  A, j = snf(abelian_group(fK))
  B = reduce(vcat, [j(i).coeff for i in gens(A)])
  Frame = lattice(V, B*basis_matrix(K))

  R = root_sublattice(Frame)
  mwl_rank = rank(NS) - 2 - rank(R)
  ade_type = root_lattice_recognition(R)[1]
  barR = primitive_closure(Frame, R)
  torsion = abelian_group(change_base_ring(ZZ,solve_left(basis_matrix(barR), basis_matrix(R))))
  return mwl_rank, snf(torsion)[1], ade_type
end


@doc Markdown.doc"""
    find_section(L::ZLat, f) ->

Given an isotropic ``f \in L`` return a vector ``x \in L`` with ``x^2 = -2``
and ``x.f=1`` if it exists.

We try to be clever and return a small ``x`` by solving a closest vector problem.
"""
function find_section(L::ZLat, f::fmpq_mat)
  V = ambient_space(L)
  @req inner_product(V, f, f)==0 "f must be isotropic"
  g = [abs(i) for i in vec(collect(inner_product(ambient_space(L),f,basis_matrix(L))))]
  if 1 in g
    i = findfirst(x->x==1,g)
    s = basis_matrix(L)[i,:]
    s = sign(inner_product(ambient_space(L),f,s)[1,1])*s
  else
    # search a smallish section using a cvp
    A = change_base_ring(ZZ,basis_matrix(L)*gram_matrix(V)*transpose(f))
    ss = solve_left(A,identity_matrix(ZZ,1))
    s = ss*basis_matrix(L)
    k, K = left_kernel(A)
    Kl = Zlattice(gram=K*transpose(K))
    # project ss to K
    sK = solve(change_base_ring(QQ,K*transpose(K)),change_base_ring(QQ,K*transpose(ss)))
    a = QQ(1)
    cv = []
    while true
      cv = Hecke.close_vectors(Kl, vec(collect(-sK)), a, check=false)
      if length(cv)>0
        break
      end
      a = a+2
    end
    sK = transpose(sK)
    v0 = 0
    for (v,_) in cv
      v = matrix(QQ, 1, rank(Kl), v)
      v1 = v+sK
      aa = (v1*gram_matrix(Kl)*transpose(v1))[1,1]
      if aa < a
        a = aa
        v0 = v
      end
    end
    s = (v0*K + ss)*basis_matrix(L)
  end

  s = s - (inner_product(V,s,s)[1,1]//2+1) * f
  # confirm the output
  @hassert :Lattice 1 inner_product(V,s,f)[1,1]==1
  @hassert :Lattice 1 inner_product(V,s,s)[1,1]==-2
  return s
end

find_section(L::ZLat, f::Vector{fmpq}) = vec(collect(find_section(L, matrix(QQ, 1, degree(L), f))))

fibration_type(NS::ZLat, f) = fibration_type(NS, matrix(QQ, 1, degree(NS), f))

################################################################################
# entropy
################################################################################

function _common_invariant(Gamma)
  return left_kernel(reduce(hcat,[g-1 for g in Gamma]))
end


function has_zero_entropy(S; rank_unimod=26)
  L, S, weyl = borcherds_method_preprocessing(S, rank_unimod)
  @vprint :K3Auto 1 "Weyl vector: $(weyl)\n"
  data, K3Autgrp, chambers, rational_curves, _ = borcherds_method(L, S, weyl, entropy_abort=true)
  C = lattice(rational_span(S), _common_invariant(K3Autgrp)[2])
  d = diagonal(rational_span(C))

  return maximum(push!([sign(i) for i in d],-1)), data, K3Autgrp, chambers, rational_curves
end
