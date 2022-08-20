export root_lattice,highest_root,coxeter_number,weyl_vector,has_zero_entropy,check_zero_entropy,common_invariant

Hecke.add_assert_scope(:K3Auto)
#set_assert_level(:K3Auto, 3)
Hecke.add_verbose_scope(:K3Auto)
#set_verbose_level(:K3Auto, 3)


################################################################################
# Types
################################################################################

mutable struct BorcherdsData
  L::ZLat
  S::ZLat
  SS::ZLat
  R::ZLat
  deltaR::Vector{fmpz_mat}
  prRdelta
  membership_test
  gramL::fmpz_mat
  gramS::fmpz_mat
  prS::fmpq_mat
end

mutable struct Chamber
  weyl_vector::fmpz_mat
  walls::Vector{fmpz_mat}  # represented as gram_S*v_S #and v_S is the S^\vee primitive minimal def. vector
  parent_wall::fmpz_mat # for the spanning tree
  data::BorcherdsData
  function Chamber()
    return new()
  end
end

function Chamber(data::BorcherdsData, weyl_vector, parent_wall)
  D = Chamber()
  D.weyl_vector = weyl_vector
  D.parent_wall = parent_wall
  D.data = data
  return D
end

function hash(C::Chamber)
  return hash(C.weyl_vector[:,1:rank(C.data.S)])
end

function Base.:(==)(C::Chamber,D::Chamber)
  @req C.data===D.data "must be in the same space"
  return C.weyl_vector[:,1:rank(C.data.S)] == D.weyl_vector[:,1:rank(D.data.S)]
end

function walls(D::Chamber)
  if !isdefined(D, :walls)
    D.walls = walls_of_chamber(D.data, D.weyl_vector)
    @assert length(D.walls)>=rank(D.data.S)
  end
  return D.walls
end

function Base.show(io::IO, c::Chamber)
  print(IOContext(io, :compact => true), "Chamber  in dimension $(length(walls(c)[1])) with $(length(walls(c))) walls")
end

@doc Markdown.doc"""

"""
function fingerprint(D::Chamber)
  v = sum(walls(D))
  G = change_base_ring(ZZ,gram_matrix(D.data.SS))
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
  m5 = []
  for i in keys(V)
    vi = sum(V[i])
    push!(m5, [i,sort!([(vi*G*transpose(j))[1,1] for j in walls(D)])])
  end
  sort!(m5)
  =#
  return hash((m1, m2, m3, m4))
end


################################################################################
# close vector functions
################################################################################

@doc Markdown.doc"""
    quadratic_triple ->

Return $\{x \in Z^n : x Q x^T + 2xb^T + c <=0\}$.
Input:
- `Q` - positive definite matrix
- `b` - vector
- `c` - rational number
"""
function quadratic_triple(Q, b, c; algorithm=:short_vectors, equal=false)
  if algorithm == :short_vectors
    L, p, dist = Hecke.convert_type(Q, b, QQ(c))
    #@vprint :K3Auto 1 ambient_space(L), basis_matrix(L), p, dist
    cv = Hecke.closest_vectors(L, p, dist, check=false, equal=equal)
    #@show isone(basis_matrix(L)),gram_matrix(L),p,dist
  end
  # using :pqt seems unfeasible
  if algorithm == :pqt
    cv = Hecke.closest_vectors(Q,b,c)#, check=false)
  end
  return cv
end

@doc Markdown.doc"""
Return $\{x \in S : x^2=d, x.v=\alpha \}$.

- `v` - row vector with $v^2 > 0$

Algorithm 2.2 in [Shimada]
"""
function short_vectors_affine(S::ZLat, v::MatrixElem, alpha::fmpq, d)
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
  cv = quadratic_triple(-Q, -b,-QQ(c),equal=true)
  cv = [transpose(x)+transpose(matrix(u))*K for u in cv]
  @hassert :K3Auto 1 all((v*gram*transpose(u))[1,1]==alpha for u in cv)
  @hassert :K3Auto 1 all((u*gram*transpose(u))[1,1]== d for u in cv)
  return cv #[u for u in cv if (u*gram*transpose(u))[1,1]==d]
end


@doc Markdown.doc"""
Return ${x in S : x^2=d, x.v>0, x.h<0}$.

- `S` - a hyperbolic lattice
- `d` - a negative integer
- `v`,`h` - vectors of positive square
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

  S_W = quadratic_triple(Q, zero_matrix(QQ,n-1,1), d)
  S_W = [transpose(matrix(x))*bW for x in S_W]
  S = fmpq_mat[]
  h = change_base_ring(QQ,h)
  for rp in S_W
    rho = abs(d - (rp*gram*transpose(rp))[1,1])*ch^-1
    t,rho = issquare_with_sqrt(rho)
    if !t
      continue
    end
    r = rho*h + rp
    if denominator(r)==1 && (r*gram*transpose(h))[1,1]>0 && (r*gram*transpose(v))[1,1] < 0
      push!(S,r)
    end
  end
  return S
end

# alg317 ... is part of walls_of_chamber
@doc Markdown.doc"""

"""
function find_basis(rays::Vector, dim=Nothing)
  r = rays[1]
  n = ncols(r)
  if dim==Nothing
    dim = n
  end
  B = zero_matrix(base_ring(r), 0, n)
  d = 0
  for r in rays
    Br = vcat(B, r)
    rk = rank(Br)
    if rk > d
      d = rk
      B = Br
    end
    if rk == dim
      break
    end
  end
  return B
end


function is_in_G(S::ZLat,g::fmpz_mat)
  D = discriminant_group(S)
  imgs = [D(vec(matrix(QQ,1,rank(S),lift(d))*g)) for d in gens(D)]
  return all(imgs[i] == gens(D)[i] for i in 1:length(gens(D))) || all(imgs[i] == -gens(D)[i] for i in 1:length(gens(D)))
  # OD = orthogonal_group(D)
  # g1 = hom(D,D,[D(lift(d)*g) for d in gens(D)])
  # gg = OD(g1)
  # return isone(gg) || gg == OD(-matrix(one(OD)))
end


hom(D::Chamber, E::Chamber) = alg319(gram_matrix(D.data.SS), walls(D), walls(E), D.data.membership_test)
aut(D::Chamber) = hom(D, D)

# worker for hom and aut
function alg319(gram::MatrixElem, raysD::Vector{fmpz_mat}, raysE::Vector{fmpz_mat}, membership_test)
  n = ncols(gram)
  partial_homs = [zero_matrix(ZZ, 0, n)]
  basis = find_basis(raysD, n)
  gram_basis = basis*gram*transpose(basis)
  # breadth first search
  # we need a depth first search with early abort
  # and one with the option to use orbits
  for i in 1:n
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
  is_in_hom_D_E(f) = all(r*f in raysE for r in raysD) # could possibly be made faster using an interior point as in Remark 3.20
  for f in partial_homs
    f = basisinv*f
    if denominator(f)!=1
      continue
    end
    fz = change_base_ring(ZZ,f)
    if !membership_test(fz) || !is_in_hom_D_E(fz)
      continue
    end
    push!(homs,fz)
  end
  @hassert :K3Auto 1 all(f*gram*transpose(f)==gram for f in homs)
  return homs
end


@doc Markdown.doc"""
Compute Delta_w

Output:

Tuples (r_S, r) where r is an element of Delta_w and r_S is the
orthogonal projection of `r` to `S`.

Algorithm 5.8 in [Shi]
"""
function alg58(L::ZLat, S::ZLat, R::ZLat,prRdelta, w; is_S_nondeg=true)
  V = ambient_space(L)
  d = exponent(discriminant_group(S))
  @hassert :K3Auto 1 V == ambient_space(S)
  n_R = [QQ(i)//d for i in (-2*d+1):0 if mod(d*i,2)==0]
  Rdual = dual(R)
  Sdual = dual(S)
  rkR = rank(R)
  delta_w = Tuple{fmpq_mat,fmpq_mat}[]
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
          if myin(vec(vv),L)
            push!(delta_w,(vs, vv))
          end
        end
      end
    end
  end
  # the chamber should intersect the boundary only at the QQ-rational points
  @hassert :K3Auto 1 !is_S_nondeg || rank(S) == rank(reduce(vcat,[s[1] for s in delta_w]))
  return delta_w
end

function alg58(L::ZLat, S::ZLat, R::ZLat, w::MatrixElem; is_S_nondeg=true)
  Rdual = dual(R)
  sv = short_vectors(rescale(Rdual,-1), 2)
  # not storing the following for efficiency
  # append!(sv,[(-v[1],v[2]) for v in sv])
  # but for convenience we include zero
  T = typeof(sv).parameters[1].parameters[1].parameters[1]
  push!(sv,(zeros(T, rank(Rdual)), QQ(0)))
  rkR = rank(R)
  prRdelta = [(matrix(QQ, 1, rkR, v[1])*basis_matrix(Rdual),v[2]) for v in sv]
  return alg58(L,S,R,prRdelta,w, is_S_nondeg=is_S_nondeg)
end


@doc Markdown.doc"""
Compute Delta_w

Output:

Tuples (r_S, r) where r is an element of Delta_w and r_S is the
orthogonal projection of `r` to `S^\vee` given in the basis of S.

Algorithm 5.8 in [Shi]
"""
function alg58(data::BorcherdsData, w::fmpz_mat)
  V = ambient_space(data.L)
  S = data.S
  d = exponent(discriminant_group(S))
  @hassert :K3Auto 1 V == ambient_space(S)
  @hassert :K3Auto 2 basis_matrix(data.L)==1
  n_R = [QQ(i)//d for i in (-2*d+1):0 if mod(d*i,2)==0]
  SSdual = dual(data.SS)
  delta_w = Tuple{fmpq_mat,fmpq_mat}[]
  wS = w*data.prS
  #wS = solve_left(gram_matrix(S),w*gram_matrix(V)*transpose(basis_matrix(S)))
  Vw = data.gramL*transpose(w)
  # since we do repeated cvp in the same lattice
  # we do the preprocessing here

  # collect the cvp inputs to avoid repeated calculation of the same cvp
  cvp_inputs = Dict{Tuple{fmpq,fmpq},Vector{fmpq_mat}}()
  for c in n_R
    cm = -c
    for (vr,vsquare) in data.prRdelta
      if vsquare != cm
        continue
      end
      a = (vr*Vw)[1,1]
      if (1-a,-2-c) in keys(cvp_inputs)
        push!(cvp_inputs[(1-a,-2-c)], vr)
      else
        cvp_inputs[(1-a,-2-c)] = [vr]
      end
      if c != 0
        if (1+a,-2-c) in keys(cvp_inputs)
          push!(cvp_inputs[(1+a,-2-c)], -vr)
        else
          cvp_inputs[(1+a,-2-c)] = [-vr]
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
  w = transpose(wS)
  tmp = FakeFmpqMat(w)
  wn = numerator(tmp)
  wd = denominator(tmp)
  _, K = left_kernel(wn)
  K = lll!(K)  # perhaps doing this has no gain?
  # (x + y*K)*gram*(x + y*K) = x gram x + 2xGKy + y K G K y

  # now I want to formulate this as a cvp
  # (x +y K) gram (x+yK) ==d
  # (x
  GK = gram*transpose(K)
  Q = K * GK
  #  Qf = change_base_ring(ZZ,Q*denominator(Q))

  for (alpha,d) in keys(cvp_inputs)
    #Sdual_na = short_vectors_affine(SSdual, wS, a, b)
    b, x = can_solve_with_solution(transpose(wn), matrix(ZZ, 1, 1, [alpha*wd]))
    if !b
      continue
    end
    b = transpose(x) * GK
    b = transpose(b)
    c = (transpose(x)*gram*x)[1,1] - d
    # solve the quadratic triple
    cv = quadratic_triple(-Q, -b,-QQ(c),equal=true)
    cv = [(transpose(x)+transpose(matrix(u))*K) for u in cv]
    @hassert :K3Auto 1 all((u*w)[1,1]==alpha for u in cv)
    @hassert :K3Auto 1 all((u*gram*transpose(u))[1,1]== d for u in cv)
    Sdual_na = [u*basis_matrix(SSdual) for u in cv]
    for vr in cvp_inputs[(alpha,d)]
      for vs in Sdual_na
        vv = vs*basis_matrix(data.S) +  vr
        if denominator(vv)==1
          push!(delta_w, (vs, vv))
        end
      end
    end
  end

  # the chamber should intersect the boundary only at the QQ-rational points
  @hassert :K3Auto 2 rank(S) == rank(reduce(vcat,[s[1] for s in delta_w]))
  return delta_w
end


@doc Markdown.doc"""

Return the walls of the L|S chamber induced by `w`.

Input:

- `w` - a Weyl vector in `L`.

Calls Polymake.
Corresponds Algorithm 5.11 in [Shi]
"""
function walls_of_chamber(data::BorcherdsData, w)
  walls1 = alg58(data, w)
  if length(walls1)==rank(data.S)
    # shortcut which avoids calling Polymake
    d = rank(data.S)
    walls = Vector{fmpz_mat}(undef,d)
    for i in 1:d
      vs = numerator(FakeFmpqMat(walls1[i][1]))
      g = gcd(vec(vs))
      if g != 1
        vs = divexact(vs, g)
      end
      walls[i] = vs
    end
    return walls
  end
  i = zero_matrix(QQ, 0, degree(data.SS))
  D = reduce(vcat, (v[1] for v in walls1), init=i)
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

function is_S_nondegenerate(L::ZLat, S::ZLat, w::fmpq_mat)
  R = Hecke.orthogonal_submodule(L, S)
  Delta_w = alg58(L, S, R, w; is_S_nondeg=false)
  V = ambient_space(L)
  G = gram_matrix(V)
  prSDelta_w = [v[1]*G for v in Delta_w]
  i = zero_matrix(QQ,0,degree(S))
  D = reduce(vcat,prSDelta_w,init=i)
  P = positive_hull(D)
  return ispointed(P)
end

function inner_point_in_S(L::ZLat, S::ZLat, w::fmpq_mat)
  R = Hecke.orthogonal_submodule(L, S)
  Delta_w = alg58(L, S, R, w)
  V = ambient_space(L)
  G = gram_matrix(V)
  prSDelta_w = [v[1]*G*transpose(basis_matrix(S)) for v in Delta_w]
  i = zero_matrix(QQ,0,rank(S))
  D = reduce(vcat,prSDelta_w,init=i)
  P = positive_hull(D)
  Pd = polarize(P)
  @hassert :K3Auto 1 is_pointed(Pd)
  h = sum(rays(Pd))
  h = matrix(QQ,1,rank(S),collect(h))*basis_matrix(S)
  @hassert :K3Auto 1 all(0<inner_product(V,v[1],h)[1,1] for v in Delta_w)
  return h
end


@doc Markdown.doc"""

Return the (-2)-walls of L containing $v^{perp_S}$.

Algorithm 5.13 in [Shi]
"""
function unproject_wall(data::BorcherdsData, vS::fmpz_mat)
  d = gcd(vec(vS*data.gramS))
  v = QQ(1//d)*(vS*basis_matrix(data.S))  # primitive in Sdual
  vsq = QQ((vS*data.gramS*transpose(vS))[1,1],d^2)

  @hassert :K3Auto 1 vsq>=-2
  rkR = rank(data.R)
  Pv = copy(data.deltaR)  # TODO: do these root matter at all?
  for alpha in 1:Int64(floor(sqrt(Float64(-2//vsq))))
    c = 2 + alpha^2*vsq
    alphav = alpha*v
    for (vr,cc) in data.prRdelta
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
Return return the L|S chamber adjacent to `D` via the wall defined by `v`.
"""
function adjacent_chamber(D::Chamber, v)
  gramL = D.data.gramL
  dimL = ncols(gramL)
  Pv = unproject_wall(D.data, v)
  @hassert :K3Auto 1 length(Pv) == length(unique(Pv))
  a = 10000
  @label getu
  a = 10*a
  rep = Tuple{Int,fmpq}[]
  u = matrix(ZZ, 1, dimL, rand(-a:a, dimL))
  Vw = gramL*transpose(D.weyl_vector)
  Vu = gramL*transpose(u)
  for i in 1:length(Pv)
    r = Pv[i]
    s = (r*Vu)[1,1]//(r*Vw)[1,1]
    if any(x[2]==s for x in rep)
      @goto getu
    end
    push!(rep,(i,s))
  end
  @hassert :K3Auto 2 length(unique([r[2] for r in rep]))==length(rep)
  sort!(rep, by=x->x[2])
  w = D.weyl_vector
  for (i,s) in rep
    r = Pv[i]
    @hassert :K3Auto 3 (r*gramL*transpose(r))[1,1]==-2
    w = w + (r*gramL*transpose(w))[1,1]*r
  end
  return Chamber(D.data, w, v)
end


function complete_to_basis(B::fmpz_mat,C::fmpz_mat)
  basis = B
  for j in 1:nrows(C)-nrows(B)
    for i in 1:nrows(C)
      c = C[i,:]
      A = vcat(basis,c)
      h = snf(A)
      if h[1:end,1:nrows(h)]==1
        basis = A
        break
      end
    end
  end
  return basis
end


@doc Markdown.doc"""
Compute the automorphism group of a K3

- `w` - initial Weyl vector
"""
function K3Auto(L::ZLat, S::ZLat, w::fmpq_mat; entropy_abort=false, compute_OR=true, max_nchambers=-1)
  # transform L to have the standard basis
  # we complete a basis of S to a basis of L
  vcat(basis_matrix(S),basis_matrix(L))
  #

  L1 = Zlattice(gram=gram_matrix(L))
  V = ambient_space(L1)
  S = lattice(V,basis_matrix(S)*inverse_basis_matrix(L))

  #=
  # if we complete R to a basis, then we can throw away the last few coords.
  # the following completes S to a basis....
  basis = complete_to_basis(change_base_ring(ZZ,basis_matrix(S)),change_base_ring(ZZ,basis_matrix(L1)))
  L2 = lattice(V, basis)
  S = lattice(V, basis[1:rank(S),:])

  L3 = Zlattice(gram=gram_matrix(L2))
  V = ambient_space(L3)
  S = lattice(V,basis_matrix(S)*inverse_basis_matrix(L2))

  w = change_base_ring(ZZ,w*inverse_basis_matrix(L)*inverse_basis_matrix(L2))
  L = L3
  #return L, S, w
  =#

  w = change_base_ring(ZZ,w*inverse_basis_matrix(L))
  L = L1


  SS = Zlattice(gram=gram_matrix(S))

  # precomputations
  R = lll(Hecke.orthogonal_submodule(L, S))
  bSR = vcat(basis_matrix(S),basis_matrix(R))
  ibSR = inv(bSR)
  I = identity_matrix(QQ,degree(L))
  prS = ibSR*I[:,1:rank(S)]#*basis_matrix(S)
  @assert prS[rank(S)+1,:]==0

  m = maximum(diagonal(-gram_matrix(R)))
  if m > 6
    @vprint :K3Auto 2 "skipping orthogonal group computation since the diagonal contains $(m)\n"
    # otherwise we run out of memory ....
    # TODO: Improve orthogonal group computation using a decompositon
    compute_OR = false
  end

  if compute_OR
    @vprint :K3Auto 3 "computing orthogonal group\n"
    OR = orthogonal_group_decomp(R)
    @vprint :K3Auto 3 "done\n"
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
    membership_test(g) = is_in_G(SS,g)
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
  deltaR = [change_base_ring(ZZ,matrix(QQ, 1, rkR, v[1])*basis_matrix(R)) for v in short_vectors(rescale(R,-1),2)]

  data = BorcherdsData(L, S, SS, R, deltaR, prRdelta, membership_test,change_base_ring(ZZ,gram_matrix(L)),change_base_ring(ZZ,gram_matrix(S)),prS)
  # for G-sets
  F = FreeModule(ZZ,rank(S))
  # initialization
  chambers = Dict{UInt64,Vector{Chamber}}()
  explored = Set{Chamber}()
  D = Chamber(data, w, zero_matrix(ZZ, 1, rank(S)))
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
    fp = fingerprint(D)
    if !haskey(chambers,fp)
      chambers[fp] = Chamber[]
    end
    is_explored = false
    for E in chambers[fp]
      gg = hom(D, E)
      if length(gg) > 0
        # enough to add a single homomorphism
        if !(gg[1] in automorphisms)
          push!(automorphisms, gg[1])
          if entropy_abort
            C = lattice(rational_span(S),common_invariant(automorphisms)[2])
            d = diagonal(rational_span(C))
            if 0 > maximum(push!([sign(i) for i in d],-1))
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
    push!(explored,D)
    nchambers = nchambers+1
    @vprint :K3Auto 3 "new weyl vector $(D.weyl_vector)\n"

    autD = aut(D)
    autD = [a for a in autD if !isone(a)]
    if length(autD) > 0
      for f in autD
        push!(automorphisms, f)
      end
      @vprint :K3Auto 1 "Found a chamber with $(length(autD)) automorphisms\n"
      # compute the orbits
      @vprint :K3Auto 2 "computing orbits"
      Omega = [F(v) for v in walls(D)]
      W = gset(matrix_group(autD),Omega)
      vv = F(D.parent_wall)
      wallsDmodAutD = [representative(w).v for w in orbits(W) if !(vv in w)]
    else
      # the minus shouldnt be necessary ... but who knows?
      wallsDmodAutD = (v for v in walls(D) if !(v==D.parent_wall || -v==D.parent_wall))
    end
    # now we need the orbits of the walls only
    # compute the adjacent chambers to be explored
    # TODO: iterate over orbit representatives only
    for v in wallsDmodAutD
      # does v come from a -2 curve?
      if -2 == (v*gram_matrix(S)*transpose(v))[1,1]
        push!(rational_curves, v)
        continue
      end
      Dv = adjacent_chamber(D, v)
      push!(waiting_list, Dv)
    end
    if max_nchambers != -1 && ntry > max_nchambers
      return data, automorphisms, chambers, rational_curves, false
    end
  end
  return data, automorphisms, chambers, rational_curves, true
end

function dist(V::Hecke.QuadSpace, r::fmpq_mat, h1::fmpq_mat, h2::fmpq_mat)
  if inner_product(V,h1-h2,r)!=0
    return inner_product(V, h1, r)[1,1]//inner_product(V, h1 - h2, r)[1,1]
  else
    return PosInf()
  end
end

function chain_reflect(V::Hecke.QuadSpace, h1, h2, w, separating_walls::Vector{fmpq_mat})
  @hassert :K3Auto 1 inner_product(V,h1,h2)[1,1]>0
  @hassert :K3Auto 1 all(inner_product(V,h1,r)[1,1]>=0 for r in separating_walls)
  @hassert :K3Auto 1 all(inner_product(V,h2,r)[1,1]<=0 for r in separating_walls)
  di(r) = dist(V, r, h1, h2)
  #sort!(separating_walls, by=di)
  separating_walls0 = deepcopy(separating_walls)
  while length(separating_walls0)>0
    _,i = findmax(di(r) for r in separating_walls0)
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

# returns QQ(D(weyl)\cap S)
function span_in_S(L, S, weyl)
  R = Hecke.orthogonal_submodule(L, S)
  V = ambient_space(L)
  Delta_w = alg58(L, S, R, weyl; is_S_nondeg=false)
  G = gram_matrix(V)
  prSDelta_w = [v[1]*G for v in Delta_w]
  i = zero_matrix(QQ, 0, degree(S))
  R = Hecke.orthogonal_submodule(L, S)
  Ddual = reduce(vcat, prSDelta_w, init=i)
  Ddual = vcat(Ddual, basis_matrix(R))
  Ddual = vcat(Ddual, -basis_matrix(R))
  Ddual = positive_hull(Ddual)
  D = polarize(Ddual)
  gensN = [matrix(QQ, 1, degree(S), v) for v in vcat(rays(D),lineality_space(D))]
  gensN = reduce(vcat, gensN, init=i)
  r = Hecke.rref!(gensN)
  gensN = gensN[1:r,:]
  QQDcapS = lattice(V, gensN)
  return QQDcapS
end

function nondeg_weyl_new(L::ZLat, S::ZLat, u0::fmpq_mat, weyl::fmpq_mat, ample0::fmpq_mat, perturbation_factor=1000)
  V = ambient_space(L)
  ample = ample0
  u = u0

  separating_walls = separating_hyperplanes(L, u, ample, -2)

  u, weyl = chain_reflect(V, ample, u, weyl, separating_walls)

  QQDcapS = span_in_S(L,S,weyl)

  N = Hecke.orthogonal_submodule(L, QQDcapS)
  sv = short_vectors(N^(-1//1), 2)
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

  @assert is_S_nondegenerate(L,S,weyl)
  return weyl, u, h
end

@doc Markdown.doc"""
Embed `S` into an hyperbolic, even unimodular lattice `L` of rank $n$.

If `n` is `26`, then the orthogonal complement $R = S^\perp$ in `L` has a (-2)-vector.
Or an error is produced (does not enumerate the genus of $R$).
"""
function embed_in_unimodular(S::ZLat, n)
  @vprint :K3Auto 1 "computing embedding in L_$(n) \n"
  r = n - rank(S)
  DS = discriminant_group(S)
  DR = rescale(DS, -1)  # discriminant group of R = S^\perp in L as predicted by Nikulin
  G = genus(DR, (0, r))  # genus of R
  R = representative(G)
  R = lll(R)
  R = Zlattice(gram=gram_matrix(R))
  if n==26 && maximum(diagonal(gram_matrix(R)))<-2
    @vprint :K3Auto 2 "checking the embedding"
    @hassert :K3Auto 1 minimum(rescale(R,-1))==2
  end
  SR, iS, iR = orthogonal_sum(S, R)
  V = ambient_space(SR)
  S = lattice(V,basis_matrix(S)*iS.matrix)
  R = lattice(V,basis_matrix(R)*iR.matrix)
  DS,_ = normal_form(discriminant_group(S))
  DR = discriminant_group(R)
  DRn,_ = normal_form(rescale(DR,-1))
  gensDS = [lift(x) for x in gens(DS)]
  imgsDS = [lift(x) for x in gens(DRn)]
  glue = reduce(vcat, [matrix(QQ,1,degree(SR),gensDS[i]+imgsDS[i]) for i in 1:length(gensDS)],init=zero_matrix(QQ,0,degree(S)))
  gensL = vcat(basis_matrix(SR), glue)
  L = lattice(V, gensL, isbasis=false)
  @hassert :K3Auto 1 abs(det(L))==1
  @hassert :K3Auto 1 denominator(gram_matrix(L))==1
  return L, S, iS, R, iR
end


@doc Markdown.doc"""
weyl_vector(L::ZLat, U0::ZLat)

Return a Weyl vector of `L`.

For `L` of signature (1,25) it uses the 24 holy constructions of the Leech
lattice.

Input:

`L` - an even unimodular lattice of signature (1,9), (1,17) or (1,25)
`U0` - a sublattice of `L` with gram matrix `[0 1; 1 -2]`
"""
# we can do the 24 constructions of the leech lattice
function weyl_vector(L::ZLat, U0::ZLat)
  @vprint :K3Auto 1 "computing an initial Weyl vector \n"
  @hassert :K3Auto 1 gram_matrix(U0) == QQ[0 1; 1 -2]
  V = ambient_space(L)
  U = U0
  R = Hecke.orthogonal_submodule(L,U)
  R0 = Hecke.orthogonal_submodule(L,U)
  if rank(L)==10
    E8 = R0
    # normalize the basis
    e8 = rescale(root_lattice(:E,8), -1)
    _, T = isisometric(e8, E8, ambient_representation=false)
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
      @vprint :K3Auto 1 "starting isometry test"
      isiso, T = isisometric(e8e8, R, ambient_representation=false)
      @vprint :K3Auto 1 "done"
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

      #s = change_base_ring(ZZ, basis_matrix(R)*gram_matrix(V)*transpose(f2))
      #e2 = solve_left(s, matrix(ZZ,1,1,[1]))*basis_matrix(R)
      @hassert :K3Auto 2 inner_product(V, f2, e2)[1,1] == 1
      @hassert :K3Auto 2 inner_product(V,e2,e2)[1,1]==-2
      #e2 = e2 - (inner_product(V,e2,e2)[1,1]*(1//2) + 1)*f2
      u = vcat(f2,e2)
      U = lattice(V,u)
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
        #=
        e8 = rescale(root_lattice(:E,8), -1)
        e8e8,_,_ = orthogonal_sum(e8, e8)
        e8e8e8,_,_ = orthogonal_sum(e8e8, e8)
        # the isometry test seems to be expensive sometimes
        isiso,T = isisometric(e8e8e8, R, ambient_representation=false)
        @vprint :K3Auto 2 root_type(R)[2]
        @vprint :K3Auto 2 "\n"
        if isiso
          break
        end
        U = U0
        R = R0
        =#
        leech,v,h = leech_from_root_lattice(rescale(R,-1))
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

function find_section(L::ZLat, f)
  V = ambient_space(L)
  g = [abs(i) for i in vec(inner_product(ambient_space(L),f,basis_matrix(L)))]
  if 1 in g
    i = findfirst(x->x==1,g)
    s = basis_matrix(L)[i,:]
    s = sign(inner_product(ambient_space(L),f,s)[1,1])*s
  else
    # search a smallish section using a cvp
    @hassert :K3Auto 1 inner_product(V,f,f)==0
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
      cv = Hecke.closest_vectors(Kl,-sK, a)
      if length(cv)>0
        break
      end
      a = a+2
    end
    sK = transpose(sK)
    v0 = 0
    for v in cv
      v = matrix(QQ,1,rank(Kl),v)
      v1 = v+sK
      aa = (v1*gram_matrix(Kl)*transpose(v1))[1,1]
      if aa < a
        a = aa
        v0 = v
      end
    end
    s = (v0*K + ss)*basis_matrix(L)
  end

  @vprint :K3Auto 2 "found section of size $(s*transpose(s))\n"
  @hassert :K3Auto 1 inner_product(V,s,f)[1,1]==1
  #return s,f
  s = s - (inner_product(V,s,s)[1,1]//2+1) * f
  @hassert :K3Auto 1 inner_product(V,s,s)[1,1]==-2
  return s
end

function preprocessingK3Auto(S, n)
  # another example
  S = Zlattice(gram=gram_matrix(S))
  L,S,iS, R,iR = oscar.embed_in_unimodular(S::ZLat, n)
  V = ambient_space(L)
  # find a hyperbolic plane
  G = gram_matrix(L)
  g,u = oscar.lll_gram_indefinite(change_base_ring(ZZ,G))
  B = transpose(u)*basis_matrix(L)
  B = vcat(B[1,:],B[end,:]-B[1,:])
  U = lattice(V, B)
  @hassert :K3Auto 1 inner_product(V,B,B) == QQ[0 1; 1 -2]
  weyl, u0 = oscar.weyl_vector(L, U)

  #find a random ample vector ... or use a perturbation of the weyl vector?
  while true
    h = matrix(ZZ,1,rank(S),rand(-10:10, rank(S)))*basis_matrix(S)
    # confirm that h is in the interior of a weyl chamber,
    # i.e. check that Q does not contain any -2 vector and h^2>0
    if inner_product(V,h,h)[1,1]<=0
      continue
    end
    @hassert :K3Auto 1 0 < inner_product(V,h,h)[1,1]
    Q = Hecke.orthogonal_submodule(S, lattice(V, h))
    if length(short_vectors(rescale(Q, -1), 2)) == 0
      break
    end
  end
  if inner_product(V,weyl,h)[1,1]<0
    h = -h
  end
  weyl1,u,hh = oscar.nondeg_weyl_new(L,S,u0, weyl,h)
  return L,S,weyl1#L,S,u0, weyl,weyl1, h
end


function parse_zero_entropy(filename="/home/simon/Dropbox/Math/MyPapers/zero entropy/CandidatesZeroEntropy_elliptic")
  io = open(filename,"r")
  s = read(io, String)
  s = split(s,"\n")
  s = [a for a in s if length(a) > 0 && a[1:1]!="/"]
  s = [split(a," ") for a in s]
  s = [[ZZ(parse(Int64,i)) for i in a] for a in s]
  res = []
  u = ZZ[0 1; 1 -2]
  for g in s
    n = Int64(sqrt(ZZ(length(g))))
    m = matrix(ZZ,n,n,g)
    m = block_diagonal_matrix([u,-m])
    push!(res,m)
  end
  return res
end

function root_type(L::ZLat)
  return _connected_components(root_sublattice(L))
end

function _connected_components(L::ZLat)
  L = lll(L)
  V = ambient_space(L)
  B = basis_matrix(L)
  B = [B[i,:] for i in 1:nrows(B)]
  C = fmpq_mat[]
  SS = ZLat[]
  ADE = Tuple{Symbol,Int64}[]
  while length(B) > 0
    CC = fmpq_mat[]
    b = pop!(B)
    push!(CC, b)
    flag = true
    while flag
      flag = false
      for c in B
        if any([inner_product(V,a,c)!=0 for a in CC])
          push!(CC,c)
          deleteat!(B,findfirst(==(c),B))
          flag = true
          break
        end
      end
    end
    S = lattice(ambient_space(L),reduce(vcat,CC))
    ade, trafo = ADE_type_with_isometry(S)
    push!(ADE, ade)
    BS = trafo*basis_matrix(S)
    S = lattice(ambient_space(L),BS)
    push!(C, BS)
    push!(SS,S)
  end
  c = reduce(vcat, C)
  @hassert :K3Auto 1 nrows(c)==rank(L)
  return lattice(V, c), ADE,SS
end


function ADE_type(G)
  r = rank(G)
  d = abs(det(G))
  if r == 8 && d==1
    return (:E,8)
  end
  if r == 7 && d == 2
    return (:E,7)
  end
  if r == 6 && d ==3
    return (:E,6)
  end
  if d == r + 1
    return (:A, r)
  end
  if d == 4
    return (:D, r)
  end
  error("not a definite root lattice")
end

function ADE_type_with_isometry(L)
  ADE = ADE_type(gram_matrix(L))
  R = root_lattice(ADE...)
  e = sign(gram_matrix(L)[1,1])
  if e == -1
    R = rescale(R,-1)
  end
  t, T = is_isometric(R,L,ambient_representation=false)
  @hassert :K3Auto 1 t
  return ADE, T
end

function root_sublattice(L::ZLat)
  V = ambient_space(L)
  if is_negative_definite(L)
    L = rescale(L,-1)
  end
  sv = matrix(ZZ,reduce(hcat,[a[1] for a in short_vectors(L, 2)]))
  sv = transpose(sv)
  sv = hnf(sv)[1:rank(L),:]*basis_matrix(L)
  return lattice(V,sv)
end

# 23 constructions of the leech lattice
function coxeter_number(ADE::Symbol, n)
  if ADE == :A
    return n+1
  elseif ADE == :D
    return 2*(n-1)
  elseif ADE == :E && n == 6
    return 12
  elseif ADE == :E && n == 7
    return 18
  elseif ADE == :E && n == 8
    return 30
  end
end

function highest_root(ADE::Symbol, n)
  if ADE == :A
    w = [1 for i in 1:n]
  elseif ADE == :D
    w = vcat([1,1],[2 for i in 3:n-1])
    w = vcat(w,[1])
  elseif ADE == :E && n == 6
    w = [1,2,3,2,1,2]
  elseif ADE == :E && n == 7
    w = [2,3,4,3,2,1,2]
  elseif ADE == :E && n == 8
    w = [2,4,6,5,4,3,2,3]
  end
  w = matrix(ZZ, 1, n, w)
  g = gram_matrix(root_lattice(ADE,n))
  @hassert :K3Auto 2 all(0<=i for i in collect(w*g))
  @hassert :K3Auto 2 (w*g*transpose(w))[1,1]==2
  return w
end

function weyl_vector(R::ZLat)
  weyl = matrix(ZZ,1,rank(R),ones(1,rank(R)))*inv(gram_matrix(R))
  return weyl*basis_matrix(R)
end

function leech_from_root_lattice(N::ZLat)
  # construct the leech lattice from one of the 23 holy constructions in SPLAG
  # we follow Ebeling
  # there seem to be some signs wrong in Ebeling?
  V = ambient_space(N)
  ADE, ade, RR = root_type(N)
  F = basis_matrix(ADE)
  for i in 1:length(ade)
    F = vcat(F, -highest_root(ade[i]...)*basis_matrix(RR[i]))
  end
  rho = sum(weyl_vector(r) for r in RR)
  h = coxeter_number(ade[1]...)
  @hassert :K3Auto 1 inner_product(V,rho,rho)== 2*h*(h+1)
  @hassert :K3Auto 1 all(h==coxeter_number(i...) for i in ade)
  rhoB = solve_left(basis_matrix(N),rho)
  v = QQ(1,h)*transpose(rhoB)
  A = Zlattice(gram=gram_matrix(N))
  c = QQ(2*(1+1//h))
  sv = [matrix(QQ,1,24,vec(v)-i)*basis_matrix(N) for i in Hecke.closest_vectors(A, v ,c,equal=true)]
  @hassert :K3Auto 1 all(inner_product(V,i,i)==2*(1+1//h) for i in sv)
  @hassert :K3Auto 1 length(sv)^2 == abs(det(ADE))
  G = reduce(vcat,sv)
  FG = vcat(F,G)
  K = transpose(kernel(matrix(ZZ,ones(Int,1,nrows(FG))))[2])
  B = change_base_ring(QQ,K)*FG
  B = hnf(FakeFmpqMat(B))
  B = QQ(1,B.den)*change_base_ring(QQ,B.num[end-23:end,:])
  @hassert :K3Auto 1 rank(B)==24
  lambda = lattice(V,B)
  @hassert :K3Auto 1 denominator(gram_matrix(lambda))==1
  lambda = lll(lambda)
  @hassert :K3Auto 1 det(lambda)==1
  @hassert :K3Auto 1 minimum(lambda)==4

  T = torsion_quadratic_module(lambda,intersect(lambda,N))
  @hassert :K3Auto 1 length(gens(T))==1 "I just expect this ... but did not really prove it"
  w = transpose(matrix(lift(gens(T)[1])))

  vN = matrix(ZZ,hcat(ones(Int,1,nrows(F)),zeros(Int,1,nrows(G))))
  vleech = matrix(ZZ,hcat(zeros(Int,1,nrows(F)),ones(Int,1,nrows(G))))
  K = transpose(kernel(vcat(vN,vleech))[2])
  return lambda, h*w, h
end

function common_invariant(Gamma)
  return left_kernel(reduce(hcat,[g-1 for g in Gamma]))
end

function has_zero_entropy(S; rank_unimod=26, preprocessing_only = false)
  L,S,iS,R,iR = oscar.embed_in_unimodular(S,rank_unimod)
  V = ambient_space(L)
  U = lattice(V,basis_matrix(S)[1:2, :])
  @hassert :K3Auto 1 det(U)==-1
  weyl,u0 = oscar.weyl_vector(L, U)
  #v = matrix(QQ,ones(Int,1,rank(S)-2))*inv(gram_matrix(S)[3:end,3:end])
  #v = denominator(v)*v
  #h = hcat(QQ[0 0 ], v) *basis_matrix(S)  #an ample vector
  u = basis_matrix(U)
  h = zero_matrix(QQ,1,rank(S))
  v = 3*u[1,:] + u[2,:]
  fudge = 1
  nt = 0
  while true
    h = matrix(QQ,1,rank(S)-2,rand(-5:5,rank(S)-2))
    h = hcat(zero_matrix(QQ,1,2),h)*basis_matrix(S)
    b = inner_product(V,h,h)[1,1]
    bb = ZZ(ceil(sqrt(Float64(abs(b)))/2))+fudge
    h = h + bb*v
    @hassert :K3Auto 1 inner_product(V,h,h)[1,1]>0
    # confirm that h is in the interior of a weyl chamber,
    # i.e. check that Q does not contain any -2 vector and h^2>0
    Q = rescale(Hecke.orthogonal_submodule(S, lattice(V, h)),-1)
    Q = lll(Q)
    @vprint :K3Auto 1 "testing ampleness $(inner_product(V,h,h)[1,1])\n"
    sv = short_vectors(Q,2)
    nt = nt+1
    if nt >10
      fudge = fudge+1
      nt = 0
    end
    if length(sv)>0
      @vprint :K3Auto 1 "not ample\n"
      continue
    end
    @vprint :K3Auto 1 "found ample class $(h)\n"
    @vprint :K3Auto 1 "computing an S-non-degenerate weyl vector\n"
    weyl1,u1 = oscar.nondeg_weyl_new(L,S,u0,weyl,h)
    if is_S_nondegenerate(L,S,weyl1)
      weyl = weyl1
      break
    end
  end
  @assert is_S_nondegenerate(L,S,weyl)

  @vprint :K3Auto 1 "preprocessing completed \n"
  if preprocessing_only
    return L,S,weyl
  end
  data, K3Autgrp, chambers, rational_curves, _ = oscar.K3Auto(L,S,weyl, entropy_abort=true)
  C = lattice(rational_span(S),common_invariant(K3Autgrp)[2])
  d = diagonal(rational_span(C))

  return maximum(push!([sign(i) for i in d],-1)), data, K3Autgrp, chambers, rational_curves
end


function check_zero_entropy(candidate::ZLat, filename="")
  z, data, K3Autgrp, chambers, rational_curves = has_zero_entropy(candidate)
  io = open(filename, "w")
  println(io, gram_matrix(candidate))
  if z > 0
    println(io, "elliptic")
  elseif z == 0
    println(io, "parabolic")
  elseif z < 0
    println(io, "hyperbolic")
  end
  close(io)
end

function check_zero_entropy(candidates::Vector,postfix="",wa="a")
  ioelliptic = open("elliptic$(postfix)", wa)
  ioparabolic = open("parabolic$(postfix)", wa)
  iohyperbolic = open("hyperbolic$(postfix)", wa)
  close(ioelliptic)
  close(ioparabolic)
  close(iohyperbolic)
  for S in candidates
    e = has_zero_entropy(S)[1]
    ioelliptic = open("elliptic$(postfix)", "a")
    ioparabolic = open("parabolic$(postfix)", "a")
    iohyperbolic = open("hyperbolic$(postfix)", "a")
    if e>0
      println(ioelliptic, gram_matrix(S))
    elseif e==0
      println(ioparabolic, gram_matrix(S))
    elseif e < 0
      println(iohyperbolic, gram_matrix(S))
    end
    close(ioelliptic)
    close(ioparabolic)
    close(iohyperbolic)
  end
end

@attr function inverse_basis_matrix(L::ZLat)
  @hassert :K3Auto 1 degree(L) == rank(L)
  return inv(basis_matrix(L))::fmpq_mat
end

function myin(v, L::ZLat)
  return all(denominator(i)==1 for i in v*inverse_basis_matrix(L))
end

function glue_map(L,S,R)
  rem = Hecke.orthogonal_submodule(lattice(ambient_space(L)),S+R)
  bSR = vcat(basis_matrix(S),basis_matrix(R),basis_matrix(rem))
  ibSR = inv(bSR)
  I = identity_matrix(QQ,degree(L))
  prS = ibSR*I[:,1:rank(S)]*basis_matrix(S)
  prR = ibSR*I[:,rank(S)+1:rank(R)+rank(S)]*basis_matrix(R)
  bL = basis_matrix(L)
  DS = discriminant_group(S)
  DR = discriminant_group(R)
  gens = Hecke.TorQuadModElem[]
  imgs = Hecke.TorQuadModElem[]
  for i in 1:rank(L)
    d = bL[i,:]
    g = DS(vec(d*prS))
    if all(0==i for i in lift(g))
      continue
    end
    push!(gens,g)
    push!(imgs,DR(vec(d*prR)))
  end
  HS,iS = sub(DS,gens)
  HR,iR = sub(DR,imgs)
  glue_map = hom(HS, HR, [HR(lift(i)) for i in imgs])
  @assert overlattice(glue_map) == L
  return glue_map,iS,iR
end

function overlattice(glue_map)
  S = relations(domain(glue_map))
  R = relations(codomain(glue_map))
  glue = [lift(g)+lift(glue_map(g)) for g in gens(domain(glue_map))]
  z = zero_matrix(QQ,0,degree(S))
  glue = reduce(vcat,[matrix(QQ,1,degree(S),g) for g in glue],init=z)
  glue = vcat(basis_matrix(S+R),glue)
  glue = FakeFmpqMat(glue)
  B = hnf(glue)
  B = QQ(1,denominator(glue))*change_base_ring(QQ,numerator(B))
  return lattice(ambient_space(S),B[end-rank(S)-rank(R)+1:end,:])
end


function orthogonal_group_decomp(L)
  if gram_matrix(L)[1,1]<0
    L = rescale(L,-1)
  end
  L = lll(L)
  V = ambient_space(L)
  G = gram_matrix(L)
  mi = minimum(diagonal(G))
  ma = maximum(diagonal(G))
  sv = short_vectors(L,mi)
  h =  hnf(matrix(ZZ,transpose(reduce(hcat,(v[1] for v in sv)))))
  h = h[1:rank(h),:]*basis_matrix(L)
  M1 = lattice(V,h)
  if rank(M1)==rank(L)
    return M1
  end
  M2 = lll(Hecke.orthogonal_submodule(L,M1))
  phi,i1,i2 = glue_map(L,M1,M2)

  H1 = domain(phi)
  H2 = codomain(phi)
  O1 = orthogonal_group(M1)
  O2 = orthogonal_group(M2)

  # could also be done on the level of discriminant groups
  # this leads to too many generators
  # ... and reducing their number seems infeasible
  # first project to the discriminant_group and then lift?
  G1,_ = stabilizer(O1,cover(H1), on_sublattices)
  G2,_ = stabilizer(O2,cover(H2), on_sublattices)
  set_nice_monomorphism(M1,G1)
  set_nice_monomorphism(M2,G2)

  G1q =  _orthogonal_group(H1, fmpz_mat[hom(H1,H1,Hecke.TorQuadModElem[H1(lift(x)*matrix(g)) for x in gens(H1)]).map_ab.map for g in gens(G1)])
  G2q =  _orthogonal_group(H2, fmpz_mat[hom(H2,H2,Hecke.TorQuadModElem[H2(lift(x)*matrix(g)) for x in gens(H2)]).map_ab.map for g in gens(G2)])


  psi1 = hom(G1, G1q, gens(G1q), check=false)
  psi2 = hom(G2, G2q, gens(G2q), check=false)

  K = [matrix(g) for g in gens(kernel(psi1)[1])]

  append!(K,[matrix(g) for g in gens(kernel(psi2)[1])])
  append!([preimage(psi1,g)*preimage(psi2, G2q(inv(phi)*hom(g)*phi)) for g in gens(G1q)])
  G = matrix_group(K)
  @assert all(on_sublattices(L,g)==L for g in gens(G))
  F = free_module(QQ, degree(L))
  sv_decomp = [F(matrix(QQ,1,rank(M1),v[1])*basis_matrix(M1)) for v in short_vectors(M1, mi)]
  append!(sv_decomp, [F(matrix(QQ,1,rank(M2),v[1])*basis_matrix(M2)) for v in short_vectors(M2, ma)])
  set_nice_monomorphism(L,G, sv_decomp)
  return G
end

function on_sublattices(L, g::MatrixGroupElem{fmpq,fmpq_mat})
  V = ambient_space(L)
  return lattice(V,basis_matrix(L)*matrix(g), check=false)
end

function set_nice_monomorphism(L, G, svF=Nothing)
  if svF === Nothing
    sv = short_vectors(L,minimum(diagonal(gram_matrix(L))))
    F = free_module(QQ,degree(L))
    svF = [F(matrix(QQ,1,rank(L),x[1])*basis_matrix(L)) for x in sv]
  end
  # TODO: custom action function for Vectors?
  phi = action_homomorphism(gset(G,svF))
  GAP.Globals.SetIsHandledByNiceMonomorphism(G.X,true)
  GAP.Globals.SetNiceMonomorphism(G.X,phi.map)
end

#=
pyexec("L = []",Main)
pyexec("def callback(new_sol_coord): L.append(new_sol_coord); return True",Main)
cb = pyeval("callback",Main)
=#
