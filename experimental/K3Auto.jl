# We work with row vectors.
@doc Markdown.doc"""

Return $\{x \in Z^n : x Q x^T + 2xb^T + c <=0\}$.
Input:
- `Q` - positive definite matrix
- `b` - vector
- `c` - rational number
"""
function quadratic_triple(Q, b, c, algorithm=:short_vectors)
  if algorithm == :short_vectors
    L, p, dist = Hecke.convert_type(Q, b, QQ(c))
    #@show ambient_space(L), basis_matrix(L), p, dist
    cv = Hecke.closest_vectors(L, p, dist)
  end
  if algorithm == :pqt
    cv = Hecke.closest_vectors(Q,b,c)
  end
  return cv
end

function alg22(gram::MatrixElem, v::MatrixElem, alpha::fmpq, d)
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
  GK = gram*transpose(K)
  Q = K * GK
  b = transpose(x) * GK
  c = (transpose(x)*gram*x)[1,1] - d
  # solve the quadratic triple
  Q = change_base_ring(QQ, Q)
  b = change_base_ring(QQ, transpose(b))
  cv = quadratic_triple(-Q, -b,-QQ(c))
  cv = [transpose(x)+transpose(matrix(u))*K for u in cv]
  #@assert all((v*gram*transpose(u))[1,1]==alpha for u in cv)
  #@assert all((u*gram*transpose(u))[1,1]>= d for u in cv)
  return [u for u in cv if (u*gram*transpose(u))[1,1]==d]
end

@doc Markdown.doc"""
Return $\{x \in S : x^2=d, x.v=\alpha \}$.

- `v` - row vector with $v^2 > 0$
"""
function alg22(S::ZLat, v::MatrixElem, alpha::fmpq, d)
  gram = gram_matrix(S)
  tmp = v*gram_matrix(ambient_space(S))*transpose(basis_matrix(S))
  v_S = solve_left(gram_matrix(S),tmp)
  sol = alg22(gram, v_S, alpha, d)
  B = basis_matrix(S)
  return [s*B for s in sol]
end

function alg23(gram::fmpq_mat, v::fmpq_mat, h::fmpq_mat, d)
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

@doc Markdown.doc"""
Return ${x in S : x^2=d, x.v>0, x.h<0}$.

- `S` - a hyperbolic lattice
- `d` - a negative integer
- `v`,`h` - vectors of positive square
"""
function alg23(S::ZLat, v::fmpq_mat, h::fmpq_mat, d)
  V = ambient_space(S)
  @assert inner_product(V,v,v)[1,1]>0
  @assert inner_product(V,h,h)[1,1]>0
  gram = gram_matrix(S)
  B = basis_matrix(S)
  vS = solve_left(B,v)
  hS = solve_left(B,h)
  return [a*B for a in alg23(gram,vS,hS,d)]
end

# alg317 ... is part of alg511
@doc Markdown.doc"""

"""
function find_basis(rays::Vector{fmpq_mat}, dim=Nothing)
  r = rays[1]
  n = ncols(r)
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

# need an example of an interesting cone

#=
gram = fmpq[-2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, -2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, -2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -2]
gram = matrix(ZZ,10,10,gram)
grami = inv(gram)
rays = [grami[i,1:10] for i in 1:10]
is_in_G(x) = true
oscar.alg318(gram, rays, is_in_G)
# returns the identity matrix

gram = QQ[-2 1 0 0; 1 -2 1 1; 0 1 -2 1; 0 1 1 -2]
rays = [inv(gram)[i,:] for i in 1:nrows(gram)]
julia> oscar.alg318(G, rays)
2-element Vector{fmpq_mat}:
 [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
 [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0]
=#

# aut_G(Cone)
alg318(gram, rays) = alg319(gram, rays, rays)
@doc Markdown.doc"""
Return if the isometry `g` of `S` acts as +-1 on the discriminant group of `S`.
"""
function is_in_G(S::ZLat,g::fmpq_mat)
  D = discriminant_group(S)
  OD = orthogonal_group(D)
  g1 = hom(D,D,[D(lift(d)*g) for d in gens(D)])
  gg = OD(g1)
  return isone(gg) || gg == OD(-matrix(one(OD)))
end

@doc Markdown.doc"""
Return $hom_G(D,E)$ for two chambers D and E.

For D = E we get $hom_G(D,D) = Aut_G(D)$.
"""
function alg319(gram::fmpq_mat, raysD::Vector{fmpq_mat}, raysE::Vector{fmpq_mat}, membership_test)
  n = ncols(gram)
  i = 0
  imgs = [zero_matrix(QQ, 0, n)]
  basis = find_basis(raysD, n)
  gram_basis = basis*gram*transpose(basis)
  # breadth first search
  # we need a depth first search with early abort
  # and one with the option to use orbits
  while i < n
    imgs_new = fmpq_mat[]
    for img in imgs
      append!(imgs_new, _one_up(gram, gram_basis, raysE, img))
    end
    imgs = imgs_new
    i = i+1
  end
  basisinv = inv(basis)
  imgs = [basisinv*f for f in imgs]
  is_in_hom_D_E(f) = all(r*f in raysD for r in raysE) # can be made faster using an interior point as in Remark 3.20
  imgs = [f for f in imgs if denominator(f)==1 && membership_test(f) && is_in_hom_D_E(f)]
  return imgs
end

function alg319(L::ZLat, S::ZLat, raysD::Vector{fmpq_mat}, raysE::Vector{fmpq_mat})
  B = basis_matrix(S)
  V = ambient_space(S)
  R = Hecke.orthogonal_submodule(L, S)
  BR = vcat(B,basis_matrix(R))
  rD = [solve_left(B,r) for r in raysD]
  rE = [solve_left(B,r) for r in raysE]
  gram = gram_matrix(S)
  SS = Zlattice(gram=gram)
  membership_test(g) = is_in_G(SS, g)
  F = alg319(gram, rD, rE, membership_test)
  i = identity_matrix(QQ,rank(R))
  result = [inv(BR)*diagonal_matrix([f,i])*BR for f in F]
  @assert all(gram_matrix(V)==f*gram_matrix(V)*transpose(f) for f in result)
  return result
end

@doc Markdown.doc"""

"""
function _one_up(gram::fmpq_mat, gram_basis::fmpq_mat, rays::Vector{fmpq_mat}, img)
  extensions = []
  k = nrows(img)
  gi = gram*transpose(img)
  for r in rays
    if k > 0
    end
    if (r*gram*transpose(r))[1,1] != gram_basis[k+1,k+1] || (k>0 && r*gi != gram_basis[k+1,1:k])
      continue
    end
    # now r has the correct inner products with what we need
    push!(extensions, vcat(img,r))
  end
  return extensions
end

@doc Markdown.doc"""
Compute Delta_w

Output:

tuples (r_S, r) where r is an element of Delta_w and r_S is the
orthogonal projection of `r` to `S`.
"""
function alg58(L::ZLat, S::ZLat, w)
  V = ambient_space(L)
  d = exponent(discriminant_group(S))
  @assert V == ambient_space(S)
  n_R = [QQ(i)//d for i in (-2*d+1):0 if mod(d*i,2)==0]
  R = Hecke.orthogonal_submodule(L,S)
  Rdual = dual(R)
  Sdual = dual(S)
  rkR = rank(R)
  delta_w = Tuple{fmpq_mat,fmpq_mat}[]
  for c in n_R
    sv = short_vectors(rescale(Rdual,-1), -c)
    append!(sv,[(-v[1],v[2]) for v in sv])
    T = typeof(sv).parameters[1].parameters[1].parameters[1]
    push!(sv,(zeros(T,rank(Rdual)),QQ(0)))
    Rc = [matrix(ZZ, 1, rkR, v[1])*basis_matrix(Rdual) for v in sv if v[2]==-c]
    for vr in Rc
      a = inner_product(V,w,vr)[1,1]
      Sdual_na = alg22(Sdual, w, 1 - a, -2 - c)
      for vs in Sdual_na
        vv = vs +  vr
        if vec(vv) in L
          push!(delta_w,(vs, vv))
        end
      end
    end
  end
  # the chamber should intersect the boundary only at the QQ-rational points
  @assert rank(S) == rank(reduce(vcat,[s[1] for s in delta_w]))
  return delta_w
end

@doc Markdown.doc"""

Return the walls of the L|S chamber induced by `v`.

Input:

- `w` - a Weyl vector in `L`.
"""
function alg511(L::ZLat, S::ZLat, w)
  Delta_w = alg58(L, S, w)
  G = gram_matrix(ambient_space(L))
  prSDelta_w = [v[1]*G for v in Delta_w]
  i = zero_matrix(QQ,0,degree(S))
  Ginv = inv(G)
  D = reduce(vcat,prSDelta_w,init=i)
  P = positive_hull(D)
  r = rays(P)
  @assert ispointed(P) "the chamber is degenerate"
  r = [matrix(QQ,1,degree(S),v)*Ginv for v in r]
  r = [v*denominator(solve_left(basis_matrix(dual(S)),v)) for v in r]
  return r
end

function is_S_nondegenerate(L::ZLat, S::ZLat, w::fmpq_mat)
  Delta_w = alg58(L, S, w)
  V = ambient_space(L)
  G = gram_matrix(V)
  prSDelta_w = [v[1]*G for v in Delta_w]
  i = zero_matrix(QQ,0,degree(S))
  Ginv = inv(G)
  D = reduce(vcat,prSDelta_w,init=i)
  P = positive_hull(D)
  r = rays(P)
  return ispointed(P)
end

@doc Markdown.doc"""

Return the walls of L which contain v.
"""
function alg513(L::ZLat, S::ZLat ,v::fmpq_mat)
  P = []
  V = ambient_space(L)
  @assert V == ambient_space(S)
  Sd = dual(S)
  Rv = lattice(V,v)
  v = basis_matrix(intersect(Sd , Rv)) # primitive in Sd
  @assert nrows(v)==1 "v not in S"
  d = inner_product(V,v,v)[1,1]
  R = Hecke.orthogonal_submodule(L,S)
  Rdual = dual(R)
  Sdual = dual(S)
  rkR = rank(R)
  Pv = fmpq_mat[]
  for alpha in 0:Int64(floor(sqrt(Float64(-2//d))))
    c = -2 - alpha^2*d
    sv = short_vectors(rescale(Rdual,-1), -c)
    append!(sv,[(-v[1],v[2]) for v in sv])
    T = typeof(sv).parameters[1].parameters[1].parameters[1]
    push!(sv,(zeros(T,rank(Rdual)),QQ(0)))
    Rc = [matrix(ZZ, 1, rkR, v[1])*basis_matrix(Rdual) for v in sv if v[2]==-c]
    for vr in Rc
      vv = alpha*v + vr
      if vv in L
        push!(Pv, vv)
      end
    end
  end
  return Pv
end

@doc Markdown.doc"""
Return a Weyl vector of the L|S chamber adjacent to `D(w)` via the wall defined by `v`.
"""
function alg514(L::ZLat, S::ZLat, w::fmpq_mat, v::fmpq_mat)
  V = ambient_space(L)
  d = degree(L)
  Pv = alg513(L, S, v)
  rep = fmpq_mat[]
  for r in Pv
    if r in rep || -r in rep
      continue
    end
    push!(rep,r)
  end
  u = matrix(QQ, 1, d, rand(Int,d))
  s(x) = inner_product(V,u,x)[1,1]//inner_product(V,w,x)[1,1]
  @assert length(unique([s(r) for r in rep]))==length(rep)
  sort!(rep, by=s)
  for r in rep
    @assert inner_product(V,r,r)[1,1]==-2
    w = w + inner_product(V, r, w)*r
  end
  return w
end

@doc Markdown.doc"""
Computes the automorphism group of a K3.

- `w0` - initial Weyl vector
"""
function alg61(L, S, w0)
  V = ambient_space(L)
  W = [[(w0,zero_matrix(QQ,1,degree(S)))]] # chambers represented by Weyl vectors
  D = alg511(L, S, w0)
  DD = [[D]]
  Gamma = alg319(L, S, D, D)
  B = Set{fmpq_mat}()
  Wl = W[1]
  Dl = DD[1]
  while length(Wl)>0
    Dlnew = Vector{fmpq_mat}[]
    Wlnew = Tuple{fmpq_mat,fmpq_mat}[]
    push!(W,Wlnew)
    push!(DD,Dlnew)
    for (w,backv) in Wl
      @assert inner_product(V,w,w)[1,1] >= 0
      Delta = alg511(L,S,w)
      autD = alg319(L, S, Delta, Delta)
      autD = [a for a in autD if !isone(a)]
      if length(autD)>0
        @show length(autD)
      end
      # TODO: iterate over orbit representatives only
      for v in Delta
        if v == backv || v ==-backv
          continue
        end
        if definesminus2hyperplane(S,v)
          push!(B,v)
          continue
        end
        w_new = alg514(L, S, w, v)
        Delta_new = alg511(L, S, w_new)
        # check G-congruence
        is_G_cong = false
        for Di in DD
          for D in Di
            gg = alg319(L, S, Delta_new, D)
            if length(gg) > 0
              append!(Gamma, gg)
              is_G_cong = true
              break
            end
          end
          if is_G_cong
            break
          end
        end
        if !is_G_cong
          # not G congruent to anything before
          @show w_new
          push!(Wlnew, (w_new,v))
          push!(Dlnew, Delta_new)
          append!(Gamma, alg319(L, S, Delta_new, Delta_new))
        end
      end
    end
    Dl, Wl = Dlnew, Wlnew
    @show length(W),length(Wl)
  end
  return Gamma, W, DD,B
end

@doc Markdown.doc"""
    definesminus2hyperplane(S::ZLat, v::fmpq_mat) -> bool

Return whether the line spanned by `v` contains a (-2) vector of `S`.
"""
function definesminus2hyperplane(S::ZLat, v::fmpq_mat)
  inS, vS = can_solve_with_solution(basis_matrix(S), v, side=:left)
  if !inS
    return false
  end
  d = denominator(vS)
  if d == 1
    r = d*vS
  else
    g = gcd([numerator(a) for a in vec(vS)])
    r = 1//g * vS
  end
  # now r is primitive in S
  r = vS *  basis_matrix(S)
  return -2 == inner_product(ambient_space(S),r,r)[1,1]
end

@doc Markdown.doc"""
Return an S-nondegenerate Weyl vector.

This seems to be a random search in essence ...  no guarantee it works.

For L_10 and L_18 we can take u0=weyl
Input:
`L` - hyperbolic even unimodular lattice of rank 10,18 or 26
`S` - hyperbolic sublattice of `L`
`u0` - interior point of the chamber defined by `weyl` in `L`
`weyl0` - a weyl vector
`ample` - ample vector in S
"""
function nondeg_weyl_shimada(L::ZLat, S::ZLat, u0::fmpq_mat, weyl0::fmpq_mat, ample::fmpq_mat;max_trys=200)
  V = ambient_space(L)
  weyl = weyl0
  ntry = 0
  while !is_S_nondegenerate(L, S, weyl)
    weyl = weyl0
    h = 20*ample+matrix(ZZ,1,rank(S),rand(-10:10,rank(S)))*basis_matrix(S)
    @show h
    if ntry > max_trys
      error("did not find an S-nondegenerate weyl vector. Since this is randomized you can retry")
    end
    ntry = ntry + 1
    R = alg23(L, u0, h, -2)
    f(r) = inner_product(V, h, r)[1,1]//inner_product(V, h-u0, r)[1,1]
    if length(Set([f(r) for r in R]))<length(R)
      continue
    end
    sort!(R, by=f)
    for r in R
      weyl = weyl + inner_product(V, weyl, r)*r
    end
    @show weyl
  end
  return weyl
  # does not check S-nondegenerateness ... so the output may be randomly wrong sometimes ... then retry again.
end


function dist(V::Hecke.QuadSpace, r::fmpq_mat, h1::fmpq_mat, h2::fmpq_mat)
  if inner_product(V,h1-h2,r)!=0
    return inner_product(V, h1, r)[1,1]//inner_product(V, h1 - h2, r)[1,1]
  else
    return PosInf()
  end
end

function chain_reflect(V::Hecke.QuadSpace, h1, h2, w, separating_walls::Vector{fmpq_mat})
  @assert inner_product(V,h1,h2)[1,1]>0
  @assert all(inner_product(V,h1,r)[1,1]>=0 for r in separating_walls)
  @assert all(inner_product(V,h2,r)[1,1]<=0 for r in separating_walls)
  di(r) = dist(V, r, h1, h2)
  sort!(separating_walls, by=di)
  separating_walls0 = copy(separating_walls)
  for k in 1:length(separating_walls)
    _,i = findmax(di(r) for r in separating_walls)
    r = separating_walls[i]
    deleteat!(separating_walls,i)
    h2 = h2 + inner_product(V, h2, r)*r
    w = w + inner_product(V, w, r)*r
    # should be decreasing
    # @show length([s for s in separating_walls0 if 0>sign(inner_product(V,h2,s)[1,1])])
  end
  # confirm output
  @assert all(inner_product(V,h2,r)[1,1]>=0 for r in separating_walls0)
  return h2, w
end

# returns QQ(D(weyl)\cap S)
function span_in_S(L, S, weyl)
  V = ambient_space(L)
  Delta_w = alg58(L, S, weyl)
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
  separating_walls = alg23(L, u, ample, -2)


  u, weyl = chain_reflect(V, ample, u, weyl, separating_walls)

  QQDcapS = span_in_S(L,S,weyl)

  N = Hecke.orthogonal_submodule(L, QQDcapS)
  sv = short_vectors(N^(-1//1), 2)
  relevant_roots = [matrix(QQ,1,rank(N),a[1])*basis_matrix(N) for a in sv]
  T = Hecke.orthogonal_submodule(S, QQDcapS)
  if rank(T)==0
    return weyl,u,u
  end
  @show rank(T)
  h = perturbation_factor*ample + matrix(QQ,1,rank(T),rand(-2:2,rank(T)))*basis_matrix(T)
  separating = fmpq_mat[r for r in relevant_roots if sign(inner_product(V, h, r)[1,1])*sign(inner_product(V, u, r)[1,1])<0]
  # fix signs
  for i in 1:length(separating)
    r = separating[i]
    if inner_product(V, u, r)[1,1]>0
      separating[i] = -r
    end
    r = separating[i]
  end
  @assert all(inner_product(V,h,r)[1,1] > 0 for r in separating)
  @assert all(inner_product(V,u,r)[1,1] < 0 for r in separating)

  u, weyl = chain_reflect(V, h, u, weyl, separating)
  @assert all(inner_product(V,u,r)[1,1] < 0 for r in separating)


  return weyl, u, h
end

isless(::PosInf, ::fmpq) = false
isless(::fmpq, ::PosInf) = true

@doc Markdown.doc"""
Embed `S` into an hyperbolic, even unimodular lattice `L` of rank $n$.

If `n` is `26`, then the orthogonal complement $R = S^\perp$ in `L` has a (-2)-vector.
Or an error is produced (does not enumerate the genus of $R$).
"""
function embed_in_unimodular(S::ZLat, n)
  r = n - rank(S)
  DS = discriminant_group(S)
  DR = rescale(DS, -1)  # discriminant group of R = S^\perp in L as predicted by Nikulin
  G = genus(DR, (0, r))  # genus of R
  R = representative(G)
  if n==26
    @assert minimum(rescale(R,-1))==2
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
  @assert abs(det(L))==1
  @assert denominator(gram_matrix(L))==1
  return L, S, iS, R, iR
end


# instead of doing neighbor steps to the one lattice we know
# we can do the 24 constructions of the leech lattice
# so basically hardcoding the neighbor steps
function weyl_vector(L, U0)
  @assert gram_matrix(U0) == QQ[0 1; 1 -2]
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
    @assert inner_product(V, weyl, weyl)[1,1] == 1240
    return weyl, weyl
  elseif rank(L) == 18
    # normalize the basis
    e8 = rescale(root_lattice(:E,8), -1)
    e8e8,_,_ = orthogonal_sum(e8, e8)
    while true
      R = Hecke.orthogonal_submodule(L,U)
      @show "starting isometry test"
      isiso, T = isisometric(e8e8, R, ambient_representation=false)
      @show "done"
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
      @assert mod(inner_product(V,v,v)[1,1],8)==0
      u = basis_matrix(U)
      f1 = u[1,:]
      e1 = u[2,:] + u[1,:]
      f2 = -inner_product(V, v, v)*1//4*f1 + 2*e1 + v
      @assert inner_product(V, f2, f2)==0

      e2 = find_section(L, f2)

      #s = change_base_ring(ZZ, basis_matrix(R)*gram_matrix(V)*transpose(f2))
      #e2 = solve_left(s, matrix(ZZ,1,1,[1]))*basis_matrix(R)
      #@assert inner_product(V, f2, e2)[1,1] == 1
      #e2 = e2 - (inner_product(V,e2,e2)[1,1]*(1//2) + 1)*f2
      u = vcat(f2,e2)
      U = lattice(V,u)
      @assert gram_matrix(U) == QQ[0 1; 1 -2]
    end
    E8E8 = lattice(V, T * basis_matrix(E8E8))
    B = vcat(basis_matrix(U), basis_matrix(E8E8))
    Bdual = inv(gram_matrix(V) * transpose(B))
    # this one does not have ample projection
    weyl = QQ[30 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1] * Bdual
    @assert inner_product(V, weyl, weyl)[1,1] == 620
    return weyl, weyl
  elseif rank(L)==26
    while true
        R = Hecke.orthogonal_submodule(L,U)
        m = minimum(rescale(R,-1))
        @show m
        if m==4
          # R is isomorphic to the Leech lattice
          # how to get u0?
          return basis_matrix(U)[1,:]
        end
        e8 = rescale(root_lattice(:E,8), -1)
        e8e8,_,_ = orthogonal_sum(e8, e8)
        e8e8e8,_,_ = orthogonal_sum(e8e8, e8)
        isiso,T = isisometric(e8e8e8, R, ambient_representation=false)
        if isiso
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
        @assert mod(inner_product(V,v,v)[1,1],8)==0
        u = basis_matrix(U)
        f1 = u[1,:]
        e1 = u[2,:] + u[1,:]
        f2 = -inner_product(V, v, v)*1//4*f1 + 2*e1 + v
        @assert inner_product(V, f2, f2)==0

        s = change_base_ring(ZZ, basis_matrix(R)*gram_matrix(V)*transpose(f2))
        e2 = solve_left(s, matrix(ZZ,1,1,[1]))*basis_matrix(R)
        @assert inner_product(V, f2, e2)[1,1] == 1
        e2 = e2 - (inner_product(V,e2,e2)[1,1]*(1//2) + 1)*f2
        u = vcat(f2,e2)
        U = lattice(V,u)
        @show u
        @assert gram_matrix(U) == QQ[0 1; 1 -2]
      end
      E8E8E8 = lattice(V, T * basis_matrix(R))
      B = vcat(basis_matrix(U), basis_matrix(E8E8E8))
      Bdual = inv(gram_matrix(V) * transpose(B))
      # this one does not have ample projection
      weyl = QQ[30 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1] * Bdual
      @assert inner_product(V, weyl, weyl)[1,1] == 0
      s = find_section(L,weyl)
      u0 = 3*weyl + s
      @assert inner_product(V, u0,u0)[1,1] == 4
      @assert inner_product(V, u0,weyl)[1,1] == 1
      return weyl, u0
    end
  error("L must be even, hyperbolic unimodular of rank 10,18,26")
end


# TODO: find a section with small coefficients (this is a cvp) ... as otherwise
# the close_vector algorithms searching for separating hyperplanes die.
function find_section(L,f)
  V = ambient_space(L)
  @assert inner_product(V,f,f)==0
  A = change_base_ring(ZZ,basis_matrix(L)*gram_matrix(V)*transpose(f))
  s = solve_left(A,identity_matrix(ZZ,1))
  @show s
  s = s*basis_matrix(L)
  k, K = left_kernel(A)

  @assert inner_product(V,s,f)[1,1]==1
  #return s,f
  s = s - (inner_product(V,s,s)[1,1]//2+1) * f
  @assert inner_product(V,s,s)[1,1]==-2
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
  @assert inner_product(V,B,B) == QQ[0 1; 1 -2]
  weyl, u0 = oscar.weyl_vector(L, U)

  #find a random ample vector ... or use a perturbation of the weyl vector?
  while true
    h = matrix(ZZ,1,rank(S),rand(-10:10, rank(S)))*basis_matrix(S)
    # confirm that h is in the interior of a weyl chamber,
    # i.e. check that Q does not contain any -2 vector and h^2>0
    if inner_product(V,h,h)[1,1]<=0
      continue
    end
    @assert 0 < inner_product(V,h,h)[1,1]
    Q = Hecke.orthogonal_submodule(S, lattice(V, h))
    if minimum(rescale(Q, -1)) > 2
      break
    end
  end
  if inner_product(V,weyl,h)[1,1]<0
    h = -h
  end
  weyl1,u,hh = oscar.nondeg_weyl_new(L,S,u0, weyl,h)
  return L,S,weyl1#L,S,u0, weyl,weyl1, h
end

# might be worth to introduce this data structure
mutable struct InducedChamber
  weyl_vector::fmpq_mat
  walls::Vector{fmpq_mat}
  aut::Vector{fmpq_mat}
  stairs::fmpq_mat  # wall to the previous level
end

function find_isotropic(L)
  V = ambient_space(L)
  while true
    v = matrix(QQ, 1, rank(L), rand(-10:10,rank(L)))*basis_matrix(L)
    if inner_product(V,v,v)==0 && v!=0
      return v
    end
  end
end

function parse_zero_entropy()
  io = open("/home/simon/Dropbox/Math/code/ZeroEntropy/giacomo lattices","r")
  s = read(io, String)
  s = split(s,"\n")
  s = [a for a in s if length(a) > 0]
  s = [split(a," ") for a in s]
  s = [[ZZ(parse(Int64,i)) for i in a] for a in s]
  res = []
  for g in s
    n = Int64(sqrt(ZZ(length(g))))
    m = matrix(ZZ,n,n,g)
    m[2,2] = -2
    push!(res,m)
  end
  return res
end
