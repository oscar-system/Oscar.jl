# We work with row vectors.
@doc Markdown.doc"""

Return $\{x \in Z^n : x Q x^T + 2xb^T + c <=0\}$.
Input:
- `Q` - positive definite matrix
- `b` - vector
- `c` - rational number
"""
function quadratic_triple(Q, b, c, algorithm=:short_vector)
  L, p, dist = Hecke.convert_type(Q, b, QQ(c))
  cv = Hecke.closest_vectors(L, p, dist)
  return cv
end

function alg22(gram::MatrixElem, v::MatrixElem, alpha::fmpq, d)
  # find a solution <x,v> = alpha with x in L
  w = gram*transpose(v)
  tmp = FakeFmpqMat(w)
  wn = numerator(tmp)
  wd = denominator(tmp)
  x = solve(transpose(wn), matrix(ZZ, 1, 1, [alpha*wd]))
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
  S = []
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
Return if the isometry `g` acts as +-1 on the discriminant group `S`.
"""
function is_in_G(S::ZLat,g::fmpq_mat)
  D = discriminant_group(S)
  OD = orthogonal_group(D)
  g1 = hom(D,D,[D(g*lift(d)) for d in gens(D)])
  gg = OD(g1)
  return isone(gg) || gg == OD(-matrix(one(OD)))
end

@doc Markdown.doc"""
Return hom_G(D,E) for two chambers D and E.

For D=E we get hom(D,D) = Aut_G(D).
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
  imgs = [f for f in imgs if membership_test(f) && is_in_hom_D_E(f)]
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
    Dlnew = []
    Wlnew = []
    for (w,backv) in Wl
      @assert inner_product(V,w,w)[1,1] > 0
      Delta = alg511(L,S,w)
      autD = alg319(L, S, Delta, Delta)
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
        @show w_new
        Delta_new = alg511(L, S, w_new)
        # check G-congruence
        flag = false
        for Di in DD
          for D in Di
            gg = alg319(L, S, Delta_new, D)
            if length(gg) > 0
              append!(Gamma, gg)
              flag = true
              break
            end
          end
          if flag
            break
          end
        end
        if !flag
          # not G congruent to anything before
          push!(Wlnew, (w_new,v))
          push!(Dlnew, Delta_new)
          append!(Gamma, alg319(L, S, Delta_new, Delta_new))
        end
      end
    end
    Dl, Wl = Dlnew, Wlnew
    push!(W,Wl)
    push!(DD,Dl)
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
function nondeg_weyl(L::ZLat, S::ZLat, u0::fmpq_mat, weyl0::fmpq_mat, ample::fmpq_mat;max_trys=200)
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

@doc Markdown.doc"""
Embed `S` into an hyperbolic, even unimodular lattice of rank $n$.
"""
function embed_in_unimodular(S::ZLat, n)
  r = n - rank(S)
  DS = discriminant_group(S)
  DR = rescale(DS, -1)  # discriminant group of R = S^\perp in L as predicted by Nikulin
  G = genus(DR, (0, r))  # genus of R
  R = representative(G)
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
  return L,iS,iR
end


# might be worth to introduce this data structure
mutable struct InducedChamber
  weyl_vector::fmpq_mat
  walls::Vector{fmpq_mat}
  aut::Vector{fmpq_mat}
  stairs::fmpq_mat  # wall to the previous level
end
