
################################################################################
#
# Computation of Automorphism groups of Enriques surfaces following
# https://doi.org/10.1307/mmj/20195769
# and
# http://www.arxiv.org/abs/1909.10813
#
################################################################################
solve_left(A,B)=solve(A,B;side=:left)

mutable struct EnriquesBorcherdsCtx
  # SY < SX < L26
  L26::ZZLat
  SX::ZZLat
  SY::ZZLat
  # the following are given with respect to the basis of SY
  initial_walls::Vector{ZZMatrix}
  initial_rays::Vector{ZZMatrix}
  initial_isotropic_rays::Vector{ZZMatrix}
  # automorphisms of the initial chamber
  # written w.r.t. the basis of SY
  initial_automorphisms::Vector{ZZMatrix}
  initial_automorphisms_mod2::Vector{FqMatrix}  # same order as above but mod 2
  initial_chamber::K3Chamber   # a non-degenerate L26/SX chamber
  # whether a -2 vector defines an outer wall is a mod 2 condition
  roots_mod2::Set{FqMatrix}  # mod 2 images of the rational curves
  gramSY::ZZMatrix  # to avoid coercions
  gramSX::ZZMatrix  # to avoid coercions
  membership_test
  # needed for the membership test
  # the following are in SY/2SY w.r.t the basis of SY mod 2
  Dplus_perp::FqMatrix
  Dplus::FqMatrix
  imgs_mod2::Set{FqMatrix}
  # w.r.t the basis given by Dplus
  Gplus_mat::MatrixGroup{FqFieldElem, FqMatrix}
  function EnriquesBorcherdsCtx()
    return new()
  end

  volume_index::ZZRingElem
  orderGbar::ZZRingElem
end

function Base.show(io::IO, dat::EnriquesBorcherdsCtx)
  print(io, "Enriques Borcherds context with det(SX) = $(det(dat.SX)).")
end

function deltaYbarplus(SY,SX)
  Sm = orthogonal_submodule(SX, SY)
  phi, inc_Dminus, inc_Dplus = glue_map(SX, Sm, SY)
  # H_Sm = pi_Sm(SX) note that H_Sm/Sm is
  H_Sm = cover(domain(phi))
  sv2 = [1//2*(i[1]*basis_matrix(Sm)) for i in short_vectors(rescale(Sm,-1),4) if i[2]==4]
  sv2 = [domain(phi)(i) for i in sv2 if i in H_Sm]
  return roots_mod2 = Set([change_base_ring(GF(2),solve_left(basis_matrix(SY),matrix(QQ,1,26,2*lift(phi(i))))) for i in sv2])
end
  
function tau_type(SY, SX)
  Sm = orthogonal_submodule(SX, SY)
  phi, inc_Dminus, inc_Dplus = glue_map(SX, Sm, SY)
  # H_Sm = pi_Sm(SX) note that H_Sm/Sm is
  H_Sm = cover(domain(phi))
  sv2 = [1//2*(i[1]*basis_matrix(Sm)) for i in short_vectors(rescale(Sm,-1),4) if i[2]==4]
  V = ambient_space(Sm)
  R =  rescale(lattice(V, 2*matrix(QQ,length(sv2),dim(V),transpose(reduce(hcat,sv2)));isbasis=false),1//2)
  return R, root_lattice_recognition(R)
end

function EnriquesBorcherdsCtx(SY::ZZLat, SX::ZZLat, L26::ZZLat, weyl::ZZMatrix)
  # X K3 ---> Y Enriques
  ECtx = EnriquesBorcherdsCtx()
  ECtx.L26 = L26
  ECtx.SX = SX
  ECtx.SY = SY
  ECtx.gramSY = change_base_ring(ZZ, gram_matrix(SY))
  ECtx.gramSX = change_base_ring(ZZ, gram_matrix(SX))

  @vprint :K3Auto 2 "computing Borcherds context\n"
  dataY,_ = BorcherdsCtx(L26, SY, weyl; compute_OR=false)
  dataY.membership_test = (x -> true)
  ECtx.initial_chamber = chamber(dataY, dataY.weyl_vector; check=true)

  # SY + Sm < SX is a primitive extension with glue map phi: D(Sm) -> D(SY)
  Sm = orthogonal_submodule(SX, SY)
  phi, inc_Dminus, inc_Dplus = glue_map(SX, Sm, SY)
  # H_Sm = pi_Sm(SX) note that H_Sm/Sm is
  H_Sm = cover(domain(phi))
  sv2 = [1//2*(i[1]*basis_matrix(Sm)) for i in short_vectors(rescale(Sm,-1),4) if i[2]==4]
  sv2 = [domain(phi)(i) for i in sv2 if i in H_Sm]
  ECtx.roots_mod2 = Set([change_base_ring(GF(2),solve_left(basis_matrix(SY),matrix(QQ,1,26,2*lift(phi(i))))) for i in sv2])

  # Cook up the membership test
  @vprint :K3Auto 2 "computing orthogonal group\n"
  OSm = orthogonal_group(Sm)
  GSm,_= image_in_Oq(Sm)
  phiSm = hom(OSm, GSm, GSm.(gens(OSm)); check=false)
  DSm = discriminant_group(Sm)
  Dminus = domain(phi)
  @vprint :K3Auto 2 "computing stabilizer "
  stab_Dminus, inc_stab = stabilizer(GSm, inc_Dminus)
  @vprint :K3Auto 2 "done\n"
  Dminus_perp, inc_Dminus_perp = orthogonal_submodule(DSm, Dminus)
  res, inc = restrict_automorphism_group(stab_Dminus, inc_Dminus_perp)
  res2, inc2 = restrict_automorphism_group(stab_Dminus, inc_Dminus)
  @vprint :K3Auto 2 "computing kernel "
  Gminus, inc_Gminus = kernel(inc)
  @vprintln :K3Auto 2 "done"
  ECtx.orderGbar = order(Gminus)

  # compute Gplus
  Q = orthogonal_submodule(L26,SY)
  DQ = discriminant_group(Q)
  glue_SY_Q,inc1,inc2 = glue_map(L26,SY, Q)
  DO = codomain(glue_SY_Q)

  @vprint :K3Auto 2 "preparing membership test\n"
  DSY = discriminant_group(SY)
  ODSY = orthogonal_group(DSY)
  tmp = [preimage(phiSm, i) for i in small_generating_set(Gminus)]
  tmp_hom = [hom(DO,DO, [DO(lift(j)*i) for j in gens(DO)]) for i in tmp]
  gens_Gplus = elem_type(ODSY)[ODSY(inv(inc1)*glue_SY_Q*i*inv(glue_SY_Q)*inc1) for i in tmp_hom]
  Gplus,_ = sub(ODSY, gens_Gplus)
  @assert order(Gplus) == order(Gminus)
  # iso_Gplus_minus = hom(Gplus, Gminus, gens(Gplus), small_generating_set(Gminus))

  Dplus = codomain(phi)
  ODplus = orthogonal_group(Dplus)
  # turn this into a mod 2 matrix group
  # with respect to the basis of SY mod 2
  B = 1//2*basis_matrix(SY)
  gens_Gplus_mat = FqMatrix[]
  for g in gens(Gplus)
    g2 = reduce(vcat,[change_base_ring(GF(2),solve_left(B,matrix(QQ,1,26, lift(g(DSY(vec(B[i,:])))))))  for i in 1:10])
    push!(gens_Gplus_mat, g2)
  end
  Gplus_mat = matrix_group(GF(2),10, gens_Gplus_mat)
  ECtx.Gplus_mat = Gplus_mat

  snf_Dplus, i_snf = snf(Dplus)
  gens_Dplus = [change_base_ring(GF(2),solve_left(basis_matrix(SY),matrix(QQ,1,26,2*lift(inc_Dplus(i_snf(i)))))) for i in gens(snf_Dplus)]

  gens_Dplus = reduce(vcat, gens_Dplus)
  ECtx.Dplus = gens_Dplus


  Dplus_perp, _ = orthogonal_submodule(DSY, Dplus)
  Dplus_perp,_ = snf(Dplus_perp)

  # needed for hash
  ECtx.Dplus_perp = reduce(vcat,[change_base_ring(GF(2),2*solve_left(basis_matrix(SY),matrix(QQ,1,26,lift(i)))) for i in gens(Dplus_perp)])


  @vprint :K3Auto 2 "done\n"
  ECtx.membership_test = membership_test_as_group
  #upper_bound = mass(ECtx)*length(initial_automorphisms(ECtx))
  return ECtx

end

@doc raw"""
Let $V_0$ be a complete set of representatives of the $L_{26}|S_Y$ chambers
making up $Nef(Y)/aut(Y)$. Then the mass is defined as
$\sum_{D \in V_0} \# Aut_G(D)$.
It can be computed as $volume_index(D0)/\#\bar G_{X-}$.
"""
function mass(ECtx::EnriquesBorcherdsCtx)
  return volindex(ECtx)[1]//ECtx.orderGbar
end

function membership_test_as_group(data::EnriquesBorcherdsCtx, f::FqMatrix)
  # test if f lies in Gplus
  # f must
  # act as identity on Dplus^perp
  b =  f in data.Gplus_mat
  return b
end

function initial_rays(ECtx::EnriquesBorcherdsCtx)
  if !isdefined(ECtx, :initial_rays)
    ECtx.initial_rays = rays(ECtx.initial_chamber)
  end
  return ECtx.initial_rays
end

function initial_isotropic_rays(ECtx::EnriquesBorcherdsCtx)
  if !isdefined(ECtx, :initial_isotropic_rays)
    iso_rays = filter(r->0==r*ECtx.gramSY*transpose(r), rays(ECtx.initial_chamber))
    ECtx.initial_isotropic_rays = iso_rays
  end
  return ECtx.initial_isotropic_rays
end


function _assure_membership_test_as_set(ECtx::EnriquesBorcherdsCtx)
  if isdefined(ECtx,:imgs_mod2)
    return
  end
  @vprint :K3Auto 2 "computing $(order(ECtx.Gplus_mat)) images mod 2\n"
  imgs_mod2 = Set{FqMatrix}()
  for g in ECtx.Gplus_mat
    h = matrix(g)
    push!(imgs_mod2, ECtx.Dplus*h)
  end
  ECtx.imgs_mod2 = imgs_mod2
end

function membership_test_set(data::EnriquesBorcherdsCtx, f::FqMatrix)
  # test if f lies in Gplus
  # f must
  # act as identity on Dplus^perp
  data.Dplus_perp*f == data.Dplus_perp || return false
  # its restriction to Dplus
  # must lie in in Gplus_res
  f_of_Dplus = data.Dplus*f
  return f_of_Dplus in data.imgs_mod2
end


gens(L::ZZLat) = [vec(basis_matrix(L)[i,:]) for i in 1:rank(L)]
basis(L::ZZLat) = gens(L)
# in principle it would be enough to just store
# tau and parent wall
# then tau can be recomputed via the spanning tree
# that would reduce storage by a factor of about 10
# ... but increase computation time
struct EnriquesChamber  # no need to make this mutable.
  data::EnriquesBorcherdsCtx
  tau::ZZMatrix
  parent_wall::ZZMatrix # for the spanning tree, needed?
end

function ==(x::EnriquesChamber, y::EnriquesChamber)
  x.data === y.data || error("bad comparison")
  return x.tau == y.tau
end

function Base.hash(D::EnriquesChamber, h::UInt)
  return hash(fingerprint(D), h)
end

function Base.in(x::ZZMatrix, D::EnriquesChamber)
  G = D.data.gramSY
  Gx = G*transpose(x)
  return all((v*Gx)[1,1]>=0 for v in walls(D))
end

# inv(g)*h in G => fingerprint(g) == fingerprint(h)
function fingerprint(D::EnriquesChamber)
  #return 0
  # the computation of this hash will be a bit expensive
  # t.f.a.e.
  # D0*g G-cong D0*g
  # g^-1 f h in G for some f in aut(D0)
  # B g^-1 f h cong B mod 2 for some f in aut(D0)
  # B g^-1 f cong B h^-1 for some f in aut(D0)
  # B g^-1 aut(D0) = Bh^-1 aut(D0)
  # We can hash this set, or something canonically computed from it
  B = D.data.Dplus_perp
  T = B*inv(D.tau)
  function _lt(x::FqMatrix,y::FqMatrix)
    n = nrows(x)
    m = ncols(x)
    for i in 1:n
      for j in 1:m
        x0 = iszero(x[i,j])
        y0 = iszero(y[i,j])
        if x0 && !y0
          return true
        elseif !x0 && y0
          return false
        end
      end
    end
    # now x == y
    return false
  end
  S = sort!([T*f for f in initial_automorphisms(D.data)],lt=_lt)
  return S
  # this hash will be cheaper if we precompute the sum
  S = sum(change_base_ring(GF(2),i) for i in initial_automorphisms(D.data))
  return B*inv(D.tau)*S
end


function initial_automorphisms(Y::EnriquesBorcherdsCtx)
  if !isdefined(Y, :initial_automorphisms)
    @vprint :K3Auto 2 "computing automorphisms"
    Y.initial_automorphisms = aut(Y.initial_chamber)
    @vprintln :K3Auto 2 " done found $(length(Y.initial_automorphisms)) automorphisms"
  end
  return Y.initial_automorphisms
end

function initial_automorphisms_mod2(Y::EnriquesBorcherdsCtx)
  if !isdefined(Y, :initial_automorphisms_mod2)
    Y.initial_automorphisms_mod2 = [change_base_ring(GF(2), i) for i in initial_automorphisms(Y)]
  end
  return Y.initial_automorphisms_mod2
end

function initial_walls(Y::EnriquesBorcherdsCtx)
  if !isdefined(Y, :initial_walls)
    @vprint :K3Auto 2 "computing walls\n"
    Y.initial_walls = walls(Y.initial_chamber)
  end
  return Y.initial_walls
end

function Base.show(io::IO, ::MIME"text/plain", D::EnriquesChamber)
  println(io, "EnriquesChamber")
  show(io, MIME("text/plain"), change_base_ring(GF(2), D.tau))
end

function Base.show(io::IO, D::EnriquesChamber)
  print(io, "EnriquesChamber")
end

function reflection(gram::ZZMatrix, v::ZZMatrix)
  n = ncols(gram)
  # E = identity_matrix(base_ring(gram), n)
  gramv = gram * transpose(v)
  c =  (v * gramv)[1,1]
  ref = -2 * gramv * v
  divexact!(ref,ref, c)
  for k in 1:n
    ref[k,k] = ref[k,k] + 1
  end
  return ref
end

@doc raw"""
    adjacent_chamber(D::K3Chamber, v::ZZMatrix) -> K3Chamber

Return return the chamber adjacent to `D` via the wall defined by `v`.
"""
function adjacent_chamber(D::EnriquesChamber, v::ZZMatrix)
  gramS = D.data.gramSY
  ref = reflection(gramS, v)
  tau2 = D.tau*ref
  Dnew = EnriquesChamber(D.data, tau2, v)
  return Dnew
end

function rays(D::EnriquesChamber)
  r = initial_rays(D.data)
  return [i*D.tau for i in r]
end

function isotropic_rays(D::EnriquesChamber)
  r = initial_isotropic_rays(D.data)
  return [i*D.tau for i in r]
end

function walls(D::EnriquesChamber)
  gramS = D.data.gramSY
  walls0 = initial_walls(D.data)
  return [r*D.tau for r in walls0]
end

function is_inner_wall(D::EnriquesChamber, v::ZZMatrix)
  v2 = change_base_ring(GF(2), v)
  return v2 in D.data.roots_mod2
end

function hom(D1::EnriquesChamber, D2::EnriquesChamber)
  result = ZZMatrix[]
  tau1inv = inv(D1.tau)
  n = 0
  #=
  for g in initial_automorphisms(D1.data)
    n = n+1
    if n==10
      break
    end
    @show n
    if D1.data.membership_test(g)[1]
      h = tau1inv*g*D2.tau
      push!(result, h)
    end
  end
  =#
  t1 = change_base_ring(GF(2),tau1inv)
  t2 = change_base_ring(GF(2), D2.tau)
  for (i,g2) in enumerate(initial_automorphisms_mod2(D1.data))
    h2 = t1*g2*t2
    if D1.data.membership_test(D1.data, h2)
      g = initial_automorphisms(D1.data)[i]
      h = tau1inv*g*D2.tau
      push!(result, h)
    end
  end
  return result
end

function hom_first(D1::EnriquesChamber, D2::EnriquesChamber)
  result = ZZMatrix[]
  tau1inv = inv(D1.tau)
  t1 = change_base_ring(GF(2),tau1inv)
  t2 = change_base_ring(GF(2), D2.tau)
  for (i,g2) in enumerate(initial_automorphisms_mod2(D1.data))
    h2 = t1*g2*t2
    if D1.data.membership_test(D1.data, h2)
      g = initial_automorphisms(D1.data)[i]
      h = tau1inv*g*D2.tau
      push!(result, h)
      break
    end
  end
  return result
end

aut(D::EnriquesChamber) = hom(D, D)

function is_outer_wall(D::EnriquesChamber, v::ZZMatrix)
  @assert v*D.data.gramSY*transpose(v) == -4
  return change_base_ring(GF(2),v) in D.data.roots_mod2
end

function borcherds_method(data::EnriquesBorcherdsCtx; max_nchambers=-1)
  S = data.SY
  # for G-sets
  F = FreeModule(ZZ,rank(S), cached=false)
  # initialization
  n = rank(S)
  D = EnriquesChamber(data, identity_matrix(ZZ, n), zero_matrix(ZZ, 1, n))
  waiting_list = EnriquesChamber[D]

  chambers = Dict{UInt64,Vector{EnriquesChamber}}()
  chambers[hash(fingerprint(D))] = EnriquesChamber[D]
  waiting_list = [D]

  automorphisms = Set{ZZMatrix}()
  rational_curves = Set{ZZMatrix}()

  # apparently computing the small generating set is very slow
  autD = aut(D)
  autD_grp = matrix_group(autD)
  @vprint :K3Auto 4 "computing small generating set "
  autD_mod2 = matrix_group([change_base_ring(GF(2),i) for i in autD])
  iso = hom(autD_grp, autD_mod2, gens(autD_mod2),check=false)
  autD = [matrix(preimage(iso,i)) for i in small_generating_set(autD_mod2)]
  if order(autD_mod2) != length(autD)
    K,i = kernel(iso)
    append!(autD, matrix.(gens(K)))
  end
  @vprintln :K3Auto 4 "done"
  # the following was too slow
  #order(autD_grp)  # somehow computing the order makes it faster
  #autD = [matrix(g) for g in small_generating_set(autD_grp)]
  for f in autD
    push!(automorphisms, f)
  end

  massY = mass(data)
  mass_explored = QQ(0)
  ntry = 0
  nchambers = 1
  while length(waiting_list) > 0
    ntry = ntry + 1
    if mod(ntry, 5)==0
      @vprint :K3Auto 2 "largest bucket: $(maximum(length(i) for i in values(chambers))) "
      @vprint :K3Auto 1 "buckets: $(length(chambers)) explored: $(nchambers) unexplored: $(length(waiting_list)) gens: $(length(automorphisms)); rat. curv. $(length(rational_curves))\n"
      @vprint :K3Auto 1 "mass left $(massY - mass_explored)"
    end
    D = popfirst!(waiting_list)

    autD = aut(D)
    mass_explored = mass_explored + inv(QQ(length(autD)))
    # we need the orbits of the walls only
    if length(autD) > 1
      # compute the orbits
      @vprint :K3Auto 3 "computing orbits"
      Omega = [F(v) for v in walls(D)]
      W = gset(matrix_group(autD),Omega)
      vv = F(D.parent_wall)
      wallsDmodAutD = [representative(w).v for w in orbits(W) if !(vv in w)]
      @vprintln :K3Auto 3 " done"
    else
      # the minus shouldn't be necessary ... but who knows?
      wallsDmodAutD = (v for v in walls(D) if !(v==D.parent_wall || -v==D.parent_wall))
    end
    # compute the adjacent chambers to be explored
    for v in wallsDmodAutD
      if is_outer_wall(D, v)
        # v comes from a rational curve
        push!(rational_curves, v)
        continue
      end
      Dv = adjacent_chamber(D, v)
      # check G-congruence
      @vprint :K3Auto 5 "checking G-congruence "
      fp = hash(fingerprint(Dv))
      if !haskey(chambers, fp)
        chambers[fp] = EnriquesChamber[]
      end
      is_explored = false
      for E in chambers[fp]
        gg = hom_first(Dv, E)
        if length(gg) > 0
          # enough to add a single homomorphism
          push!(automorphisms, gg[1])
          is_explored = true
          break
        end
      end
      @vprintln :K3Auto 5 "done"
      if is_explored
        continue
      end
      push!(chambers[fp], Dv)
      push!(waiting_list, Dv)
      nchambers = nchambers+1
    end
    if max_nchambers != -1 && ntry > max_nchambers
      return data, collect(automorphisms), reduce(append!,values(chambers), init=EnriquesChamber[]), collect(rational_curves), false
    end
  end
  @vprint :K3Auto 1 "$(length(automorphisms)) automorphism group generators\n"
  @vprint :K3Auto 1 "$(nchambers) congruence classes of chambers \n"
  @vprint :K3Auto 1 "$(length(rational_curves)) orbits of rational curves\n"
  @vprint :K3Auto 1 "$mass $(massY)\n"
  @vprint :K3Auto 1 "$mass explored $(mass_explored)\n"
  @assert massY == mass_explored
  return data, collect(automorphisms), reduce(append!,values(chambers), init=EnriquesChamber[]), collect(rational_curves), true
end


function frame_lattice(L, v)
  @assert v*gram_matrix(ambient_space(L))*transpose(v)==0
  F = orthogonal_submodule(L,lattice(ambient_space(L),v))
  vF = change_base_ring(ZZ,solve_left(basis_matrix(F), v))
  @assert gcd(vec(vF))==1
  b = Hecke._complete_to_basis(vF)
  return lll(lattice(ambient_space(F),b[1:end-1,:]*basis_matrix(F)))
end


function tau_taubar_general(R::ZZLat)

end

function volindex(enr::EnriquesBorcherdsCtx)
  Sm = orthogonal_submodule(enr.L26, enr.SY)
  rt = root_lattice_recognition(Sm)[1]
  sort!(rt)
  d = _emb_types()
  j = findfirst(i->rt == i[3], d)
  return d[j]
end

function _emb_types()
  d = [
  [2^12 * 3^5 * 5^2 * 7, "12A", [(:D, 8)]],
  [2^12 * 3^3 * 5 * 7,"12B", [(:A, 7)]],
  [2^8 * 3^4 * 5 * 7,"20A", [(:D, 4), (:D, 5)]],
  [2^10 * 3^2 * 5 * 7, "20B", [(:D, 4), (:D, 4)]],
  [2^6 * 3^3 * 5 * 7, "20C", vcat([(:A, 1) for i in 1:10], [(:D, 6)])],
  [2^6 * 3^3 * 5 * 7, "20D", [(:A, 3), (:A, 4)]],
  [2^7 * 3^4 * 5, "20E", vcat([(:A, 1) for i in 1:5], [(:A, 5)])],
  [2^9 * 3^2 * 5, "20F", [(:A, 3), (:A, 3)]],
  [2^7 * 3^2 * 5, "40A", vcat([(:A, 1) for i in 1:4], [(:A, 3),(:A, 3)])],
  [2^3 * 3^2 * 5 * 7, "40B", vcat([(:A, 1) for i in 1:8], [(:D, 4),(:D, 4)])],
  [2^3 * 3^2 * 5 * 7, "40C", vcat([(:A, 1) for i in 1:6], [(:A, 3)])],
  [2^5 * 3^2 * 5, "40D", vcat([(:A, 1) for i in 1:12], [(:D, 4)])],
  [2^5 * 3^2 * 5, "40E", vcat([(:A, 1) for i in 1:2], [(:A, 2),(:A, 2)])],
  [2^5 * 3^2,"96A", [(:A, 1) for i in 1:8]],
  [2^3 * 3^2,"96B", [(:A, 1) for i in 1:16]],
  [2^3 * 3^2,"96C", [(:A, 1) for i in 1:4]]]
  for i in d
    sort!(i[3])
  end
  return d
end

function can_extend_with_extension(S, L, fS)
  R = orthogonal_submodule(L, S)
  @req is_definite(R) "R must be definite"
  OR = orthogonal_group(R)
  DR = discriminant_group(R)
  ODR,_= orthogonal_group(DR)
  phiR = hom(OR, ODR, ODR.(gens(OR)); check=false)
  GR = phiR(OR)[1]
  phi,incDS,incDR = glue_map(S,R,L)
  DS = domain(phi)
  @vprint :K3Auto 2 "computing stabilizer "
  stab_DR, inc_stabDR = stabilizer(GR, incDR)
  resDR, inc_resDR = restrict_automorphism_group(stab_DR, inc_DR)
  # to be continued
end

function compute_Gplus(SY,SX,L26)

  # SY + Sm < SX is a primitive extension with glue map phi: D(Sm) -> D(SY)
  Sm = orthogonal_submodule(SX, SY)
  phi, inc_Dminus, inc_Dplus = glue_map(SX, Sm, SY)

  # Cook up the membership test
  @vprint :K3Auto 2 "computing orthogonal group\n"
  OSm = orthogonal_group(Sm)
  GSm,_= image_in_Oq(Sm)
  phiSm = hom(OSm, GSm, GSm.(gens(OSm)); check=false)
  DSm = discriminant_group(Sm)
  Dminus = domain(phi)
  @vprint :K3Auto 2 "computing stabilizer "
  stab_Dminus, inc_stab = stabilizer(GSm, inc_Dminus)
  @vprint :K3Auto 2 "done\n"
  Dminus_perp, inc_Dminus_perp = orthogonal_submodule(DSm, Dminus)
  res, inc = restrict_automorphism_group(stab_Dminus, inc_Dminus_perp)
  res2, inc2 = restrict_automorphism_group(stab_Dminus, inc_Dminus)
  @vprint :K3Auto 2 "computing kernel "
  Gminus, inc_Gminus = kernel(inc)
  @vprintln :K3Auto 2 "done"

  # compute Gplus
  Q = orthogonal_submodule(L26,SY)
  DQ = discriminant_group(Q)
  glue_SY_Q,inc1,inc2 = glue_map(L26,SY, Q)
  DO = codomain(glue_SY_Q)

  @vprint :K3Auto 2 "preparing membership test\n"
  DSY = discriminant_group(SY)
  ODSY = orthogonal_group(DSY)
  tmp = [preimage(phiSm, i) for i in small_generating_set(Gminus)]
  tmp_hom = [hom(DO,DO, [DO(lift(j)*i) for j in gens(DO)]) for i in tmp]
  gens_Gplus = elem_type(ODSY)[ODSY(inv(inc1)*glue_SY_Q*i*inv(glue_SY_Q)*inc1) for i in tmp_hom]
  Gplus,_ = sub(ODSY, gens_Gplus)
  @assert order(Gplus) == order(Gminus)
  # iso_Gplus_minus = hom(Gplus, Gminus, gens(Gplus), small_generating_set(Gminus))
  return Gplus
end

function ellfib_number(SY, SX, L26)
  Gplus = compute_Gplus(SY,SX,L26)
  DSY = discriminant_group(SY)
  isoY = [i for i in DSY if quadratic_product(i)==0 && !is_zero(i)]
  X = gset(Gplus, (x,g)->g(x),isoY)
  orbX = orbits(X)
  return length(orbX), [representative(i) for i in orbX]
end

function fibration_types(fbar::TorQuadModuleElem, SY, SX, L26)
  DY = discriminant_group(SY)
  Sm = orthogonal_submodule(SX, SY)
  phi, inc_Dplus, _ = glue_map(SX, SY, Sm)
  Dplus = domain(phi)
  sv2 = [1//2*(i[1]*basis_matrix(Sm)) for i in short_vectors(rescale(Sm,-1),4) if i[2]==4]
  sv2 = [codomain(phi)(i) for i in sv2 if i in cover(codomain(phi))]
  DeltabarY = Set([inc_Dplus(inv(phi)(i)) for i in sv2])
  ebar = 0*DY[1]
  for x in DY
    if inner_product(x,fbar)!=0
      ebar = x
      break
    end
  end
  @assert !iszero(ebar)


  Delta_fbar = [x for x in DeltabarY if inner_product(fbar,x) ==0]
  Delta1 = [ x for x in Delta_fbar if !(x+fbar in Delta_fbar)]
  Delta2 = [x for x in Delta_fbar if x+fbar in Delta_fbar && inner_product(x,ebar)==0]
  return _identify_fibers(Delta1,fbar),_identify_fibers(Delta2,fbar)
end

function _identify_fibers(D::Vector{TorQuadModElem},fbar)
  n = length(D)
  g = zero_matrix(GF(2),n,n)
  for i in 1:n
    for j in i+1:n
      if !iszero(inner_product(D[i],D[j]))
        g[i,j] = 1
        g[j,i] = 1
      end
    end
  end
  return _identify_ADE(g)
end

function _identify_ADE(g::MatElem)
  G = graph_from_adjacency_matrix(Undirected,g)
  # We identify an ADE lattice by
  # ((number of roots) // 2, rank of gram mod 2)
  An = [(1, 0), (3, 2), (6, 2), (10, 4), (15, 4), (21, 6), (28, 6), (36, 8), (45, 8), (55, 10)]
  Dn = [(12, 2), (20, 4), (30, 4), (42, 6), (56, 6), (72, 8), (90, 8)]
  En = [(36, 6), (63, 6), (120, 8)]
  fiber_types = Tuple{Symbol,Int}[]
  for c in connected_components(G)
    l = (length(c),rank(g[c,c]))
    if l in An
      i = findfirst(==(l),An)
      t = (:A,i)
    elseif l in Dn
      i = 3+findfirst(==(l),Dn)
      t = (:D,i)
    elseif l in En
      i = 5+findfirst(==(l),En)
      t = (:E,i)
    end
    push!(fiber_types,t)
  end
  return fiber_types
end

function contraction_type(hbar, enr::EnriquesBorcherdsCtx)
  phi = enr.roots_mod2
  k = base_ring(hbar)
  SY = change_base_ring(k, 1//2*gram_matrix(enr.SY))
  SYh = SY*transpose(hbar)
  phi_h = reduce(vcat,[r for r in phi if r*SYh ==0], init=zero_matrix(k,0,10))
  g = phi_h*SY*transpose(phi_h)
  return _identify_ADE(g)
end 


