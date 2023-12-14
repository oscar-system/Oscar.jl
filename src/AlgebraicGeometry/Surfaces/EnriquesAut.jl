
################################################################################
#
# Computation of Automorphism groups of Enriques surfaces following
# https://doi.org/10.1307/mmj/20195769
# and
# http://www.arxiv.org/abs/1909.10813
#
################################################################################
mutable struct EnriquesBorcherdsCtx
  # SY < SX < L26
  L26::ZZLat
  SX::ZZLat
  SY::ZZLat
  # the following are given with respect to the basis of SY
  initial_walls::Vector{ZZMatrix}
  # automorphisms of the initial chamber
  # written w.r.t. the basis of SY
  initial_automorphisms::Vector{ZZMatrix}
  initial_automorphisms_mod2::Vector{fpMatrix}  # same order as above but mod 2
  initial_chamber::K3Chamber   # a non-degenerate L26/SX chamber
  # whether a -2 vector defines an outer wall is a mod 2 condition
  roots_mod2::Set{fpMatrix}  # mod 2 images of the rational curves
  gramSY::ZZMatrix  # to avoid coercions
  gramSX::ZZMatrix  # to avoid coercions
  membership_test
  # needed for the membership test
  # the following are in SY/2SY w.r.t the basis of SY mod 2
  Dplus_perp::fpMatrix
  Dplus::fpMatrix
  imgs_mod2::Set{fpMatrix}
  # w.r.t the basis given by Dplus
  Gplus_mat::MatrixGroup{fpFieldElem, fpMatrix}
  function EnriquesBorcherdsCtx()
    return new()
  end

  volume_index::ZZRingElem
  orderGbar::ZZRingElem
end

function Base.show(io::IO, dat::EnriquesBorcherdsCtx)
  print(io, "Enriques Borcherds context with det(SX) = $(det(dat.SX)).")
end


function EnriquesBorcherdsCtx(SY::ZLat, SX::ZLat, L26::ZLat, weyl::ZZMatrix)
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
  @vprint :K3Auto 2 "computing walls\n"
  ECtx.initial_walls = walls(ECtx.initial_chamber)
  @vprint :K3Auto 2 "computing automorphisms\n"
  ECtx.initial_automorphisms = aut(ECtx.initial_chamber)
  ECtx.initial_automorphisms_mod2 = [change_base_ring(GF(2), i) for i in ECtx.initial_automorphisms]

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
  stab_Dminus, inc_stab = stabilizer(GSm, Dminus, Oscar._on_subgroups)
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
  gens_Gplus = [ODSY(inv(inc1)*glue_SY_Q*i*inv(glue_SY_Q)*inc1) for i in tmp_hom]
  Gplus,_ = sub(ODSY, gens_Gplus)
  @assert order(Gplus) == order(Gminus)
  # iso_Gplus_minus = hom(Gplus, Gminus, gens(Gplus), small_generating_set(Gminus))

  Dplus = codomain(phi)
  ODplus = orthogonal_group(Dplus)
  # turn this into a mod 2 matrix group
  # with respect to the basis of SY mod 2
  B = 1//2*basis_matrix(SY)
  gens_Gplus_mat = fpMatrix[]
  for g in gens(Gplus)
    g2 = reduce(vcat,[change_base_ring(GF(2),solve_left(B,matrix(QQ,1,26, lift(g(DSY(vec(B[i,:])))))))  for i in 1:10])
    push!(gens_Gplus_mat, g2)
  end
  Gplus_mat = matrix_group(gens_Gplus_mat)
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
  upper_bound = mass(ECtx)*length(ECtx.initial_automorphisms)
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

function membership_test_as_group(data::EnriquesBorcherdsCtx, f::fpMatrix)
  # test if f lies in Gplus
  # f must
  # act as identity on Dplus^perp
  b =  f in data.Gplus_mat
  return b
end

function _assure_membership_test_as_set(ECtx::EnriquesBorcherdsCtx)
  if isdefined(ECtx,:imgs_mod2)
    return
  end
  @vprint :K3Auto 2 "computing $(order(ECtx.Gplus_mat)) images mod 2\n"
  imgs_mod2 = Set{fpMatrix}()
  for g in ECtx.Gplus_mat
    h = matrix(g)
    push!(imgs_mod2, h*ECtx.Dplus)
  end
  ECtx.imgs_mod2 = imgs_mod2
end

function membership_test_set(data::EnriquesBorcherdsCtx, f::fpMatrix)
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
  return hash(0, h)
  # the computation of this hash will be a bit expensive
  # t.f.a.e.
  # D0*g G-cong D0*g
  # g^-1 f h in G for some f in aut(D0)
  # B g^-1 f h cong B mod 2 for some f in aut(D0)
  # B g^-1 f cong B h^-1 for some f in aut(D0)
  # B g^-1 aut(D0) = Bh^-1 aut(D0)
  # We can hash this set, or something canonically computed from it
  B = D.data.Dplus_perp
  S = Set([B*inv(D.tau)*f for f in D.data.initial_automorphisms])
  return hash(S, h)
  # this hash will be cheaper if we precompute the sum
  S = sum(change_base_ring(GF(2),i) for i in D.data.initial_automorphisms)
  return hash(B*inv(D.tau)*S, h)
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

# inv(g)*h in G => fingerprint(g) == fingerprint(h)
# So far no good idea?
function fingerprint(D::EnriquesChamber)
  return 0
  return hnf(D.data.Dplus_perp*change_base_ring(GF(2),D.tau))
end

function walls(D::EnriquesChamber)
  gramS = D.data.gramSY
  walls0 = D.data.initial_walls
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
  for g in D1.data.initial_automorphisms
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
  for (i,g2) in enumerate(D1.data.initial_automorphisms_mod2)
    h2 = t1*g2*t2
    if D1.data.membership_test(D1.data, h2)
      g = D1.data.initial_automorphisms[i]
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
  for (i,g2) in enumerate(D1.data.initial_automorphisms_mod2)
    h2 = t1*g2*t2
    if D1.data.membership_test(D1.data, h2)
      g = D1.data.initial_automorphisms[i]
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
  waiting_list = [D]

  automorphisms = Set{ZZMatrix}()
  rational_curves = Set{ZZMatrix}()

  massY = mass(data)
  mass_explored = QQ(0)
  ntry = 0
  nchambers = 0
  while length(waiting_list) > 0
    ntry = ntry + 1
    if mod(ntry, 5)==0
      @vprint :K3Auto 2 "largest bucket: $(maximum(length(i) for i in values(chambers))) "
      @vprint :K3Auto 1 "buckets: $(length(chambers)) explored: $(nchambers) unexplored: $(length(waiting_list)) gens: $(length(automorphisms)); rat. curv. $(length(rational_curves))\n"
      @vprint :K3Auto 1 "mass left $(massY - mass_explored)"
    end
    D = popfirst!(waiting_list)
    # check G-congruence
    fp = hash(fingerprint(D))
    if !haskey(chambers, fp)
      chambers[fp] = EnriquesChamber[]
    end
    is_explored = false
    for E in chambers[fp]
      gg = hom_first(D, E)
      if length(gg) > 0
        # enough to add a single homomorphism
        push!(automorphisms, gg[1])
        is_explored = true
        break
      end
    end
    if is_explored
      continue
    end
    push!(chambers[fp], D)
    nchambers = nchambers+1

    autD = aut(D)
    mass_explored = mass_explored + inv(QQ(length(autD)))
    # we need the orbits of the walls only
    if length(autD) > 1
      # apparently computing the small generating set is very slow
      autD_grp = matrix_group(autD)
      @vprint :K3Auto 4 "computing small generating set "
      autD_mod2 = matrix_group([change_base_ring(GF(2),i) for i in autD])
      iso = hom(autD_grp, autD_mod2, gens(autD_mod2),check=false)
      autD = [matrix(preimage(iso,i)) for i in small_generating_set(autD_mod2)]
      @vprintln :K3Auto 4 "done"
      # the following was too slow
      #order(autD_grp)  # somehow computing the order makes it faster
      #autD = [matrix(g) for g in small_generating_set(autD_grp)]
      for f in autD
        push!(automorphisms, f)
      end
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
  @assert massY == mass_explored
  return data, collect(automorphisms), reduce(append!,values(chambers), init=EnriquesChamber[]), collect(rational_curves), true
end


function frame_lattice(L, v)
  @assert v*gram_matrix(ambient_space(L))*transpose(v)==0
  F = orthogonal_submodule(L,lattice(ambient_space(L),v))
  vF = solve_left(change_base_ring(ZZ,basis_matrix(F)), v)
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

