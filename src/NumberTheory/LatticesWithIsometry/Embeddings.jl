export admissible_equivariant_primitive_extensions,
       equivariant_primitive_extensions,
       primitive_embeddings_in_primary_lattice,
       primitive_embeddings_of_elementary_lattice

GG = GAP.Globals

################################################################################
#
# Miscellaneous
#
################################################################################

function inner_orthogonal_sum(T::TorQuadMod, U::TorQuadMod)
  @assert modulus_bilinear_form(T) == modulus_bilinear_form(U)
  @assert modulus_quadratic_form(T) == modulus_quadratic_form(U)
  @assert ambient_space(cover(T)) === ambient_space(cover(U))
  cS = cover(T)+cover(U)
  rS = relations(T)+relations(U)
  geneT = [lift(a) for a in gens(T)]
  geneU = [lift(a) for a in gens(U)]
  S = torsion_quadratic_module(cS, rS, gens = unique([g for g in vcat(geneT, geneU) if !is_zero(g)]), modulus = modulus_bilinear_form(T), modulus_qf = modulus_quadratic_form(T))
  TinS = hom(T, S, S.(geneT))
  UinS = hom(U, S, S.(geneU))
  if order(S) == 1
    S.gram_matrix_quadratic = matrix(QQ,0,0,[])
    S.gram_matrix_bilinear = matrix(QQ,0,0,[])
    set_attribute!(S, :is_degenerate, false)
    S.ab_grp = abelian_group()
  end
  return S, TinS, UinS
end

function embedding_orthogonal_group(i::TorQuadModMor, j::TorQuadModMor)
  A = domain(i)
  B = domain(j)
  D = codomain(i)
  Dorth = direct_sum(A, B)[1]
  ok, phi = is_isometric_with_isometry(Dorth, D)
  @assert ok
  OD = orthogonal_group(D)
  OA = orthogonal_group(A)
  OB = orthogonal_group(B)

  geneOAinDorth = TorQuadModMor[]
  for f in gens(OA)
    m = block_diagonal_matrix([matrix(f), identity_matrix(ZZ, ngens(B))])
    m = hom(Dorth, Dorth, m)
    push!(geneOAinDorth, m)
  end

  geneOBinDorth = TorQuadModMor[]
  for f in gens(OB)
    m = block_diagonal_matrix([identity_matrix(ZZ, ngens(A)), matrix(f)])
    m = hom(Dorth, Dorth, m)
    push!(geneOBinDorth, m)
  end
  geneOAinOD = [OD(compose(inv(phi), compose(g, phi)), check = false) for g in geneOAinDorth]
  OAtoOD = hom(OA, OD, geneOAinOD, check = false)
  geneOBinOD = [OD(compose(inv(phi), compose(g, phi)), check = false) for g in geneOBinDorth]
  OBtoOD = hom(OB, OD, geneOBinOD, check = false)
  return OAtoOD, OBtoOD
end

function _embedding_orthogonal_group(i::TorQuadModMor)
  @req is_injective(i) "i must be injective"
  A = domain(i)
  D = codomain(i)
  OA = orthogonal_group(A)
  OD = orthogonal_group(D)
  if order(OA) == 1
    return hom(OA, OD, [one(OD)], check = false)
  end
  B = orthogonal_submodule(D, A)[1]
  Dorth = inner_orthogonal_sum(A, B)[1]
  ok, phi = is_isometric_with_isometry(Dorth, D)
  @assert ok

  geneOAinDorth = TorQuadModMor[]
  for f in gens(OA)
    m = block_diagonal_matrix([matrix(f), identity_matrix(ZZ, ngens(B))])
    m = hom(Dorth, Dorth, m)
    push!(geneOAinDorth, m)
  end
  geneOAinOD = [OD(compose(inv(phi), compose(g, phi)), check = false) for g in geneOAinDorth]
  OAtoOD = hom(OA, OD, geneOAinOD, check = false)
  return OAtoOD::GAPGroupHomomorphism{AutomorphismGroup{TorQuadMod}, AutomorphismGroup{TorQuadMod}}
end

function _restrict(f::TorQuadModMor, i::TorQuadModMor)
  imgs = TorQuadModElem[]
  V = domain(i)
  for g in gens(V)
    h = f(i(g))
    hV = preimage(i, h)
    push!(imgs, hV)
  end
  return hom(V, V, imgs)
end

function _restrict(f::AutomorphismGroupElem{TorQuadMod}, i::TorQuadModMor)
  return _restrict(hom(f), i)
end

function _is_invariant(f::TorQuadModMor, i::TorQuadModMor)
  fab = f.map_ab
  V = domain(i)
  for a in gens(V)
    b = f(i(a))
    haspreimage(fab, data(b))[1] || return false
  end
  return true
end

function _is_invariant(f::AutomorphismGroupElem{TorQuadMod}, i::TorQuadModMor)
  return _is_invariant(hom(f), i)
end

function _is_invariant(aut::AutomorphismGroup{TorQuadMod}, i::TorQuadModMor)
  return all(g -> _is_invariant(g, i), gens(aut))
end

function _get_V(L, q, f, fq, mu, p)
  f_mu = mu(f)
  if !is_zero(f_mu)
    L_sub = intersect(lattice_in_same_ambient_space(L, inv(f_mu)*basis_matrix(L)), dual(L))
    B = basis_matrix(L_sub)
    V, Vinq = sub(q, q.([vec(collect(B[i,:])) for i in 1:nrows(B)]))
  else 
    V, Vinq = q, id_hom(q)
  end  
  pV, pVinV = primary_part(V, p)
  pV, pVinV = sub(V, pVinV.([divexact(order(g), p)*g for g in gens(pV) if !(order(g)==1)]))
  pVinq = compose(pVinV, Vinq)
  fpV = _restrict(fq, pVinq)
  return pV, pVinq, fpV
end

function _rho_functor(q::TorQuadMod, p, l::Union{Integer, fmpz})
  N = relations(q)
  if l == 0 
    Gl = N
    Gm = intersect(1//p*N, dual(N))
    rholN = torsion_quadratic_module(Gl, p*Gm, modulus = QQ(1), modulus_qf = QQ(2))
    gene = [q(lift(g)) for g in gens(rholN)]
    return sub(q, gene)
  end  
  k = l-1
  m = l+1
  Gk = intersect(1//p^k*N, dual(N))
  Gl = intersect(1//p^l*N, dual(N))
  Gm = intersect(1//p^m*N, dual(N))
  rholN = torsion_quadratic_module(Gl, Gk+p*Gm, modulus = QQ(1), modulus_qf = QQ(2))
  gene = [q(lift(g)) for g in gens(rholN)]
  return sub(q, gene)
end

##############################################################################
#
#  Orbits and stabilizers of discriminant subgroups
#
##############################################################################

function on_subgroup(H::Oscar.GAPGroup, g::AutomorphismGroupElem)
  G = domain(parent(g))
  return sub(G, g.(gens(H)))[1]
end

function stabilizer(O::AutomorphismGroup{TorQuadMod}, i::TorQuadModMor)
  @req domain(O) === codomain(i) "Incompatible arguments"
  q = domain(O)
  togap = get_attribute(O, :to_gap)
  tooscar = get_attribute(O, :to_oscar)
  A = codomain(togap)
  OA = automorphism_group(A)
  OinOA, _ = sub(OA, OA.([g.X for g in gens(O)]))
  N, _ = sub(A, togap.(i.(gens(domain(i)))))
  stab, _ = stabilizer(OinOA, N, on_subgroup)
  return sub(O, O.([h.X for h in gens(stab)]))
end

function _as_Fp_vector_space_quotient(HinV, p, f)
  if f isa AutomorphismGroupElem
    f = hom(f)
  end
  i = HinV.map_ab
  H = domain(HinV)
  Hab = domain(i)
  Hs, HstoHab = snf(Hab)
  V = codomain(HinV)
  Vab = codomain(i)
  Vs, VstoVab = snf(Vab)

  function _VtoVs(x::TorQuadModElem)
    return inv(VstoVab)(data(x))
  end

  function _VstoV(x::GrpAbFinGenElem)
    return V(VstoVab(x))
  end

  VtoVs = Hecke.MapFromFunc(_VtoVs, _VstoV, V, Vs)

  n = ngens(Vs)
  F = GF(p)
  MVs = matrix(compose(VstoVab, compose(f.map_ab, inv(VstoVab))))
  Vp = VectorSpace(F, n)

  function _VstoVp(x::GrpAbFinGenElem)
    v = x.coeff
    return Vp(vec(collect(v)))
  end

  function _VptoVs(v::ModuleElem{gfp_fmpz_elem})
    x = lift.(v.v)
    return Vs(vec(collect(x)))
  end

  VstoVp = Hecke.MapFromFunc(_VstoVp, _VptoVs, Vs, Vp)
  VtoVp = compose(VtoVs, VstoVp)
  gene = gfp_fmpz_mat[matrix(F.((i(HstoHab(a))).coeff)) for a in gens(Hs)]
  subgene = [Vp(vec(collect(transpose(v)))) for v in gene]
  Hp, _ = sub(Vp, subgene)
  Qp, VptoQp = quo(Vp, Hp)
  fVp = change_base_ring(F, MVs)
  ok, fQp = can_solve_with_solution(VptoQp.matrix, fVp*VptoQp.matrix)
  @assert ok


  return Qp, VtoVp, VptoQp, fQp
end

function subgroups_orbit_representatives_and_stabilizers_elementary(Vinq::TorQuadModMor,
                                                                    G::AutomorphismGroup{TorQuadMod},
                                                                    ord::Hecke.IntegerUnion,
                                                                    f::Union{TorQuadModMor, AutomorphismGroupElem{TorQuadMod}} = id_hom(codomain(Vinq)),
                                                                    l::Hecke.IntegerUnion = 0)
  if ord == 1
    if l != 0
      return Tuple{TorQuadModMor, AutomorphismGroup{TorQuadMod}}[]
    end
    _, triv = sub(codomain(Vinq), TorQuadModElem[])
    return Tuple{TorQuadModMor, AutomorphismGroup{TorQuadMod}}[(triv, G)]
  end
  V = domain(Vinq)
  q = codomain(Vinq)
  p = elementary_divisors(V)[1]
  g = valuation(ord, p)
  @req all(a -> haspreimage(Vinq.map_ab, data(l*a))[1], gens(q)) "lq is not contained in V"
  H0, H0inV = sub(V, [preimage(Vinq, (l*a)) for a in gens(q)])
  Qp, VtoVp, VptoQp, fQp = _as_Fp_vector_space_quotient(H0inV, p, _restrict(f, Vinq))
  Vp = codomain(VtoVp)

  if dim(Qp) == 0
    if order(V) == ord
      return Tuple{TorQuadModMor, AutomorphismGroup{TorQuadMod}}[(Vinq, G)]
    end
    return Tuple{TorQuadModMor, AutomorphismGroup{TorQuadMod}}[]
  end

  OV = orthogonal_group(V)
  gene_GV = elem_type(OV)[OV(_restrict(g, Vinq), check = false) for g in gens(G)]
  GV, _ = sub(OV, gene_GV)
  GVinG = hom(GV, G, gens(GV), gens(G), check = false)

  act_GV_Vp = gfp_fmpz_mat[change_base_ring(base_ring(Qp), matrix(gg)) for gg in gens(GV)]
  act_GV_Qp = gfp_fmpz_mat[solve(VptoQp.matrix, g*VptoQp.matrix) for g in act_GV_Vp]
  MGp = matrix_group(act_GV_Qp)
  @assert fQp in MGp
  MGptoGV = hom(MGp, GV, gens(GV), check = false)
  MGptoG = compose(MGptoGV, GVinG)

  res = Tuple{TorQuadModMor, AutomorphismGroup{TorQuadMod}}[]

  if g-ngens(snf(abelian_group(H0))[1]) > dim(Qp)
    return res
  end

  F = base_ring(Qp)
  k, K = kernel(VptoQp.matrix, side = :left)
  gene_H0p = ModuleElem{gfp_elem}[Vp(vec(collect(K[i,:]))) for i in 1:k]
  orb_and_stab = orbit_representatives_and_stabilizers(MGp, g-k)
  for (orb, stab) in orb_and_stab
    i = orb.map
    _V = codomain(i)
    for v in gens(orb)
      vv = _V(i(v)*fQp)
      if !can_solve(i.matrix, vv.v, side = :left)
        @goto non_fixed
      end
    end
    gene_orb_in_Qp = AbstractAlgebra.Generic.QuotientModuleElem{gfp_fmpz_elem}[Qp(vec(collect(i(v).v))) for v in gens(orb)]
    gene_orb_in_Vp = AbstractAlgebra.Generic.ModuleElem{gfp_fmpz_elem}[preimage(VptoQp, v) for v in gene_orb_in_Qp]
    gene_submod_in_Vp = vcat(gene_orb_in_Vp, gene_H0p)
    gene_submod_in_V = TorQuadModElem[preimage(VtoVp, Vp(v)) for v in gene_submod_in_Vp]
    gene_submod_in_q = TorQuadModElem[image(Vinq, v) for v in gene_submod_in_V]
    orbq, orbqinq = sub(q, gene_submod_in_q)
    @assert order(orbq) == ord
    stabq, _ = image(MGptoG, stab)
    push!(res, (orbqinq, stabq))
    @label non_fixed
  end
  return res
end

function subgroups_orbit_representatives_and_stabilizers(O::AutomorphismGroup{TorQuadMod}; order::Hecke.IntegerUnion = -1)
  togap = get_attribute(O, :to_gap)
  tooscar = get_attribute(O, :to_oscar)
  q = domain(O)
  A = abelian_group(q)
  it = order == -1 ? subgroups(A) : subgroups(A, order=order)
  Agap = codomain(togap)
  coAgap = [sub(Agap, togap.(q.(j[2].(gens(j[1])))))[1] for j in it]
  OAgap = automorphism_group(Agap)
  OinOAgap, j = sub(OAgap, OAgap.([g.X for g in gens(O)]))
  m = gset(OinOAgap, on_subgroup, coAgap)
  orbs = orbits(m)
  res = Tuple{TorQuadModMor, AutomorphismGroup{TorQuadMod}}[]
  for orb in orbs
    rep = representative(orb)
    stab, _ = stabilizer(OinOAgap, rep, on_subgroup)
    _, rep = sub(q, TorQuadModElem[tooscar(Agap(g)) for g in gens(rep)])
    stab, _ = sub(O, O.([h.X for h in gens(stab)]))
    push!(res, (rep, stab))
  end 
  return res
end

function classes_conjugate_subgroups(O::AutomorphismGroup{TorQuadMod}, q::TorQuadMod)                             
  sors = subgroups_orbit_representatives_and_stabilizers(O, order=order(q))
  return filter(d -> is_isometric_with_isometry(domain(d[1]), q)[1], sors)
end

#################################################################################
#
# Embeddings in elementary lattices
#
#################################################################################

# here for convenience, we choose in entry N and M to be of full rank and
# with basis matrix equal to the identity matrix

function _isomorphism_classes_primitive_extensions(N::ZLat, M::ZLat, H::TorQuadMod)
  @assert is_one(basis_matrix(N))
  @assert is_one(basis_matrix(M))
  results = Tuple{ZLat, ZLat, ZLat}[]
  GN, _ = image_in_Oq(N)
  GM, _ = image_in_Oq(M)
  qN = domain(GN)
  qM = domain(GM)

  D, qNinD, qMinD = orthogonal_sum(qN, qM)
  OD = orthogonal_group(D)

  subsN = classes_conjugate_subgroups(GN, rescale(H, -1))
  @assert !isempty(subsN)
  subsM = classes_conjugate_subgroups(GM, H)
  @assert !isempty(subsM)

  for H1 in subsN, H2 in subsM
    ok, phi = is_anti_isometric_with_anti_isometry(domain(H1[1]), domain(H2[1]))
    @assert ok

    HNinqN, stabN = H1
    OHN = orthogonal_group(HN)

    HMinqM, stabM = H2
    OHM = orthogonal_group(HM)

    actN = hom(stabN, OHN, [OHN(_restrict(x, HNinqN)) for x in gens(stabN)])

    actM = hom(stabM, OHM, [OHM(_restrict(x, HMinqM)) for x in gens(stabM)])
    imM, _ = image(actM)

    stabNphi = AutomorphismGroupElem{TorQuadMod}[OHM(compose(inv(phi), compose(hom(actN(g)), phi))) for g in gens(stabN)]
    stabNphi, _ = sub(OHM, stabNphi)
    reps = double_cosets(OHM, stabNphi, imM)
    @info "$(length(reps)) isomorphism classe(s) of primitive extensions"
    for g in reps
      g = representative(g)
      phig = compose(phi, hom(g))
      _glue = Vector{fmpq}[lift(qNinD(HNinqN(g))) + lift(qMinD(HMinqM(phig(g)))) for g in gens(domain(phig))]
      z = zero_matrix(QQ, 0, degree(N)+degree(M))
      glue = reduce(vcat, [matrix(QQ, 1, degree(N)+degree(M), g) for g in _glue], init = z)
      glue = vcat(identity_matrix(QQ, rank(N)+rank(M)), glue)
      glue = FakeFmpqMat(glue)
      _B = hnf(glue)
      _B = QQ(1, denominator(glue))*change_base_ring(QQ, numerator(_B))
      L = lattice(ambient_space(cover(D)), _B[end-rank(N)-rank(M)+1:end, :])
      N2 = lattice_in_same_ambient_space(L, identity_matrix(QQ, rank(L))[1:rank(N),:])
      @assert genus(N) == genus(N2)
      M2 = lattice_in_same_ambient_space(L, identity_matrix(QQ,rank(L))[rank(N)+1:end, :])
      @assert genus(M) == genus(M2)
      push!(results, (L, M2, N2))
      @info "Gluing done"
      GC.gc()
    end
  end
  return results
end

function _isomorphism_classes_primitive_extensions_along_elementary(N::ZLat, M::ZLat, H::TorQuadMod)
  @assert is_one(basis_matrix(N))
  @assert is_one(basis_matrix(M))
  bool, p = is_elementary_with_prime(H)
  @assert bool
  @assert is_elementary(M, p)
  results = Tuple{ZLat, ZLat, ZLat}[]
  GN, _ = image_in_Oq(N)
  GM, _ = image_in_Oq(M)
  qN = domain(GN)
  qM = domain(GM)

  D, qNinD, qMinD = orthogonal_sum(qN, qM)
  OD = orthogonal_group(D)
  VN, VNinqN, _ = _get_V(N, qN, identity_matrix(QQ, rank(N)), id_hom(qN), minpoly(identity_matrix(QQ,1)), p)
  subsN = subgroups_orbit_representatives_and_stabilizers_elementary(VNinqN, GN, p)
  filter!(HN -> is_anti_isometric_with_anti_isometry(domain(HN[1]), H), subsN)
  @assert !isempty(subsN)
  subsM = subgroups_orbit_representatives_and_stabilizers_elementary(id_hom(qM), GM, p)
  filter!(HM -> is_isometric_with_isometry(domain(HM[1]), H), subsM)
  @assert !isempty(subsM)

  for H1 in subsN, H2 in subsM
    ok, phi = is_anti_isometric_with_anti_isometry(domain(H1[1]), domain(H2[1]))
    @assert ok

    HNinqN, stabN = H1
    OHN = orthogonal_group(HN)

    HMinqM, stabM = H2
    OHM = orthogonal_group(HM)

    actN = hom(stabN, OHN, [OHN(_restrict(x, HNinqN)) for x in gens(stabN)])

    actM = hom(stabM, OHM, [OHM(_restrict(x, HMinqM)) for x in gens(stabM)])
    imM, _ = image(actM)

    stabNphi = AutomorphismGroupElem{TorQuadMod}[OHM(compose(inv(phi), compose(hom(actN(g)), phi))) for g in gens(stabN)]
    stabNphi, _ = sub(OHM, stabNphi)
    reps = double_cosets(OHM, stabNphi, imM)
    @info "$(length(reps)) isomorphism classe(s) of primitive extensions"
    for g in reps
      g = representative(g)
      phig = compose(phi, hom(g))
      _glue = Vector{fmpq}[lift(qNinD(HNinqN(g))) + lift(qMinD(HMinqM(phig(g)))) for g in gens(domain(phig))]
      z = zero_matrix(QQ, 0, degree(N)+degree(M))
      glue = reduce(vcat, [matrix(QQ, 1, degree(N)+degree(M), g) for g in _glue], init = z)
      glue = vcat(identity_matrix(QQ, rank(N)+rank(M)), glue)
      glue = FakeFmpqMat(glue)
      _B = hnf(glue)
      _B = QQ(1, denominator(glue))*change_base_ring(QQ, numerator(_B))
      L = lattice(ambient_space(cover(D)), _B[end-rank(N)-rank(M)+1:end, :])
      N2 = lattice_in_same_ambient_space(L, identity_matrix(QQ, rank(L))[1:rank(N),:])
      @assert genus(N) == genus(N2)
      M2 = lattice_in_same_ambient_space(L, identity_matrix(QQ,rank(L))[rank(N)+1:end, :])
      @assert genus(M) == genus(M2)
      push!(results, (L, M2, N2))
      @info "Gluing done"
      GC.gc()
    end
  end
  return results
end


@doc Markdown.doc"""
    primitive_embeddings_in_primary_lattice(L::ZLat, M::ZLat)
                                     -> Vector{Tuple{ZLat, ZLat, ZLat}}

Given a `p`-primary lattice `L`, unique in its genus, and a lattice `M`,
compute representatives for all isomorphism classes of primitive embeddings
of `M` in `L` up to the actions of $\bar{O}(M)$ and $O(L)$. Here
$\bar{O}(M)$ denotes the image of $O(M)\to O(q_M)$.

The output is given in terms of triples `(L', M', N')` where `L'` is
isometric to `L`, `M'` is a sublattice of `L'` isometric to `M` and
`N'` is the orthogonal complement of `M'` in `L'`.
"""
function primitive_embeddings_in_primary_lattice(L::ZLat, M::ZLat, GL::AutomorphismGroup{TorQuadMod} = orthogonal_group(discriminant_group(L)); check::Bool = false)
  bool, p = is_primary_with_prime(L)
  @req bool "L must be unimodular or elementary"
  el = is_elementary(L, p)
  if check
    @req length(genus_representatives(L)) == 1 "L must be unique in its genus"
  end
  M = Zlattice(gram = gram_matrix(M))
  results = Tuple{ZLat, ZLat, ZLat}[]
  qL = rescale(discriminant_group(L), -1)
  GL = Oscar._orthogonal_group(qL, matrix.(gens(GL)))
  pL, _, nL = signature_tuple(L)
  pM, _, nM = signature_tuple(M)
  @req (pL-pM >= 0 && nL-nM >= 0) "Impossible embedding"
  if ((pL, nL) == (pM, nM))
    genus(M) != genus(L) && return results
    return [(M, M, lattice_in_same_ambient_space(M, zero_matrix(QQ, 0, degree(M))))]
  end
  qM = discriminant_group(M)
  D, (qMinD, qLinD), proj = Hecke._orthogonal_sum_with_injections_and_projections([qM, qL])
  if el
    VM, VMinqM, _ = _get_V(M, qM, identity_matrix(QQ, rank(M)), id_hom(qM), minpoly(identity_matrix(QQ,1)), p)
  else
    VM, VMinqM = primary_part(qM, p)
  end
  for k  in divisors(gcd(order(VM), order(qL)))
    @info "Glue order: $(k)"
    if el
      subsL = subgroups_orbit_representatives_and_stabilizers_elementary(id_hom(GL), GL, k)
    else
      subsL = subgroups_orbit_representatives_and_stabilizers(GL, order = k)
    end
    @info "$(length(subsL)) subgroup(s)"
    for H in subsL
      HL = domain(H[1])
      it = subgroups(abelian_group(VM), order = order(HL))
      subsM = [sub(qM, VMinqM.(VM.(j[2].(gens(j[1])))))[2] for j in it]
      filter!(HM -> is_anti_isometric_with_anti_isometry(domain(HM), HL)[1], subsM)
      isempty(subsM)  && continue
      @info "Possible gluing"
      HM = subsM[1]
      ok, phi = is_anti_isometric_with_anti_isometry(domain(HM), HL)
      @assert ok
      _glue = [lift(qMinD(HM(g))) + lift(qLinD(H[1](phi(g)))) for g in gens(domain(HM))]
      ext, _ = sub(D, D.(_glue))
      perp, j = orthogonal_submodule(D, ext)
      disc = torsion_quadratic_module(cover(perp), cover(ext), modulus = modulus_bilinear_form(perp),
                                                               modulus_qf = modulus_quadratic_form(perp))
      disc = rescale(disc, -1)
      !is_genus(disc, (pL-pM, nL-nM))  && continue
      G = genus(disc, (pL-pM, nL-nM))
      @info "We can glue: $G"
      Ns = representatives(G)
      @info "$(length(Ns)) possible orthogonal complement(s)"
      Ns = lll.(Ns)
      Ns = ZLat[Zlattice(gram=gram_matrix(N)) for N in Ns]
      ext2, _ = sub(D, [qMinD(HM(g)) for g in gens(domain(HM))])
      perp2, j = orthogonal_submodule(D, ext2)
      qM2, _ = sub(qM, [proj[1](j(g)) for g in gens(perp2)])
      for N in Ns
        append!(results, _isomorphism_classes_primitive_extensions(N, M, qM2))
        GC.gc()
      end
      GC.gc()
    end
  end
  @assert all(triple -> genus(triple[1]) == genus(L), results)
  return results
end

function primitive_embeddings_of_elementary_lattice(L::ZLat, M::ZLat, GL::AutomorphismGroup{TorQuadMod} = orthogonal_group(discriminant_group(L)); classification::Bool = false, check::Bool = false)
  bool, p = is_elementary_with_prime(M)
  @req bool "M must be elementary"
  if check
    @req length(genus_representatives(L)) == 1 "L must be unique in its genus"
  end
  results = Tuple{ZLat, ZLat, ZLat}[]
  qL = rescale(discriminant_group(L), -1)
  GL = Oscar._orthogonal_group(qL, matrix.(gens(GL)))
  pL, _, nL = signature_tuple(L)
  pM, _, nM = signature_tuple(M)
  @req (pL-pM >= 0 && nL-nM >= 0) "Impossible embedding"
  if ((pL, nL) == (pM, nM))
    if genus(M) != genus(L) 
      return false, results
    else
      return true, [(M, M, lattice_in_same_ambient_space(M, zero_matrix(QQ, 0, degree(M))))]
    end
  end
  qM = discriminant_group(M)
  D, (qMinD, qLinD), proj = Hecke._orthogonal_sum_with_injections_and_projections([qM, qL])
  VL, VLinqL, _ = _get_V(L, qL, identity_matrix(QQ, rank(L)), id_hom(qL), minpoly(identity_matrix(QQ,1)), p)
  for k  in divisors(gcd(order(qM), order(VL)))
    @info "Glue order: $(k)"
    subsL = subgroups_orbit_representatives_and_stabilizers_elementary(VLinqL, GL, k)
    @info "$(length(subsL)) subgroup(s)"
    subsM = subgroups_orbit_representatives_and_stabilizers_elementary(id_hom(qM), orthogonal_group(qM), k)
    for H in subsL
      HL = domain(H[1])
      _subsM = filter(HM -> is_anti_isometric_with_anti_isometry(domain(HM[1]), HL)[1], subsM)
      isempty(_subsM)  && continue
      @info "Possible gluing"
      HM = _subsM[1][1]
      ok, phi = is_anti_isometric_with_anti_isometry(domain(HM), HL)
      @assert ok
      _glue = [lift(qMinD(HM(g))) + lift(qLinD(H[1](phi(g)))) for g in gens(domain(HM))]
      ext, _ = sub(D, D.(_glue))
      perp, j = orthogonal_submodule(D, ext)
      disc = torsion_quadratic_module(cover(perp), cover(ext), modulus = modulus_bilinear_form(perp),
                                                               modulus_qf = modulus_quadratic_form(perp))
      disc = rescale(disc, -1)
      !is_genus(disc, (pL-pM, nL-nM))  && continue
      G = genus(disc, (pL-pM, nL-nM))
      !classification && return true, G
      @info "We can glue: $G"
      Ns = representatives(G)
      @info "$(length(Ns)) possible orthogonal complement(s)"
      Ns = lll.(Ns)
      Ns = ZLat[Zlattice(gram=gram_matrix(N)) for N in Ns]
      ext2, _ = sub(D, [qMinD(HM(g)) for g in gens(domain(HM))])
      perp2, j = orthogonal_submodule(D, ext2)
      qM2, _ = sub(qM, [proj[1](j(g)) for g in gens(perp2)])
      for N in Ns
        append!(results, _isomorphism_classes_primitive_extensions_along_elementary(N, M, qM2))
        GC.gc()
      end
      GC.gc()
    end
  end
  @assert all(triple -> genus(triple[1]) == genus(L), results)
  return (length(results) >0), results
end

####################################################################################
#
# Admissible equivariant primitive extensions
#
####################################################################################

function equivariant_primitive_extensions(A::LatticeWithIsometry, GA::AutomorphismGroup{TorQuadMod},
                                         B::LatticeWithIsometry, GB::AutomorphismGroup{TorQuadMod},
                                         L::ZLat; check::Bool = false)
  
  if check
    pA, nA = signature_tuple(A)
    pB, nB = signature_tuple(B)
    pL, nL = signature_tuple(L)
    @req pA+pB == pL "Incompatible signatures"
    @req nA+nB == nL "Incompatible signatures"
    @req ambient_space(A) === ambient_space(B) === ambient_space(L) "Lattices must all live in the same ambient space"
  end

  results = LatticeWithIsometry[]

  qA, fqA = discriminant_group(A)
  qB, fqB = discriminant_group(B)
  
  if check
    @req all(g -> compose(g, compose(fqA, inv(g))) == fqA, GA) "GA does not centralize fqA"
    @req all(g -> compose(g, compose(fqB, inv(g))) == fqB, GB) "GB does not centralize fqB"
  end
  D, qAinD, qBinD = inner_orthogonal_sum(qA, qB)
  OD = orthogonal_group(D)
  OqAinOD = embedding_orthogonal_group(qAinD)
  OqBinOD = embedding_orthogonal_group(qBinD)
  OqA = domain(OqAinOD)
  OqB = domain(OqBinOD)

  gamma, _, _ = glue_map(L, lattice(A), lattice(B))
  subsA = classes_conjugate_subgroups(GA, domain(gamma))
  @assert !isempty(subsN)
  subsB = classes_conjugate_subgroups(GB, codomain(gamma))
  @assert !isempty(subsM)

  for (H1, H2) in Hecke.cartesian_product_iterator([subsA, subsB], inplace=true)
    ok, phi = is_anti_isometric_with_anti_isometry(H1, H2)
    @assert ok
    SAinqA, stabA = H1
    OSAinOqA = _embedding_orthogonal_group(SAinqA)
    OSA = domain(OSAinOqA)
    OSAinOD = compose(OSAinOqA, OqAinOD)

    SBinqB, stabB = H2
    SB = domain(SBinqB)
    OSBinOqB = _embedding_orthogonal_group(SBinqB)
    OSB = domain(OSBinOqB)
    OSBinOD = compose(OSBinOqB, OqBinOD)

    # we compute the image of the stabalizers in the respective OS* and we keep track
    # of the elements of the stabilizers acting trivially in the respective S*
    actA = hom(stabA, OSA, [OSA(Oscar._restrict(x, SAinqA)) for x in gens(stabA)])
    imA, _ = image(actA)
    kerA = AutomorphismGroupElem{TorQuadMod}[OqAinOD(x) for x in gens(kernel(actA)[1])]
    fSA = OSA(_restrict(fqA, SAinqA))

    actB = hom(stabB, OSB, [OSB(Oscar._restrict(x, SBinqB)) for x in gens(stabB)])
    imB, _ = image(actB)
    kerB = AutomorphismGroupElem{TorQuadMod}[OqBinOD(x) for x in gens(kernel(actB)[1])]
    fSB = OSB(_restrict(fqB, SBinqB))


    fSAinOSB = OSB(compose(inv(phi), compose(hom(fSA), phi)))
    bool, g0 = representative_action(OSB, fSAinOSB, fSB)
    bool || continue
    phi = compose(phi, hom(OSB(g0)))
    fSAinOSB = OSB(compose(inv(phi), compose(hom(fSA), phi)))
    # Now the new phi is "sending" the restriction of fA to this of fB.
    # So we can glue SA and SB.
    @assert fSAinOSB == fSB

    # Now it is time to compute generators for O(SB, rho_l(qB), fB), and the induced
    # images of stabA|stabB for taking the double cosets next
    center, _ = centralizer(OSB, fSB)
    center, _ = sub(OSB, [OSB(c) for c in gens(center)])
    stabSAphi = AutomorphismGroupElem{TorQuadMod}[OSB(compose(inv(phi), compose(hom(actA(g)), phi))) for g in gens(stabA)]
    stabSAphi, _ = sub(center, stabSAphi)
    stabSB, _ = sub(center, [actB(s) for s in gens(stabB)])
    reps = double_cosets(center, stabSB, stabSAphi)

    for g in reps
      g = representative(g)
      phig = compose(phi, hom(g))

      _glue = Vector{fmpq}[lift(g) + lift(phig(g)) for g in gens(domain(phig))]
      z = zero_matrix(QQ, 0, degree(A))
      glue = reduce(vcat, [matrix(QQ, 1, degree(A), g) for g in _glue], init=z)
      glue = vcat(basis_matrix(lattice(A)+lattice(B)), glue)
      glue = FakeFmpqMat(glue)
      _B = hnf(glue)
      _B = QQ(1, denominator(glue))*change_base_ring(QQ, numerator(_B))
      C2 = lattice(ambient_space(A), _B[end-rank(A)-rank(B)+1:end, :])
      fC2 = block_diagonal_matrix([isometry(A), isometry(B)])
      _B = solve_left(reduce(vcat, basis_matrix.([A,B])), basis_matrix(C2))
      fC2 = _B*fC2*inv(_B)
      L2 = lattice_with_isometry(C2, fC2, ambient_representation = false)
      @assert genus(L2) == genus(L)
      push!(results, L2)
    end
  end
  return results
end

@doc Markdown.doc"""
    admissible_equivariant_primitive_extensions(Afa::LatticeWithIsometry,
                                                Bfb::LatticeWithIsometry,
                                                Cfc::LatticeWithIsometry,
                                                p::Int; check=true)
                                                     -> Vector{LatticeWithIsometry}

Given a triple of lattices with isometry `(A, fa)`, `(B, fb)` and `(C, fc)` and a
prime number `p`, such that `(A, B, C)` is `p`-admissible, return a set of
representatives of the double coset $G_B\backslash S\slash/G_A$ where:

  - $G_A$ and $G_B$ are the respective images of the morphisms
    $O(A, fa) -> O(q_A, \bar{fa})$ and $O(B, fb) -> O(q_B, \bar{fb})$;
  - $S$ is the set of all primitive extensions $A \perp B \subseteq C'$ with isometry
    $fc'$ where $p\cdot C' \subseteq A\perpB$ and the type of $(C', fc'^p)$ is equal
    to the type of $(C, fc)$.

If `check == true` the input triple is checked to a `p`-admissible triple of
integral lattices (with isometry) with `fA` and `fB` having relatively coprime
irreducible minimal polynomials and imposing that `A` and `B` are orthogonal
if `A`, `B` and `C` lie in the same ambient quadratic space.

Note that `Afa` and `Bfb` must be of pure type, i.e. the minimal polynomials
of the associated isometries must be irreducible (and relatively coprime).

See Algorithm 2 of [BH22].
"""
function admissible_equivariant_primitive_extensions(A::LatticeWithIsometry,
                                                     B::LatticeWithIsometry,
                                                     C::LatticeWithIsometry,
                                                     p::Integer; check=true)
  # requirement for the algorithm of BH22
  @req is_prime(p) "p must be a prime number" 

  if check 
    @req all(L -> is_integral(L), [A, B, C]) "Underlying lattices must be integral"
    @req is_of_hermitian_type(A) && is_of_hermitian_type(B) "Afa and Bfb must be of hermitian type"
    @req gcd(minpoly(A), minpoly(B)) == 1 "Minimal irreducible polynomials must be relatively coprime"
    @req is_admissible_triple(A, B, C, p) "Entries, in this order, do not define an admissible triple"
    if ambient_space(A) === ambient_space(B) === ambient_space(C)
      G = gram_matrix(ambient_space(C))
      @req iszero(basis_matrix(A)*G*transpose(basis_matrix(B))) "Lattice in same ambient space must be orthogonal"
    end 
  end 

  results = LatticeWithIsometry[]

  # this is the glue valuation: it is well-defined because the triple in input is admissible
  g = div(valuation(divexact(det(A)*det(B), det(C)), p), 2)

  qA, fqA = discriminant_group(A)
  qB, fqB = discriminant_group(B)
  GA = image_centralizer_in_Oq(A)
  GB = image_centralizer_in_Oq(B)

  # this is where we will perform the glueing
  if ambient_space(A) === ambient_space(B)
    D, qAinD, qBinD = inner_orthogonal_sum(qA, qB)
    OD = orthogonal_group(D)
    OqAinOD = embedding_orthogonal_group(qAinD)
    OqBinOD = embedding_orthogonal_group(qBinD)
  else
    D, qAinD, qBinD = orthogonal_sum(qA, qB)
    OD = orthogonal_group(D)
    OqAinOD, OqBinOD = embedding_orthogonal_group(qAinD, qBinD)
  end 

  OqA = domain(OqAinOD)
  OqB = domain(OqBinOD)

    # if the glue valuation is zero, then we glue along the trivial group and we don't
  # have much more to do. Since the triple is p-admissible, A+B = C
  if g == 0
    geneA = AutomorphismGroupElem{TorQuadMod}[OqAinOD(OqA(a.X)) for a in gens(GA)]
    geneB = AutomorphismGroupElem{TorQuadMod}[OqBinOD(OqB(b.X)) for b in gens(GB)]
    gene = vcat(geneA, geneB)
    GC2, _ = sub(OD, gene)
    if ambient_space(A) === ambient_space(B) === ambient_space(C)
      C2 = lattice(A)+lattice(B)
      fC2 = block_diagonal_matrix([isometry(A), isometry(B)])
      _B = solve_left(reduce(vcat, basis_matrix.([A,B])), basis_matrix(C2))
      fC2 = _B*fC2*inv(_B)
      @assert fC2*gram_matrix(C2)*transpose(fC2) == gram_matrix(C2)
    else
      C2 = direct_sum(lattice(A), lattice(B))[1]
      fC2 = block_diagonal_matrix([isometry(A), isometry(B)])
    end
    if is_of_type(lattice_with_isometry(C2, fC2^p, ambient_representation = false), type(C))
      C2fc2 = lattice_with_isometry(C2, fC2, ambient_representation=false)
      set_attribute!(C2fc2, :image_centralizer_in_Oq, GC2)
      push!(results, C2fc2)
    end
    return results
  end

  # these are GA|GB-invariant, fA|fB-stable, and should contain the kernels of any glue map
  VA, VAinqA, fVA = _get_V(lattice(A), qA, isometry(A), fqA, minpoly(B), p)
  VB, VBinqB, fVB = _get_V(lattice(B), qB, isometry(B), fqB, minpoly(A), p)

  # since the glue kernels must have order p^g, in this condition, we have nothing
  if min(order(VA), order(VB)) < p^g
    return results
  end

  # scale of the dual: any glue kernel must contain the multiples of l of the respective
  # discriminant groups
  l = level(genus(C))
  

  # We look for the GA|GB-invariant and fA|fB-stable subgroups of VA|VB which respectively
  # contained lqA|lqB. This is done by computing orbits and stabilisers of VA/lqA (resp VB/lqB)
  # seen as a F_p-vector space under the action of GA (resp. GB). Then we check which ones
  # are fA-stable (resp. fB-stable)

  subsA = subgroups_orbit_representatives_and_stabilizers_elementary(VAinqA, GA, p^g, fVA, ZZ(l))
  subsB = subgroups_orbit_representatives_and_stabilizers_elementary(VBinqB, GB, p^g, fVB, ZZ(l))
  # once we have the potential kernels, we create pairs of anti-isometric groups since glue
  # maps are anti-isometry
  R = Tuple{eltype(subsA), eltype(subsB), TorQuadModMor}[]
  for H1 in subsA, H2 in subsB
    ok, phi = is_anti_isometric_with_anti_isometry(domain(H1[1]), domain(H2[1]))
    !ok && continue
    push!(R, (H1, H2, phi))
  end

  # now, for each pair of anti-isometric potential kernels, we need to see whether
  # it is (fA,fB)-equivariant, up to conjugacy. For each working pair, we compute the
  # corresponding overlattice and check whether it satisfies the type condition
  for (H1, H2, phi) in R
    SAinqA, stabA = H1
    OSAinOqA = _embedding_orthogonal_group(SAinqA)
    OSA = domain(OSAinOqA)
    OSAinOD = compose(OSAinOqA, OqAinOD)

    SBinqB, stabB = H2
    SB = domain(SBinqB)
    OSBinOqB = _embedding_orthogonal_group(SBinqB)
    OSB = domain(OSBinOqB)
    OSBinOD = compose(OSBinOqB, OqBinOD)

    # we compute the image of the stabalizers in the respective OS* and we keep track
    # of the elements of the stabilizers acting trivially in the respective S*
    actA = hom(stabA, OSA, [OSA(Oscar._restrict(x, SAinqA)) for x in gens(stabA)])
    imA, _ = image(actA)
    kerA = AutomorphismGroupElem{TorQuadMod}[OqAinOD(x) for x in gens(kernel(actA)[1])]
    fSA = OSA(_restrict(fqA, SAinqA))

    actB = hom(stabB, OSB, [OSB(Oscar._restrict(x, SBinqB)) for x in gens(stabB)])
    imB, _ = image(actB)
    kerB = AutomorphismGroupElem{TorQuadMod}[OqBinOD(x) for x in gens(kernel(actB)[1])]
    fSB = OSB(_restrict(fqB, SBinqB))

    # we get all the elements of qB of order exactly p^{l+1}, which are not mutiple of an
    # element of order p^{l+2}. In theory, glue maps are classified by the orbit of phi
    # under the action of O(SB, rho_{l+1}(qB), fB)
    rB, rBinqB = _rho_functor(qB, p, valuation(l, p)+1)
    @assert Oscar._is_invariant(stabB, rBinqB)
    rBinSB = hom(rB, SB, TorQuadModElem[SBinqB\(rBinqB(k)) for k in gens(rB)])
    @assert is_injective(rBinSB) # we indeed have rho_{l+1}(qB) which is a subgroup of SB

    # We compute the generators of O(SB, rho_{l+1}(qB))
    OSBrB, _ = stabilizer(OSB, rBinSB)
    @assert fSB in OSBrB

    # phi might not "send" the restriction of fA to this of fB, but at least phi*fA*phi^-1
    # should be conjugate to fB inside O(SB, rho_l(qB)) for the glueing.
    # If not, we try the next potential pair.
    fSAinOSBrB = OSBrB(compose(inv(phi), compose(hom(fSA), phi)))
    bool, g0 = representative_action(OSBrB, fSAinOSBrB, OSBrB(fSB))
    bool || continue
    phi = compose(phi, hom(OSB(g0)))
    fSAinOSBrB = OSB(compose(inv(phi), compose(hom(fSA), phi)))
    # Now the new phi is "sending" the restriction of fA to this of fB.
    # So we can glue SA and SB.
    @assert fSAinOSBrB == fSB

    # Now it is time to compute generators for O(SB, rho_l(qB), fB), and the induced
    # images of stabA|stabB for taking the double cosets next
    center, _ = centralizer(OSBrB, OSBrB(fSB))
    center, _ = sub(OSB, [OSB(c) for c in gens(center)])
    stabSAphi = AutomorphismGroupElem{TorQuadMod}[OSB(compose(inv(phi), compose(hom(actA(g)), phi))) for g in gens(stabA)]
    stabSAphi, _ = sub(center, stabSAphi)
    stabSB, _ = sub(center, [actB(s) for s in gens(stabB)])
    reps = double_cosets(center, stabSB, stabSAphi)

    # now we iterate over all double cosets and for each representative, we compute the
    # corresponding overlattice in the glueing. If it has the wanted type, we compute
    # the image of the centralizer in OD from the stabA and stabB.

    for g in reps
      g = representative(g)
      phig = compose(phi, hom(g))

      if ambient_space(A) === ambient_space(B)
        _glue = Vector{fmpq}[lift(g) + lift(phig(g)) for g in gens(domain(phig))]
        z = zero_matrix(QQ, 0, degree(A))
        glue = reduce(vcat, [matrix(QQ, 1, degree(A), g) for g in _glue], init=z)
        glue = vcat(basis_matrix(lattice(A)+lattice(B)), glue)
        glue = FakeFmpqMat(glue)
        _B = hnf(glue)
        _B = QQ(1, denominator(glue))*change_base_ring(QQ, numerator(_B))
        C2 = lattice(ambient_space(A), _B[end-rank(A)-rank(B)+1:end, :])
        fC2 = block_diagonal_matrix([isometry(A), isometry(B)])
        _B = solve_left(reduce(vcat, basis_matrix.([A,B])), basis_matrix(C2))
        fC2 = _B*fC2*inv(_B)
      else
        _glue = Vector{fmpq}[lift(qAinD(SAinqA(g))) + lift(qBinD(SBinqB(phig(g)))) for g in gens(domain(phig))]
        z = zero_matrix(QQ,0,degree(A)+degree(B))
        glue = reduce(vcat, [matrix(QQ,1,degree(A)+degree(B),g) for g in _glue], init=z)
        glue = vcat(block_diagonal_matrix([basis_matrix(A), basis_matrix(B)]), glue)
        glue = FakeFmpqMat(glue)
        _B = hnf(glue)
        _B = QQ(1, denominator(glue))*change_base_ring(QQ, numerator(_B))
        C2 = lattice(ambient_space(cover(D)), _B[end-rank(A)-rank(B)+1:end, :])
        fC2 = block_diagonal_matrix([isometry(A), isometry(B)])
        __B = solve_left(block_diagonal_matrix(basis_matrix.([A,B])), basis_matrix(C2))
        fC2 = __B*fC2*inv(__B)
        @assert fC2*gram_matrix(C2)*transpose(fC2) == gram_matrix(C2)
      end

      if !is_of_type(lattice_with_isometry(C2, fC2^p, ambient_representation=false), type(C))
        continue
      end

      ext, _ = sub(D, D.(_glue))
      perp, j = orthogonal_submodule(D, ext)
      disc = torsion_quadratic_module(cover(perp), cover(ext), modulus = modulus_bilinear_form(perp),
                                                               modulus_qf = modulus_quadratic_form(perp))

      qC2 = discriminant_group(C2)
      ok, phi2 = is_isometric_with_isometry(qC2, disc)
      @assert ok

      C2 = lattice_with_isometry(C2, fC2, ambient_representation=false)
      im2_phi, _ = sub(OSA, OSA.([compose(phig, compose(hom(g1), inv(phig))) for g1 in gens(imB)]))
      im3, _ = sub(imA, gens(intersect(imA, im2_phi)[1]))
      stab = Tuple{AutomorphismGroupElem{TorQuadMod}, AutomorphismGroupElem{TorQuadMod}}[(x, imB(compose(inv(phig), compose(hom(x), phig)))) for x in gens(im3)]
      stab = AutomorphismGroupElem{TorQuadMod}[OSAinOD(x[1])*OSBinOD(x[2]) for x in stab]
      stab = union(stab, kerA)
      stab = union(stab, kerB)
      stab = TorQuadModMor[_restrict(g, j) for g in stab]
      stab = TorQuadModMor[hom(disc, disc, [disc(lift(g(perp(lift(l))))) for l in gens(disc)]) for g in stab]
      stab = Oscar._orthogonal_group(discriminant_group(C2), [compose(phi2, compose(g, inv(phi2))).map_ab.map for g in stab])
      set_attribute!(C2, :image_centralizer_in_Oq, stab)
      push!(results, C2)
    end
  end

  return results
end

