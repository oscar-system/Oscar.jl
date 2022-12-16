GG = GAP.Globals

################################################################################
#
# Inner orthogonal sum
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

##############################################################################
#
#  Group functionalities
#
##############################################################################

function embedding_orthogonal_group(i::TorQuadModMor)
  ok, j = has_complement(i)
  @req ok "domain(i) needs to have a complement in codomain(i)"
  A = domain(i)
  B = domain(j)
  D = codomain(i)
  gene = vcat(i.(gens(A)), j.(gens(B)))
  n = ngens(A)
  OD = orthogonal_group(D)
  OA = is_degenerate(A) ? orthogonal_group_degenerate(A) : orthogonal_group(A)

  geneOA = elem_type(OD)[]
  for f in gens(OA)
    imgs = [i(f(a)) for a in gens(A)]
    imgs = vcat(imgs, gene[n+1:end])
    if can_solve(reduce(vcat, [data(a).coeff for a in gene]), reduce(vcat, [data(a).coeff for a in imgs]))
      _f = solve(reduce(vcat, [data(a).coeff for a in gene]), reduce(vcat, [data(a).coeff for a in imgs]))
      _f = hom(D.ab_grp, D.ab_grp, _f)
    else
      _f = solve(reduce(vcat, [data(a).coeff for a in imgs]), reduce(vcat, [data(a).coeff for a in gene]))
      _f = inv(hom(D.ab_grp, D.ab_grp, _f))
    end
    f = TorQuadModMor(D, D, _f)
    push!(geneOA, OD(f))
  end
  OAtoOD = hom(OA, OD, gens(OA), geneOA)
  return OAtoOD::GAPGroupHomomorphism
end

function embedding_orthogonal_group(i1, i2)
  D = codomain(i1)
  A = domain(i1)
  B = domain(i2)
  gene = vcat(i1.(gens(A)), i2.(gens(B)))
  n = ngens(A)
  OD, OA, OB = orthogonal_group.([D, A, B])

  geneOA = elem_type(OD)[]
  for f in gens(OA)
    imgs = [i1(f(a)) for a in gens(A)]
    imgs = vcat(imgs, gene[n+1:end])
    _f = map_entries(ZZ, solve(reduce(vcat, [data(a).coeff for a in gene]), reduce(vcat, [data(a).coeff for a in imgs])))
    _f = hom(D.ab_grp, D.ab_grp, _f)
    f = TorQuadModMor(D, D, _f)                                                                    
    push!(geneOA, OD(f))
  end
  geneOB = elem_type(OD)[]
  for f in gens(OB)
    imgs = [i2(f(a)) for a in gens(B)]
    imgs = vcat(gene[1:n], imgs)
    _f = map_entries(ZZ, solve(reduce(vcat, [data(a).coeff for a in gene]), reduce(vcat, [data(a).coeff for a in imgs])))
    _f = hom(D.ab_grp, D.ab_grp, _f)
    f = TorQuadModMor(D, D, _f)
    push!(geneOB, OD(f))
  end

  OAtoOD = hom(OA, OD, gens(OA), geneOA)
  OBtoOD = hom(OB, OD, gens(OB), geneOB)
  return (OAtoOD, OBtoOD)::Tuple{GAPGroupHomomorphism, GAPGroupHomomorphism}
 end

function _as_Fp_vector_space_quotient(HinV, p, f)
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

function _subgroups_representatives(Vinq::TorQuadModMor, G::AutomorphismGroup{TorQuadMod}, g::Int, f::Union{TorQuadModMor, AutomorphismGroupElem{TorQuadMod}} = one(G), l::Union{Int, fmpz} = 0)
  V = domain(Vinq)
  q = codomain(Vinq)
  p = elementary_divisors(V)[1]
  @req all(a -> haspreimage(Vinq.map_ab, data(l*a))[1], gens(q)) "lq is not contained in V"
  H0, H0inV = sub(V, [preimage(Vinq, (l*a)) for a in gens(q)])
  Qp, VtoVp, VptoQp, fQp = _as_Fp_vector_space_quotient(H0inV, p, f)
  Vp = codomain(VtoVp)

  if dim(Qp) == 0
    if order(V) == p^g
      return Tuple{TorQuadModMor, AutomorphismGroup{TorQuadMod}}[(Vinq, G)]
    end
    return Tuple{TorQuadModMor, AutomorphismGroup{TorQuadMod}}[]
  end

  OV = orthogonal_group(V)
  gene_GV = TorQuadModMor[OV(_restrict(g, Vinq)) for g in gens(G)]
  GV, _ = sub(OV, gene_GV)
  @assert GV(f) in GV
  GVinG = hom(GV, G, gens(GV), gens(G))

  act_GV_Vp = gfp_fmpz_mat[change_base_ring(base_ring(Qp), matrix(gg)) for gg in gens(GV)]
  act_GV_Qp = gfp_fmpz_mat[solve(VptoQp.matrix, g*VptoQp.matrix) for g in act_GV_Vp]
  MGp = matrix_group(act_GV_Qp)
  @assert fQp in MGp
  MGptoGV = hom(MGp, GV, gens(GV))
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
      if !can_solve_with_solution(i.matrix, vv.v, side = :left)[1]
        @goto non_fixed
      end
    end
    gene_orb_in_Qp = AbstractAlgebra.Generic.QuotientModuleElem{gfp_fmpz_elem}[Qp(vec(collect(i(v).v))) for v in gens(orb)]
    gene_orb_in_Vp = AbstractAlgebra.Generic.ModuleElem{gfp_fmpz_elem}[preimage(VptoQp, v) for v in gene_orb_in_Qp]
    gene_submod_in_Vp = vcat(gene_orb_in_Vp, gene_H0p)
    gene_submod_in_V = TorQuadModElem[preimage(VtoVp, Vp(v)) for v in gene_submod_in_Vp]
    gene_submod_in_q = TorQuadModElem[image(Vinq, v) for v in gene_submod_in_V]
    orbq, orbqinq = sub(q, gene_submod_in_q)
    @assert order(orbq) == p^g
    stabq, _ = image(MGptoG, stab)
    push!(res, (orbqinq, stabq))
    @label non_fixed
  end
  return res
end

