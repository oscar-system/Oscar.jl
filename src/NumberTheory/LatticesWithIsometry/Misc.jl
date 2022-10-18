GG = GAP.Globals

###############################################################################
#
#  Cyclotomic polynomials
#
###############################################################################

function _cyclotomic_polynomial(n::Int64)
  @assert n > 0
  _, x = QQ["x"]
  return Hecke.cyclotomic(n, x)
end

function _is_cyclotomic_polynomial(p::Union{fmpz_poly, fmpq_poly})
  n = degree(p)
  R = parent(p)
  x = gen(R)
  list_cyc = union(Int64[k for k in euler_phi_inv(n)], [1])
  list_poly = [Hecke.cyclotomic(k, x) for k in list_cyc]
  return any(q -> R(collect(coefficients(q))) == p, list_poly)
end

###############################################################################
#
#  Exponent of fmpq/fmpz_mat
#
###############################################################################

function _is_of_finite_exponent(f::Union{fmpq_mat, fmpz_mat})
  !Hecke.is_squarefree(minpoly(f)) && return false
  chi = charpoly(f)
  fact = collect(factor(chi))
  return all(p -> _is_cyclotomic_polynomial(p[1]), fact)
end

function _exponent(f::Union{fmpq_mat, fmpz_mat})
  !_is_of_finite_exponent(f) && return -1
  degs = unique(degree.([p[1] for p in collect(factor(minpoly(f)))]))
  exps = euler_phi_inv(degs[1])::Vector{fmpz}
  for i in 2:length(degs)
    union!(exps, euler_phi_inv(degs[i]))
  end
  maxdeg = lcm(exps)
  divmd = divisors(maxdeg)
  n = findfirst(k -> isone(f^k), divmd)
  @assert n !== nothing
  return divmd[n]
end

function _image(p::Union{fmpz_poly, fmpq_poly}, f::TorQuadModMor)
  p = change_base_ring(ZZ, p)
  M = f.map_ab.map
  M = p(M)
  pf = hom(domain(f), codomain(f), M)
  return pf
end

function _image(p::Union{fmpz_poly, fmpz_poly}, f::AutomorphismGroupElem{TorQuadMod})
  D = domain(f)
  M = matrix(f)
  return _image(p, hom(D, D, M))
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
  D = domain(f)
  M = matrix(f)
  return _restrict(hom(D, D, M), i)
end

##############################################################################
#
#  Group functionalities
#
##############################################################################

function embedding_orthogonal_group(i1, i2)
   D = codomain(i1)
   A = domain(i1)
   B = codomain(i2)
   gene = vcat(i1.(gens(A)), i2.(gens(B)))
   n = ngens(A)
   OD, OA, OB = orthogonal_group.([D, A, B])

   geneOA = elem_type(OD)[]
   for f in gens(OA)
     imgs = [i1(f(a)) for a in gens(A)]
     imgs = vcat(imgs, gene[n+1:end])
     _f = hom(lift.(gene), lift.(imgs))
     f = TorQuadModMor(D, D, _f)                                                                    
     push!(geneOA, f)
   end
   geneOB = elem_type(OD)[]
   for f in gens(OB)
     imgs = [i2(f(a)) for a in gens(B)]
     imgs = vcat(gene[1:n], imgs)
     _f = hom(lift.(gene), lift.(imgs))
     f = TorQuadModMor(D, D, _f)
     push!(geneOB, f)
   end

   OAtoOD = hom(OA, OD, gens(OA), geneOA)
   OBtoOD = hom(OB, OD, gens(OB), geneOB)
   return OAtoOD, OBtoOD
 end

function _as_Fp_vector_space_quotient(HinV, p, f)
  i = HinV.map_ab
  H = domain(HinV)
  Hab = domain(i)
  Hs, HstoHab = snf(Hab)
  f = _restrict(f, HinV)
  V = codmain(HinV)
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
  MVs = matrix(compose(inv(VstoVab), compose(f.map_ab, VstoVab)))
  Vp = VectorSpace(F, n)

  function _VstoVp(x::GrpAbFinGenElem)
    v = x.coeff
    return Vp(vec(collect(v)))
  end

  function _VptoVs(v::ModuleElem{gfp_elem})
    x = lift.(v.v)
    return Vs(vec(collect(x)))
  end

  VstoVp = Hecke.MapFromFunc(_VstoVp, _VptoVs, Vs, Vp)
  VtoVp = compose(VtoVs, VstoVp)
  M = reduce(vcat, [(i(HstoHab(a))).coeff for a in gens(Hs)])
  subgene = [Vp(M*v.v) for v in gens(Vp)]
  Hp, _ = sub(Vp, subgene)
  Qp, VptoQp = quo(Vp, Hp)
  fVp = change_base_ring(F, matrix(MVs))
  ok, fQp = can_solve_with_solution(transpose(VptoQp.matrix), transpose(Vp.toQp.matrix)*fVp)
  @assert ok


  return Qp, VtoVp, VptoQp, fQp
end

function _subgroups_representatives(Vinq::TorQuadModMor, G, g, f = one(G), l::Int = 0)
  V = domain(Vinq)
  q = codomain(Vinq)
  p = elementary_divisors(V)[1]
  @req all(a -> haspreimage(Vtoq.map_ab, data(l*a))[1], gens(q)) "l*q is not contained in V"
  H0, H0inV = sub(V, [preimage(Vtoq, (l*a)) for a in gens(q)])
  Qp, VtoVp, VptoQp, fQp = _as_Fp_vector_space_quotient(H0inV, p, f)

  gene_GV = [_restrict(g, Vinq) for g in gens(G)]
  GV = sub(orthogonal_group(V), gene_GV)
  @assert f in GV
  GVinG = hom(GV, G, gens(GV), gens(G))

  act_GV_Vp = [change_base_ring(base_ring(Qp), matrix(gg)) for gg in gens(GV)]
  act_GV_Qp = [solve_left(transpose(VptoQp).matrix, transpose(VptoQp).matrix*g) for g in act_GV_Vp]
  MGp = matrix_group(act_GV_Qp)
  @assert fQp in MGp
  MGptoGV = hom(MGp, GV, gens(GV))
  MGptoG = compose(MGptoGV, GVtoG)

  res = []
  if g-ngens(snf(abelian_group(H10))) > dim(Qp)
    return res
  end
  F = base_ring(Qp)
  k, K = kernel(VptoQp.matrix, side = :left)
  gene_H0p = [Vp(vec(collect(K[i,:]))) for i in 1:k]
  orb_and_stab = orbit_representatives_and_stabilizers(MGp, g-k)
  for (orb, stab) in orb_and_stab
    i = orb.map
    _V = codomain(i)
    for v in gens(ord)
      vv = _V(transpose(fQp*i(v)))
      try
        preimage(i, vv)
      catch
        @goto non_fixed
      end
    end
    gene_orb_in_Qp = [Qp(vec(collect(i(v).v))) for v in gens(orb)]
    gene_orb_in_Vp = [preimage(VptoQp, v) for v in gene_orb_in_Qp]
    gene_submod_in_Vp = vcat(gene_ord_in_Vp, gene_H0p)
    gene_submod_in_V = [preimage(VtoVp, v) for v in gene_submod_in_Vp]
    gene_submod_in_q = [image(Vinq, v) for v in gene_submod_in_V]
    orbq, _ = sub(q, gene_submod_in_q)
    @assert order(orb) == p^g
    stabq = image(MGptoG, stab)
    push!(res, (orbq, stabq))
    @label non_fixed
  end
  return res
end

