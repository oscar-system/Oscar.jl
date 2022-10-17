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

function _image(p::fmpz_poly, f::TorQuadModMor)
  M = f.map_ab.map
  M = p(M)
  fab = hom(domain(f.map_ab), codomain(f.map_ab), M)
  pf = hom(domain(f), codomain(f), fab)
  return pf
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

function __orbits_and_stabilizers_elementary_abelian_subgroups(G, aut, gens_aut, gens_act, min_order, max_order)
  Gap = G.X
  pcgs = GG.Pcgs(Ggap)
  p = GG.Order(pcgs[1])
  if p == 2
    act = GG.OnSubspacesByCanonicalBasisGF2
  else
    act = GG.OnSubspacesByCanonicalBasis
  end
  F = GF(p)
  Fgap = GG.GF(p)
  function matrix_representation(f, pcgs)
      return  [[GG.ExponentsOfPcElements(pcgs, GG.Image(f,i))[1]*GG.One(Fgap)] for i in pcgs]
  end
  
  gens_mats = GAP.julia_to_gap([matrix_representation(f, pcgs) for f in gens_act], recursive=true)

  n = length(pcgs)
  V = VectorSpace(F, n)
  valmin = valuation(min_order, p)
  valmax = valuation(max_order, p)
  valmax = max(valmax, n)

  if minval == 0
    subs = [(sub(G, elem_type(G)[])[1], aut)]
    minval = 1
  else
    subs = []
  end

  for k in valmin:valmax
    it = Hecke.subsets(n, k)
    b = basis(V)
    _subs = GAP.julia_to_gap([GG.LinearCombination(b, GAP.julia_to_gap(vec(collect(lift.(b[l].v))), recursive=true)) for l in it], recursive=true)
    orbs = GG.OrbitsDomain(aut, _subs, gens_aut, gens_mats, act)

    for orb in orbs
      orb = orb[1]
      stab = GG.Stabilizer(aut, orb, gens_aut, gens_mats, act)
      gene = [GG.PcElementByExponents(pcgs, i) for i in orb]
      gene = G.(gene)
      subgrp = sub(G, gene)[1]
      push!(subs, (subgrp, stab))
    end
  end
  return subs
end

function __orbits_and_stabilizers_elementary_abelian_subgroups_equivariant(G, aut, gens_aut, gens_act, min_order, max_order, g)
  @assert g.X in GG.Center(aut.X)

  Gap = G.X
  pcgs = GG.Pcgs(Ggap)
  p = GG.Order(pcgs[1])
  if p == 2
    act = GG.OnSubspacesByCanonicalBasisGF2
  else
    act = GG.OnSubspacesByCanonicalBasis
  end
  F = GF(p)
  Fgap = GG.GF(p)
  function matrix_representation(f, pcgs)
      return  [[GG.ExponentsOfPcElements(pcgs, GG.Image(f,i))[1]*GG.One(Fgap)] for i in pcgs]
  end
  
  gens_mats = GAP.julia_to_gap([GAPmatrix_representation(f, pcgs) for f in gens_act], recursive=true)

  n = length(pcgs)
  V = VectorSpace(F, n)
  valmin = valuation(min_order, p)
  valmax = valuation(max_order, p)
  valmax = max(valmax, n)

  if minval == 0
    subs = [(sub(G, elem_type(G)[])[1], aut)]
    minval = 1
  else
    subs = []
  end

  for k in valmin:valmax
    it = Hecke.subsets(n, k)
    b = basis(V)
    _subs = [GG.LinearCombination(b, GAP.julia_to_gap(vec(collect(lift.(b[l].v))), recursive=true)) for l in it]
    _subs = GAP.julia_to_gap(filter!(i -> act(i, g) == i, _subs), recursive=true)
    orbs = GG.OrbitsDomain(aut, _subs, gens_aut, gens_mats, act)

    for orb in orbs
      orb = orb[1]
      stab = GG.Stabilizer(aut, orb, gens_aut, gens_mats, act)
      gene = [GG.PcElementByExponents(pcgs, i) for i in orb]
      gene = G.(gene)
      subgrp = sub(G, gene)[1]
      push!(subs, (subgrp, stab))
    end
  end
  return subs
end

function __subgroup_representatives1(epi, aut, gens_aut, gens_act, max_order)
  G0 = domain(epi)
  G1 = codomain(epi)
  N = GG.InvariantElementaryAbelianSeries(G1.X, gens_act, G1.X, true)
  N = N[end-1]
  invs = GG.InvariantSubgroupsElementaryAbelianGroup(N, gens_act)
  filter!(i -> GG.IsNonTrivial(i), invs)
  m = minimum([GG.Size(i) for i in invs])
  for N in invs
    if GG.Size(N) == m
      break
    end
  end
  subgrp_reps = __orbits_and_stabilizers_elementary_abelian_subgroups(N, aut, gens_aut, gens_act, 1, max_order)
  
  if GG.Size(N) == order(G1)
    return subgrp_reps
  end

  epi_mod = GG.NaturalHomomorphismByNormalSubgroupNC(G1.X, N)
  GmodN = GG.Image(epi_mod)
  aut_on_GmodN = GAP.julia_to_gap([GG.InducedAutomorphism(epi_mod, i) for i in gens_act], recursive=true)
  epi_new = GG.CompositionMapping(epi_mod, epi)
  Alist = __subgroup_representatives1(epi_new, aut, gens_aut, aut_on_GmodN, max_order)
  for j in 1:length(Alist)
    A = Alist[j]
    repr = GG.PreImage(epi_mod, A[1].X)
    repr = Oscar._oscar_group(repr, A[1])
    Alist[j] = (repr, A[1])
  end
  @assert order(Alist[1][1]) == GG.Size(N)
  popfirst!(Alist)
  subgrp_reps = vcat(subgrp_reps, Alist)

  for A in Alist
    stab = A[2]
    A = A[1]
    act = GAP.julia_to_gap([GG.InducedAutomorphism(epi, i) for i in GG.GeneratorsOfGroup(stab)], recursive=true)
    Blist_A = __orbits_and_stabilizers_elementary_abelian_subgroups(N, stab, GG.GeneratorsOfGroup(stab), act, 1, max_order)
    if order(Blist_A[end]) == order(N)
      pop!(Blist_A)
    end
    for B in Blist_A
      C_AB = B[2]
      B = B[1]
      epiB = GG.NaturalHomomorphismByNormalSubgroupNC(A.X, B.X)
      AmodB = GG.Image(epiB)
      NmodB = GG.Image(epiB, N.X)
      UmodBList = GG.ComplementClassesRepresentatives(AmodB, NmodB)
      gens_C_AB = GG.GeneratorsOfGroup(C_AB)
      act_C_AB = GAP.julia_to_gap([GG.InducedAutomorphism(epi, i) for i in gens_C_AB], recursive=true)
      act_C_AB = GAP.julia_to_gap([GG.InducedAutomorphism(epiB, i) for i in act_C_AB], recursive=true)
      complement_reps = GG.ExternalOrbitsStabilisers(C_AB, UmodBList, gens_C_AB, act_C_AB, GAP.evalstr("function(x,g); return Image(g,x);end;"))
      complement_reps = [(Oscar._oscar_group(GG.PreImage(epiB, GG.Representative(i)), A), GG.StabilizerOfExternalSet(i)) for i in complement_reps]
      subgrp_reps = vcat(subgrp_reps, complement_reps)
    end
  end
  sort!(subgrp_reps, s -> order(s[1]))
  return subgrp_reps
end

function __subgroup_representatives(G, aut, max_order)
  gens_aut = GG.GeneratorsOfGroup(aut)
  epi = GG.GroupHomomorphismByImages(G.X, G.X, GG.GeneratorsOfGroup(G.X), GG.GeneratorsOfGroup(G.X))
  return __subgroup_representatives1(epi, aut, gens_aut, gens_aut, max_order)
end

function __subgroup_representatives_elementary(G, aut, order)
  gens_aut = GG.GeneratorsOfGroup(aut)
  return __orbits_and_stabilizers_elementary_abelian_subgroups(G, aut, gens_aut, gens_aut, order, order)
end

function __subgroup_representatives_elementary_equivariant(G, aut, order, g)
  gens_aut = GG.GeneratorsOfGroup(aut)
  return __orbits_and_stabilizers_elementary_abelian_subgroups_equivariant(G, aut, gens_aut, gens_aut, order, order, g)
end

function _subgroups_representatives(H, G, order, g = 1)                                             
  g = G(g)
  if order(H) == 1
    return [(H, G.X)]
  end

  if !is_prime(elementary_divisors(H)[1])
    subgrp_reps = __subgroup_representatives(H, G.X, order)
    subgrp_reps = [S for S in subgrp_reps if order(S[1]) == order]
  elseif order(H) < ZZ(2^8)
    subgrp_reps = __subgroup_representatives_elementary_equivariant(H, G.X, order, g.X)
    subgrp_reps = [S for S in subgrp_reps if order(S[1]) == order]
  else
    a=1
  end
 end

