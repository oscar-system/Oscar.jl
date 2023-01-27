export possible_ideals_for_cubics,
       possible_ideals_for_k3_models

function _subgroup_trivial_action_on_H20(RR::RepRing, chi1::Oscar.GAPGroupClassFunction, chi2::Oscar.GAPGroupClassFunction)
  h =  conjugacy_classes(chi1.table)
  dchi1 = determinant(chi1)
  dchi2 = determinant(chi2)
  ker = elem_type(underlying_group(RR))[]
  for i in 1:length(h)
    if dchi1[i] == dchi2[i]
        append!(ker, collect(h[i]))
    end
  end
  ker, _ = sub(underlying_group(RR), ker)
  return ker
end

function _id_symplectic_subgroup_k3(RR::RepRing, chi1::Oscar.GAPGroupClassFunction, chi2::Oscar.GAPGroupClassFunction, p::GAPGroupHomomorphism)
  ker = _subgroup_trivial_action_on_H20(RR, chi1, chi2)
  return p(ker)[1].X == codomain(p).X
end

function possible_ideals_for_k3_models(G::Tuple{Int, Int}, Gs::Tuple{Int, Int}, n::Int, d::Int, t::Int)
  bool, RR, sum_index, p = _has_pfr(small_group(G[1], G[2]), n)
  if bool == false
    return Tuple{ProjRep, Vector{SymmetricCompleteIntersections}}[]
  end
  @info "$(length(sum_index)) possible action(s) to consider"
  F = splitting_field(RR)
  S, _ = grade(PolynomialRing(F, "x" => 0:n-1)[1])
  _, j = homogeneous_component(S, d)
  E = underlying_group(RR)
  Irr = irreducible_characters_underlying_group(RR)
  res = Tuple{ProjRep, Vector{SymmetricCompleteIntersections}}[]
  for l in sum_index
    @info "Test a character"
    chi = sum(Irr[l])
    detchi = determinant(chi)
    chid = symmetric_power(conj(chi), d)
    ct = constituents(chid, t)
    chis = Oscar.GAPGroupClassFunction[]
    for chi2 in ct
      detchi2 = determinant(chi2)
      if detchi == detchi2
        push!(chis, chi2)
      end
    end
    length(chis) == 0 ? continue : nothing
    @info "$(length(chis)) possibility.ies"
    rep = affording_representation(RR, chi)
    poss = CharacterGrassmannian[character_grassmannian(homogeneous_polynomial_representation(rep, d), nu) for nu in chis]
    prep = ProjRep(rep, p)
    poss = SymmetricCompleteIntersections[SymmetricCompleteIntersections(prep, M, j) for M in poss]
    push!(res, (prep, poss))
  end
  return res
end

function possible_ideals_for_cubics(G)
  bool, RR, sum_index, p = _has_pfr(G, 6)
  !bool && return Tuple{ProjRep, Vector{SymmetricCompleteIntersections}}[]
  @info "$(length(sum_index)) possible actions to consider"
  F = splitting_field(RR)
  S, _ = grade(PolynomialRing(F, "x" => 0:5)[1])
  _, j = homogeneous_component(S, 3)
  E = underlying_group(RR)
  Irr = irreducible_characters_underlying_group(RR)
  res = Tuple{ProjRep, Vector{SymmetricCompleteIntersections}}[]
  for l in sum_index
    @info "Test a character"
    chi = sum(Irr[l])
    detchi = determinant(chi)
    chid = symmetric_power(conj(chi), 3)
    ct = constituents(chid, 1)
    chis = Oscar.GAPGroupClassFunction[]
    for chi2 in ct
      detchi2 = determinant(chi2)
      if detchi == detchi2*detchi2
        push!(chis, chi2)
      end
    end
    length(chis) == 0 && continue
    @info "$(length(chis)) possibility.ies"
    rep = affording_representation(RR, chi)
    poss = CharacterGrassmannian[character_grassmannian(homogeneous_polynomial_representation(rep, 3), nu) for nu in chis]
    prep = ProjRep(rep, p)
    poss = SymmetricCompleteIntersections[SymmetricCompleteIntersections(prep, M, j) for M in poss]
    push!(res, (prep, poss))
  end
  return res
end

