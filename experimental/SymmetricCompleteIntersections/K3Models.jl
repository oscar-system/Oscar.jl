export possible_ideals_for_cubics,
       possible_ideals_for_k3_models

function _subgroup_trivial_action_on_H20(RR::RepRing, chi1::Oscar.GAPGroupClassFunction, chi2::Oscar.GAPGroupClassFunction)
  h =  conjugacy_classes(underlying_group(RR))
  dchi1 = determinant(chi1)
  dchi2 = determinant(chi2)
  ker = elem_type(underlying_group(RR))[]
  for e in h
    if dchi1(representative(e)) == dchi2(representative(e))
      append!(ker, collect(e))
    end
  end
  ker, _ = sub(underlying_group(RR), ker)
  return ker
end

function _subgroup_trivial_action_on_H31(RR::RepRing, chi1::Oscar.GAPGroupClassFunction, chi2::Oscar.GAPGroupClassFunction)
  h = conjugacy_classes(underlying_group(RR))
  dchi1 = determinant(chi1)
  dchi2 = determinant(chi2)
  ker = elem_type(underlying_group(RR))[]
  for e in h
    if dchi1(representative(e)) == dchi2(representative(e))^2
      append!(ker, collect(e))
    end
  end
  ker, _ = sub(underlying_group(RR), ker)
  return ker
end

function _id_symplectic_subgroup_k3(RR::RepRing, chi1::Oscar.GAPGroupClassFunction, chi2::Oscar.GAPGroupClassFunction, p::GAPGroupHomomorphism)
  ker = _subgroup_trivial_action_on_H20(RR, chi1, chi2)
  return small_group_identification(p(ker)[1])
end

function _id_symplectic_subgroup_cubics(RR::RepRing, chi1::Oscar.GAPGroupClassFunction, chi2::GAPGroupClassFunction, p::GAPGroupHomomorphism)
  ker = _subgroup_trivial_action_on_H31(RR, chi1, chi2)
  return small_group_identification(p(ker)[1])
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
    chid = symmetric_power(conj(chi), d)
    ct = constituents(chid, t)
    chis = Oscar.GAPGroupClassFunction[]
    for chi2 in ct
      if _id_symplectic_subgroup_k3(RR, chi, chi2, p) == Gs
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

function possible_ideals_for_cubics(G::Tuple{Int, Int})
  bool, RR, sum_index, p = _has_pfr(small_group(G[1], G[2]), 6)
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
    chid = symmetric_power(conj(chi), 3)
    ct = constituents(chid, 1)
    chis = Oscar.GAPGroupClassFunction[]
    for chi2 in ct
      if _id_symplectic_subgroup_cubics(RR, chi, chi2, p) == G
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

