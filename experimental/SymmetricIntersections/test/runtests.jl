using Test
using Oscar

@testset "Elevators" begin
  W = sort(rand(1:5, 10))
  S, x = graded_polynomial_ring(QQ, 10; weights=W, cached=false)

  function weight(p)
    return degree(p).coeff[1]
  end

  for i in 1:5
    el = Oscar.elevator(x, weight, i)
    Oscar.number_of_elevations(el) == 0 && continue
    Si, SitoS = homogeneous_component(S, i)
    @test Oscar.number_of_elevations(el) == vector_space_dim(Si)
    @test Set([prod(x[l]) for l in el]) == Set(SitoS.(gens(Si)))
  end

  S, x = graded_polynomial_ring(QQ, 10; cached=false)
  el = Oscar.elevator(x, weight, 5, lbs=[2,0,0,0,0,0,0,0,0, 0])
  @test Oscar.number_of_elevations(el) == binomial(12, 3)
  @test Oscar.underlying_list(el) === x
  @test Oscar.degree_of_elevations(el) == 5
  @test Oscar.associated_function(el) === weight
  @test typeof(Oscar.underlying_iterator(el)) == SubObjectIterator{PointVector{ZZRingElem}}
end

@testset "Linear representations" begin
  E = small_group(8, 4)
  RR = @inferred representation_ring(E)
  @test base_field(RR) isa AbsSimpleNumField
  @test underlying_group(RR) === E
  @test order(Oscar.character_table_underlying_group(RR)) == 8
  @test all(is_irreducible, Oscar.irreducible_characters_underlying_group(RR))
  @test all(in(E), Oscar.generators_underlying_group(RR))

  chis = @inferred all_characters(RR, 2)
  @test length(chis) == 11

  RR, reps = @inferred all_irreducible_representations(E)
  @test all(is_irreducible, reps)
  chis = Oscar.irreducible_characters_underlying_group(RR)
  @test length(reps) == 5
  @test all(i -> character_representation(RR, representation_mapping(reps[i])) == chis[i], 1:length(chis))

  co = rand(1:2, 3)
  chi = sum([co[j]*chis[j] for j in 1:3])

  chid = symmetric_power(conj(chi), 4)
  cd = @inferred character_decomposition(chid)
  @test sum([c[1]*Int(degree(c[2])) for c in cd]) == Int(degree(chid))
  @test all(c -> is_isotypical(c[1]*c[2]), cd)
  cs = constituents(chid, 3)
  @test allunique(cs)
  @test all(nu -> is_constituent(chid, nu), cs)
  @test constituents(chid, Int(degree(chid))) == [chid]
  chidt = exterior_power(chid, 3)
  @test all(nu -> is_constituent(chidt, det(nu)), cs)

  rep = @inferred affording_representation(RR, chi)
  @test representation_ring(rep) === RR
  @test underlying_group(rep) === E
  @test character_representation(rep) == chi

  f = representation_mapping(rep)
  @test domain(f) === E
  @test codomain(f) isa MatrixGroup
  @test character_representation(RR, f) == chi

  n = @inferred dimension_representation(rep)
  mr = @inferred matrix_representation(rep)
  @test all(m -> isone(m^8), mr)
  @test all(m -> size(m) == (n,n), mr)
  @test !is_faithful(rep) || is_isomorphic(matrix_group(mr), E)

  repv = @inferred dual_representation(rep)
  @test dimension_representation(repv) == dimension_representation(rep)
  @test character_representation(RR, representation_mapping(repv)) == conj(chi)

  repv3 = @inferred symmetric_power_representation(repv, 2)
  @test dimension_representation(repv3) == binomial(1+dimension_representation(repv), 2)
  @test character_representation(RR, representation_mapping(repv3)) == symmetric_power(character_representation(repv), 2)
  cds = character_decomposition(repv3)
  B = @inferred basis_isotypical_component(repv3, cds[1][2])
  @test length(B) == cds[1][1]
  B = reduce(vcat, B)
  @test is_submodule(repv3, B)
  rB = @inferred action_on_submodule(repv3, B)
  @test character_representation(rB) == cds[1][1]*cds[1][2]
  @test is_isotypical_component(repv3, rB)
  @test is_equivalent(repv3, rB) || (length(cds) > 1)

  reph = @inferred homogeneous_polynomial_representation(rep, 2)
  @test dimension_representation(reph) == binomial(dimension_representation(rep) +1,2)
  @test character_representation(RR, representation_mapping(reph)) == symmetric_power(conj(chi), 2)
  ic = @inferred isotypical_components(reph)
  inj = [v[1] for v in values(ic)]
  proj = [v[2] for v in values(ic)]
  @test all(ii -> rank(ii) == nrows(ii), inj)
  @test all(pp -> rank(pp) == ncols(pp), proj)
  @test all(j -> isone(inj[j]*proj[j]), 1:length(inj))
  @test all(j -> all(k -> (j==k)||iszero(inj[j]*proj[k]), 1:length(inj)), 1:length(inj))
  B = reduce(vcat, inj[2,:])
  if !(rank(B) == 0)
    B2 = @inferred complement_submodule(reph, B)
    repQ, proj = @inferred quotient_representation(reph, B2)
    @test is_isotypical(repQ)
    @test is_constituent(character_representation(reph), character_representation(repQ))
  else
    @test is_isotypical(reph)
  end
end

@testset "Projective representations" begin
  fpr = faithful_projective_representations(small_group(16, 12), 4)
  p = Oscar.associated_schur_cover(fpr[1])
  RR = representation_ring_linear_lift(fpr[1])
  G = underlying_group(fpr[1])
  @test allunique(character_linear_lift, fpr)

  prep = rand(fpr)
  @test underlying_group(prep) === G
  @test Oscar.associated_schur_cover(prep) === p
  @test representation_ring_linear_lift(prep) === RR
  @test dimension_representation(prep) == 4

  chi = @inferred character_linear_lift(prep)
  @test is_projective(chi, p)
  @test is_irreducible(prep) == is_irreducible(chi)
  @test is_isotypical(prep) == is_isotypical(chi)

  rep = @inferred linear_lift(prep)
  @test is_faithful(prep)
  @test is_faithful(rep, p)

  @test count(prep2 -> is_similar(prep, prep2), fpr) == 1

  f = representation_mapping_linear_lift(prep)
  @test order(domain(f)) == 64
  @test codomain(f) isa MatrixGroup
  @test character_representation(RR, f) == chi

  n = @inferred dimension_representation(prep)
  mr = @inferred matrix_representation(prep)
  @test all(m -> size(m) == (n,n), mr)
  mg = matrix_group(mr)
  Z, _ = center(mg)
  els = filter(m -> is_diagonal(matrix(m)) && (length(eigenvalues(matrix(m))) == 1), collect(Z))
  Q, _ = quo(mg, els)
  @test is_isomorphic(G, Q)

  prepv = @inferred dual_representation(prep)
  @test Oscar.associated_schur_cover(prepv) === p
  @test representation_ring_linear_lift(prepv) === RR
  @test underlying_group(prepv) === G
  @test dimension_representation(prepv) == dimension_representation(prep)
  @test character_linear_lift(prepv) == conj(chi)

  prepsq = @inferred prep + prep
  @test Oscar.associated_schur_cover(prepsq) === p
  @test representation_ring_linear_lift(prepsq) === RR
  @test underlying_group(prepsq) === G
  @test dimension_representation(prepsq) == 2*dimension_representation(prep)
  @test character_linear_lift(prepsq) == 2*chi

  prepv3 = @inferred symmetric_power_representation(prepv, 2)
  @test Oscar.associated_schur_cover(prepv3) === p
  @test representation_ring_linear_lift(prepv3) === RR
  @test underlying_group(prepsq) === G
  @test dimension_representation(prepv3) == binomial(5, 2)
  @test character_linear_lift(prepv3) == symmetric_power(conj(chi), 2)
  cds = @inferred character_decomposition(prepv3)
  B = @inferred basis_isotypical_component(linear_lift(prepv3), cds[1][2])
  @test length(B) == cds[1][1]
  B = reduce(vcat, B)
  @test is_submodule(linear_lift(prepv3), B)
  rB = @inferred action_on_submodule(linear_lift(prepv3), B)
  prepB = @inferred projective_representation(RR, representation_mapping(rB), p)
  @test character_linear_lift(prepB) == cds[1][1]*cds[1][2]
  @test is_subrepresentation(prepv3, prepB)
  @test is_isotypical_component(prepv3, prepB)

  prep2 = @inferred exterior_power_representation(prep, 2)
  @test Oscar.associated_schur_cover(prep2) === p
  @test representation_ring_linear_lift(prep2) === RR
  @test underlying_group(prep2) === G
  @test dimension_representation(prep2) == 6
  @test character_linear_lift(prep2) == exterior_power(chi, 2)
  pd, B = @inferred to_equivalent_block_representation(prep2)
  @test is_similar(prep2, pd)
end

@testset "Symmetric Grassmannians" begin
  fpr = faithful_projective_representations(small_group(16, 8), 4)
  RR = representation_ring_linear_lift(fpr[1])

  prep = rand(fpr)
  preph = homogeneous_polynomial_representation(prep, 2)
  chi = character_linear_lift(preph)
  cd = character_decomposition(chi)
  nu = cd[1][1]*cd[1][2]
  @test is_isotypical(nu)
  isog = Oscar.isotypical_grassmannian(preph, nu)
  @test Oscar.submodule_character(isog) == nu
  @test Oscar.submodule_dimension(isog) == Int(degree(nu))
  @test Oscar.module_representation(isog) === linear_lift(preph)
  @test is_irreducible(isog)
  @test !is_empty(isog)
  Iisog = defining_ideal(isog)
  @test projective_dimension(isog) == dim(Iisog)-1
  std_el = standard_element(isog)
  std_el = reduce(vcat, std_el)
  @test is_submodule(linear_lift(preph), std_el)
  r = action_on_submodule(linear_lift(preph), std_el)
  @test character_representation(r) == nu
  cs = constituents(chi, 3)
  nu = rand(cs)
  cg = Oscar.character_grassmannian(preph, nu)
  cd = character_decomposition(cg)
  @test Oscar.submodule_character(cg) == nu
  @test Oscar.submodule_dimension(cg) == Int(degree(nu))
  @test Oscar.module_representation(cg) === linear_lift(preph)
  @test is_irreducible(cg)
  @test !is_empty(cg)
  @test_throws ErrorException defining_ideal(cg)
  d = 0
  for nu2 in cd
    isog = Oscar.isotypical_factor(cg, nu2)
    d += dim(defining_ideal(isog)) - 1
  end
  @test projective_dimension(cg) == d
  std_el = standard_element(cg)
  std_el = reduce(vcat, [reduce(vcat, b) for b in std_el])
  @test is_submodule(linear_lift(preph), std_el)
  r = action_on_submodule(linear_lift(preph), std_el)
  @test character_representation(r) == nu
  ig = Oscar.invariant_grassmannian(prep, 2)
  @test Oscar.submodule_dimension(ig) == 2
  @test Oscar.module_representation(ig) == linear_lift(prep)
  @test !is_irreducible(ig)
  @test !is_empty(ig)
  Iig = defining_ideal(ig)
  @test projective_dimension(ig) == dim(Iig) - 1
  for nu in constituents(character_linear_lift(prep), 2)
    cg = Oscar.irreducible_component(ig, nu)
    @test projective_dimension(cg) <= projective_dimension(ig)
  end

  ac = all_characters(RR, 1)
  prep2 = prep+prep
  l = rand(ac)
  dg = Oscar.determinant_grassmannian(prep2, l, 2)
  while is_empty(dg)
    l = rand(ac)
    dg = Oscar.determinant_grassmannian(prep2, l, 2)
  end
  @test Oscar.submodule_dimension(dg) == 2
  @test Oscar.module_representation(dg) == linear_lift(prep2)
  @test Oscar.submodule_determinant_character(dg) == l
  Idg = defining_ideal(dg)
  @test projective_dimension(dg) == dim(Idg) - 1
  @test is_irreducible(dg) == (length(primary_decomposition(Idg)) == 1)

  cg = rand(irreducible_components(dg))
  b = reduce(vcat, [reduce(vcat,b) for b in standard_element(cg)])
  nu = character_representation(action_on_submodule(linear_lift(prep2), b))
  @test det(nu) == l
end

@testset "Symmetric intersections" begin
  si = symmetric_intersections(small_group(21, 1), 4, 4, 1)
  (prep, symci) = rand(si)
  @test all(S -> projective_group_action(S) === prep, symci)
  S = symci[1]
  M = underlying_space_of_modules(S)
  chi = Oscar.submodule_character(M)
  fs, n = parametrization_data(S)[1]
  @test n == 1
  @test length(fs)-1 == projective_dimension(M)
  @test all(f -> f[1] isa MPolyDecRingElem, fs)
  @test all(f -> is_semi_invariant_polynomial(linear_lift(prep), f[1]), fs)
  f = rand(fs)[1]
  r = linear_representation(linear_lift(prep), f)
  @test character_representation(r) == chi
  std_el = standard_element(S)
  @test is_invariant_ideal(linear_lift(prep), std_el)

  si = symmetric_intersections(small_group(96, 204), 4, 2, 3)
  (prep, symci) = rand(si)
  Is = standard_element.(symci)
  Ms = underlying_space_of_modules.(symci)
  chis = Oscar.submodule_character.(Ms)
  @test all(I -> is_invariant_ideal(linear_lift(prep), I), Is)
  i = rand(1:length(Is))
  I = Is[i]
  r = linear_representation(linear_lift(prep), I)
  @test character_representation(r) == chis[i]
end

