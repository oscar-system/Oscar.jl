using Oscar.SymInt

@testset "Elevators" begin
  R, _ = PolynomialRing(QQ, 10, cached=false)
  W = sort(rand(1:5, 10))
  S, x = grade(R, W)
  
  function weight(p)
    return degree(p).coeff[1]
  end

  for i in 1:10
    el = elevator(x, weight, i)
    number_of_elevations(el) == 0 && continue
    Si, SitoS = homogeneous_component(S, i)
    @test number_of_elevations(el) == dim(Si)
    @test Set([prod(x[l]) for l in el]) == Set(SitoS.(gens(Si))) 
  end

  S, x = grade(R)
  el = elevator(x, weight, 5, lbs=[2,0,0,0,0,0,0,0,0, 0])
  @test number_of_elevations(el) == binomial(12, 3)
  @test underlying_list(el) === x
  @test degree_of_elevations(el) == 5
  @test associated_function(el) === weight
  @test typeof(underlying_iterator(el)) == SubObjectIterator{PointVector{fmpz}}
end

@testset "Linear representations" begin
  E = small_group(8, 4)
  RR = @inferred representation_ring(E)
  @test base_field(RR) isa AnticNumberField
  @test underlying_group(RR) === E
  @test order(character_table_underlying_group(RR)) == 8
  @test all(chi -> is_irreducible(chi), irreducible_characters_underlying_group(RR))
  @test all(h -> h in E, generators_underlying_group(RR))

  chis = @inferred all_characters(RR, 2)
  @test length(chis) == 11

  RR, reps = @inferred all_irreducible_representations(E)
  @test all(r -> is_irreducible(r), reps)
  chis = irreducible_characters_underlying_group(RR)
  @test length(reps) == 5
  @test all(i -> character_representation(RR, representation_mapping(reps[i])) == chis[i], 1:length(chis))

  for i in 1:5
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
    @test all(nu -> is_constituent(chidt, determinant(nu)), cs)
    Z, _ = center_of_character(chid)
    @test is_subgroup(E, Z)[1]

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

    repv3 = @inferred symmetric_power_representation(repv, 3)
    @test dimension_representation(repv3) == binomial(2+dimension_representation(repv), 3)
    @test character_representation(RR, representation_mapping(repv3)) == symmetric_power(character_representation(repv), 3)
    cds = character_decomposition(repv3)
    B = @inferred basis_isotypical_component(repv3, cds[1][2])
    @test length(B) == cds[1][1]
    B = reduce(vcat, B)
    @test is_submodule(repv3, B)
    rB = @inferred action_on_submodule(repv3, B)
    @test character_representation(rB) == cds[1][1]*cds[1][2]
    @test is_isotypical_component(repv3, rB)
    @test is_equivalent(repv3, rB) || (length(cds) > 1)

    rep2 = @inferred exterior_power_representation(rep, 2)
    @test dimension_representation(rep2) == binomial(dimension_representation(rep), 2)
    @test character_representation(RR, representation_mapping(rep2)) == exterior_power(character_representation(rep), 2)
    rd, B = @inferred to_equivalent_block_representation(rep2)
    @test is_equivalent(rep2, rd)

    reph = @inferred homogeneous_polynomial_representation(rep, 2)
    @test dimension_representation(reph) == binomial(dimension_representation(rep) +1,2)
    @test character_representation(RR, representation_mapping(reph)) == symmetric_power(conj(chi), 2)
    ic = @inferred isotypical_components(reph)
    inj = [v[1] for v in collect(values(ic))]
    proj = [v[2] for v in collect(values(ic))]
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
end

@testset "Projective representations" begin
  fpr = faithful_projective_representations(32, 49, 4)
  p = associated_schur_cover(fpr[1])
  RR = representation_ring_linear_lift(fpr[1])
  G = underlying_group(fpr[1])
  @test allunique(character_linear_lift.(fpr))

  for i in 1:5
    prep = fpr[i]
    @test underlying_group(prep) === G
    @test associated_schur_cover(prep) === p
    @test representation_ring_linear_lift(prep) === RR
    @test dimension_linear_lift(prep) == 4

    chi = @inferred character_linear_lift(prep)
    @test is_projective(chi, p)
    @test is_irreducible(prep) == is_irreducible(chi)
    @test is_isotypical(prep) == is_isotypical(chi)

    rep = @inferred linear_lift(prep)
    @test is_faithful(prep)
    @test is_faithful(rep, p)

    @test count(prep2 -> is_similar(prep, prep2), fpr) == 1
   
    f = representation_mapping_linear_lift(prep)
    @test order(domain(f)) == 1024
    @test codomain(f) isa MatrixGroup
    @test character_representation(RR, f) == chi

    n = @inferred dimension_linear_lift(prep)
    mr = @inferred matrix_representation_linear_lift(prep)
    @test all(m -> isone(m^4), mr)
    @test all(m -> size(m) == (n,n), mr)
    mg = matrix_group(mr)
    Z, _ = center(mg)
    els = filter(m -> is_diagonal(matrix(m)) && (length(eigvals(matrix(m))) == 1), collect(Z))
    Q, _ = quo(mg, els)
    @test is_isomorphic(G, Q)

    prepv = @inferred dual_representation(prep)
    @test associated_schur_cover(prepv) === p
    @test representation_ring_linear_lift(prepv) === RR
    @test underlying_group(prepv) === G
    @test dimension_linear_lift(prepv) == dimension_linear_lift(prep)
    @test character_linear_lift(prepv) == conj(chi)

    prepsq = @inferred prep + prep
    @test associated_schur_cover(prepsq) === p
    @test representation_ring_linear_lift(prepsq) === RR
    @test underlying_group(prepsq) === G
    @test dimension_linear_lift(prepsq) == 2*dimension_linear_lift(prep)
    @test character_linear_lift(prepsq) == 2*chi

    prepv3 = @inferred symmetric_power_representation(prepv, 3)
    @test associated_schur_cover(prepv3) === p
    @test representation_ring_linear_lift(prepv3) === RR
    @test underlying_group(prepsq) === G
    @test dimension_linear_lift(prepv3) == binomial(6, 3)
    @test character_linear_lift(prepv3) == symmetric_power(conj(chi), 3)
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
    @test associated_schur_cover(prep2) === p
    @test representation_ring_linear_lift(prep2) === RR
    @test underlying_group(prep2) === G
    @test dimension_linear_lift(prep2) == 6
    @test character_linear_lift(prep2) == exterior_power(chi, 2)
    pd, B = @inferred to_equivalent_block_representation(prep2)
    @test is_similar(prep2, pd)

    preph = @inferred homogeneous_polynomial_representation(prep, 2)
    @test associated_schur_cover(preph) === p
    @test representation_ring_linear_lift(preph) === RR
    @test underlying_group(preph) === G
    @test dimension_linear_lift(preph) == 10
    @test character_linear_lift(preph) == symmetric_power(conj(chi), 2)
  end
end

@testset "Symmetric Grassmannians" begin
  fpr = faithful_projective_representations(16, 8, 4)
  RR = representation_ring_linear_lift(fpr[1])
  for i in 1:4
    prep = fpr[i]
    preph = homogeneous_polynomial_representation(prep, 4)
    chi = character_linear_lift(preph)
    cd = character_decomposition(chi)
    nu = cd[1][1]*cd[1][2]
    @test is_isotypical(nu)
    isog = @inferred isotypical_grassmannian(preph, nu)
    @test submodule_character(isog) == nu
    @test submodule_dimension(isog) == Int(degree(nu))
    @test module_representation(isog) === linear_lift(preph)
    @test is_irreducible(isog)
    @test !is_empty(isog)
    Iisog = @inferred defining_ideal(isog)
    @test projective_dimension(isog) == dim(Iisog)-1
    @test sprint(show, describe(isog)) isa String
    std_el = @inferred standard_element(isog)
    std_el = reduce(vcat, std_el)
    @test is_submodule(linear_lift(preph), std_el)
    r = action_on_submodule(linear_lift(preph), std_el)
    @test character_representation(r) == nu

    cs = constituents(chi, 8)
    nu = rand(cs)
    cg = @inferred character_grassmannian(preph, nu)
    cd = character_decomposition(cg)
    @test submodule_character(cg) == nu
    @test submodule_dimension(cg) == Int(degree(nu))
    @test module_representation(cg) === linear_lift(preph)
    @test is_irreducible(cg)
    @test !is_empty(cg)
    @test_throws ErrorException defining_ideal(cg)
    @test sprint(show, describe(cg)) isa String
    d = 0
    for nu2 in cd
      isog = @inferred isotypical_factor(cg, nu2)
      d += dim(defining_ideal(isog)) - 1
    end
    @test projective_dimension(cg) == d
    std_el = standard_element(cg)
    std_el = reduce(vcat, [reduce(vcat, b) for b in std_el])
    @test is_submodule(linear_lift(preph), std_el)
    r = action_on_submodule(linear_lift(preph), std_el)
    @test character_representation(r) == nu

    ig = @inferred invariant_grassmannian(prep, 2)
    @test submodule_dimension(ig) == 2
    @test module_representation(ig) == linear_lift(prep)
    @test !is_irreducible(ig)
    @test !is_empty(ig)
    Iig = @inferred defining_ideal(ig)
    @test projective_dimension(ig) == dim(Iig) - 1
    for nu in constituents(character_linear_lift(prep), 2)
      cg = @inferred irreducible_component(ig, nu)
      @test projective_dimension(cg) <= projective_dimension(ig)
    end

    ac = all_characters(RR, 1)
    prep2 = prep+prep
    dgs = DeterminantGrassmannian[determinant_grassmannian(prep2, l, 2) for l in ac]
    for j in 1:length(dgs)
      dg = dgs[j]
      @test submodule_dimension(dg) == 2
      @test module_representation(dg) == linear_lift(prep2)
      @test submodule_determinant_character(dg) == ac[j]
      is_empty(dg) && continue
      Idg = @inferred defining_ideal(dg)
      @test projective_dimension(dg) == dim(Idg) - 1
      @test is_irreducible(dg) == (length(primary_decomposition(Idg)) == 1)
      for cg in irreducible_components(dg)
        b = reduce(vcat, [reduce(vcat,b) for b in standard_element(cg)])
        nu = character_representation(action_on_submodule(linear_lift(prep2), b))
        @test determinant(nu) == ac[j]
      end
    end
  end
end

@testset "Symmetric intersections" begin
  si = symmetric_intersections((21, 1), 4, 4, 1)
  for (prep, symci) in si
    @test all(S -> projective_group_action(S) === prep, symci)
    S = symci[1]
    M = @inferred underlying_moduli_space_of_modules(S)
    chi = submodule_character(M)
    fs, n = parametrization_data(S)[1]
    @test n == 1
    @test length(fs)-1 == projective_dimension(M)
    @test all(f -> f[1] isa MPolyElem_dec, fs)
    @test all(f -> is_semi_invariant_polynomial(linear_lift(prep), f[1]), fs)
    for f in fs
      f = f[1]
      r = @inferred action_on_polynomial(linear_lift(prep), f)
      @test character_representation(r) == chi
    end
    std_el = @inferred standard_element(S)
    @test is_invariant_ideal(linear_lift(prep), std_el)
  end

  si = symmetric_intersections((96, 204), 4, 3, 4)
  for (prep, symci) in si
    Is = standard_element.(symci)
    Ms = underlying_moduli_space_of_modules.(symci)
    chis = submodule_character.(Ms)
    @test all(I -> is_invariant_ideal(linear_lift(prep), I), Is)
    for i in 1:length(Is)
      I = Is[i]
      r = @inferred action_on_ideal(linear_lift(prep), I)
      @test character_representation(r) == chis[i]
    end
  end
end


