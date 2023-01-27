@testset "Elevators" begin
  R, _ = PolynomialRing(QQ, 10, cached=false)
  W = sort(rand(1:5, 10))
  S, x = grade(R, W)
  
  function weight(p)
    return degree(p).coeff[1]
  end

  for i in 1:10
    el = elevator(x, weight, i)
    Si, SitoS = homogeneous_component(S, i)
    @test number_of_elevations(el) == dim(Si)
    @test Set([prod(x[l]) for l in el]) == Set(SitoS.(gens(Si))) 
  end

  S, x = grade(R)
  el = elevator(x, weight, 5, lbs=[2,0,0,0,0,0,0,0,0])
  @test number_of_elevations(el) == binomial(12, 3)
  @test underlying_list(el) === x
  @test degree_of_elevations(el) == 5
  @test associated_function(el) === weight
  @test typeof(underlying_iterator(el)) == SubObjectIterator{PointVector{fmpz}}
end

@testset "Representations" begin
  # representation rings
  E = small_group(8, 4)
  RR = @inferred representation_ring(E, small=false)
  @test splitting_field(RR) isa QQAbField
  @test underlying_group(RR) === E
  @test order(character_table_underlying_group(RR)) == 8
  @test all(chi -> is_irreducible(chi), irreducible_characters_underlying_group(RR))
  @test all(h -> h in E, generators_underlying_group(RR))

  RR = @inferred representation_ring(E)
  chis = @inferred all_characters(RR, 2)
  @test length(chis) == 11

  reps = @inferred all_irreducible_representations(RR)
  chis = irreducible_characters_underlying_group(RR)
  @test length(reps) == 5
  @test all(i -> character_representation(RR, representation_mapping(reps[i])) == chis[i], 1:length(chis))

  # characters
  chi = chis[end] + chis[1] + chis[2]
  chid = symmetric_power(conj(chi), 4)
  cd = @inferred character_decomposition(chid)
  @test sum([c[1]*Int(degree(c[2])) for c in cd]) == 35
  @test !is_isotypical(chid)
  @test all(c -> is_isotypical(c[1]*c[2]), cd)
  cs = @inferred constituents(chid, 3)
  @test allunique(cs)
  @test all(nu -> is_constituent(chid, nu), cs)
  @test constituents(chid, 35) == [chid]
  chidt = exterior_power(chid, 3)
  @test all(nu -> is_constituent(chidt, determinant(nu)), cs)
  Z, _ = @inferred center_of_character(chid)
  @test order(Z) == 1

  # linear representations
  rep = @inferred affording_representation(RR, chi)
  @test representation_ring(rep) === RR
  @test underlying_group(rep) === E
  @test character_representation(rep) == chi

  f = representation_mapping(rep)
  @test domain(f) === E
  @test codomain(f) isa MatrixGroup
  @test character_representation(RR, f) == chi
  @test all(r -> is_irreducible(r), reps)
  
  mr = @inferred matrix_representation(rep)
  @test all(m -> isone(m^8), mr)
  @test all(m -> size(m) == (4,4), mr)
  @test is_isomorphic(matrix_group(mr), E)
  @test is_faithful(rep)


end
