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

  RR = @inferred small_group(E)
  chis = @inferred all_characters(RR, 2)
  @test length(chis) == 11

  reps = @inferred all_irreducible_representations(RR)
  chis = irreducible_characters_underlying_group(RR)
  @test length(reps) == 4
  @test all(i -> character_representation(RR, representation_mapping(reps[i])) == chis[i], 1:length(chis))

  # characters
end
