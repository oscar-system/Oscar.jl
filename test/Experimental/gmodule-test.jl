@testset "Experimental.gmodule" begin

  G = small_group(7*3, 1)
  z = Oscar.RepPc.reps(abelian_closure(QQ)[1], G)
  @test length(z) == 5

  z = irreducible_modules(G)
  @test length(z) == 5

  z = irreducible_modules(ZZ, G)
  @test length(z) == 5

  l = irreducible_modules(AnticNumberField, small_group(48, 17), minimal_degree = true)
  ds = degree.(base_ring.(l))
  @test length(l) == 12
  @test count(isequal(1), ds) == 8
  @test count(isequal(2), ds) == 4

  l = irreducible_modules(AnticNumberField, small_group(48, 29), minimal_degree = true)
  ds = degree.(base_ring.(l))
  @test length(l) == 8
  @test count(isequal(1), ds) == 6
  @test count(isequal(2), ds) == 2

  G = SL(2, 3)
  @test length(Oscar.RepPc.reps(abelian_closure(QQ)[1], G)) == 7
end


@testset "Experimental LocalH2" begin
  Qx, x = QQ["x"]
  k, a = number_field(x^6+108)
  
  G, mG = automorphism_group(k)
  for g = G
    for h = G
      @test mG(g) * mG(h) == mG(g*h)
    end
  end

  l2 = prime_decomposition(maximal_order(k), 2)
  k2, _ = Hecke.generic_completion(k, l2[1][1], 120)
  z = Hecke.local_fundamental_class_serre(k2, prime_field(k2))
  C, mG, mU = Oscar.GrpCoh.gmodule(k2, prime_field(k2))
  G = domain(mG)

  for g = G
    for h = G
      @test mG(g) * mG(h) == mG(g*h)
    end
  end

  q = Oscar.GrpCoh.cohomology_group(C, 2)
  @test order(q[1]) == 6 && is_cyclic(q[1])

  c = Oscar.GrpCoh.CoChain{2,elem_type(C.G),elem_type(C.M)}(C,
     Dict{NTuple{2, elem_type(C.G)}, elem_type(C.M)}(
       ((h), (g)) => preimage(mU, z(mG(g), mG(h))) for g = G for h = G))

  @test order(preimage(q[2], c)) == 6
end

