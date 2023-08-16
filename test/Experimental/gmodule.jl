@testset "Experimental.gmodule" begin
  function testrels(G::Oscar.GAPGroup)
    ggens = gens(G)
    rels = relations(G)
    return all(p -> map_word(p[1], ggens) == map_word(p[2], ggens), rels)
  end

  G = dihedral_group(8)
  @test testrels(G)             # full pc group
  @test testrels(center(G)[1])  # subgroup of full pc group
  F = FPGroup(G)
  @test testrels(F)             # full f.p. group
  @test testrels(center(F)[1])  # subgroup of full f.p. group
  @test testrels(PermGroup(G))  # nothing of the above
end


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
  k, a = number_field(x^6+108, cached = false)
  
  G, mG = automorphism_group(k)
  for g = G
    for h = G
      @test mG(g) * mG(h) == mG(g*h)
    end
  end

  l2 = prime_decomposition(maximal_order(k), 2)
  k2, _ = Hecke.completion(k, l2[1][1], 120)
  z = Hecke.local_fundamental_class_serre(k2, prime_field(k2))
  C, mG, mU = Oscar.GrpCoh.gmodule(k2, prime_field(k2))
  G = domain(mG)

  pe = gen(k2)
  for g = G
    for h = G
      mu = (mG(g) * mG(h))(pe) - mG(g*h)(pe)
      @test mG(g) * mG(h) == mG(g*h)
    end
  end

  q = Oscar.GrpCoh.cohomology_group(C, 2)
  @test order(q[1]) == 6 && is_cyclic(q[1])

  c = Oscar.GrpCoh.CoChain{2,elem_type(C.G),elem_type(C.M)}(C,
     Dict{NTuple{2, elem_type(C.G)}, elem_type(C.M)}(
       ((g), (h)) => preimage(mU, z(mG(g), mG(h))) for g = G for h = G))

  @test Oscar.GrpCoh.istwo_cocycle(c)         

  @test order(preimage(q[2], c)) == 6
end

@testset "GlobalFundClass" begin
  k, a = quadratic_field(-101)
  H = hilbert_class_field(k)
  G, _ = galois_group(H, QQ)
  @test describe(G) == "D28"

  r = ray_class_field(49*maximal_order(k), n_quo = 7)
  G, _ = galois_group(r, QQ)
  @test describe(G) == "C7 x D14"

  lp = collect(keys(factor(7*maximal_order(k))))
  s = ray_class_field(lp[1])
  G, _ = galois_group(s, QQ)
  @test describe(G) == "C6 x D42"

  @test degree(normal_closure(s)) == 126
end

@testset "BrauerGroup" begin
  G = transitive_group(24, 201)
  T = character_table(G)
  R = gmodule(T[9])
  S = gmodule(CyclotomicField, R)

  @test schur_index(T[9]) == 2

  B, mB = relative_brauer_group(base_ring(S), character_field(S))
  b = B(S)

  C = grunwald_wang(Dict(2*ZZ => 2), Dict(complex_embeddings(QQ)[1] => 2))
  @test degree(C) == 2
  K, m1, m2 = compositum(base_ring(S), absolute_simple_field(number_field(C))[1])

  SS = Oscar.GModuleFromGap.gmodule_over(m2, gmodule(K, S))

  @test degree(base_ring(SS)) == 2
end
