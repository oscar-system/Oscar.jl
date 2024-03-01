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

  l = irreducible_modules(AbsSimpleNumField, small_group(48, 17), minimal_degree = true)
  ds = degree.(base_ring.(l))
  @test length(l) == 12
  @test count(isequal(1), ds) == 8
  @test count(isequal(2), ds) == 4

  l = irreducible_modules(AbsSimpleNumField, small_group(48, 29), minimal_degree = true)
  ds = degree.(base_ring.(l))
  @test length(l) == 8
  @test count(isequal(1), ds) == 6
  @test count(isequal(2), ds) == 2

  G = SL(2, 3)
  @test length(Oscar.RepPc.reps(abelian_closure(QQ)[1], G)) == 7
  @test length(Oscar.RepPc.reps(GF(7, 6), G)) == 7
  @test length(Oscar.RepPc.reps(GF(2, 6), G)) == 3

  G = SL(2,3)
  F = GF(3,2)
  L = [
          matrix(F, [0 0 1; 1 0 0; 0 1 0]),
          matrix(F, [2 0 0; 0 1 0; 0 0 2]),
          matrix(F, [2 0 0; 0 2 0; 0 0 1]),
          matrix(F, [1 0 0; 0 1 0; 0 0 1])
      ]
  m = free_module(F, 3)
  M = GModule(m, G, [hom(m, m, a) for a in L])

  E = GF(3,6)
  phi = embed(F, E)
  LE = [
          matrix(E, [0 0 1; 1 0 0; 0 1 0]),
          matrix(E, [2 0 0; 0 1 0; 0 0 2]),
          matrix(E, [2 0 0; 0 2 0; 0 0 1]),
          matrix(E, [1 0 0; 0 1 0; 0 0 1])
      ]
  mE = free_module(E, 3)

  @test extension_of_scalars(M, phi) == GModule(mE, G, [hom(mE, mE, a) for a in LE])

  G = Oscar.GrpCoh.fp_group_with_isomorphism(gens(G))[1]
  q, mq = maximal_abelian_quotient(PcGroup, G)
  @test length(Oscar.RepPc.brueckner(mq)) == 24
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

@testset "Various" begin
  @test length(all_extensions(abelian_group(2), cyclic_group(PermGroup, 2))) == 2
  k, a = cyclotomic_field(5)
  zk = maximal_order(k)

  U, mU = Oscar.GaloisCohomology_Mod.units_mod_ideal(5*zk)
  A, mA = automorphism_group(PermGroup, k)

  C = gmodule(A, mU)
  U, mU = unit_group_fac_elem(zk)
  C = gmodule(A, mU)
  C = C ⊗ C
  @test elementary_divisors(C.M) == ZZRingElem[10, 10, 10, 0]

  C = C ⊕ C

  D, _ = Oscar.GModuleFromGap.ghom(C, C)
  @test dim(D) == 4

  C = gmodule(GF(5), C)
  i = indecomposition(C)
  @test length(i) == 8

  G = dihedral_group(8)
  z = irreducible_modules(G)
  @test dim((z[1] ⊕ z[2]) ⊗ z[3]) == 2

  k, a = quadratic_field(3)
  r, mr = ray_class_group(7*5*maximal_order(k), n_quo = 2)
  z = gmodule(automorphism_group(PermGroup, k)[1], mr)
end

@testset "H^3" begin
  # This is C = pi_2(X_P) for the standard presentation of Q_8
  # We know that H^3(G, C) = Z/|G|Z
  genss = [@perm((1,2,6,3)(4,8,5,7)), @perm((1,4,6,5)(2,7,3,8))];
  G, = sub(symmetric_group(8), genss);
  @test small_group_identification(G) == (8, 4)
  mats = [[0, -1, 0, 0, 0, 0, 1, 0, -1, 1, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 1, -1, 1, 0, -1, 0, 0, 0, -1, 0, 1, -1, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1]];
  F = free_abelian_group(7);
  M1, M2 = matrix(ZZ, 7, 7, mats[1]), matrix(ZZ, 7, 7, mats[2]);
  C = gmodule(G, [hom(F, F, M1), hom(F, F, M2)]);
  q = cohomology_group(C, 3)[1]
  @test order(q) == 8
  @test is_cyclic(q)
  C = gmodule(GF(5), C)
  i = indecomposition(C)
  @test length(i) == 5

  C, _ = Oscar.GModuleFromGap.ghom(C, C)
  @test dim(C) == 49
end
