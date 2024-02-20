@testset "elliptic fibrations" begin
  B = matrix(FlintQQ, 16, 16 ,[2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3//2, 1//2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1//2, 3//2, 3//2, 1//2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1//2, 3//2, 0, 1//2, 1//2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1//2, 1//2, 1//2, 0, 1//2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1//2, 1, 1//2, 0, 1//2, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1//2, 0, 1//2, 0, 3//5, 1//10]);
  G = matrix(FlintQQ, 16, 16 ,[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -170]);
  NS = integer_lattice(B, gram = G);
  V = ambient_space(NS)
  f =  QQFieldElem[2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
  r, torsion, ade = Oscar.fibration_type(NS,f)
  @test r == 4
  @test order(torsion)==1
  @test sort(ade)==[(:A, 1), (:A, 1), (:E, 8)]

  s = Oscar.find_section(NS, f)
  @test inner_product(V, s, s) == -2
  @test inner_product(V, s, f) == 1

  h = Oscar.ample_class(NS)
  @test inner_product(V, h, h )[1,1]> 0
  K = orthogonal_submodule(NS, lattice(V, h))
  @test length(short_vectors(rescale(K,-1), 2))==0
end


@testset "walls of chamber" begin
  S = integer_lattice(gram=QQ[-2 1 0 0; 1 -2 1 1; 0 1 -2 1; 0 1 1 -2])
  # fix an embedding
  B = matrix(FlintQQ, 10, 10 ,[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1//3, 2//3, 1//3, 2//3, 2//3, 2//3, 1//3, 1//3, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]);
  G = matrix(FlintQQ, 10, 10 ,[-2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, -2, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, -1, 1, -1, -1, -1, 0, 0, 0, 0, -1, -2, 1, -1, 0, -1, 0, 0, 0, 0, 1, 1, -2, 0, 0, 1, 0, 0, 0, 0, -1, -1, 0, -2, -1, -1, 0, 0, 0, 0, -1, 0, 0, -1, -2, -1, 0, 0, 0, 0, -1, -1, 1, -1, -1, -2]);
  L = integer_lattice(B, gram = G);

  B = matrix(FlintQQ, 4, 10 ,[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]);
  G = matrix(FlintQQ, 10, 10 ,[-2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, -2, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, -1, 1, -1, -1, -1, 0, 0, 0, 0, -1, -2, 1, -1, 0, -1, 0, 0, 0, 0, 1, 1, -2, 0, 0, 1, 0, 0, 0, 0, -1, -1, 0, -2, -1, -1, 0, 0, 0, 0, -1, 0, 0, -1, -2, -1, 0, 0, 0, 0, -1, -1, 1, -1, -1, -2]);
  S = integer_lattice(B, gram = G);

  weyl = QQ[31   61   52   71   5   -6   5   -2   -7   8]
  weylk3 = change_base_ring(ZZ,Oscar.solve(basis_matrix(L), weyl; side = :left))
  k3,_ = BorcherdsCtx(L, S, weylk3; compute_OR=false)
  walls = Oscar._walls_of_chamber(k3, weylk3)
  @test length(walls)==4
  walls1 =  [
  ZZ[0   0   2   1],
  ZZ[1   1   1   2],
  ZZ[0   0   -1   -1],
  ZZ[-1   0   0   0]]
  @test issetequal(walls, walls1)
end

@testset "K3 surface automorphism groups" begin
  S = integer_lattice(gram=QQ[-2 1 0 0; 1 -2 1 1; 0 1 -2 1; 0 1 1 -2])
  _, k3aut, chambers, rational_mod_aut = borcherds_method(S, 10, compute_OR=true)
  @test order(matrix_group(k3aut))==2
  @test length(chambers) == 1
  @test length(rational_mod_aut) == 3

  C = chambers[1]
  p = Oscar.inner_point(C)
  @test all((p*C.data.gramS*transpose(v))[1,1]>0 for v in walls(C))
  QQC = Oscar.span_in_S(C.data.L, C.data.S, change_base_ring(QQ, weyl_vector(C)))
  @test rank(QQC)==rank(S)
  wall1 = Oscar._walls_of_chamber(C.data,weyl_vector(C),:short)
  wall2 = Oscar._walls_of_chamber(C.data,weyl_vector(C),:close)
  @test Set(wall1)==Set(wall2)

  #=
  _, k3aut, chambers, rational_mod_aut = borcherds_method(S, 18, compute_OR=true)
  @test order(matrix_group(k3aut))==2
  @test length(chambers) == 1
  @test length(rational_mod_aut) == 3

  _, k3aut, chambers, rational_mod_aut = borcherds_method(S, 26, compute_OR=true)
  @test order(matrix_group(k3aut))==2
  @test length(chambers) == 1
  @test length(rational_mod_aut) == 3
  =#
  # Another example with finite automorphism group
  S,_ = direct_sum(integer_lattice(gram=ZZ[0 1; 1 -2]),rescale(root_lattice(:D,4),-1))
  _, k3aut, chambers, rational_mod_aut = borcherds_method(S, 10, compute_OR=true)
  @test order(matrix_group(k3aut))==6
  @test length(chambers) == 1
  @test length(rational_mod_aut) == 4

  #_, k3aut, chambers, rational_mod_aut = borcherds_method(S, 18, compute_OR=true)
  #@test order(matrix_group(k3aut))==6
  #@test length(chambers) == 1
  #@test length(rational_mod_aut) == 4

  _, k3aut, chambers, rational_mod_aut = borcherds_method(S, 10, compute_OR=false)
  @test length(k3aut)==0
  @test length(chambers) == 6
  @test length(rational_mod_aut) == 6

  #_, k3aut, chambers, rational_mod_aut = borcherds_method(S, 26, compute_OR=true)
  #@test order(matrix_group(k3aut))==6
  #@test length(chambers) == 1
  #@test length(rational_mod_aut) == 4

  S,_ = direct_sum(integer_lattice(gram=ZZ[0 1; 1 -2]),rescale(root_lattice(:D,4),-1))
  # The preprocessing of the following involves some random choices which can sometimes be suboptimal
  #  we test them separately
  S = integer_lattice(gram=gram_matrix(S))
  n = 26
  L, S, iS = embed_in_unimodular(S::ZZLat, 1, n-1,primitive=true,even=true)
  V = ambient_space(L)
  U = lattice(V,basis_matrix(S)[1:2,:])
  weyl, u0 = weyl_vector(L, U)

  @test iszero(weyl*gram_matrix(V)*transpose(weyl))
  # This is the a little expensive bit ... we leave it to the lower dimensional tests
  # weyl1, u, hh = Oscar.weyl_vector_non_degenerate(L, S, u0, weyl, h)
  weyl1 = QQ[80 30 -4 -14 -27 11 2 -12 -14 1 9 7 -1 -2 16 -12 4 7 11 6 -1 1 0 3 1 0]
  weyl2 = change_base_ring(ZZ, Oscar.solve(basis_matrix(L), weyl1; side = :left))
  _, k3aut, chambers, rational_mod_aut =borcherds_method(L, S, weyl2, compute_OR=true)

  @test order(matrix_group(k3aut))==6
  @test length(chambers) == 1
  @test length(rational_mod_aut) == 4



  # one with parabolic automorphism group

  # we hardcode the embedding of the following lattice
  # because length(chambers) depends on the embedding
  #=
  S, _ = direct_sum(integer_lattice(gram=ZZ[0 1; 1 -2]),integer_lattice(gram=ZZ[-50;]))
  k3aut, chambers, rational_mod_aut = borcherds_method(S, 10, compute_OR=true)
  @test length(k3aut)==2
  @test length(chambers) == 74
  @test length(rational_mod_aut) == 4
  =#

  #=
  t = [QQFieldElem[0 1 0 0 0 0 0 0 0 0; 1 -2 0 0 0 0 0 0 0 0; 0 0 -50 0 0 0 0 0 0 0; 0 0 0 -2 1 0 0 0 0 1; 0 0 0 1 -2 0 0 0 0 -1; 0 0 0 0 0 -2 -1 1 0 1; 0 0 0 0 0 -1 -2 0 0 0; 0 0 0 0 0 1 0 -2 0 -1; 0 0 0 0 0 0 0 0 -2 -1; 0 0 0 1 -1 1 0 -1 -1 -4],
   QQFieldElem[1 0 0 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0 0 0; 0 0 1//50 21//25 4//25 19//25 31//50 31//50 37//50 13//25; 0 0 0 1 0 0 0 0 0 0; 0 0 0 0 1 0 0 0 0 0; 0 0 0 0 0 1 0 0 0 0; 0 0 0 0 0 0 1 0 0 0; 0 0 0 0 0 0 0 1 0 0; 0 0 0 0 0 0 0 0 1 0; 0 0 0 0 0 0 0 0 0 1],
   QQFieldElem[1 0 0 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0 0 0; 0 0 1 0 0 0 0 0 0 0]]
  t = [matrix(QQ,i) for i in t]
  gramV, basisL, basisS = t
  V = quadratic_space(QQ, gramV)
  L = lattice(V, basisL)
  S = lattice(V, basisS)

  weyl = ZZ[76   30   275   -231   -45   -208   -172   -171   -203   -143]
  k3 = Oscar.BorcherdsCtx(L,S, weyl, false)

  _, k3aut, chambers, rational_mod_aut = borcherds_method(k3,entropy_abort=false)

  @test length(k3aut)==2
  @test length(chambers) == 127
  @test length(rational_mod_aut) == 4

  SS = integer_lattice(gram=gram_matrix(S))
  C = lattice(ambient_space(SS),Oscar._common_invariant(k3aut)[2])
  d = diagonal(rational_span(C))
  @test d[1] == 0 # a common invariant isotropic ray.
  =#
end

