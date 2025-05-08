@testset "elliptic fibrations" begin
  B = matrix(QQ, 16, 16 ,[2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3//2, 1//2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1//2, 3//2, 3//2, 1//2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1//2, 3//2, 0, 1//2, 1//2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1//2, 1//2, 1//2, 0, 1//2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1//2, 1, 1//2, 0, 1//2, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1//2, 0, 1//2, 0, 3//5, 1//10]);
  G = matrix(QQ, 16, 16 ,[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -170]);
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
  B = matrix(QQ, 10, 10 ,[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1//3, 2//3, 1//3, 2//3, 2//3, 2//3, 1//3, 1//3, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]);
  G = matrix(QQ, 10, 10 ,[-2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, -2, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, -1, 1, -1, -1, -1, 0, 0, 0, 0, -1, -2, 1, -1, 0, -1, 0, 0, 0, 0, 1, 1, -2, 0, 0, 1, 0, 0, 0, 0, -1, -1, 0, -2, -1, -1, 0, 0, 0, 0, -1, 0, 0, -1, -2, -1, 0, 0, 0, 0, -1, -1, 1, -1, -1, -2]);
  L = integer_lattice(B, gram = G);

  B = matrix(QQ, 4, 10 ,[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]);
  G = matrix(QQ, 10, 10 ,[-2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, -2, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, -1, 1, -1, -1, -1, 0, 0, 0, 0, -1, -2, 1, -1, 0, -1, 0, 0, 0, 0, 1, 1, -2, 0, 0, 1, 0, 0, 0, 0, -1, -1, 0, -2, -1, -1, 0, 0, 0, 0, -1, 0, 0, -1, -2, -1, 0, 0, 0, 0, -1, -1, 1, -1, -1, -2]);
  S = integer_lattice(B, gram = G);

  weyl = QQ[31   61   52   71   5   -6   5   -2   -7   8]
  weylk3 = change_base_ring(ZZ,Oscar.solve(basis_matrix(L), weyl; side = :left))
  k3,_ = BorcherdsCtx(L, S, weylk3; compute_OR=false)
  walls1 = Oscar._walls_of_chamber(k3, weylk3)
  @test length(walls1)==4
  walls2 =  [
  ZZ[0   0   2   1],
  ZZ[1   1   1   2],
  ZZ[0   0   -1   -1],
  ZZ[-1   0   0   0]]
  @test issetequal(walls1, walls2)
  
  c = chamber(k3, weylk3, zero_matrix(ZZ,1,rank(S)), walls2)
  @test length(walls(c)) == 4
  @test Oscar.alg319(c.data.gramS, walls(c), walls(c), c.data.membership_test) == aut(c)

end

@testset "K3 surface automorphism groups" begin
  S = integer_lattice(gram=QQ[-2 1 0 0; 1 -2 1 1; 0 1 -2 1; 0 1 1 -2])
  @test 1 == Oscar.has_zero_entropy(S; rank_unimod=10)[1]
  _, k3aut, chambers, rational_mod_aut = borcherds_method(S, 10, compute_OR=true)
  @test order(matrix_group(k3aut))==2
  @test length(chambers) == 1
  @test length(rational_mod_aut) == 3
  @test length(rays(chambers[1])) == 4

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
  _, k3aut, chambers, rational_mod_aut = borcherds_method(S, 10; compute_OR=true)
  @test order(matrix_group(k3aut))==6
  @test length(chambers) == 1
  @test length(rational_mod_aut) == 4

  _, k3aut, chambers, rational_mod_aut = borcherds_method(S, 18, compute_OR=true)
  @test order(matrix_group(k3aut))==6
  @test length(chambers) == 1
  @test length(rational_mod_aut) == 4

  _, k3aut, chambers, rational_mod_aut = borcherds_method(S, 10; compute_OR=false)
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

@testset "find_section" begin
  L = integer_lattice(gram=ZZ[0 2 3;  2 -2 1; 3 1 -2])
  f = QQ[1 0 0;]
  s = Oscar.find_section(L,f)
  V = ambient_space(L)
  @test inner_product(V,f,s)==1
  @test inner_product(V,s,s)==-2
end

@testset "weyl_vector_non_degenerate" begin 
  B = matrix(QQ, 10, 10 ,[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1//3, 2//3, 1//3, 2//3, 2//3, 2//3, 1//3, 1//3]);
  G = matrix(QQ, 10, 10 ,[-2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, -2, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, -1, 1, -1, -1, -1, 0, 0, 0, 0, -1, -2, 1, -1, 0, -1, 0, 0, 0, 0, 1, 1, -2, 0, 0, 1, 0, 0, 0, 0, -1, -1, 0, -2, -1, -1, 0, 0, 0, 0, -1, 0, 0, -1, -2, -1, 0, 0, 0, 0, -1, -1, 1, -1, -1, -2]);
  L = integer_lattice(B, gram = G);
  B = matrix(QQ, 4, 10 ,[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]);
  G = matrix(QQ, 10, 10 ,[-2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, -2, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, -1, 1, -1, -1, -1, 0, 0, 0, 0, -1, -2, 1, -1, 0, -1, 0, 0, 0, 0, 1, 1, -2, 0, 0, 1, 0, 0, 0, 0, -1, -1, 0, -2, -1, -1, 0, 0, 0, 0, -1, 0, 0, -1, -2, -1, 0, 0, 0, 0, -1, -1, 1, -1, -1, -2]);
  S = integer_lattice(B, gram = G);
  u = QQ[90 157 218//3 346//3 -7//3 13//3 16//3 -11//3 -1//3 32//3]
  weyl = QQ[90 157 218//3 346//3 -7//3 13//3 16//3 -11//3 -1//3 32//3]
  ample0 = QQ[15 28 13 21 0 0 0 0 0 0]
  @test !Oscar.is_S_nondegenerate(L,S,weyl)
  weyl1, _,_ = Oscar.weyl_vector_non_degenerate(L,S,u,weyl,ample0)
  @test Oscar.is_S_nondegenerate(L,S,weyl1)

end 

@testset "preprocessing" begin
  g = -gram_matrix(root_lattice(:A,9))
  g[8,9] = g[9,8] = 0
  g[9,5] = g[5,9] = 1
  UE7 = g
  @assert det(UE7) == 2
  UE7 = integer_lattice(gram=UE7)
  L,S,w = Oscar.borcherds_method_preprocessing(rescale(UE7,2),18)
  @test (w*gram_matrix(rational_span(L))*transpose(w))[1,1] == 620
end 

@testset "#4804 is fixed" begin
  gramL = QQ[-1702038 -1020375 -1347502 698 392 392 -655 698 0 43 -612 -306 -306 349 -961 -698 1396 -961 -1616 1788 -263 741 698 1133 349 349; -1020375 -611722 -807837 418 234 234 -393 418 0 25 -368 -184 -184 209 -577 -418 836 -577 -970 1070 -159 443 418 677 209 209; -1347502 -807837 -1066826 552 309 309 -519 552 0 33 -486 -243 -243 276 -762 -552 1104 -762 -1281 1413 -210 585 552 894 276 276; 698 418 552 -2 -1 -1 1 -1 0 0 0 0 0 0 1 0 0 0 0 -1 1 1 -1 -1 0 1; 392 234 309 -1 -2 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 1 0 -1 -1 0 0; 392 234 309 -1 -1 -2 1 -1 0 0 0 0 0 0 0 0 0 0 0 -1 0 1 0 -1 0 0; -655 -393 -519 1 0 1 -2 1 0 0 0 0 0 0 -1 0 0 0 0 1 0 -1 0 0 0 0; 698 418 552 -1 0 -1 1 -2 0 0 0 0 0 0 1 0 0 0 0 -1 0 1 0 0 0 1; 0 0 0 0 0 0 0 0 -2 1 1 -1 1 1 -1 0 1 -1 -1 0 0 0 0 -1 1 1; 43 25 33 0 0 0 0 0 1 -2 -1 1 -1 -1 1 0 -1 0 1 0 0 0 0 0 0 -1; -612 -368 -486 0 0 0 0 0 1 -1 -2 0 -1 0 0 0 0 0 0 0 0 0 0 1 0 0; -306 -184 -243 0 0 0 0 0 -1 1 0 -2 1 1 -1 0 1 -1 -1 0 0 0 0 0 1 1; -306 -184 -243 0 0 0 0 0 1 -1 -1 1 -2 -1 1 0 0 0 0 0 0 0 0 0 0 -1; 349 209 276 0 0 0 0 0 1 -1 0 1 -1 -2 1 0 -1 0 1 0 0 0 0 0 0 -1; -961 -577 -762 1 0 0 -1 1 -1 1 0 -1 1 1 -4 0 0 0 -1 1 -1 0 1 0 1 0; -698 -418 -552 0 0 0 0 0 0 0 0 0 0 0 0 -2 1 -1 -1 1 1 1 -1 1 1 1; 1396 836 1104 0 0 0 0 0 1 -1 0 1 0 -1 0 1 -4 2 2 -1 0 -1 0 -1 -1 -2; -961 -577 -762 0 0 0 0 0 -1 0 0 -1 0 0 0 -1 2 -4 -2 0 1 1 0 1 2 2; -1616 -970 -1281 0 0 0 0 0 -1 1 0 -1 0 1 -1 -1 2 -2 -4 1 0 1 0 0 1 1; 1788 1070 1413 -1 -1 -1 1 -1 0 0 0 0 0 0 1 1 -1 0 1 -4 0 -1 -1 -1 -1 0; -263 -159 -210 1 1 0 0 0 0 0 0 0 0 0 -1 1 0 1 0 0 -4 -1 2 -1 -1 -1; 741 443 585 1 0 1 -1 1 0 0 0 0 0 0 0 1 -1 1 1 -1 -1 -4 0 -1 -1 -1; 698 418 552 -1 -1 0 0 0 0 0 0 0 0 0 1 -1 0 0 0 -1 2 0 -4 -1 0 0; 1133 677 894 -1 -1 -1 0 0 -1 0 1 0 0 0 0 1 -1 1 0 -1 -1 -1 -1 -6 0 -1; 349 209 276 0 0 0 0 0 1 0 0 1 0 0 1 1 -1 2 1 -1 -1 -1 0 0 -4 -1; 349 209 276 1 0 0 0 1 1 -1 0 1 -1 -1 0 1 -2 2 1 0 -1 -1 0 -1 -1 -6]
  L = integer_lattice(gram=gramL)
  S = lattice_in_same_ambient_space(L,QQ[27 0 0 -718 2461 -5151 -7004 8148 -2023 -1993 -5916 -1891 -2244 353 115 -2550 2457 -1126 -6715 5126 -1772 3613 2155 2478 -3825 -17; 6 27 0 -599 2019 -4242 -5768 6687 -1666 -1649 -4872 -1565 -1848 283 87 -2100 2008 -935 -5530 4206 -1467 2960 1767 2033 -3150 -14; 0 0 2 -43 144 -303 -412 477 -119 -118 -348 -112 -132 20 6 -150 143 -67 -395 300 -105 211 126 145 -225 -1])
  weyl = ZZ[-9299   -2023   -1649   315664   -1076618   2255939   3067477   -3564886   885992   874062   2590974   829387   982777   -153391   -49161   1116798   -1073654   494350   2940911   -2242571   777276   -1579935   -942594   -1084060   1675196   7446]
  k3 = Oscar.BorcherdsCtx(L, S, weyl;compute_OR=false)[1]
  c = chamber(k3, weyl)
  W = walls(c)
  @test !(ZZ[7096 1602 17739] in W) # not a wall 
  @test length(W)==4
  @test all( v in W for v in [ZZ[-5573 -1287 -14094], ZZ[-214 -54 -567], ZZ[-22 0 -27], ZZ[4511 1044 11421]])
end
