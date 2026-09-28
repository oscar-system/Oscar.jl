@testset "Lattices with isometry" begin
  A3 = root_lattice(:A, 3)
  agg = QQMatrix[
                 matrix(QQ, 3, 3, [-1 0 0; 0 -1 0; 0 0 -1]),
                 matrix(QQ, 3, 3, [1 1 1; 0 -1 -1; 0 1 0]),
                 matrix(QQ, 3, 3, [0 1 1; -1 -1 -1; 1 1 0]),
                 matrix(QQ, 3, 3, [1 0 0; -1 -1 -1; 0 0 1]),
                 matrix(QQ, 3, 3, [1 0 0; 0 1 1; 0 0 -1])
                ]
  OA3 = matrix_group(agg)
  set_attribute!(A3, :isometry_group, OA3)
  f = agg[2]
  g = agg[4]

  L = integer_lattice(gram = matrix(QQ, 0, 0, []))
  Lf = integer_lattice_with_isometry(L; neg=true)
  @test order_of_isometry(Lf) == -1

  F, C = invariant_coinvariant_pair(Lf)
  @test rank(F) + rank(C) == rank(Lf)
  @test order_of_isometry(C) == order_of_isometry(Lf)
  @test is_one(isometry(F))

  @test is_primary_with_prime(integer_lattice_with_isometry(root_lattice(:E, 6)))[1]
  @test is_elementary_with_prime(integer_lattice_with_isometry(root_lattice(:E, 7)))[1]
  @test is_unimodular(integer_lattice_with_isometry(hyperbolic_plane_lattice()))

  L = @inferred integer_lattice_with_isometry(A3)
  @test is_primary(L, 2)
  @test !is_elementary(L, 2)
  @test length(unique([L, L, L])) == 1
  @test ambient_space(L) isa QuadSpaceWithIsom
  @test isone(isometry(L))
  @test isone(ambient_isometry(L))
  @test isone(order_of_isometry(L))
  @test order(image_centralizer_in_Oq(L)[1]) == 2
  @test rational_spinor_norm(L) > 0

  for func in [rank, genus, basis_matrix, is_positive_definite,
               gram_matrix, det, scale, norm, is_integral, is_negative_definite,
               degree, is_even, discriminant, signature_tuple, is_definite]
    k = @inferred func(L)
    @test k == func(A3)
  end

  LfQ = @inferred rational_span(L)
  @test LfQ isa QuadSpaceWithIsom
  @test evaluate(minimal_polynomial(L), 1) == 0
  @test evaluate(characteristic_polynomial(L), 0) == -1

  @test minimum(L) == 2
  @test is_positive_definite(L)
  @test is_definite(L)

  nf = multiplicative_order(f)
  @test_throws ArgumentError integer_lattice_with_isometry(A3, zero_matrix(QQ, 0, 0))

  L2 = @inferred integer_lattice_with_isometry(A3, f; ambient_representation=false)
  @test order_of_isometry(L2) == nf
  L2v = @inferred dual(L2)
  @test order_of_isometry(L2v) == nf
  @test ambient_isometry(L2v) == ambient_isometry(L2)
  
  L3 = @inferred integer_lattice_with_isometry(A3, g; ambient_representation=true)
  @test order_of_isometry(L3) == multiplicative_order(g)
  @test L3^(order_of_isometry(L3)+1) == L3
  @test genus(lll(L3; same_ambient=false)) == genus(L3)

  L4 = @inferred rescale(L3, QQ(1//4))
  @test !is_integral(L4)
  @test order_of_isometry(L4) == order_of_isometry(L3)
  @test_throws ArgumentError dual(L4)
  @test ambient_isometry(lll(L4)) == ambient_isometry(L4)

  @test rank(direct_sum(L2, L3)[1]) == rank(L2) + rank(L3)

  L5 = @inferred lattice(ambient_space(L2))
  @test (L2 == L5)

  B = matrix(QQ, 8, 8, [1 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0; 0 0 1 0 0 0 0 0; 0 0 0 1 0 0 0 0; 0 0 0 0 1 0 0 0; 0 0 0 0 0 1 0 0; 0 0 0 0 0 0 1 0; 0 0 0 0 0 0 0 1]);
  G = matrix(QQ, 8, 8, [-4 2 0 0 0 0 0 0; 2 -4 2 0 0 0 0 0; 0 2 -4 2 0 0 0 2; 0 0 2 -4 2 0 0 0; 0 0 0 2 -4 2 0 0; 0 0 0 0 2 -4 2 0; 0 0 0 0 0 2 -4 0; 0 0 2 0 0 0 0 -4]);
  L = integer_lattice(B; gram=G);
  f = matrix(QQ, 8, 8, [1 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0; 0 0 1 0 0 0 0 0; -2 -4 -6 -4 -3 -2 -1 -3; 2 4 6 5 4 3 2 3; -1 -2 -3 -3 -3 -2 -1 -1; 0 0 0 0 1 0 0 0; 1 2 3 3 2 1 0 2]);
  Lf = integer_lattice_with_isometry(L, f);

  GLf, _ = @inferred image_centralizer_in_Oq(Lf)
  @test order(GLf) == 600

  M = @inferred coinvariant_lattice(Lf)
  @test is_of_hermitian_type(M)
  H = hermitian_structure(M)
  @test H isa HermLat

  qL, fqL = @inferred discriminant_group(Lf)
  @test divides(ZZ(order_of_isometry(M)), order(fqL))[1]
  @test is_elementary(qL, 2)

  S = @inferred collect(values(signatures(M)))
  @test S[1] .+ S[2] == signature_pair(genus(M))

  @test rank(invariant_lattice(M)) == 0
  @test rank(invariant_lattice(Lf)) == rank(Lf) - rank(M)

  t = type(Lf)
  @test length(collect(keys(t))) == 2
  @test is_of_type(Lf, t)
  @test !is_of_same_type(Lf, M)
  @test is_hermitian(type(M))

  B = matrix(QQ, 4, 8, [0 0 0 0 3 0 0 0; 0 0 0 0 1 1 0 0; 0 0 0 0 1 0 1 0; 0 0 0 0 2 0 0 1]);
  G = matrix(QQ, 8, 8, [-2 1 0 0 0 0 0 0; 1 -2 0 0 0 0 0 0; 0 0 2 -1 0 0 0 0; 0 0 -1 2 0 0 0 0; 0 0 0 0 -2 -1 0 0; 0 0 0 0 -1 -2 0 0; 0 0 0 0 0 0 2 1; 0 0 0 0 0 0 1 2]);
  L = integer_lattice(B; gram=G);
  f = matrix(QQ, 8, 8, [1 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0; 0 0 1 0 0 0 0 0; 0 0 0 1 0 0 0 0; 0 0 0 0 0 1 0 0; 0 0 0 0 -1 1 0 0; 0 0 0 0 0 0 0 1; 0 0 0 0 0 0 -1 1]);
  Lf = integer_lattice_with_isometry(L, f);
  GL = image_centralizer_in_Oq(Lf)[1]
  @test order(GL) == 72

  B = matrix(QQ, 4, 6, [0 0 0 0 -2 1; 0 0 0 0 3 -4; 0 0 1 0 -1 0; 0 0 0 1 0 -1]);
  G = matrix(QQ, 6, 6, [2 1 0 0 0 0; 1 -2 0 0 0 0; 0 0 2//5 4//5 2//5 -1//5; 0 0 4//5 -2//5 -1//5 3//5; 0 0 2//5 -1//5 2//5 4//5; 0 0 -1//5 3//5 4//5 -2//5]);
  L = integer_lattice(B; gram=G);
  f = matrix(QQ, 6, 6, [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1; 0 0 -1 0 0 1; 0 0 0 -1 1 -1]);
  Lf = integer_lattice_with_isometry(L, f);
  GL = image_centralizer_in_Oq(Lf)[1]
  @test order(GL) == 2

  F, C, _ = invariant_coinvariant_pair(A3, OA3)
  @test rank(F) == 0
  @test C == A3
  _, _, G = invariant_coinvariant_pair(A3, OA3; ambient_representation=false)
  @test order(G) == order(OA3)
  C, _ = coinvariant_lattice(A3, sub(OA3, elem_type(OA3)[OA3(agg[2]), OA3(agg[4])])[1]; ambient_representation=false)
  @test is_sublattice(A3, C)

  B = matrix(QQ, 3, 3, [1 0 0; 0 1 0; 0 0 1]);
  G = matrix(QQ, 3, 3, [2 -1 0; -1 2 0; 0 0 -4]);
  L = integer_lattice(B, gram = G);
  f = matrix(QQ, 3, 3, [0 -1 0; 1 1 0; 0 0 1]);
  Lf = integer_lattice_with_isometry(L, f);
  @test is_bijective(image_centralizer_in_Oq(Lf)[2])
end
