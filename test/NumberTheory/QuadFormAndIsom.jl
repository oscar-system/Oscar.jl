Oscar.set_lwi_level(2)
set_verbosity_level(:ZZLatWithIsom, -1)

@testset "Printings" begin
  function _show_details(io::IO, X::Union{ZZLatWithIsom, QuadSpaceWithIsom})
    return show(io, MIME"text/plain"(), X)
  end
  # Finite order
  L = root_lattice(:A, 2)
  Lf = integer_lattice_with_isometry(L)
  Vf = ambient_space(Lf)
  for X in [Lf, Vf]
    @test sprint(_show_details, X) isa String
    @test sprint(Oscar.to_oscar, X) isa String
    @test sprint(show, X) isa String
    @test sprint(show, X; context=:supercompact => true) isa String
  end
  # Infinite order
  L = integer_lattice(; gram = QQ[1 2; 2 1])
  f = QQ[4 -1; 1 0]
  Lf = integer_lattice_with_isometry(L, f)
  Vf = ambient_space(Lf)
  for X in [Lf, Vf]
    @test sprint(_show_details, X) isa String
    @test sprint(Oscar.to_oscar, X) isa String
    @test sprint(show, X) isa String
    @test sprint(show, X; context=:supercompact => true) isa String
  end
end

@testset "Spaces with isometry" begin
  D5 = root_lattice(:D, 5)
  V = ambient_space(D5)
  Vf = @inferred quadratic_space_with_isometry(V; neg=true)
  @test is_one(-isometry(Vf))
  @test order_of_isometry(Vf) == 2
  @test space(Vf) === V

  for func in [rank, dim, gram_matrix, det, discriminant, is_positive_definite,
              is_negative_definite, is_definite, diagonal, signature_tuple]
    k = func(Vf)
    @test k == func(V)
  end

  @test evaluate(minimal_polynomial(Vf), -1) == 0
  @test evaluate(characteristic_polynomial(Vf), 0) == 1

  G = matrix(QQ, 6, 6 ,[3, 1, -1, 1, 0, 0, 1, 3, 1, 1, 1, 1, -1, 1, 3, 0, 0, 1, 1, 1, 0, 4, 2, 2, 0, 1, 0, 2, 4, 2, 0, 1, 1, 2, 2, 4])
  V = quadratic_space(QQ, G)
  f = matrix(QQ, 6, 6, [1 0 0 0 0 0; 0 0 -1 0 0 0; -1 1 -1 0 0 0; 0 0 0 1 0 -1; 0 0 0 0 0 -1; 0 0 0 0 1 -1])
  Vf = @inferred quadratic_space_with_isometry(V, f)
  @test order_of_isometry(Vf) == 3
  @test order_of_isometry(rescale(Vf, -3)) == 3
  L = lattice(rescale(Vf, -2))
  @test is_even(L)
  @test is_negative_definite(L) == is_positive_definite(Vf)
  @test rational_spinor_norm(Vf) > 0

  @test order_of_isometry(direct_sum(Vf, Vf, Vf)[1]) == 3

  @test Vf != quadratic_space_with_isometry(V; neg=true)
  @test length(unique([Vf, quadratic_space_with_isometry(V, isometry(Vf))])) == 1
  @test Vf^(order_of_isometry(Vf)+1) == Vf

  V = quadratic_space(QQ, matrix(QQ, 0, 0, []))
  Vf = @inferred quadratic_space_with_isometry(V)
  @test order_of_isometry(Vf) == -1
end

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

@testset "Enumeration of lattices with finite isometries" begin
  A4 = root_lattice(:A, 4)
  OA4 = matrix_group([
                      matrix(QQ, 4, 4, [-1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 -1]),
                      matrix(QQ, 4, 4, [1 1 1 1; 0 -1 -1 -1; 0 1 0 0; 0 0 1 0]),
                      matrix(QQ, 4, 4, [0 1 1 1; -1 -1 -1 -1; 1 1 0 0; 0 0 1 0]),
                      matrix(QQ, 4, 4, [1 0 0 0; -1 -1 -1 -1; 0 0 0 1; 0 0 1 0]),
                      matrix(QQ, 4, 4, [1 0 0 0; 0 1 1 1; 0 0 0 -1; 0 0 -1 0]),
                      matrix(QQ, 4, 4, [1 0 0 0; -1 -1 -1 0; 0 0 0 -1; 0 0 1 1]),
                      matrix(QQ, 4, 4, [1 0 0 0; 0 1 0 0; 0 0 1 1; 0 0 0 -1])
                    ])
  cc = conjugacy_classes(OA4)

  D = Oscar._test_isometry_enumeration(A4, 6)
  for n in keys(D)
    @test length(D[n]) == length(filter(c -> order(representative(c)) == n, cc))
  end
  
  N = rand(D[6])
  ONf = image_centralizer_in_Oq(integer_lattice_with_isometry(lattice(N), ambient_isometry(N)))[1]
    # for N, the image in OqN of the centralizer of fN in ON is directly
    # computing during the construction of the admissible primitive extension.
    # We compare if at least we obtain the same orders (we can't directly
    # compare the groups since they do not act exactly on the same module...
    # and isomorphism test might be slow)
  @test order(ONf) == order(image_centralizer_in_Oq(N)[1])

  E6 = root_lattice(:E, 6)
  @test length(enumerate_classes_of_lattices_with_isometry(E6, 20)) == 0
  @test length(enumerate_classes_of_lattices_with_isometry(E6, 9)) == 1
  @test length(enumerate_classes_of_lattices_with_isometry(genus(E6), 1)) == 1

  @test length(admissible_triples(E6, 2; IpA=[2])) == 2
  @test length(admissible_triples(rescale(E6, 2), 2; IpB=[4])) == 2
  @test length(admissible_triples(E6, 3; IpA=[2], IpB=[4])) == 1

  Zx, x = ZZ["x"]
  U = hyperbolic_plane_lattice()
  E8 = root_lattice(:E, 8)
  L = direct_sum(U, E8)[1]
  r = enumerate_classes_of_lattices_with_isometry(L, 30; char_poly=(x-1)^2*cyclotomic(30, x))
  @test length(r) == 1
  @test det(invariant_lattice(r[1])) == -1
  
  B = matrix(QQ, 10, 10 ,[2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3//2, 1//2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1//2, 1//2, 1//2, 1//2, 1, 0, 0, 0, 0, 0, 3//2, 1, 0, 1//2, 1//2, 1//2, 0, 0, 0, 0, 1//2, 1, 0, 1//2, 0, 0, 1//2, 0, 0, 0, 1, 1//2, 1//2, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 10, 10 ,[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -10]);
  M = integer_lattice(B, gram = G);
  r_local = enumerate_classes_of_lattices_with_isometry(M, 2; char_poly=(x-1)^4*(x+1)^6, _local=true)
  r_global = enumerate_classes_of_lattices_with_isometry(M, 2; char_poly=(x-1)^4*(x+1)^6, _local=false)
  @test length(r_local)==9
  @test length(r_global)==11

  M = direct_sum(U, U, U, U, U)[1]
  r = enumerate_classes_of_lattices_with_isometry(M, 4; min_poly=(x^2-1)*cyclotomic(4,x), pos_sigs=[(1,2), (2,1), (4,2)], neg_sigs=[(2,0)], fix_root=4)
  @test length(r) == 5
  @test length(filter(N -> det(kernel_lattice(N, 2)) == 4, r)) == 3
end

@testset "Enumeration of lattices with isometry of hermitian type" begin
  # Infinite isometry: chi is a Salem polynomial
  G = genus(torsion_quadratic_module(QQ[0;]), (9, 1))
  _, x = QQ[:x]
  chi = x^(10)+x^9-x^7-x^6-x^5-x^4-x^3+x+1
  rht = @inferred representatives_of_hermitian_type(G, chi)
  @test !isempty(rht)
  @test all(N -> !is_finite(order_of_isometry(N)), rht)

  # Galois orbits
  U = hyperbolic_plane_lattice()
  L = direct_sum(U, U, U, U)[1]
  @test length(representatives_of_hermitian_type(L, 5)) == 3
  @test length(representatives_of_hermitian_type(L, 5, 5)) == 2

  L = integer_lattice(; gram=matrix(QQ, 2, 2, [2 0; 0 2]))
  I = integer_lattice_with_isometry(L; neg=true)
  r = splitting(I, 2)
  @test is_of_hermitian_type(r[1])
end

@testset "Primitive extensions and embeddings" begin
  # Compute orbits of short vectors
  k = integer_lattice(; gram=matrix(QQ,1,1,[4]))
  E8 = root_lattice(:E, 8)
  ok, sv = primitive_embeddings(E8, k; classification=:sub)
  @test ok
  @test length(sv) == 1
  ok, sv = primitive_embeddings(rescale(E8, 2), rescale(k, QQ(1//2)); check=false)
  @test !ok
  @test is_empty(sv)
  @test isempty(primitive_embeddings(rescale(E8, -1), k; check=false)[2])

  k = integer_lattice(; gram=matrix(QQ,1,1,[6]))
  E7 = root_lattice(:E, 7)
  ok, sv = primitive_embeddings(E7, k; classification=:emb)
  @test ok
  @test length(sv) == 2
  q = discriminant_group(E7)
  p, z, n = signature_tuple(E7)
  ok, _ = primitive_embeddings(q, (p,n), E7; classification=:none)
  @test ok
  A5 = root_lattice(:A, 5)
  ok, sv = primitive_embeddings(A5, k)
  @test ok
  @test length(sv) == 2

  k = integer_lattice(; gram=matrix(QQ,1,1,[2]))
  ok, sv = primitive_embeddings(E7, k)
  @test ok
  @test length(sv) == 1

  ok, _ = primitive_embeddings(q, (p,n), E7; classification=:none)
  @test ok

  @test !primitive_embeddings(rescale(E7, 2), k; classification=:none, check=false)[1]

  B = matrix(QQ, 8, 8, [1 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0; 0 0 1 0 0 0 0 0; 0 0 0 1 0 0 0 0; 0 0 0 0 1 0 0 0; 0 0 0 0 0 1 0 0; 0 0 0 0 0 0 1 0; 0 0 0 0 0 0 0 1]);
  G = matrix(QQ, 8, 8, [-4 2 0 0 0 0 0 0; 2 -4 2 0 0 0 0 0; 0 2 -4 2 0 0 0 2; 0 0 2 -4 2 0 0 0; 0 0 0 2 -4 2 0 0; 0 0 0 0 2 -4 2 0; 0 0 0 0 0 2 -4 0; 0 0 2 0 0 0 0 -4]);
  L = integer_lattice(B; gram=G);
  f = matrix(QQ, 8, 8, [1 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0; 0 0 1 0 0 0 0 0; -2 -4 -6 -4 -3 -2 -1 -3; 2 4 6 5 4 3 2 3; -1 -2 -3 -3 -3 -2 -1 -1; 0 0 0 0 1 0 0 0; 1 2 3 3 2 1 0 2]);
  Lf = integer_lattice_with_isometry(L, f);
  F = invariant_lattice(Lf)
  C = coinvariant_lattice(Lf)
  reps = @inferred admissible_equivariant_primitive_extensions(F, C, Lf^0, 5)
  @test length(reps) == 1
  @test is_of_same_type(Lf, reps[1])
  _, reps = @inferred primitive_extensions(lattice(F), lattice(C); form_over=[discriminant_group(L)], classification=:embemb)
  @test length(reps) == 2
  @test reps[1] == (L, lattice(F), lattice(C))

  C = integer_lattice(; gram = QQ[-4 -2 -2  0  2 -2 -1 -3 -2 -1 -3 -2 -4  0 -2 -1;
                                  -2 -4 -1 -3  1 -1  1 -3 -1 -5 -3  2 -1  2 -2 -1;
                                  -2 -1 -4  0  2 -2 -1 -3 -2 -2 -3 -3 -4  0 -2  0;
                                   0 -3  0 -8  0 -2  4 -2  0 -7 -1  4  1  2  1  0;
                                   2  1  2  0 -4  0  0  2  2  1  1  2  3  0  1  0;
                                  -2 -1 -2 -2  0 -4  0 -4 -1 -2 -2 -2 -2  0  0 -1;
                                  -1  1 -1  4  0  0 -4  0 -1  3 -1 -3 -2 -1 -3 -1;
                                  -3 -3 -3 -2  2 -4  0 -8 -2 -4 -3 -2 -2  1 -1 -3;
                                  -2 -1 -2  0  2 -1 -1 -2 -4 -2 -3 -3 -4 -1 -3 -1;
                                  -1 -5 -2 -7  1 -2  3 -4 -2 -10 -4 3 -1  2 -2 -1;
                                  -3 -3 -3 -1  1 -2 -1 -3 -3 -4 -6 -1 -4  1 -5 -1;
                                  -2  2 -3  4  2 -2 -3 -2 -3  3 -1 -8 -5 -3 -1 -1;
                                  -4 -1 -4  1  3 -2 -2 -2 -4 -1 -4 -5 -8 -2 -4  0;
                                   0  2  0  2  0  0 -1  1 -1  2  1 -3 -2 -4  0 -1;
                                  -2 -2 -2  1  1  0 -3 -1 -3 -2 -5 -1 -4  0 -8 -2;
                                  -1 -1  0  0  0 -1 -1 -3 -1 -1 -1 -1  0 -1 -2 -4])

  L = integer_lattice(; gram = QQ[0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
                                1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
                                0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
                                0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
                                0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
                                0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
                                0 0 0 0 0 0 -2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
                                0 0 0 0 0 0 1 -2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
                                0 0 0 0 0 0 0 1 -2 1 0 0 0 1 0 0 0 0 0 0 0 0 0;
                                0 0 0 0 0 0 0 0 1 -2 1 0 0 0 0 0 0 0 0 0 0 0 0;
                                0 0 0 0 0 0 0 0 0 1 -2 1 0 0 0 0 0 0 0 0 0 0 0;
                                0 0 0 0 0 0 0 0 0 0 1 -2 1 0 0 0 0 0 0 0 0 0 0;
                                0 0 0 0 0 0 0 0 0 0 0 1 -2 0 0 0 0 0 0 0 0 0 0;
                                0 0 0 0 0 0 0 0 1 0 0 0 0 -2 0 0 0 0 0 0 0 0 0;
                                0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 1 0 0 0 0 0 0 0;
                                0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 -2 1 0 0 0 0 0 0;
                                0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 -2 1 0 0 0 1 0;
                                0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 -2 1 0 0 0 0;
                                0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 -2 1 0 0 0;
                                0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 -2 1 0 0;
                                0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 -2 0 0;
                                0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 -2 0;
                                0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4])

  ok, reps = primitive_embeddings(L, C; check=false)
  @test length(reps) == 1

  ## Odd case
  B = matrix(QQ, 5, 5 ,[1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1]);
  G = matrix(QQ, 5, 5 ,[3, 1, 0, 0, 0, 1, 3, 1, 1, -1, 0, 1, 3, 0, 0, 0, 1, 0, 3, 0, 0, -1, 0, 0, 3]);
  L = integer_lattice(B, gram = G);
  k = lattice_in_same_ambient_space(L, B[2:2, :])
  N = orthogonal_submodule(L, k)
  ok, reps = primitive_extensions(k, N; glue_order=[3])
  @test ok
  @test length(reps) == 2

  GL = genus(torsion_quadratic_module(QQ[1;]), (3, 3))
  A2 = root_lattice(:A, 2)
  ok, reps = primitive_embeddings(GL, A2)
  @test ok
  @test !is_even(reps[1][3])

  k = integer_lattice(; gram=QQ[4;])
  I = direct_sum(hyperbolic_plane_lattice(), k)[1]
  @test primitive_embeddings(I, k; classification=:none)[1]
end

@testset "Equivariant primitive extensions" begin
  A2 = root_lattice(:A, 2)
  Fs = integer_lattice_with_isometry(root_lattice(:E, 6); neg=false)
  q = torsion_quadratic_module(QQ[0;])

  ok, reps = equivariant_primitive_extensions(A2, Fs; form_over=[q], classification=:embemb)
  @test ok
  @test length(reps) == 6

  ok, reps = equivariant_primitive_extensions(reps[1][2], reps[1][3]; form_over=[q], classification=:embemb)
  @test ok
  @test length(reps) == 2

  Fs2 = integer_lattice_with_isometry(integer_lattice(;gram = QQ[3;]); neg=false)
  ok, reps = equivariant_primitive_extensions(A2, Fs2; even=false)

  @test ok
  @test length(reps) == 9
end

@testset "Fix hermitian miranda-morrison" begin
  B = matrix(QQ, 4, 4, [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]);
  G = matrix(QQ, 4, 4, [0 0 0 4; 0 0 4 0; 0 4 0 0; 4 0 0 0]);
  L = integer_lattice(B, gram=G);
  f = matrix(QQ, 4, 4, [1 0 2 0; 0 1 0 -2; -1 0 -1 0; 0 1 0 -1]);
  Lf = integer_lattice_with_isometry(L, f);
  GLf, _ = image_centralizer_in_Oq(Lf)
  @test order(GLf) == 96

  B = matrix(QQ, 6, 6, [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1]);
  G = matrix(QQ, 6, 6, [0 0 0 0 0 3; 0 0 0 0 3 0; 0 0 -6 0 0 0; 0 0 0 -6 0 0; 0 3 0 0 0 0; 3 0 0 0 0 0]);
  L = integer_lattice(B, gram=G);
  f = matrix(QQ, 6, 6, [0 1 -1 1 2 0; -1 -2 1 -1 0 -2; 0 0 0 1 2 -2; 2 2 -1 0 -2 2; -1 0 0 1 2 -1; 0 1 0 1 1 0]);
  Lf = integer_lattice_with_isometry(L, f);
  GLf, _ = image_centralizer_in_Oq(Lf)
  @test order(GLf) == 24192
end

@testset "Fix Galois action" begin
  U = hyperbolic_plane_lattice()
  L, _ = direct_sum(U, U)
  reps = representatives_of_hermitian_type(L, 4, 4)
  @test length(reps) == 1
  # The rest is for code coverage
  Lf = first(reps)
  Lf = rescale(Lf, 5)
  G, _ = image_centralizer_in_Oq(Lf)
  _, qLf = discriminant_group(Lf)
  @test all(g -> g*G(qLf) == G(qLf)*g, gens(G))
end

@testset "Fix type condition" begin
  B = matrix(QQ, 6, 6 ,[-1//2, 0, 1//2, 0, 0, 0, 0, 1//2, 0, 0, -1//2, 0, 0, 0, 0, -1//2, 0, 1//2, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 1, 0, -1, 0, -1, 0, 0, 0])
  G = matrix(QQ, 6, 6 ,[1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, -13, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, -91])
  L = integer_lattice(B, gram = G)
  r = enumerate_classes_of_lattices_with_isometry(L, 14; char_poly=cyclotomic_polynomial(14), fix_root=14)
  @test length(r) == 3
end

@testset "New functionalities" begin
  r = representatives_of_hermitian_type((2, 4), 3^5, cyclotomic_polynomial(9), 9)
  @test length(r) == 1
  Lf = only(r)
  @test order_of_isometry(Lf) == 9
  @test is_of_hermitian_type(Lf)
  G = genus(Lf)
  @test signature_pair(G) == (2, 4)
  @test det(G) == 3^5

  r = representatives_of_hermitian_type(Tuple{Int, Int}[(2,2), (2, 10), (2,18)], 1, cyclotomic_polynomial(12), 12)
  @test length(r) == 3
  @test all(Lf -> order_of_isometry(Lf) == 12, r)
  @test all(is_of_hermitian_type, r)
  @test all(is_unimodular, r)
  @test all(Lf -> signature_tuple(Lf)[1] == 2, r)

  r = representatives_of_hermitian_type(Tuple{Int, Int}[(2,2), (2,6)], 1:16, cyclotomic_polynomial(4), 4)
  @test length(r) == 8
  @test all(Lf -> det(Lf) in 1:16, r)
  @test all(Lf -> order_of_isometry(Lf) == 4, r)
  @test all(is_of_hermitian_type, r)
  @test count(Lf -> is_one(order(discriminant_group(Lf)[2])), r) == 4
end
