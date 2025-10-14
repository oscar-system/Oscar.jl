
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
  
  r = enumerate_classes_of_lattices_with_isometry(E8, 12; char_poly=x^8 - x^6 + x^2 - 1)
  @test length(r) == 1
  r = enumerate_classes_of_lattices_with_isometry(E8, 12; char_poly=x^8 + x^6 + x^2 + 1)
  @test length(r) == 1

  
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
