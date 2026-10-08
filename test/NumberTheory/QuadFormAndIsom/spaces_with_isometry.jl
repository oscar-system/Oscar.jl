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
