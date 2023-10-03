@testset "exterior powers of modules" begin
  I = Oscar.ordered_multi_index([1, 2, 5], 5)
  I1 = Oscar.ordered_multi_index([1], 5)
  I2 = Oscar.ordered_multi_index([2], 5)
  I5 = Oscar.ordered_multi_index([5], 5)
  J = Oscar.ordered_multi_index([3, 4], 5)
  K = Oscar.ordered_multi_index([2, 4], 5)
  Oscar.wedge(I, J)
  sign, ind = Oscar.wedge([I2, I5, I1])
  @test ind == I
  @test sign == 1
  sign, ind = Oscar.wedge([I2, I1, I5])
  @test ind == I
  @test sign == -1

  sign, ind = Oscar.wedge(J, K)
  @test sign == 0
  sign, ind = Oscar.wedge(I, K)
  @test sign == 0

  c = [ind for ind in Oscar.OrderedMultiIndexSet(3, 5)]
  @test length(c) == binomial(5, 3)

  @test all(x->c[Oscar.linear_index(x)] == x, Oscar.OrderedMultiIndexSet(3, 5))

  @test [Oscar.ordered_multi_index(k, 3, 5) for k in 1:binomial(5, 3)] == collect(Oscar.OrderedMultiIndexSet(3, 5))

  R, (x, y, u, v, w) = QQ[:x, :y, :u, :v, :w]
  F = FreeMod(R, 5)
  Fwedge3 = Oscar.exterior_power(F, 3)
  Fwedge1 = Oscar.exterior_power(F, 1)
  Fwedge2 = Oscar.exterior_power(F, 2)

  success, orig_mod, q = Oscar.is_exterior_power(Fwedge3)
  @test success
  @test orig_mod === F
  @test q == 3

  v = sum(gens(R)[i]*F[i] for i in 1:5)
  phi0 = Oscar.wedge_multiplication_map(exterior_power(F, 0), F, v)
  phi1 = Oscar.wedge_multiplication_map(F, exterior_power(F, 2), v)
  @test iszero(compose(phi0, phi1))
  @test image(phi0)[1] == kernel(phi1)[1]
  phi2 = Oscar.wedge_multiplication_map(Fwedge1, Fwedge2, v)
  phi3 = Oscar.wedge_multiplication_map(Fwedge2, Fwedge3, v)
  @test !iszero(phi2)
  @test !iszero(phi3)
  img, _ = image(phi2)
  ker, _ = kernel(phi3)
  @test img == ker
  @test iszero(compose(phi2, phi3))

  @test Oscar.wedge([F[1], F[2], F[4]]) == - Oscar.wedge([F[2], F[1], F[4]]) == Oscar.wedge([F[4], F[1], F[2]])

  K = koszul_complex(v)
  @test all(i->iszero(homology(K, i)), 1:5)

  for i in 1:5
    @test iszero(koszul_homology(v, i))
    @test iszero(koszul_homology(v, F, i))
  end

  @test !iszero(koszul_homology(v, 0))
  @test !iszero(koszul_homology(v, F, 0))
end

@testset "exterior powers of graded modules" begin
  R, (x, y, u, v, w) = QQ[:x, :y, :u, :v, :w]
  S, (x, y, u, v, w) = grade(R)
  F = graded_free_module(S, [1, 1, 1, 1, -2])
  Fwedge1 = Oscar.exterior_power(F, 1)
  Fwedge2 = Oscar.exterior_power(F, 2)
  Fwedge3 = Oscar.exterior_power(F, 3)
  @test is_graded(Fwedge3)

  S1 = graded_free_module(S, 1)
  I, inc_I = sub(S1, [f^3*S1[1] for f in gens(S)])

  Oscar.koszul_dual(Fwedge2[3])

  dual_basis = Oscar.koszul_dual(gens(Fwedge1))
  tmp = [Oscar.wedge(u, v) for (u, v) in zip(dual_basis, gens(Fwedge1))]
  Fwedge5 = Oscar.exterior_power(F, 5)
  @test all(x->x==Fwedge5[1], tmp)

  dual_basis = Oscar.koszul_dual(gens(Fwedge2))
  tmp = [Oscar.wedge(u, v) for (u, v) in zip(dual_basis, gens(Fwedge2))]
  @test all(x->x==Fwedge5[1], tmp)
end
