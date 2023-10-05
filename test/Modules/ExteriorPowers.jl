@testset "exterior powers of modules" begin
  I = Oscar.ordered_multi_index([1, 2, 5], 5)
  I1 = Oscar.ordered_multi_index([1], 5)
  I2 = Oscar.ordered_multi_index([2], 5)
  I5 = Oscar.ordered_multi_index([5], 5)
  J = Oscar.ordered_multi_index([3, 4], 5)
  K = Oscar.ordered_multi_index([2, 4], 5)
  @test "$(Oscar._wedge(I, J)[2])" == "0 < 1 < 2 < 3 < 4 < 5 <= 5"
  sign, ind = Oscar._wedge([I2, I5, I1])
  @test ind == I
  @test sign == 1
  sign, ind = Oscar._wedge([I2, I1, I5])
  @test ind == I
  @test sign == -1

  sign, ind = Oscar._wedge(J, K)
  @test sign == 0
  sign, ind = Oscar._wedge(I, K)
  @test sign == 0

  c = [ind for ind in Oscar.OrderedMultiIndexSet(3, 5)]
  @test length(c) == binomial(5, 3)

  @test all(x->c[Oscar.linear_index(x)] == x, Oscar.OrderedMultiIndexSet(3, 5))

  @test [Oscar.ordered_multi_index(k, 3, 5) for k in 1:binomial(5, 3)] == collect(Oscar.OrderedMultiIndexSet(3, 5))

  R, (x, y, u, v, w) = QQ[:x, :y, :u, :v, :w]
  F = FreeMod(R, 5)
  Fwedge3 = Oscar.exterior_power(F, 3)
  tmp = Oscar.exterior_power(F, 3, cached=false)
  @test tmp !== Fwedge3
  @test Fwedge3 === Oscar.exterior_power(F, 3)
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
  phi2_alt = Oscar.wedge_multiplication_map(Oscar.exterior_power(F, 1, cached=false), Oscar.exterior_power(F, 2, cached=false), v)
  @test domain(phi2) !== domain(phi2_alt)
  @test codomain(phi2) !== codomain(phi2_alt)
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

  dual_basis = Oscar.koszul_duals(gens(Fwedge1))
  tmp = [Oscar.wedge(u, v) for (u, v) in zip(dual_basis, gens(Fwedge1))]
  Fwedge5 = Oscar.exterior_power(F, 5)
  @test all(x->x==Fwedge5[1], tmp)

  dual_basis = Oscar.koszul_duals(gens(Fwedge2))
  tmp = [Oscar.wedge(u, v) for (u, v) in zip(dual_basis, gens(Fwedge2))]
  @test all(x->x==Fwedge5[1], tmp)
end

@testset "induced maps on exterior powers" begin
  R, (x, y, u, v, w) = QQ[:x, :y, :u, :v, :w]

  R5 = FreeMod(R, 5)
  R4 = FreeMod(R, 4)

  A = R[x y u v; 7*y u v w; u 5*v w x; v w x y; w 4*x y u]
  phi = hom(R5, R4, A)

  phi_2 = Oscar.induced_map_on_exterior_power(phi, 2)

  for ind in Oscar.OrderedMultiIndexSet(2, 5)
    imgs = [phi(R5[i]) for i in Oscar.indices(ind)]
    img = Oscar.wedge(imgs)
    @test img == phi_2(domain(phi_2)[Oscar.linear_index(ind)])
  end

  phi_3 = Oscar.induced_map_on_exterior_power(phi, 3)

  A3 = matrix(phi_3)
  for ind1 in Oscar.OrderedMultiIndexSet(3, 5)
    for ind2 in Oscar.OrderedMultiIndexSet(3, 4)
      @test A3[Oscar.linear_index(ind1), Oscar.linear_index(ind2)] == det(A[Oscar.indices(ind1), Oscar.indices(ind2)])
    end
  end

  psi = hom(R4, R5, transpose(A))
  psi_2 = Oscar.induced_map_on_exterior_power(psi, 2)
  @test compose(phi_2, psi_2) == Oscar.induced_map_on_exterior_power(compose(phi, psi), 2)
  psi_3 = Oscar.induced_map_on_exterior_power(psi, 3)
  @test compose(phi_3, psi_3) == Oscar.induced_map_on_exterior_power(compose(phi, psi), 3)
  psi_3_alt = Oscar.induced_map_on_exterior_power(psi, 3, domain = Oscar.exterior_power(R4, 3, cached=false), codomain = Oscar.exterior_power(R5, 3, cached=false))
  @test matrix(psi_3) == matrix(psi_3_alt)
  @test domain(psi_3) !== domain(psi_3_alt)
  @test codomain(psi_3) !== codomain(psi_3_alt)
end
