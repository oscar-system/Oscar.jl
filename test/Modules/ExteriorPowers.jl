@testset "exterior powers of modules" begin
  R, (x, y, u, v, w) = QQ[:x, :y, :u, :v, :w]
  F = FreeMod(R, 5)
  Fwedge3, _ = Oscar.exterior_power(F, 3)
  tmp, _ = Oscar.exterior_power(F, 3, cached=false)
  @test tmp !== Fwedge3
  @test Fwedge3 === Oscar.exterior_power(F, 3)[1]
  Fwedge1, _ = Oscar.exterior_power(F, 1)
  Fwedge2, _ = Oscar.exterior_power(F, 2)

  success, orig_mod, q = Oscar._is_exterior_power(Fwedge3)
  @test success
  @test orig_mod === F
  @test q == 3

  v = sum(gens(R)[i]*F[i] for i in 1:5)
  phi0 = Oscar.wedge_multiplication_map(exterior_power(F, 0)[1], F, v)
  phi1 = Oscar.wedge_multiplication_map(F, exterior_power(F, 2)[1], v)
  @test iszero(compose(phi0, phi1))
  @test image(phi0)[1] == kernel(phi1)[1]
  phi2 = Oscar.wedge_multiplication_map(Fwedge1, Fwedge2, v)
  phi2_alt = Oscar.wedge_multiplication_map(Oscar.exterior_power(F, 1, cached=false)[1], Oscar.exterior_power(F, 2, cached=false)[1], v)
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
  S, _ = graded_polynomial_ring(QQ, 5)
  F = graded_free_module(S, [1, 1, 1, 1, -2])
  Fwedge1, _ = Oscar.exterior_power(F, 1)
  Fwedge2, _ = Oscar.exterior_power(F, 2)
  Fwedge3, _ = Oscar.exterior_power(F, 3)
  @test is_graded(Fwedge3)

  S1 = graded_free_module(S, 1)
  I, inc_I = sub(S1, [f^3*S1[1] for f in gens(S)])

  Oscar.koszul_dual(Fwedge2[3])

  dual_basis = Oscar.koszul_duals(gens(Fwedge1))
  tmp = [Oscar.wedge(u, v) for (u, v) in zip(dual_basis, gens(Fwedge1))]
  Fwedge5, _ = Oscar.exterior_power(F, 5)
  @test all(==(Fwedge5[1]), tmp)

  dual_basis = Oscar.koszul_duals(gens(Fwedge2))
  tmp = [Oscar.wedge(u, v) for (u, v) in zip(dual_basis, gens(Fwedge2))]
  @test all(==(Fwedge5[1]), tmp)
end

@testset "induced maps on exterior powers" begin
  R, (x, y, u, v, w) = QQ[:x, :y, :u, :v, :w]

  R5 = FreeMod(R, 5)
  R4 = FreeMod(R, 4)

  A = R[x y u v; 7*y u v w; u 5*v w x; v w x y; w 4*x y u]
  phi = hom(R5, R4, A)

  phi_2 = Oscar.induced_map_on_exterior_power(phi, 2)

  for ind in combinations(5, 2)
    imgs = [phi(R5[i]) for i in ind]
    img = Oscar.wedge(imgs)
    @test img == phi_2(domain(phi_2)[Oscar.linear_index(ind, 5)])
  end

  phi_3 = Oscar.induced_map_on_exterior_power(phi, 3)

  A3 = matrix(phi_3)
  for ind1 in combinations(5, 3)
    for ind2 in combinations(4, 3)
      @test A3[Oscar.linear_index(ind1, 5), Oscar.linear_index(ind2, 4)] == det(A[data(ind1), data(ind2)])
    end
  end

  psi = hom(R4, R5, transpose(A))
  psi_2 = Oscar.induced_map_on_exterior_power(psi, 2)
  @test compose(phi_2, psi_2) == Oscar.induced_map_on_exterior_power(compose(phi, psi), 2)
  psi_3 = Oscar.induced_map_on_exterior_power(psi, 3)
  @test compose(phi_3, psi_3) == Oscar.induced_map_on_exterior_power(compose(phi, psi), 3)
  psi_3_alt = Oscar.induced_map_on_exterior_power(psi, 3, domain = Oscar.exterior_power(R4, 3, cached=false)[1], codomain = Oscar.exterior_power(R5, 3, cached=false)[1])
  @test matrix(psi_3) == matrix(psi_3_alt)
  @test domain(psi_3) !== domain(psi_3_alt)
  @test codomain(psi_3) !== codomain(psi_3_alt)

  psi_3_alt = hom(domain(psi_3), codomain(psi_3), psi)
  @test psi_3_alt == psi_3
end

@testset "multiplication map" begin
  R, (x, y, u, v, w) = QQ[:x, :y, :u, :v, :w]

  F = FreeMod(R, 5)
  F3, mm = Oscar.exterior_power(F, 3)
  v = (F[1], F[3], F[4])
  u = (F[1], F[4], F[3])
  @test mm(v) == -mm(u)
  w = mm(v)
  @test Oscar.wedge_pure_function(F3)(v) == w
  @test Oscar.wedge_generator_decompose_function(F3)(w) == v
  @test preimage(mm, w) == v
end

@testset "printing" begin
  R, (x, y, u, v, w) = QQ[:x, :y, :u, :v, :w]

  F = FreeMod(R, 5)
  F3, mm = Oscar.exterior_power(F, 3)
  v = (F[1], F[3], F[4])
  u = (F[1], F[4], F[3])

  @test "$(mm(v))" == "e[1]^e[3]^e[4]"
  #@test "$(F3)" == "⋀^3(Free module of rank 5 over multivariate polynomial ring in 5 variables over QQ)"
  @test "$(F3)" == "3rd exterior power of Free module of rank 5 over multivariate polynomial ring in 5 variables over QQ"

  eu = sum(f*e for (f, e) in zip(gens(R), gens(F)))
  K = koszul_complex(eu)
  for i in 1:5
    #@test "$(K[i])" == "⋀^$(5-i)($F)"
    @test "$(K[i])" == "$(Oscar.ordinal_number_string(5-i)) exterior power of $F"
  end
end

@testset "exterior powers of subquos" begin
  R, (x, y, z, w) = QQ[:x, :y, :z, :w]

  M = R[x y z; y-1 z-2 w]

  I = ideal(R, minors(M, 2))

  A, pr = quo(R, I)

  phi = hom(FreeMod(A, 2), FreeMod(A, 3), change_base_ring(A, M))

  MM = cokernel(phi)

  MM2, mm2 = exterior_power(MM, 2)
  MM1, mm1 = exterior_power(MM, 1)
  MM0, mm0 = exterior_power(MM, 0)
  MM3, mm3 = exterior_power(MM, 3)

  @test iszero(MM3)

  (u, v, w) = gens(MM1)
  uv = wedge(u, v)
  (U, V, W) = gens(MM)
  @test uv == mm2((U, V))
  @test (U, V) == preimage(mm2, uv)
  uv = wedge(u, v)
  uw = wedge(u, w)
  vw = wedge(v, w)
  psi = hom(FreeMod(A, 3), MM2, [uv, uw, vw])
  @test !iszero(kernel(psi)[1])
  @test parent(wedge(u, v)) === MM2

  d02 = Oscar.wedge_multiplication_map(MM0, MM2, uv + vw)
  d01 = Oscar.wedge_multiplication_map(MM0, MM1, u + 2*v - w)
  d12 = Oscar.wedge_multiplication_map(MM1, MM2, u + 2*v - w)
  @test iszero(compose(d01, d12))
  @test iszero(wedge([u, v, w]))
end
