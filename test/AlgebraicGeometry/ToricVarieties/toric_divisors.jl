@testset "Torus-invariant divisors" begin
  F5 = hirzebruch_surface(NormalToricVariety, 5)
  dP3 = del_pezzo_surface(NormalToricVariety, 3)
  P2 = projective_space(NormalToricVariety, 2)

  D=toric_divisor(F5, [0, 0, 0, 0])
  D2 = divisor_of_character(F5, [1, 2])
  D3 = toric_divisor(dP3, [1, 0, 0, 0, 0, 0])
  D4 = canonical_divisor(dP3)
  D5 = anticanonical_divisor(dP3)
  D6 = trivial_divisor(dP3)
  D7 = toric_divisor(P2, [1, 1, 0])
  D8 = toric_divisor(P2, [1, -1, 0])
  D9 = toric_divisor(P2, [0, -1, 0])

  @testset "Should fail" begin
    @test_throws ArgumentError toric_divisor(F5, [0, 0, 0])
    @test_throws ArgumentError D+D3
    @test_throws ArgumentError D-D3
    @test_throws ArgumentError D==D3
  end

  @testset "Basic properties" begin
    @test is_prime(D) == false
    @test is_prime(D2) == false
    @test is_prime(D3) == true
    @test is_cartier(D) == true
    @test is_cartier(D2) == true
    @test is_principal(D) == true
    @test is_principal(D2) == true
    @test is_trivial(D) == true
    @test is_trivial(D2) == false
    @test is_basepoint_free(D) == true
    @test is_basepoint_free(D2) == true
    @test is_ample(D) == false
    @test is_ample(D2) == false
    @test is_very_ample(D) == false
    @test is_very_ample(D2) == false
    @test is_nef(D) == true
    @test is_nef(D2) == true
    @test is_integral(D) == true
    @test is_integral(D2) == true
    @test is_q_cartier(D) == true
    @test is_q_cartier(D2) == true
    @test is_prime(D) == false
    @test is_prime(D2) == false
    @test is_effective(D7) == true
    @test is_effective(D8) == false
  end

  @testset "Basic attributes" begin
    @test coefficients(D) == [0, 0, 0, 0]
    @test coefficients(D2) == [1, 2, 9, -2]
    @test dim(toric_variety(D)) == 2
    @test dim(polyhedron(D)) == 0
    @test ambient_dim(polyhedron(D)) == 2
  end

  @testset "Arithmetic" begin
    @test (D == D2) == false
    @test (D4 + D5 == D6) == true
    @test is_principal(ZZRingElem(2)*D+D2) == true
    @test is_principal(2*D-D2) == true
    @test coefficients(D2+D2) == coefficients(2*D2)
    @test coefficients(D2-D2) == [0, 0, 0, 0]
  end

  @testset "Error handling" begin
    X = normal_toric_variety(
      [[1, 2], [1, 3, 4], [2, 3, 4]], [[-1, -1, 0], [-1, -1, -2], [3, 0, 2], [0, 3, 2]]
    )
    d = canonical_divisor(X)
    @test_throws ArgumentError is_ample(d)
    @test_throws ArgumentError is_cartier(d)
  end
end
