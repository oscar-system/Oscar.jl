@testset "lifting morphisms through free resolutions" begin
  R, (x, y, z) = QQ[:x, :y, :z]

  I = ideal(R, [x, y, z])
  J = ideal(R, [x^2, y^2, z^2])

  F = FreeMod(R, 1)
  IF, _ = I*F
  JF, _ = J*F

  M1, _ = quo(F, IF)
  M2, _ = quo(F, JF)

  phi = hom(M2, M1, gens(M1))

  res1 = free_resolution(M1).C
  res2 = free_resolution(M2).C

  psi = Oscar.lift_morphism_through_free_resolutions(phi, domain_resolution=res2, codomain_resolution=res1);

  @test Oscar.extends_right(psi)
  @test !Oscar.extends_left(psi)
  @test_throws ErrorException psi[-5]
  @test Oscar.left_bound(psi) == -2
  @test Oscar.right_bound(psi) == -1
  @test matrix(psi[-2]) == R[;]
  @test matrix(psi[-1]) == R[1;]
  @test matrix(psi[0]) == R[1;]
  @test Oscar.right_bound(psi) == 0
  @test matrix(psi[2]) == R[y*z 0 0; 0 x*y 0; 0 0 x*z]
  @test Oscar.right_bound(psi) == 2
  @test matrix(psi[1]) == R[z 0 0; 0 y 0; 0 0 x]
  @test Oscar.right_bound(psi) == 2
  @test matrix(psi[3]) == R[x*y*z;]
  @test Oscar.right_bound(psi) == 3
end



