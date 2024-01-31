@testset "Projective rational points" begin
  P2 = projective_space(QQ, 2)

  @test P2([0,1,0]) == P2([0,2,0])
  @test P2([0,1,2]) == P2([0,2,4])
  (s0, s1, s2) = homogeneous_coordinates(P2)
  X = subscheme(P2, s2*s1-s0^2)

  @test_throws ErrorException P2([0,0,0])
  @test_throws ErrorException X([1,2,1])

  ideal(X([1,1,1])) == ideal(P2([1,1,1]))
  A = algebraic_set(X)
  @test A([1,1,1]) == X([1,1,1])

  @test is_subscheme(scheme(A([1,1,1])), A)


  P2Z = projective_space(ZZ, 2)
  @test_throws ErrorException P2Z([2,2,2])
  @test P2([0,1,2]) == P2([0,-1,-2])


end
