@testset "affine rational points" begin
  A2 = affine_space(GF(2), [:x, :y]);
  (x, y) = coordinates(A2);
  X = algebraic_set(x*y);

  pX = X([1,0])
  pA = A2([1,0])
  @test pX == pA

  A2a = affine_space(GF(2), [:x, :y]);
  @test A2a([1,0]) == pA

  @test pX in X
  @test pA in X
  @test parent(pX)===X

  @test is_prime(ideal(pX))
  @test px[1] == 1
  @test px[2] == 0

  @test ideal(pa) == ideal(px)

  # parametrize a cuspidal cubic rational curve
  A2 = affine_space(QQ, [:x, :y])
  (x,y) = coordinates(A2)
  At = affine_space(QQ, [:t])
  (t,) = coordinates(At)
  C = algebraic_set(x^3-y^2)
  f = SpecMor(At, C, [t^2,t^3])
  @test f(At([1,])) == C([1,1])
  @test f(At([2,])) == C([4,8])
end
