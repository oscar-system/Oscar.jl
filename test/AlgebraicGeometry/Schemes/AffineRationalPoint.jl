@testset "affine rational points" begin
  A2 = affine_space(GF(2), [:x, :y]);
  (x, y) = coordinates(A2);
  X = algebraic_set(x*y);

  pX = X([1,0])
  pA = A2(pX)
  @test pX == pA
  @test A2([1,0]) == pA
  # ring around the rosie once
  @test pX == X(rational_point_coordinates(defining_ideal(scheme(pX))))
  Oscar.closed_embedding(pA)
  Oscar.closed_embedding(pX)
  @test dim(tangent_space(pX))==1
  @test dim(tangent_space(pA))==2

  @test pX in X
  @test pA in X
  @test codomain(pX)===X
  @test ambient_scheme(pX)===X
  
  A3 = affine_space(QQ,3)
  x,y,z = coordinates(A3)
  X = algebraic_set(x^2+y^2+z^3)
  @test is_du_val_singularity(X([0,0,0]))

  @test is_prime(ideal(pX))
  @test pX[1] == 1
  @test pX[2] == 0

  @test ideal(pA) == ideal(pX)

  # parametrize a cuspidal cubic rational curve
  A2 = affine_space(QQ, [:x, :y])
  (x,y) = coordinates(A2)
  At = affine_space(QQ, [:t])
  (t,) = coordinates(At)
  C = algebraic_set(x^3-y^2)
  f = morphism(At, C, [t^2,t^3])
  @test f(At([1,])) == C([1,1])
  @test f(At([2,])) == C([4,8])
  @test tangent_space(C,C([0,0])) == algebraic_set(0*x)
  @test !is_smooth(C([0,0]))
  D = algebraic_set(x^2-y)
  @test is_smooth(D([0,0]))
  @test tangent_space(D,D([0,0])) == algebraic_set(y)

  # disjoint union of a line and a point, smooth
  E = algebraic_set(ideal([x*(y-1),y^2-y]))
  @test is_smooth(E([0,0]))
  @test is_smooth(E([0,1]))

  # union of two lines meeting in 0
  F = algebraic_set(ideal([x*y]))
  @test !is_smooth(F([0,0]))
  @test is_smooth(F([0,1]))

  P2cov = covered_scheme(projective_space(QQ,2))
  pt = P2cov[1][1]([1,0])
  @test pt in P2cov

  k = GF(2)
  L = GF(2,2)
  A2 = affine_space(k ,2)
  A2k = Oscar.RationalPointSet(spec(k), A2)
  A2L = Oscar.RationalPointSet(spec(L), A2)
  pk = A2k([1,1])
  pL = A2L(pk) # conversion
  # @test_throws ErrorException pk == pL # Does not throw anymore. Why?
  #@test !(pk == pL) # no automatic coercion -> consistent with hom interpretation since the domains differ

  # reduction mod p
  k = ZZ
  L = GF(2,2)
  A2 = affine_space(k ,2)
  A2k = Oscar.RationalPointSet(spec(k), A2)
  A2L = Oscar.RationalPointSet(spec(L), A2)
  pk = A2k([1,1])
  pL = A2L(pk) # conversion
  @test !(pk == pL) # no automatic coercion -> consistent with hom interpretation since the domains differ
end
