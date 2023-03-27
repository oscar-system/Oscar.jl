@testset "affine algebraic sets" begin
  P,(x,y) = polynomial_ring(QQ, [:x, :y])
  p = x^2+y^2+1
  X = affine_variety(p)
  Y = affine_variety(Spec(P, ideal([p^2])))
  @test X == Y
  q = x^2+y^2 # irreducible but not geometrically irreducible.
  @test_throws ErrorException affine_variety(q)
  Z = vanishing_locus(p) # works

  U = affine_variety(hypersurface_complement(X,y))
  @test closure(U) == affine_variety(p)

  A = affine_variety(x^2-y)
  B = affine_variety(y)
  ABset = set_theoretic_intersection(A,B)
  ABsch = intersect(A,B)

  @test ABset == vanishing_locus(ideal([x,y])) # reduced no mutiplicity
  @test ABset != ABsch # reduced no mutiplicity
end
