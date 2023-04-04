@testset "Projective Varieties" begin
  P,(x,y,z) = polynomial_ring(QQ, [:x, :y])
  G = grade(P)
  p = x^2+y^2+z^2
  X = projective_variety(p)
  Y = projective_scheme(P, ideal([p^2]))
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
