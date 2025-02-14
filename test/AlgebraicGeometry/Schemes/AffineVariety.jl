@testset "affine variety" begin
  P, (x,y) = polynomial_ring(QQ, [:x,:y])
  X = variety(x, check=false)
  @test is_geometrically_integral(X)

  P, (x,y) = polynomial_ring(QQ, [:x, :y])
  p = x^2 + y^2 + 1
  X =  variety(p)
  Y = variety(spec(P, ideal([p^2])))
  @test X == Y
  q = x^2 + y^2 # irreducible but not geometrically irreducible.
  @test_throws ErrorException variety(q)

  U = variety(hypersurface_complement(X, y))
  @test closure(U) == variety(p)

  A = variety(x^2 - y)
  B = variety(y)
  ABset = set_theoretic_intersection(A, B)
  ABsch = intersect(A, B)

  @test ABset == variety(ideal([x,y])) # reduced no multiplicity
  @test ABset != ABsch # reduced no multiplicity
end
