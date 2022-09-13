@testset "PBWAlgebra.constructor" begin
  r, (x, y, z) = QQ["x", "y", "z"]
  R, (x, y, z) = pbw_algebra(r, [0 x*y x*z; 0 0 y*z + 1; 0 0 0], deglex(gens(r)))

  @test elem_type(R) == PBWAlgElem{fmpq, Singular.n_Q}
  @test parent_type(x) == PBWAlgRing{fmpq, Singular.n_Q}
  @test coefficient_ring(R) == QQ
  @test coefficient_ring(x) == QQ
  @test symbols(R) == [:x, :y, :z]
  @test gens(R) == [x, y, z]
  @test ngens(R) == 3
  @test gen(R, 2) == y
  @test R[2] == y

  3*x^2 + y*z == R([3, 1], [[2, 0, 0], [0, 1, 1]])
end

@testset "PBWAlgebra.printing" begin
  r, (x, y, z) = QQ["x", "y", "z"]
  R, (x, y, z) = pbw_algebra(r, [0 x*y x*z; 0 0 y*z + 1; 0 0 0], deglex(gens(r)))

  @test length(string(R)) > 2
  @test length(string(x + y)) > 2
  @test string(one(R)) == "1"
end

@testset "PBWAlgebra.iteration" begin
  r, (x, y, z) = QQ["x", "y", "z"]
  R, (x, y, z) = pbw_algebra(r, [0 x*y x*z; 0 0 y*z + 1; 0 0 0], deglex(gens(r)))

  p = -((x*z*y)^6 - 1)

  @test iszero(constant_coefficient(p - constant_coefficient(p)))
  @test length(p) == length(collect(terms(p)))
  @test length(p) == length(collect(monomials(p)))

  s = zero(R)
  for t in terms(p)
    s += t
  end
  @test p == s
  @test p == leading_term(p) + tail(p)

  s = zero(R)
  for (c, m) in zip(coefficients(p), monomials(p))
    s += c*m
  end
  @test p == s
  @test p == leading_coefficient(p)*leading_monomial(p) + tail(p)

  s = build_ctx(R)
  for (c, e) in zip(coefficients(p), exponent_vectors(p))
    push_term!(s, c, e)
  end
  @test p == finish(s)

  @test p == R(collect(coefficients(p)), collect(exponent_vectors(p)))
end

@testset "PBWAlgebra.ideals" begin
  r, (x, y, z) = QQ["x", "y", "z"]
  R, (x, y, z) = pbw_algebra(r, [0 x*y x*z; 0 0 y*z + 1; 0 0 0], deglex(gens(r)))

  I = left_ideal([x^2, y^2])
  @test length(string(I)) > 2
  @test ngens(I) > 1
  @test !iszero(I)
  @test !isone(I)
  @test x^2 - y^2 in I
  @test !(x + 1 in I)
  @test isone(I + left_ideal([z^2]))
end

@testset "PBWAlgebra.weyl_algebra" begin
  R, (x, dx) = weyl_algebra(QQ, ["x"])
  @test dx*x == 1 + x*dx

  R, (x, y, dx, dy) = weyl_algebra(QQ, ["x", "y"])
  @test dx*x == 1 + x*dx
  @test dy*y == 1 + y*dy
  @test dx*y == y*dx
  @test dy*x == x*dy
  @test x*y == y*x
end
