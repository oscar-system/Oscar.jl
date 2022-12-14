@testset "PBWAlgebraQuo.constructor" begin
  r, (x, y, z) = QQ["x", "y", "z"]
  R, (x, y, z) = pbw_algebra(r, [0 x*y x*z; 0 0 y*z + 1; 0 0 0], deglex(r))
  Q, M = quo(R, two_sided_ideal(R, [x]))
  (X, Y, Z) = gens(Q)

  @test elem_type(Q) == PBWAlgQuoElem{fmpq, Singular.n_Q}
  @test parent_type(X) == PBWAlgQuo{fmpq, Singular.n_Q}
  @test coefficient_ring(R) == QQ
  @test coefficient_ring(X) == QQ
  @test base_ring(Q) == R
  @test base_ring(X) == R
  @test gens(Q) == [X, Y, Z]
  @test ngens(Q) == 3
  @test gen(Q, 2) == M(y)
  @test Q[2] == M(y)
end

@testset "PBWAlgebraQuo.printing" begin
  r, (x, y, z) = QQ["x", "y", "z"]
  R, (x, y, z) = pbw_algebra(r, [0 x*y x*z; 0 0 y*z + 1; 0 0 0], deglex(r))
  Q, M = quo(R, two_sided_ideal(R, [x]))
  @test length(string(Q)) > 2
  @test length(string(M(x) + M(y))) > 2
  @test string(one(R)) == "1"
end

function test_elem(Q::PBWAlgQuo{fmpq})
  R = base_ring(Q)
  return Q(R(rand(base_ring(R), 1:4, 1:4, 1:4)))
end

@testset "PBWAlgebraQuo.conformance" begin
  r, (a, h, f, e) = QQ["a", "h", "f", "e"]
  rel = @pbw_relations(e*f == f*e-h, e*h == h*e+2*e, f*h == h*f-2*f)
  R, (a, h, f, e) = pbw_algebra(r, rel, revlex(r))
  Q, _ = quo(R, two_sided_ideal(R, [h]))
  test_NCRing_interface(Q; reps = 1)
end

