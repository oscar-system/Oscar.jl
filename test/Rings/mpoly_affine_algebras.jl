@testset "mpoly_affine_algebras.normalization" begin
  R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
  Q, _ = quo(R, ideal(R, [z - x^4, z - y^6]))
  for (S, M, FI) in normalization(Q)
    @test parent(FI[1]) == Q
    @test isa(FI[2], Oscar.Ideal)
    @test parent(M(Q(x+y))) == S
  end

  stuff, deltas, tot_delta = normalization_with_delta(Q)
  @test isa(tot_delta, Int)
  @test length(stuff) == length(deltas)
  for ((S, M, FI), delta) in zip(stuff, deltas)
    @test parent(FI[1]) == Q
    @test isa(FI[2], Oscar.Ideal)
    @test parent(M(Q(x+y))) == S
    @test isa(delta, Int)
  end

  for algorithm in (:equidimDec, :primeDec)
    R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
    Q, _ = quo(R, ideal(R, [(x^2 - y^3)*(x^2 + y^2)*x]))
    for (S, M, FI) in normalization(Q, algorithm=algorithm)
      @test parent(FI[1]) == Q
      @test isa(FI[2], Oscar.Ideal)
      @test parent(M(Q(x+y))) == S
    end

    stuff, deltas, tot_delta = normalization_with_delta(Q, algorithm=algorithm)
    @test isa(tot_delta, Int)
    @test length(stuff) == length(deltas)
    for ((S, M, FI), delta) in zip(stuff, deltas)
      @test parent(FI[1]) == Q
      @test isa(FI[2], Oscar.Ideal)
      @test parent(M(Q(x+y))) == S
      @test isa(delta, Int)
    end
  end

  @test !is_normal(quo(R, ideal(R, [x^2 - y^3]))[1])
  @test is_normal(quo(R, ideal(R, [x - y^3]))[1])

  R, (x, y, z) = polynomial_ring(ZZ, ["x", "y", "z"])
  @test_throws ArgumentError is_normal(quo(R, ideal(R, [x - y^3]))[1])
end

@testset "mpoly_affine_algebras.integral_basis" begin
  R, (x, y) = polynomial_ring(QQ, ["x", "y"])
  (den, nums) = integral_basis(y^5-x^3*(x+1)^4, 2; algorithm = :hensel)
  @test all(p -> base_ring(parent(p)) == R, nums)
  @test base_ring(parent(den)) == R
  @test_throws ArgumentError integral_basis(x*y^5-x^3*(x+1)^4, 2)
  @test_throws ArgumentError integral_basis(y^5-x^3*(x+1)^4, 3)
  @test_throws ArgumentError integral_basis((x+y)*(x+y^2), 1)

  R, (x, y) = polynomial_ring(GF(2), ["x", "y"])
  @test_throws ArgumentError integral_basis(y^5-x^3*(x+1)^4, 2; algorithm = :what)
  (den, nums) = integral_basis(y^5-x^3*(x+1)^4, 2; algorithm = :normal_global)
  (den, nums) = integral_basis(y^5-x^3*(x+1)^4, 2; algorithm = :normal_local)
  @test all(p -> base_ring(parent(p)) == R, nums)
  @test base_ring(parent(den)) == R

  R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
  @test_throws ArgumentError integral_basis(y^5-x^3*(x+1)^4, 2)

  R, (x, y) = polynomial_ring(GF(2), ["x", "y"])
  @test_throws ArgumentError integral_basis(x*y^5-x^3*(x+1)^4, 2)

  R, (x, y) = polynomial_ring(RationalFunctionField(QQ, ["s", "t"])[1], ["x", "y"])
  @test_throws NotImplementedError integral_basis(y^5-x^3*(x+1)^4, 2)
end

@testset "Noether normalization" begin
  R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
  A, _ = quo(R, ideal(R, [x*y, x*z]))
  L = noether_normalization(A)
  @test length(L) == 3 # just a smoke test
end

@testset "mpoly_affine_algebra.vector_space_dimension" begin
  r, (x, y) = polynomial_ring(QQ, [:x, :y])
  @test vector_space_dimension(quo(r, ideal(r, [x^2+y^2]))[1]) == -1
  @test vector_space_dimension(quo(r, ideal(r, [x^2+y^2, x^2-y^2]))[1]) == 4

  r, (x, y) = polynomial_ring(ZZ, [:x, :y])
  @test_throws ErrorException vector_space_dimension(quo(r, ideal(r, [x, y]))[1])
end

@testset "Subalgebra membership" begin
  R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
  Q, _ = quo(R, ideal(R, [z - y^3]))
  f = Q(z + x^2 - y^6)
  v = [ Q(x), Q(y^3) ]
  fl, t = subalgebra_membership(f, v)
  @test fl
  @test t(v...) == f
  g = Q(y)
  fl, t = subalgebra_membership(g, v)
  @test !fl

  R, (x, y, z) = graded_polynomial_ring(QQ, [ "x", "y", "z" ], [ 3, 1, 3 ])
  f = x^2 - y^6 + z^2
  v = [ x, y^3, z - y^3 ]
  fl, t = subalgebra_membership_homogeneous(f, v)
  @test fl
  @test t(v...) == f
  g = y
  fl, t = subalgebra_membership_homogeneous(g, v)
  @test !fl

  I = ideal(R, [ z - y^3 ])
  Q, RtoQ = quo(R, I)
  f = Q(x)^2 + Q(y)^6 + Q(z)^2
  v = [ Q(x), Q(y)^3 ]
  fl, t = subalgebra_membership_homogeneous(f, v)
  @test fl
  @test t(v...) == f
  g = Q(y)
  fl, t = subalgebra_membership_homogeneous(g, v)
  @test !fl

  R, (x, y) = graded_polynomial_ring(QQ, ["x", "y"], [ 1, 2 ])
  V = [ x^4 + y^2, y, x ]
  V1 = minimal_subalgebra_generators(V)
  @test V1 == [ x, y ]
  V2, rels = Oscar.minimal_subalgebra_generators_with_relations(V)
  @test V2 == [ x, y ]
  for i = 1:length(V)
    @test V[i] == rels[i](V2...)
  end
end
