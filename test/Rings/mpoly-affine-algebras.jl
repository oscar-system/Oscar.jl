@testset "mpoly_affine_algebras.normalization" begin
  R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
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
    R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
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

  @test !is_normal(quo(R, ideal(R, [x^2 - y^3]))[1]; check=false)
  @test is_normal(quo(R, ideal(R, [x - y^3]))[1]; check=false)
  @test !is_normal(quo(R, ideal(R, [x^2]))[1])
  @test is_normal(quo(R, ideal(R, [z^2 - x*y]))[1])

  R, (x, y, z) = polynomial_ring(ZZ, [:x, :y, :z])
  @test_throws ArgumentError is_normal(quo(R, ideal(R, [x - y^3]))[1])

  R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
  Q, _ = quo(R, ideal(R, [(x^2 - y^3)]))
  for (S, M, FI) in normalization(Q, algorithm=:isPrime)
    @test parent(FI[1]) == Q
    @test isa(FI[2], Oscar.Ideal)
    @test parent(M(Q(x+y))) == S
  end
end

@testset "mpoly_affine_algebras.integral_basis" begin
  R, (x, y) = polynomial_ring(QQ, [:x, :y])
  (den, nums) = integral_basis(y^5-x^3*(x+1)^4, 2; algorithm = :hensel)
  @test all(p -> base_ring(parent(p)) == R, nums)
  @test base_ring(parent(den)) == R
  @test_throws ArgumentError integral_basis(x*y^5-x^3*(x+1)^4, 2)
  @test_throws ArgumentError integral_basis(y^5-x^3*(x+1)^4, 3)
  @test_throws ArgumentError integral_basis((x+y)*(x+y^2), 1)

  R, (x, y) = polynomial_ring(GF(2), [:x, :y])
  @test_throws ArgumentError integral_basis(y^5-x^3*(x+1)^4, 2; algorithm = :what)
  (den, nums) = integral_basis(y^5-x^3*(x+1)^4, 2; algorithm = :normal_global)
  (den, nums) = integral_basis(y^5-x^3*(x+1)^4, 2; algorithm = :normal_local)
  @test all(p -> base_ring(parent(p)) == R, nums)
  @test base_ring(parent(den)) == R

  R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
  @test_throws ArgumentError integral_basis(y^5-x^3*(x+1)^4, 2)

  R, (x, y) = polynomial_ring(GF(2), [:x, :y])
  @test_throws ArgumentError integral_basis(x*y^5-x^3*(x+1)^4, 2)

  R, (x, y) = polynomial_ring(rational_function_field(QQ, [:s, :t])[1], [:x, :y])
  @test_throws NotImplementedError integral_basis(y^5-x^3*(x+1)^4, 2)
end

@testset "Noether normalization" begin
  R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
  A, _ = quo(R, ideal(R, [x*y, x*z]))
  L = noether_normalization(A)
  @test length(L) == 3 # just a smoke test
end

@testset "mpoly_affine_algebra.vector_space_dim" begin
  r, (x, y) = polynomial_ring(QQ, [:x, :y])
  @test !is_finite_dimensional_vector_space(quo(r, ideal(r, [x^2+y^2]))[1])
  @test_throws AbstractAlgebra.InfiniteDimensionError vector_space_dim(quo(r, ideal(r, [x^2+y^2]))[1])
  @test vector_space_dim(quo(r, ideal(r, [x^2+y^2, x^2-y^2]))[1]) == 4
  @test is_finite_dimensional_vector_space(quo(r, ideal(r, [x^2+y^2, x^2-y^2]))[1])
  @test vector_space_dim(quo(r, ideal(r, [one(r)]))[1]) == 0
  @test is_finite_dimensional_vector_space(quo(r, ideal(r, [one(r)]))[1])

  @test !is_finite_dimensional_vector_space(r)

  r, (x, y) = polynomial_ring(ZZ, [:x, :y])
  @test_throws ErrorException vector_space_dim(quo(r, ideal(r, [x, y]))[1])
end

@testset "mpoly_affine_algebra.monomial_basis" begin
  R, (x, y) = graded_polynomial_ring(QQ, [:x, :y])
  A, _ = quo(R, ideal(R, [x^2, y^3]))
  mons = monomial_basis(A)
  @test mons == [x*y^2, y^2, x*y, y, x, 1]
  @test length(mons) == vector_space_dim(A)

  A, _ = quo(R, ideal(R, [one(R)]))
  @test isempty(monomial_basis(A))
end

@testset "Subalgebra membership" begin
  R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
  Q, _ = quo(R, ideal(R, [z - y^3]))
  f = Q(z + x^2 - y^6)
  v = [ Q(x), Q(y^3) ]
  fl, t = subalgebra_membership(f, v)
  @test fl
  @test t(v...) == f
  g = Q(y)
  fl, t = subalgebra_membership(g, v)
  @test !fl

  R, (x, y, z) = graded_polynomial_ring(QQ, [ :x, :y, :z ], [ 3, 1, 3 ])
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

  R, (x, y) = graded_polynomial_ring(QQ, [:x, :y], [ 1, 2 ])
  Q, _ = quo(R, ideal(R, [x^10 - y^5]))
  for S in [R, Q]
    u = gen(S, 1)
    v = gen(S, 2)
    V = [ u^4 + v^2, v, u ]
    V1 = @inferred minimal_subalgebra_generators(V)
    @test all(f -> parent(f) === S, V1)
    @test V1 == [ u, v ]
    V2, rels = @inferred Oscar.minimal_subalgebra_generators_with_relations(V)
    @test all(f -> parent(f) === S, V2)
    @test V2 == [ u, v ]
    for i = 1:length(V)
      @test V[i] == rels[i](V2...)
    end

    # An example, where the minimal generators aren't just the first elements
    V = [u^2, v^2, u^2*v^2, v^3]
    V1 = @inferred minimal_subalgebra_generators(V)
    @test V1 == [u^2, v^2, v^3]
    @test all(f -> parent(f) === S, V1)
    V2, rels = @inferred Oscar.minimal_subalgebra_generators_with_relations(V)
    @test V2 == [u^2, v^2, v^3 ]
    @test all(f -> parent(f) === S, V2)
    for i = 1:length(V)
      @test V[i] == rels[i](V2...)
    end
  end
end

@testset "algebraic independence" begin
  R, (x, y) = polynomial_ring(QQ, [:x, :y])
  V = [x, y]
  fl, I = is_algebraically_independent_with_relations(V)
  @test fl
  @test is_zero(I)
  @test fl == is_algebraically_independent(V)

  V = [2x, 3y+x^3, y]
  fl, I = is_algebraically_independent_with_relations(V)
  @test !fl
  t1, t2, t3 = gens(base_ring(I))
  @test I == ideal([t2 - 1//8*t1^3 - 3t3])
  @test fl == is_algebraically_independent(V)
end
