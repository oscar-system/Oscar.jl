@testset "construct affine semigroups" begin
    #Q = ZZ^2
    Q = affine_semigroup([1 0;0 1])
    @test is_pointed(Q)
    C = cone(Q)
    @test ambient_dim(C) == ambient_dim(Q)
    H = hyperplanes(Q)
    @test length(H) == 2
    @test dim(zonotope(Q)[1]) == rank(Q) 

    #non-normal example
    Q = affine_semigroup([2 3 0 0;0 0 2 3])
    @test is_pointed(Q)

    #= TODO: monoid algebras must go elsewhere!
    kQ = monoid_algebra(Q,QQ)
    @test !is_normal(kQ)

    #non-pointed example
    Q = affine_semigroup([[-1],[2]])
    @test !is_pointed(Q)
    @test_throws ArgumentError monoid_algebra(Q, QQ)
    @test_throws ArgumentError monoid_algebra([[-1], [2]], QQ)
    @test_throws ArgumentError monoid_algebra([-1 2], QQ)
    =#
end

@testset "Affine semigroups" begin
  # Test different constructors
  A = ZZ[1 0; 0 1]
  Q = Oscar.AffineSemigroup(A)
  M = [1 0; 0 1]
  Q = Oscar.affine_semigroup(M)
  M = [[2, 0], [5, 1], [3, 0]]
  Q = Oscar.affine_semigroup(M)
  C = Oscar.InjectiveResolutions.polyhedral_cone(Q)
  Q = Oscar.AffineSemigroup(C)
  @test C == Oscar.InjectiveResolutions.polyhedral_cone(Q)
  M = [1 0; 0 1; -1 0]
  Q = Oscar.affine_semigroup(M)
  C = Oscar.InjectiveResolutions.polyhedral_cone(Q)
  @test_throws AssertionError Oscar.AffineSemigroup(C) # cone is not pointed, which is not supported at the moment
  @test gens(Vector, Q) isa Vector
  @test gen(Vector, Q, 3) isa Vector
  F = Oscar.InjectiveResolutions.ambient_group(Q)
  @test repres(Q(gen(Q, 1)) + 2*Q(gen(Q, 2))) == F([1, 2])
end

@testset "test if a semigroup is saturated" begin
  A = ZZ[1 0;0 1]
  Q = Oscar.AffineSemigroup(A)
  @test Oscar.InjectiveResolutions.is_saturated(Q)
  @test [1,3] in Q

  A = ZZ[2 0;3 0]
  Q = Oscar.AffineSemigroup(A)
  @test !Oscar.InjectiveResolutions.is_saturated(Q)
  @test !([1,0] in Q)
  @test [5,0] in Q
end