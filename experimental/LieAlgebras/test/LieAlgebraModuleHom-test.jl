@testset "LieAlgebras.LieAlgebraModuleHom" begin
  @testset "Constructor and basic properties" begin
    L = special_orthogonal_lie_algebra(QQ, 3)

    V1 = standard_module(L)
    V2 = direct_sum(V1, trivial_module(L, 3))

    h = hom(V1, V2, [basis(V2, 1), basis(V2, 2), basis(V2, 3)]) # embed sl_2 into gl_2

    @test domain(h) == V1
    @test codomain(h) == V2
    @test matrix(h) == matrix(
      QQ,
      [
        1 0 0 0 0 0
        0 1 0 0 0 0
        0 0 1 0 0 0
      ],
    )
    for x in basis(L)
      for v in basis(V1)
        @test h(x * v) == x * h(v)
      end
    end
  end

  @testset "Image and kernel" begin
    L = special_orthogonal_lie_algebra(QQ, 3)
    V = direct_sum(standard_module(L), standard_module(L))
    v1, v2, v3, v4, v5, v6 = basis(V)
    h = hom(V, V, [zero(V), zero(V), zero(V), v1, v2, v3])
    @test image(h, v1) == h(v1) == zero(V)
    @test image(h, v2) == h(v2) == zero(V)
    @test image(h, v3) == h(v3) == zero(V)
    @test image(h, v4) == h(v4) == v1
    @test image(h, v5) == h(v5) == v2
    @test image(h, v6) == h(v6) == v3
  end

  @testset "Composition" begin
    L = special_orthogonal_lie_algebra(QQ, 3)
    V1 = standard_module(L)
    V2 = direct_sum(V1, V1)
    V3 = direct_sum(V1, V1, V1)
    f = hom(V1, V2, [V2([v, zero(V1)]) for v in basis(V1)])
    g = hom(
      V2,
      V3,
      [basis(V3, 1), basis(V3, 2), basis(V3, 3), basis(V3, 7), basis(V3, 8), basis(V3, 9)],
    )

    h = compose(f, g)
    @test domain(h) == V1
    @test codomain(h) == V3
    for x in basis(V1)
      @test h(x) == g(f(x))
    end
  end

  @testset "Inverses" begin
    L = special_orthogonal_lie_algebra(QQ, 3)
    V = direct_sum(standard_module(L), standard_module(L))
    v1, v2, v3, v4, v5, v6 = basis(V)
    h = hom(V, V, [zero(V), zero(V), zero(V), v1, v2, v3])
    @test !is_isomorphism(h)

    @test is_isomorphism(identity_map(V))
    @test identity_map(V) == inv(identity_map(V))

    h = hom(V, V, [v4, v5, v6, v1, v2, v3])
    @test is_isomorphism(h)
    @test h == inv(h)
    @test identity_map(V) == compose(h, inv(h))
    @test identity_map(V) == compose(inv(h), h)
  end
end
