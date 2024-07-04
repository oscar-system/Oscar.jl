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
    @test is_welldefined(h)
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

  @testset "Identity and zero map" begin
    L = special_linear_lie_algebra(QQ, 3)
    stdV = standard_module(L)
    V1 = symmetric_power(stdV, 3)[1]
    V2 = exterior_power(stdV, 2)[1]

    h = identity_map(V1)
    @test domain(h) == V1
    @test codomain(h) == V1
    @test matrix(h) == identity_matrix(QQ, dim(V1))
    @test is_welldefined(h)
    for v in basis(V1)
      @test h(v) == v
    end
    # @test image(h) == ...
    # @test kernel(h) == ...

    h = zero_map(V1, V2)
    @test domain(h) == V1
    @test codomain(h) == V2
    @test matrix(h) == zero_matrix(QQ, dim(V1), dim(V2))
    @test is_welldefined(h)
    for v in basis(V1)
      @test h(v) == zero(V2)
    end
    # @test image(h) == ...
    # @test kernel(h) == ...
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
    @test is_welldefined(h)
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
    @test is_welldefined(inv(h))
    @test h == inv(h)
    @test identity_map(V) == compose(h, inv(h))
    @test identity_map(V) == compose(inv(h), h)
  end

  @testset "Direct sum constructions" begin
    L = special_orthogonal_lie_algebra(QQ, 4)
    Vs = [
      standard_module(L), dual(standard_module(L)), exterior_power(standard_module(L), 2)[1]
    ]
    V = direct_sum(Vs...)

    @test canonical_injections(V) == [canonical_injection(V, i) for i in 1:length(Vs)]
    @test canonical_projections(V) == [canonical_projection(V, i) for i in 1:length(Vs)]

    @test all(is_welldefined, canonical_injections(V))
    @test all(is_welldefined, canonical_projections(V))

    # direct sum universal properties
    for i in 1:length(Vs)
      @test canonical_injection(V, i) * canonical_projection(V, i) == identity_map(Vs[i])
    end
    @test sum(
      proj * inj for (proj, inj) in zip(canonical_projections(V), canonical_injections(V))
    ) == identity_map(V)
  end

  @testset "hom_direct_sum" begin
    L = special_orthogonal_lie_algebra(QQ, 4)
    V11 = standard_module(L)
    V12 = dual(standard_module(L))
    V21 = exterior_power(standard_module(L), 2)[1]
    V22 = dual(dual(standard_module(L)))
    V1 = direct_sum(V11, V12)
    V2 = direct_sum(V21, V22)
    h11 = hom(V11, V21, [zero(V21) for _ in basis(V11)])
    h12 = hom(V11, V22, basis(V22))
    h21 = hom(V12, V21, [zero(V21) for _ in basis(V12)])
    h22 = hom(V12, V22, [zero(V22) for _ in basis(V12)])

    h = hom_direct_sum(V1, V2, [h11 h12; h21 h22])
    @test domain(h) == V1
    @test codomain(h) == V2
    @test is_welldefined(h)
    #@test image(h) == image(canonical_injection(V2, 2))
    #@test kernel(h) == kernel(canonical_injjection(V1, 1))
  end

  @testset "hom_tensor" begin
    L = special_orthogonal_lie_algebra(QQ, 4)
    V11 = standard_module(L)
    V12 = dual(standard_module(L))
    V21 = dual(dual(standard_module(L)))
    V22 = exterior_power(standard_module(L), 2)[1]
    V1 = tensor_product(V11, V12)
    V2 = tensor_product(V21, V22)
    h1 = hom(V11, V21, basis(V21))
    h2 = hom(V12, V22, [zero(V22) for _ in basis(V12)])

    h = hom_tensor(V1, V2, [h1, h2])
    @test domain(h) == V1
    @test codomain(h) == V2
    @test is_welldefined(h)

    V11 = standard_module(L)
    V12 = dual(standard_module(L))
    V2i = dual(dual(standard_module(L)))
    V1 = tensor_product(V11, V12)
    V2 = tensor_power(V2i, 2)[1]
    h1 = hom(V11, V2i, basis(V2i))
    h2 = hom(V12, V2i, [zero(V2i) for _ in basis(V12)])

    h = hom_tensor(V1, V2, [h1, h2])
    @test domain(h) == V1
    @test codomain(h) == V2
    @test is_welldefined(h)
  end

  @testset "hom (lift $f_power)" for f_power in
                                     [exterior_power, symmetric_power, tensor_power]
    for k in 0:3
      L = special_orthogonal_lie_algebra(QQ, 4)
      Vb = standard_module(L)
      Wb = dual(dual(standard_module(L)))
      V, _ = f_power(Vb, k)
      W, _ = f_power(Wb, k)
      hb = hom(Vb, Wb, [2 * b for b in basis(Wb)])

      h = hom(V, W, hb)
      @test domain(h) == V
      @test codomain(h) == W
      @test is_welldefined(h)
      vs = elem_type(Vb)[Vb(rand(-10:10, dim(Vb))) for _ in 1:k]
      @test h(V(vs)) == W(hb.(vs))
    end
  end
end
