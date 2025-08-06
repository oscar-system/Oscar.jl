@testset "LieAlgebras.LieAlgebraHom" begin
  @testset "Constructor and basic properties" begin
    L1 = special_linear_lie_algebra(QQ, 2)
    L2 = general_linear_lie_algebra(QQ, 2)
    h = hom(L1, L2, matrix(QQ, [0 1 0 0; 0 0 1 0; 1 0 0 -1])) # embed sl_2 into gl_2

    @test domain(h) == L1
    @test codomain(h) == L2
    @test matrix(h) == matrix(QQ, [0 1 0 0; 0 0 1 0; 1 0 0 -1])
    @test is_welldefined(h)
  end

  @testset "Image and kernel" begin
    L = abelian_lie_algebra(QQ, 6)
    x1, x2, x3, x4, x5, x6 = basis(L)
    h = hom(L, L, [zero(L), zero(L), zero(L), x3, x4, x5])
    @test image(h, x1) == h(x1) == zero(L)
    @test image(h, x2) == h(x2) == zero(L)
    @test image(h, x3) == h(x3) == zero(L)
    @test image(h, x4) == h(x4) == x3
    @test image(h, x5) == h(x5) == x4
    @test image(h, x6) == h(x6) == x5

    @test image(h) == sub(L, [x3, x4, x5])
    @test image(h, sub(L)) == sub(L, [x3, x4, x5])
    @test image(h, ideal(L)) == sub(L, [x3, x4, x5])

    @test kernel(h) == ideal(L, [x1, x2, x3])

    @test image(h, sub(L, [x2, x3, x4])) == sub(L, [x3])
    @test h(sub(L, [x2, x3, x4])) == sub(L, [x3])
    @test image(h, ideal(L, [x2, x3, x4])) == sub(L, [x3])
    @test h(ideal(L, [x2, x3, x4])) == sub(L, [x3])
  end

  @testset "Identity and zero map" begin
    L1 = special_linear_lie_algebra(QQ, 3)
    L2 = general_linear_lie_algebra(QQ, 3)

    h = identity_map(L1)
    @test domain(h) == L1
    @test codomain(h) == L1
    @test matrix(h) == identity_matrix(QQ, dim(L1))
    @test is_welldefined(h)
    for x in basis(L1)
      @test h(x) == x
    end
    @test image(h) == sub(L1)
    @test kernel(h) == ideal(L1, [zero(L1)])

    h = zero_map(L1, L2)
    @test domain(h) == L1
    @test codomain(h) == L2
    @test matrix(h) == zero_matrix(QQ, dim(L1), dim(L2))
    @test is_welldefined(h)
    for x in basis(L1)
      @test h(x) == zero(L2)
    end
    @test image(h) == sub(L2, [zero(L2)])
    @test kernel(h) == ideal(L1)
  end

  @testset "Composition" begin
    L1 = special_linear_lie_algebra(QQ, 2)
    L2 = special_linear_lie_algebra(QQ, 3)
    L3 = general_linear_lie_algebra(QQ, 3)
    f = hom(L1, L2, [basis(L2, 1), basis(L2, 4), basis(L2, 7)]) # embedding sl_2 into sl_3
    g = hom(
      L2,
      L3,
      matrix(
        QQ,
        [
          0 1 0 0 0 0 0 0 0
          0 0 1 0 0 0 0 0 0
          0 0 0 0 0 1 0 0 0
          0 0 0 1 0 0 0 0 0
          0 0 0 0 0 0 1 0 0
          0 0 0 0 0 0 0 1 0
          1 0 0 0 -1 0 0 0 0
          0 0 0 0 1 0 0 0 -1
        ],
      ),
    ) # embedding sl_3 into gl_3
    h = compose(f, g)
    @test domain(h) == L1
    @test codomain(h) == L3
    for x in basis(L1)
      @test h(x) == g(f(x))
    end
  end

  @testset "Inverses" begin
    L1 = special_linear_lie_algebra(QQ, 2)
    L2 = general_linear_lie_algebra(QQ, 2)
    h = hom(L1, L2, matrix(QQ, [0 1 0 0; 0 0 1 0; 1 0 0 -1])) # embed sl_2 into gl_2

    @test !is_isomorphism(h)

    @test is_isomorphism(identity_map(L1))
    @test identity_map(L1) == inv(identity_map(L1))

    @test !is_isomorphism(zero_map(L1))
  end

  @testset "Hum72, Exercise 2.10" begin
    L = special_linear_lie_algebra(QQ, 2)
    e, f, h = basis(L)
    h = hom(L, L, [-f, -e, -h])
    @test is_isomorphism(h)
    @test h == inv(h)
    @test identity_map(L) == compose(h, inv(h))
    @test identity_map(L) == compose(inv(h), h)
  end
end
