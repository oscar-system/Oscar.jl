@testset "WeylGroup" begin
  @testset "weyl_group(cartan_matrix::ZZMatrix)" begin
    W = weyl_group(cartan_matrix(:A, 2))
    @test isfinite(W) == true
    @test ngens(W) == 2

    W = weyl_group(cartan_matrix(:A, 3))
    @test isfinite(W) == true
    @test ngens(W) == 3

    W = weyl_group(cartan_matrix(:B, 2))
    @test isfinite(W) == true
    @test ngens(W) == 2

    W = weyl_group(ZZ[2 -2; -2 2]) # TODO: replace with cartan_matrix(A_1^(1)), once functionality for affine type is added
    @test isfinite(W) == false
  end

  @testset "accessors" begin
    W = weyl_group(:A, 2)
    @test isfinite(W) === W.finite
    @test root_system(W) === W.root_system
  end

  @testset "longest_element(W::WeylGroup)" begin
    # A1
    W = weyl_group(:A, 1)
    @test longest_element(W) == gen(W, 1)

    # A2
    W = weyl_group(:A, 2)
    @test word(longest_element(W)) == UInt8[1, 2, 1]

    # B2
    W = weyl_group(:B, 2)
    @test word(longest_element(W)) == UInt8[2, 1, 2, 1]

    # B3
    W = weyl_group(:B, 3)
    @test word(longest_element(W)) == UInt8[3, 2, 3, 1, 2, 3, 1, 2, 1]

    # F4
    W = weyl_group(:F, 4)
    @test word(longest_element(W)) ==
      UInt8[4, 3, 2, 3, 1, 2, 3, 4, 3, 2, 3, 1, 2, 3, 4, 3, 2, 3, 1, 2, 3, 1, 2, 1]

    # G2
    W = weyl_group(:G, 2)
    @test word(longest_element(W)) == UInt8[2, 1, 2, 1, 2, 1]
  end
end

@testset "Base.:(*)(x::WeylGroupElem, y::WeylGroupElem)" begin
  # test A2
  W = weyl_group(:A, 2)
  s = gens(W)

  @test parent(s[1] * s[2]) === parent(s[1]) === parent(s[2])

  @test word(s[2] * s[1]) == UInt[2, 1]
  @test word(s[1] * s[2]) == UInt[1, 2]

  @test word(s[1] * s[2] * s[1]) == UInt[1, 2, 1]
  @test word(s[2] * s[1] * s[2]) == UInt[1, 2, 1]

  # test A3
  W = weyl_group(:A, 3)
  s = gens(W)

  @test parent(s[1] * s[2]) === parent(s[1]) === parent(s[2])

  @test word(s[2] * s[1]) == UInt8[2, 1]
  @test word(s[1] * s[2]) == UInt8[1, 2]

  @test word(s[3] * s[1]) == UInt8[3, 1]
  @test word(s[1] * s[3]) == UInt8[3, 1]

  @test word(s[3] * s[2]) == UInt8[3, 2]
  @test word(s[2] * s[3]) == UInt8[2, 3]

  @test word(s[1] * s[3] * s[1]) == UInt8[3]
  @test word(s[3] * s[1] * s[3]) == UInt8[1]

  @test word(s[1] * s[2] * s[1]) == UInt8[1, 2, 1]
  @test word(s[2] * s[1] * s[2]) == UInt8[1, 2, 1]

  @test word(s[2] * s[3] * s[2]) == UInt8[2, 3, 2]
  @test word(s[3] * s[2] * s[3]) == UInt8[2, 3, 2]
end

@testset "Base.:(*)(x::WeylGroupElem, w::WeightLatticeElem)" begin
  R = root_system(:A, 2)
  W = weyl_group(R)
  
  rho = weyl_vector(R)
  @test longest_element(W)*rho == -rho
end

@testset "ReducedExpressionIterator" begin
  W = weyl_group(:A, 3)
  s = gens(W)

  # test for s1
  iter = reduced_expressions(s[1])
  @test iter.el === s[1]
  @test iter.up_to_commutation == false

  re = collect(iter)
  @test length(re) == 1
  @test re[1] == word(s[1])

  # test for w0
  w0 = longest_element(W)
  iter = reduced_expressions(w0)
  @test iter.el === w0
  @test iter.up_to_commutation == false

  re = collect(iter)
  @test length(re) == 16
  @test re[1] == word(w0)
  @test re[16] == UInt8[3, 2, 1, 3, 2, 3]

  iter = reduced_expressions(w0; up_to_commutation=true)
  @test iter.el === w0
  @test iter.up_to_commutation == true

  re = collect(iter)
  @test length(re) == 8
  @test re[1] == word(w0)
  @test re[8] == UInt8[3, 2, 3, 1, 2, 3]
end
