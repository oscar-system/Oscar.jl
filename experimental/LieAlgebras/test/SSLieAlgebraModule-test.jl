@testset verbose = true "LieAlgebras.SSLieAlgebraModule" begin
  @testset "dim_of_simple_module" begin
    # All concrete test results have been computed using the LiE CAS (http://wwwmathlabo.univ-poitiers.fr/~maavl/LiE/) v2.2.2
    let L = lie_algebra(QQ, :A, 6)
      dim = @inferred dim_of_simple_module(
        Int, L, [1, 3, 5, 0, 1, 0]
      )
      @test dim isa Int
      @test dim == 393513120
    end

    let R = root_system(:B, 7)
      dim = @inferred dim_of_simple_module(
        ZZRingElem, R, [7, 2, 5, 1, 0, 2, 6]
      )
      @test dim isa ZZRingElem
      @test dim == 307689492858882008424585750
    end

    let L = lie_algebra(QQ, :C, 3)
      dim = @inferred dim_of_simple_module(
        L, [3, 3, 3]
      )
      @test dim isa Int
      @test dim == 262144
    end

    let L = lie_algebra(QQ, :D, 5)
      dim = @inferred dim_of_simple_module(
        Int128, L, [1, 2, 3, 4, 5]
      )
      @test dim isa Int128
      @test dim == 591080490000
    end

    let L = lie_algebra(QQ, :E, 6)
      dim = @inferred dim_of_simple_module(
        ZZRingElem, L, [6, 5, 4, 3, 2, 1]
      )
      @test dim isa ZZRingElem
      @test dim == 53947263633682628459250
    end

    let L = lie_algebra(QQ, :F, 4)
      dim = @inferred dim_of_simple_module(
        BigInt, L, [2, 4, 1, 2]
      )
      @test dim isa BigInt
      @test dim == 5989283015625
    end

    let L = lie_algebra(QQ, :G, 2)
      dim = @inferred dim_of_simple_module(L, ZZ.([2, 2]))
      @test dim isa Int
      @test dim == 729
    end

    let L = special_linear_lie_algebra(QQ, 2) # type A_1 but without known root system
      dim = @inferred dim_of_simple_module(Int8, L, [15])
      @test dim isa Int8
      @test dim == 16
    end

    let L = special_orthogonal_lie_algebra(QQ, 7) # type B_3 but without known root system
      dim = @inferred dim_of_simple_module(L, ZZ.([1, 2, 1]))
      @test dim isa Int
      @test dim == 2800
    end
  end

  @testset "dominant_character" begin
    function check_dominant_character(
      LR::Union{LieAlgebra,RootSystem}, hw::Vector{<:Oscar.IntegerUnion}
    )
      domchar = @inferred dominant_character(LR, hw)
      R = LR isa RootSystem ? LR : root_system(LR)
      @test domchar isa Dict{WeightLatticeElem,Int}
      @test domchar[WeightLatticeElem(R, hw)] == 1
      @test issetequal(keys(domchar), dominant_weights(LR, hw))
      @test all(is_dominant, keys(domchar))
      @test all(>=(1), values(domchar))

      domchar_ZZ = @inferred dominant_character(ZZRingElem, LR, hw)
      @test domchar_ZZ isa Dict{WeightLatticeElem,ZZRingElem}
      @test domchar_ZZ == Dict{WeightLatticeElem,ZZRingElem}(domchar)

      return domchar
    end

    # All concrete test results have been computed using the LiE CAS (http://wwwmathlabo.univ-poitiers.fr/~maavl/LiE/) v2.2.2
    let R = root_system(Tuple{Symbol,Int}[]), hw = Int[]
      domchar = check_dominant_character(R, hw)
      @test domchar == Dict(WeightLatticeElem(R, Int[]) => 1)
    end

    let L = lie_algebra(QQ, :A, 2), hw = [0, 0]
      domchar = check_dominant_character(L, hw)
      R = root_system(L)
      @test domchar == Dict(WeightLatticeElem(R, [0, 0]) => 1)
    end

    let L = lie_algebra(QQ, :A, 3), hw = ZZ.([1, 1, 1])
      domchar = check_dominant_character(L, hw)
      R = root_system(L)
      @test domchar == Dict(
        WeightLatticeElem(R, [1, 1, 1]) => 1,
        WeightLatticeElem(R, [2, 0, 0]) => 2,
        WeightLatticeElem(R, [0, 0, 2]) => 2,
        WeightLatticeElem(R, [0, 1, 0]) => 4,
      )
    end

    let L = lie_algebra(QQ, :C, 3), hw = [2, 0, 1]
      domchar = check_dominant_character(L, hw)
      R = root_system(L)
      @test domchar == Dict(
        WeightLatticeElem(R, [2, 0, 1]) => 1,
        WeightLatticeElem(R, [0, 1, 1]) => 1,
        WeightLatticeElem(R, [3, 0, 0]) => 1,
        WeightLatticeElem(R, [1, 1, 0]) => 3,
        WeightLatticeElem(R, [0, 0, 1]) => 6,
        WeightLatticeElem(R, [1, 0, 0]) => 7,
      )
    end

    let L = lie_algebra(QQ, :D, 4), hw = [0, 3, 1, 0]
      domchar = check_dominant_character(L, hw)
      R = root_system(L)
      @test domchar == Dict(
        WeightLatticeElem(R, [0, 3, 1, 0]) => 1,
        WeightLatticeElem(R, [1, 1, 2, 1]) => 1,
        WeightLatticeElem(R, [1, 2, 0, 1]) => 2,
        WeightLatticeElem(R, [2, 0, 3, 0]) => 1,
        WeightLatticeElem(R, [2, 0, 1, 2]) => 2,
        WeightLatticeElem(R, [0, 0, 3, 2]) => 1,
        WeightLatticeElem(R, [2, 1, 1, 0]) => 3,
        WeightLatticeElem(R, [0, 1, 3, 0]) => 2,
        WeightLatticeElem(R, [0, 1, 1, 2]) => 3,
        WeightLatticeElem(R, [0, 2, 1, 0]) => 7,
        WeightLatticeElem(R, [3, 0, 0, 1]) => 4,
        WeightLatticeElem(R, [1, 0, 2, 1]) => 8,
        WeightLatticeElem(R, [1, 0, 0, 3]) => 4,
        WeightLatticeElem(R, [1, 1, 0, 1]) => 12,
        WeightLatticeElem(R, [2, 0, 1, 0]) => 16,
        WeightLatticeElem(R, [0, 0, 3, 0]) => 12,
        WeightLatticeElem(R, [0, 0, 1, 2]) => 16,
        WeightLatticeElem(R, [0, 1, 1, 0]) => 26,
        WeightLatticeElem(R, [1, 0, 0, 1]) => 36,
        WeightLatticeElem(R, [0, 0, 1, 0]) => 56,
      )
    end

    let R = root_system(:E, 6), hw = [1, 0, 1, 0, 1, 0]
      domchar = check_dominant_character(R, hw)
      @test domchar == Dict(
        WeightLatticeElem(R, [1, 0, 1, 0, 1, 0]) => 1,
        WeightLatticeElem(R, [0, 0, 0, 1, 1, 0]) => 2,
        WeightLatticeElem(R, [2, 1, 0, 0, 0, 1]) => 3,
        WeightLatticeElem(R, [0, 1, 1, 0, 0, 1]) => 6,
        WeightLatticeElem(R, [2, 0, 1, 0, 0, 0]) => 10,
        WeightLatticeElem(R, [1, 0, 0, 0, 1, 1]) => 16,
        WeightLatticeElem(R, [1, 2, 0, 0, 0, 0]) => 15,
        WeightLatticeElem(R, [0, 0, 2, 0, 0, 0]) => 20,
        WeightLatticeElem(R, [1, 0, 0, 1, 0, 0]) => 44,
        WeightLatticeElem(R, [0, 1, 0, 0, 0, 2]) => 36,
        WeightLatticeElem(R, [0, 1, 0, 0, 1, 0]) => 92,
        WeightLatticeElem(R, [2, 0, 0, 0, 0, 1]) => 104,
        WeightLatticeElem(R, [0, 0, 1, 0, 0, 1]) => 204,
        WeightLatticeElem(R, [1, 1, 0, 0, 0, 0]) => 425,
        WeightLatticeElem(R, [0, 0, 0, 0, 0, 2]) => 416,
        WeightLatticeElem(R, [0, 0, 0, 0, 1, 0]) => 836,
        WeightLatticeElem(R, [1, 0, 0, 0, 0, 0]) => 1600,
      )
    end

    let L = lie_algebra(QQ, :G, 2), hw = [1, 2]
      domchar = check_dominant_character(L, hw)
      R = root_system(L)
      @test domchar == Dict(
        WeightLatticeElem(R, [1, 2]) => 1,
        WeightLatticeElem(R, [4, 0]) => 1,
        WeightLatticeElem(R, [2, 1]) => 2,
        WeightLatticeElem(R, [0, 2]) => 2,
        WeightLatticeElem(R, [3, 0]) => 3,
        WeightLatticeElem(R, [1, 1]) => 5,
        WeightLatticeElem(R, [2, 0]) => 7,
        WeightLatticeElem(R, [0, 1]) => 7,
        WeightLatticeElem(R, [1, 0]) => 10,
        WeightLatticeElem(R, [0, 0]) => 10,
      )
    end

    let L = special_linear_lie_algebra(QQ, 2), hw = [7] # type A_1 but without known root system
      domchar = check_dominant_character(L, hw)
      R = root_system(L)
      @test domchar == Dict(
        WeightLatticeElem(R, [7]) => 1,
        WeightLatticeElem(R, [5]) => 1,
        WeightLatticeElem(R, [3]) => 1,
        WeightLatticeElem(R, [1]) => 1,
      )
    end

    let L = special_orthogonal_lie_algebra(QQ, 7), hw = ZZ.([1, 2, 0]) # type B_3 but without known root system
      domchar = check_dominant_character(L, hw)
      R = root_system(L)
      @test domchar == Dict(
        WeightLatticeElem(R, [1, 2, 0]) => 1,
        WeightLatticeElem(R, [2, 0, 2]) => 1,
        WeightLatticeElem(R, [2, 1, 0]) => 1,
        WeightLatticeElem(R, [0, 1, 2]) => 2,
        WeightLatticeElem(R, [0, 2, 0]) => 2,
        WeightLatticeElem(R, [3, 0, 0]) => 2,
        WeightLatticeElem(R, [1, 0, 2]) => 3,
        WeightLatticeElem(R, [1, 1, 0]) => 6,
        WeightLatticeElem(R, [2, 0, 0]) => 6,
        WeightLatticeElem(R, [0, 0, 2]) => 9,
        WeightLatticeElem(R, [0, 1, 0]) => 9,
        WeightLatticeElem(R, [1, 0, 0]) => 15,
        WeightLatticeElem(R, [0, 0, 0]) => 15,
      )
    end
  end

  @testset "character" begin
    function check_character(
      LR::Union{LieAlgebra,RootSystem}, hw::Vector{<:Oscar.IntegerUnion}
    )
      char = @inferred character(LR, hw)
      R = LR isa RootSystem ? LR : root_system(LR)
      @test char isa Dict{WeightLatticeElem,Int}
      @test char[WeightLatticeElem(R, hw)] == 1
      @test all(>=(1), values(char))
      @test sum(values(char)) == dim_of_simple_module(LR, hw)
      domchar = @inferred dominant_character(LR, hw)
      @test all(w -> domchar[w] == char[w], keys(domchar))

      char_ZZ = @inferred character(ZZRingElem, LR, hw)
      @test char_ZZ isa Dict{WeightLatticeElem,ZZRingElem}
      @test char_ZZ == Dict{WeightLatticeElem,ZZRingElem}(char)

      return char
    end

    # All concrete test results have been computed using the LiE CAS (http://wwwmathlabo.univ-poitiers.fr/~maavl/LiE/) v2.2.2
    let R = root_system(Tuple{Symbol,Int}[]), hw = Int[]
      domchar = check_character(R, hw)
      @test domchar == Dict(WeightLatticeElem(R, Int[]) => 1)
    end

    let L = lie_algebra(QQ, :A, 2), hw = [0, 0]
      domchar = check_character(L, hw)
      R = root_system(L)
      @test domchar == Dict(WeightLatticeElem(R, [0, 0]) => 1)
    end

    let L = lie_algebra(QQ, :A, 3), hw = ZZ.([1, 1, 0])
      char = check_character(L, hw)
      R = root_system(L)
      @test char == Dict(
        WeightLatticeElem(R, [1, 1, 0]) => 1,
        WeightLatticeElem(R, [2, -1, 1]) => 1,
        WeightLatticeElem(R, [-1, 2, 0]) => 1,
        WeightLatticeElem(R, [2, 0, -1]) => 1,
        WeightLatticeElem(R, [0, 0, 1]) => 2,
        WeightLatticeElem(R, [1, -2, 2]) => 1,
        WeightLatticeElem(R, [0, 1, -1]) => 2,
        WeightLatticeElem(R, [-2, 1, 1]) => 1,
        WeightLatticeElem(R, [1, -1, 0]) => 2,
        WeightLatticeElem(R, [-1, -1, 2]) => 1,
        WeightLatticeElem(R, [-2, 2, -1]) => 1,
        WeightLatticeElem(R, [1, 0, -2]) => 1,
        WeightLatticeElem(R, [-1, 0, 0]) => 2,
        WeightLatticeElem(R, [0, -2, 1]) => 1,
        WeightLatticeElem(R, [-1, 1, -2]) => 1,
        WeightLatticeElem(R, [0, -1, -1]) => 1,
      )
    end

    let L = lie_algebra(QQ, :C, 3), hw = [2, 0, 1]
      char = check_character(L, hw)
    end

    let R = root_system(:D, 4), hw = [0, 3, 1, 0]
      char = check_character(R, hw)
    end

    let L = lie_algebra(QQ, :E, 6), hw = [1, 0, 1, 0, 1, 0]
      char = check_character(L, hw)
    end

    let L = lie_algebra(QQ, :G, 2), hw = [3, 2]
      char = check_character(L, hw)
    end

    let L = special_linear_lie_algebra(QQ, 2), hw = [7] # type A_1 but without known root system
      char = check_character(L, hw)
      R = root_system(L)
      @test char == Dict(
        WeightLatticeElem(R, [7]) => 1,
        WeightLatticeElem(R, [5]) => 1,
        WeightLatticeElem(R, [3]) => 1,
        WeightLatticeElem(R, [1]) => 1,
        WeightLatticeElem(R, [-1]) => 1,
        WeightLatticeElem(R, [-3]) => 1,
        WeightLatticeElem(R, [-5]) => 1,
        WeightLatticeElem(R, [-7]) => 1,
      )
    end

    let L = special_orthogonal_lie_algebra(QQ, 7), hw = ZZ.([1, 2, 1]) # type B_3 but without known root system
      char = check_character(L, hw)
    end
  end

  @testset "tensor_product_decomposition" begin
    function test_tensor_product_decomposition(
      LR::Union{LieAlgebra,RootSystem},
      hw1::Vector{<:Oscar.IntegerUnion},
      hw2::Vector{<:Oscar.IntegerUnion},
    )
      dec = @inferred tensor_product_decomposition(LR, hw1, hw2)
      @test dec == @inferred tensor_product_decomposition(LR, hw2, hw1)
      @test multiplicity(dec, Int.(hw1 + hw2)) == 1
      dim_prod = dim_of_simple_module(LR, hw1) * dim_of_simple_module(LR, hw2)
      dim_dec = sum(multiplicity(dec, w) * dim_of_simple_module(LR, w) for w in unique(dec))
      @test dim_prod == dim_dec
      return dec
    end

    # All concrete test results have been computed using the LiE CAS (http://wwwmathlabo.univ-poitiers.fr/~maavl/LiE/) v2.2.2
    let R = root_system(Tuple{Symbol,Int}[]), hw = Int[]
      dec = test_tensor_product_decomposition(R, hw, hw)
      @test dec == multiset(Dict(Int[] => 1))
    end

    let L = lie_algebra(QQ, :A, 3), hw1 = [2, 1, 2], hw2 = [1, 0, 1]
      dec = test_tensor_product_decomposition(L, hw1, hw2)
      @test dec == multiset(
        Dict(
          [3, 1, 3] => 1,
          [3, 2, 1] => 1,
          [1, 2, 3] => 1,
          [4, 0, 2] => 1,
          [2, 0, 4] => 1,
          [1, 3, 1] => 1,
          [2, 1, 2] => 3,
          [2, 2, 0] => 1,
          [0, 2, 2] => 1,
          [3, 0, 1] => 1,
          [1, 0, 3] => 1,
          [1, 1, 1] => 1,
        ),
      )
    end

    let L = lie_algebra(QQ, :B, 4), hw1 = [1, 1, 0, 0], hw2 = [0, 1, 0, 1]
      dec = test_tensor_product_decomposition(L, hw1, hw2)
      @test dec == multiset(
        Dict(
          [1, 2, 0, 1] => 1,
          [2, 0, 1, 1] => 1,
          [0, 1, 1, 1] => 1,
          [2, 1, 0, 1] => 1,
          [1, 0, 0, 3] => 1,
          [0, 2, 0, 1] => 1,
          [1, 0, 1, 1] => 2,
          [3, 0, 0, 1] => 1,
          [1, 1, 0, 1] => 3,
          [0, 0, 0, 3] => 1,
          [0, 0, 1, 1] => 2,
          [2, 0, 0, 1] => 2,
          [0, 1, 0, 1] => 2,
          [1, 0, 0, 1] => 2,
          [0, 0, 0, 1] => 1,
        ),
      )
    end

    let L = lie_algebra(QQ, :C, 2), hw1 = [2, 2], hw2 = ZZ.([2, 0])
      dec = test_tensor_product_decomposition(L, hw1, hw2)
      @test dec == multiset(
        Dict(
          [4, 2] => 1,
          [2, 3] => 1,
          [4, 1] => 1,
          [0, 4] => 1,
          [2, 2] => 2,
          [4, 0] => 1,
          [0, 3] => 1,
          [2, 1] => 1,
          [0, 2] => 1,
        ),
      )
    end

    let L = lie_algebra(QQ, :D, 5), hw1 = ZZ.([1, 1, 3, 0, 2]), hw2 = [2, 1, 0, 2, 0]
      dec = test_tensor_product_decomposition(L, hw1, hw2)
    end

    let R = root_system(:E, 6), hw1 = [1, 1, 0, 0, 1, 2], hw2 = [2, 0, 1, 1, 0, 0]
      dec = test_tensor_product_decomposition(R, hw1, hw2)
    end

    let R = root_system(:G, 2), hw1 = [1, 3], hw2 = [5, 2]
      dec = test_tensor_product_decomposition(R, hw1, hw2)
    end

    let L = special_linear_lie_algebra(QQ, 2), hw1 = [7], hw2 = [2] # type A_1 but without known root system
      dec = test_tensor_product_decomposition(L, hw1, hw2)
      @test dec == multiset(Dict([9] => 1,
        [7] => 1,
        [5] => 1,
      ))
    end

    let L = special_orthogonal_lie_algebra(QQ, 7) # type B_3 but without known root system
      hw1 = ZZ.([1, 1, 0])
      hw2 = ZZ.([0, 1, 1])

      dec = test_tensor_product_decomposition(L, hw1, hw2)
      @test dec == multiset(
        Dict(
          [1, 2, 1] => 1,
          [2, 0, 3] => 1,
          [2, 1, 1] => 1,
          [0, 1, 3] => 1,
          [0, 2, 1] => 1,
          [3, 0, 1] => 1,
          [1, 0, 3] => 2,
          [1, 1, 1] => 3,
          [2, 0, 1] => 2,
          [0, 0, 3] => 2,
          [0, 1, 1] => 2,
          [1, 0, 1] => 2,
          [0, 0, 1] => 1,
        ),
      )
    end
  end

  @testset "demazure_operator" begin
    R = root_system(:B, 2)
    rho = weyl_vector(R)

    @test demazure_operator(simple_root(R, 1), rho) == Dict(
      WeightLatticeElem(R, [1, 1]) => 1,
      WeightLatticeElem(R, [-1, 3]) => 1,
    )

    @test demazure_operator(simple_root(R, 2), rho) == Dict(
      WeightLatticeElem(R, [1, 1]) => 1,
      WeightLatticeElem(R, [2, -1]) => 1,
    )

    @test demazure_operator(root(R, 3), rho) == Dict(
      WeightLatticeElem(R, [1, 1]) => 1,
      WeightLatticeElem(R, [-2, 1]) => 1,
      WeightLatticeElem(R, [0, 1]) => 1,
      WeightLatticeElem(R, [-1, 1]) => 1,
    )

    @test demazure_operator(simple_root(R, 1), -rho) == Dict()

    @test demazure_operator(simple_root(R, 2), -rho) == Dict()

    @test demazure_operator(root(R, 3), -rho) == Dict(
      WeightLatticeElem(R, [0, -1]) => -1,
      WeightLatticeElem(R, [1, -1]) => -1,
    )
  end

  @testset "demazure_character" begin
    function demazure_character_trivial_tests(R::RootSystem, w::WeightLatticeElem)
      W = weyl_group(R)

      char = demazure_character(R, w, longest_element(W))
      @test char == character(R, w)
      @test char isa Dict{WeightLatticeElem,Int}

      char = demazure_character(R, w, one(W))
      @test char == Dict(w => 1)
      @test char isa Dict{WeightLatticeElem,Int}
    end
    # The "correct" solution in the below tests has been computed using Magma V2.28-7.
    # Note that Magma considers right actions of Weyl groups on weights,
    # so to reproduce one needs to reverse all Weyl words when copying a test case to Magma.

    @testset "Type A" begin
      for n in 2:4
        R = root_system(:A, n)
        for i in 1:n
          w = fundamental_weight(R, i)
          demazure_character_trivial_tests(R, w)
        end
        demazure_character_trivial_tests(R, weyl_vector(R))
        demazure_character_trivial_tests(R, 2 * weyl_vector(R))
      end

      R = root_system(:A, 2)
      W = weyl_group(R)

      #x = epsilon
      char = demazure_character(R, fundamental_weight(R, 1), one(W))
      @test char == Dict(
        WeightLatticeElem(R, [1, 0]) => 1
      )

      char = demazure_character(R, fundamental_weight(R, 2), one(W))
      @test char == Dict(
        WeightLatticeElem(R, [0, 1]) => 1
      )

      char = demazure_character(R, weyl_vector(R), one(W))
      @test char == Dict(
        WeightLatticeElem(R, [1, 1]) => 1
      )

      #x = s[1]
      char = demazure_character(R, fundamental_weight(R, 1), W([1]))
      @test char == Dict(
        WeightLatticeElem(R, [-1, 1]) => 1,
        WeightLatticeElem(R, [1, 0]) => 1,
      )

      char = demazure_character(R, fundamental_weight(R, 2), W([1]))
      @test char == Dict(
        WeightLatticeElem(R, [0, 1]) => 1
      )

      char = demazure_character(R, weyl_vector(R), W([1]))
      @test char == Dict(
        WeightLatticeElem(R, [-1, 2]) => 1,
        WeightLatticeElem(R, [1, 1]) => 1,
      )

      #x = s[2]
      char = demazure_character(R, fundamental_weight(R, 1), W([2]))
      @test char == Dict(
        WeightLatticeElem(R, [1, 0]) => 1
      )

      char = demazure_character(R, fundamental_weight(R, 2), W([2]))
      @test char == Dict(
        WeightLatticeElem(R, [1, -1]) => 1,
        WeightLatticeElem(R, [0, 1]) => 1,
      )

      char = demazure_character(R, weyl_vector(R), W([2]))
      @test char == Dict(
        WeightLatticeElem(R, [2, -1]) => 1,
        WeightLatticeElem(R, [1, 1]) => 1,
      )

      #x = s[1]*s[2]
      char = demazure_character(R, fundamental_weight(R, 1), W([1, 2]))
      @test char == Dict(
        WeightLatticeElem(R, [0, -1]) => 1,
        WeightLatticeElem(R, [-1, 1]) => 1,
        WeightLatticeElem(R, [1, 0]) => 1,
      )

      char = demazure_character(R, fundamental_weight(R, 2), W([1, 2]))
      @test char == Dict(
        WeightLatticeElem(R, [1, -1]) => 1,
        WeightLatticeElem(R, [0, 1]) => 1,
      )

      char = demazure_character(R, weyl_vector(R), W([1, 2]))
      @test char == Dict(
        WeightLatticeElem(R, [2, -1]) => 1,
        WeightLatticeElem(R, [1, -2]) => 1,
        WeightLatticeElem(R, [-1, 2]) => 1,
        WeightLatticeElem(R, [0, 0]) => 1,
        WeightLatticeElem(R, [1, 1]) => 1,
      )

      #x = s[2]*s[1]
      char = demazure_character(R, fundamental_weight(R, 1), W([2, 1]))
      @test char == Dict(
        WeightLatticeElem(R, [-1, 1]) => 1,
        WeightLatticeElem(R, [1, 0]) => 1,
      )

      char = demazure_character(R, fundamental_weight(R, 2), W([2, 1]))
      @test char == Dict(
        WeightLatticeElem(R, [-1, 0]) => 1,
        WeightLatticeElem(R, [1, -1]) => 1,
        WeightLatticeElem(R, [0, 1]) => 1,
      )

      char = demazure_character(R, weyl_vector(R), W([2, 1]))
      @test char == Dict(
        WeightLatticeElem(R, [-1, 2]) => 1,
        WeightLatticeElem(R, [-2, 1]) => 1,
        WeightLatticeElem(R, [2, -1]) => 1,
        WeightLatticeElem(R, [0, 0]) => 1,
        WeightLatticeElem(R, [1, 1]) => 1,
      )

      #x = s[1]*s[2]*s[1]
      char = demazure_character(R, fundamental_weight(R, 1), W([1, 2, 1]))
      @test char == Dict(
        WeightLatticeElem(R, [-1, 1]) => 1,
        WeightLatticeElem(R, [0, -1]) => 1,
        WeightLatticeElem(R, [1, 0]) => 1,
      )

      char = demazure_character(R, fundamental_weight(R, 2), W([1, 2, 1]))
      @test char == Dict(
        WeightLatticeElem(R, [-1, 0]) => 1,
        WeightLatticeElem(R, [1, -1]) => 1,
        WeightLatticeElem(R, [0, 1]) => 1,
      )

      char = demazure_character(R, weyl_vector(R), W([1, 2, 1]))
      @test char == Dict(
        WeightLatticeElem(R, [-1, 2]) => 1,
        WeightLatticeElem(R, [-2, 1]) => 1,
        WeightLatticeElem(R, [2, -1]) => 1,
        WeightLatticeElem(R, [1, -2]) => 1,
        WeightLatticeElem(R, [-1, -1]) => 1,
        WeightLatticeElem(R, [0, 0]) => 2,
        WeightLatticeElem(R, [1, 1]) => 1,
      )
    end

    @testset "Type B" begin
      for n in 2:4
        R = root_system(:B, n)
        for i in 1:n
          w = fundamental_weight(R, i)
          demazure_character_trivial_tests(R, w)
        end
        demazure_character_trivial_tests(R, weyl_vector(R))
        demazure_character_trivial_tests(R, 2 * weyl_vector(R))
      end

      R = root_system(:B, 3)
      W = weyl_group(R)

      #x = epsilon
      char = demazure_character(R, fundamental_weight(R, 1), one(W))
      @test char == Dict(
        WeightLatticeElem(R, [1, 0, 0]) => 1
      )

      char = demazure_character(R, fundamental_weight(R, 2), one(W))
      @test char == Dict(
        WeightLatticeElem(R, [0, 1, 0]) => 1
      )

      char = demazure_character(R, fundamental_weight(R, 3), one(W))
      @test char == Dict(
        WeightLatticeElem(R, [0, 0, 1]) => 1
      )

      char = demazure_character(R, weyl_vector(R), one(W))
      @test char == Dict(
        WeightLatticeElem(R, [1, 1, 1]) => 1
      )

      #x = s[1]
      char = demazure_character(R, fundamental_weight(R, 1), W([1]))
      @test char == Dict(
        WeightLatticeElem(R, [-1, 1, 0]) => 1,
        WeightLatticeElem(R, [1, 0, 0]) => 1,
      )

      char = demazure_character(R, fundamental_weight(R, 2), W([1]))
      @test char == Dict(
        WeightLatticeElem(R, [0, 1, 0]) => 1
      )

      char = demazure_character(R, fundamental_weight(R, 3), W([1]))
      @test char == Dict(
        WeightLatticeElem(R, [0, 0, 1]) => 1
      )

      char = demazure_character(R, weyl_vector(R), W([1]))
      @test char == Dict(
        WeightLatticeElem(R, [-1, 2, 1]) => 1,
        WeightLatticeElem(R, [1, 1, 1]) => 1,
      )

      #x = s[2]
      char = demazure_character(R, fundamental_weight(R, 1), W([2]))
      @test char == Dict(
        WeightLatticeElem(R, [1, 0, 0]) => 1
      )

      char = demazure_character(R, fundamental_weight(R, 2), W([2]))
      @test char == Dict(
        WeightLatticeElem(R, [1, -1, 2]) => 1,
        WeightLatticeElem(R, [0, 1, 0]) => 1,
      )

      char = demazure_character(R, fundamental_weight(R, 3), W([2]))
      @test char == Dict(
        WeightLatticeElem(R, [0, 0, 1]) => 1
      )

      char = demazure_character(R, weyl_vector(R), W([2]))
      @test char == Dict(
        WeightLatticeElem(R, [2, -1, 3]) => 1,
        WeightLatticeElem(R, [1, 1, 1]) => 1,
      )

      #x = s[3]
      char = demazure_character(R, fundamental_weight(R, 1), W([3]))
      @test char == Dict(
        WeightLatticeElem(R, [1, 0, 0]) => 1
      )

      char = demazure_character(R, fundamental_weight(R, 2), W([3]))
      @test char == Dict(
        WeightLatticeElem(R, [0, 1, 0]) => 1
      )

      char = demazure_character(R, fundamental_weight(R, 3), W([3]))
      @test char == Dict(
        WeightLatticeElem(R, [0, 1, -1]) => 1,
        WeightLatticeElem(R, [0, 0, 1]) => 1,
      )

      char = demazure_character(R, weyl_vector(R), W([3]))
      @test char == Dict(
        WeightLatticeElem(R, [1, 2, -1]) => 1,
        WeightLatticeElem(R, [1, 1, 1]) => 1,
      )

      #x = s[1]*s[3]
      char = demazure_character(R, fundamental_weight(R, 1), W([1, 3]))
      @test char == Dict(
        WeightLatticeElem(R, [-1, 1, 0]) => 1,
        WeightLatticeElem(R, [1, 0, 0]) => 1,
      )

      char = demazure_character(R, fundamental_weight(R, 2), W([1, 3]))
      @test char == Dict(
        WeightLatticeElem(R, [0, 1, 0]) => 1
      )

      char = demazure_character(R, fundamental_weight(R, 3), W([1, 3]))
      @test char == Dict(
        WeightLatticeElem(R, [0, 1, -1]) => 1,
        WeightLatticeElem(R, [0, 0, 1]) => 1,
      )

      char = demazure_character(R, weyl_vector(R), W([1, 3]))
      @test char == Dict(
        WeightLatticeElem(R, [1, 1, 1]) => 1,
        WeightLatticeElem(R, [-1, 2, 1]) => 1,
        WeightLatticeElem(R, [-1, 3, -1]) => 1,
        WeightLatticeElem(R, [1, 2, -1]) => 1,
      )

      #x=s[1]*s[2]*s[3]
      char = demazure_character(R, fundamental_weight(R, 1), W([1, 2, 3]))
      @test char == Dict(
        WeightLatticeElem(R, [0, 0, 0]) => 1,
        WeightLatticeElem(R, [0, 1, -2]) => 1,
        WeightLatticeElem(R, [-1, 1, 0]) => 1,
        WeightLatticeElem(R, [1, 0, 0]) => 1,
        WeightLatticeElem(R, [0, -1, 2]) => 1,
      )

      char = demazure_character(R, fundamental_weight(R, 2), W([1, 2, 3]))
      @test char == Dict(
        WeightLatticeElem(R, [0, 1, 0]) => 1,
        WeightLatticeElem(R, [1, -1, 2]) => 1,
        WeightLatticeElem(R, [1, 0, 0]) => 1,
        WeightLatticeElem(R, [1, 1, -2]) => 1,
      )

      char = demazure_character(R, fundamental_weight(R, 3), W([1, 2, 3]))
      @test char == Dict(
        WeightLatticeElem(R, [0, 1, -1]) => 1,
        WeightLatticeElem(R, [0, 0, 1]) => 1,
      )

      char = demazure_character(R, weyl_vector(R), W([1, 2, 3]))
      @test char == Dict(
        WeightLatticeElem(R, [2, 0, 1]) => 1,
        WeightLatticeElem(R, [1, 2, -3]) => 1,
        WeightLatticeElem(R, [1, 1, -1]) => 1,
        WeightLatticeElem(R, [1, 1, 1]) => 1,
        WeightLatticeElem(R, [0, 3, -3]) => 1,
        WeightLatticeElem(R, [0, 1, 1]) => 1,
        WeightLatticeElem(R, [-1, 3, -1]) => 1,
        WeightLatticeElem(R, [2, -1, 3]) => 1,
        WeightLatticeElem(R, [0, 0, 3]) => 1,
        WeightLatticeElem(R, [1, -2, 5]) => 1,
        WeightLatticeElem(R, [0, 2, -1]) => 1,
        WeightLatticeElem(R, [-1, 2, 1]) => 1,
        WeightLatticeElem(R, [1, -1, 3]) => 1,
        WeightLatticeElem(R, [2, 2, -3]) => 1,
        WeightLatticeElem(R, [1, 2, -1]) => 1,
        WeightLatticeElem(R, [1, 0, 1]) => 1,
        WeightLatticeElem(R, [1, 3, -5]) => 1,
        WeightLatticeElem(R, [2, 1, -1]) => 1,
      )

      #x=s[3]*s[2]*s[3]*s[1]*s[2]*s[3]*s[1]*s[2]
      char = demazure_character(R, fundamental_weight(R, 1), W([3, 2, 3, 1, 2, 3, 1, 2]))
      @test char == Dict(
        WeightLatticeElem(R, [1, 0, 0]) => 1,
        WeightLatticeElem(R, [-1, 1, 0]) => 1,
        WeightLatticeElem(R, [0, -1, 2]) => 1,
        WeightLatticeElem(R, [0, 0, 0]) => 1,
        WeightLatticeElem(R, [1, -1, 0]) => 1,
        WeightLatticeElem(R, [0, 1, -2]) => 1,
      )

      char = demazure_character(R, fundamental_weight(R, 2), W([3, 2, 3, 1, 2, 3, 1, 2]))
      @test char == Dict(
        WeightLatticeElem(R, [-1, 2, -2]) => 1,
        WeightLatticeElem(R, [1, 1, -2]) => 1,
        WeightLatticeElem(R, [-2, 1, 0]) => 1,
        WeightLatticeElem(R, [2, -1, 0]) => 1,
        WeightLatticeElem(R, [0, -1, 0]) => 1,
        WeightLatticeElem(R, [-1, -1, 2]) => 1,
        WeightLatticeElem(R, [-1, 1, 0]) => 1,
        WeightLatticeElem(R, [1, -2, 2]) => 1,
        WeightLatticeElem(R, [1, 0, 0]) => 1,
        WeightLatticeElem(R, [-1, 0, 2]) => 1,
        WeightLatticeElem(R, [1, -1, 0]) => 1,
        WeightLatticeElem(R, [0, 1, 0]) => 1,
        WeightLatticeElem(R, [0, 1, -2]) => 1,
        WeightLatticeElem(R, [1, -1, 2]) => 1,
        WeightLatticeElem(R, [-1, 1, -2]) => 1,
        WeightLatticeElem(R, [0, -1, 2]) => 1,
        WeightLatticeElem(R, [0, 0, 0]) => 3,
        WeightLatticeElem(R, [1, 0, -2]) => 1,
        WeightLatticeElem(R, [-1, 0, 0]) => 1,
      )

      char = demazure_character(R, fundamental_weight(R, 3), W([3, 2, 3, 1, 2, 3, 1, 2]))
      @test char == Dict(
        WeightLatticeElem(R, [0, -1, 1]) => 1,
        WeightLatticeElem(R, [0, 0, 1]) => 1,
        WeightLatticeElem(R, [1, 0, -1]) => 1,
        WeightLatticeElem(R, [-1, 0, 1]) => 1,
        WeightLatticeElem(R, [0, 0, -1]) => 1,
        WeightLatticeElem(R, [-1, 1, -1]) => 1,
        WeightLatticeElem(R, [0, 1, -1]) => 1,
        WeightLatticeElem(R, [1, -1, 1]) => 1,
      )

      char = demazure_character(R, weyl_vector(R), W([3, 2, 3, 1, 2, 3, 1, 2]))
      @test length(char) == 124
      @test char[WeightLatticeElem(R, [-3, 3, -3])] == 1
      @test char[WeightLatticeElem(R, [0, 2, -3])] == 7
      @test char[WeightLatticeElem(R, [0, 1, -1])] == 12
      @test char[WeightLatticeElem(R, [3, 0, -3])] == 2
      @test char[WeightLatticeElem(R, [0, 1, -3])] == 6
    end

    @testset "Type C" begin
      for n in 2:4
        R = root_system(:C, n)
        for i in 1:n
          w = fundamental_weight(R, i)
          demazure_character_trivial_tests(R, w)
        end
        demazure_character_trivial_tests(R, weyl_vector(R))
        demazure_character_trivial_tests(R, 2 * weyl_vector(R))
      end
    end

    @testset "Type D" begin
      for n in 4:4
        R = root_system(:D, n)
        for i in 1:n
          w = fundamental_weight(R, i)
          demazure_character_trivial_tests(R, w)
        end
        demazure_character_trivial_tests(R, weyl_vector(R))
        demazure_character_trivial_tests(R, 2 * weyl_vector(R))
      end
    end

    @testset "Type E" begin
      for n in 6:6
        R = root_system(:E, n)
        for i in 1:n
          w = fundamental_weight(R, i)
          demazure_character_trivial_tests(R, w)
        end
        #demazure_character_trivial_tests(R, weyl_vector(R)) #takes too long for CI
        #demazure_character_trivial_tests(R, 2*weyl_vector(R)) #takes too long for CI
      end
    end

    @testset "Type F" begin
      for n in 4:4
        R = root_system(:F, n)
        for i in 1:n
          w = fundamental_weight(R, i)
          demazure_character_trivial_tests(R, w)
        end
        demazure_character_trivial_tests(R, weyl_vector(R))
        #demazure_character_trivial_tests(R, 2*weyl_vector(R)) #takes too long for CI
      end
    end

    @testset "Type G" begin
      for n in 2:2
        R = root_system(:G, n)
        for i in 1:n
          w = fundamental_weight(R, i)
          demazure_character_trivial_tests(R, w)
        end
        demazure_character_trivial_tests(R, weyl_vector(R))
        demazure_character_trivial_tests(R, 2 * weyl_vector(R))
      end
    end

    @testset "Input Conversion" begin
      R = root_system(:B, 3)
      W = weyl_group(R)
      L = lie_algebra(QQ, R)
      w_int = [0, 1, 0]
      w_weight = WeightLatticeElem(R, w_int)
      x = W([1, 2, 3])
      reduced_expression = word(x)
      result = Dict(
        WeightLatticeElem(R, [0, 1, 0]) => 1,
        WeightLatticeElem(R, [1, -1, 2]) => 1,
        WeightLatticeElem(R, [1, 0, 0]) => 1,
        WeightLatticeElem(R, [1, 1, -2]) => 1,
      )

      char = @inferred demazure_character(R, w_weight, x)
      @test char == result
      @test char isa Dict{WeightLatticeElem,Int}

      char_ZZ = @inferred demazure_character(ZZRingElem, R, w_weight, x)
      @test char_ZZ == result
      @test char_ZZ isa Dict{WeightLatticeElem,ZZRingElem}

      char = @inferred demazure_character(R, w_weight, reduced_expression)
      @test char == result
      @test char isa Dict{WeightLatticeElem,Int}

      char_ZZ = @inferred demazure_character(ZZRingElem, R, w_weight, reduced_expression)
      @test char_ZZ == result
      @test char_ZZ isa Dict{WeightLatticeElem,ZZRingElem}

      char = @inferred demazure_character(R, w_int, x)
      @test char == result
      @test char isa Dict{WeightLatticeElem,Int}

      char_ZZ = @inferred demazure_character(ZZRingElem, R, w_int, x)
      @test char_ZZ == result
      @test char_ZZ isa Dict{WeightLatticeElem,ZZRingElem}

      char = @inferred demazure_character(R, w_int, reduced_expression)
      @test char == result
      @test char isa Dict{WeightLatticeElem,Int}

      char_ZZ = @inferred demazure_character(ZZRingElem, R, w_int, reduced_expression)
      @test char_ZZ == result
      @test char_ZZ isa Dict{WeightLatticeElem,ZZRingElem}

      char = @inferred demazure_character(L, w_weight, x)
      @test char == result
      @test char isa Dict{WeightLatticeElem,Int}

      char_ZZ = @inferred demazure_character(ZZRingElem, L, w_weight, x)
      @test char_ZZ == result
      @test char_ZZ isa Dict{WeightLatticeElem,ZZRingElem}

      char = @inferred demazure_character(L, w_weight, reduced_expression)
      @test char == result
      @test char isa Dict{WeightLatticeElem,Int}

      char_ZZ = @inferred demazure_character(ZZRingElem, L, w_weight, reduced_expression)
      @test char_ZZ == result
      @test char_ZZ isa Dict{WeightLatticeElem,ZZRingElem}

      char = @inferred demazure_character(L, w_int, x)
      @test char == result
      @test char isa Dict{WeightLatticeElem,Int}

      char_ZZ = @inferred demazure_character(ZZRingElem, L, w_int, x)
      @test char_ZZ == result
      @test char_ZZ isa Dict{WeightLatticeElem,ZZRingElem}

      char = @inferred demazure_character(L, w_int, reduced_expression)
      @test char == result
      @test char isa Dict{WeightLatticeElem,Int}

      char_ZZ = @inferred demazure_character(ZZRingElem, L, w_int, reduced_expression)
      @test char_ZZ == result
      @test char_ZZ isa Dict{WeightLatticeElem,ZZRingElem}
    end
  end
end
