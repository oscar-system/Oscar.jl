
function lie_algebra_module_conformance_test(
  L::LieAlgebra{C},
  V::LieAlgebraModule{C},
  parentT::DataType=LieAlgebraModule{C},
  elemT::DataType=LieAlgebraModuleElem{C};
  num_random_tests::Int=10,
) where {C<:RingElement}
  @testset "basic manipulation" begin
    v = V(rand(-10:10, dim(V)))

    # @test parentT <: LieAlgebraModule{C}
    # @test elemT <: LieAlgebraModuleElem{C}
    @test V isa parentT
    @test v isa elemT

    @test parent_type(elemT) == parentT
    @test elem_type(parentT) == elemT

    @test parent(v) == V

    @test base_ring(v) == base_ring(V)
    @test elem_type(base_ring(V)) == C

    # this block stays only as long as `ngens` and `gens` are not specialized for Lie algebra modules
    @test dim(V) == ngens(V)
    @test basis(V) == gens(V)
    @test all(i -> basis(V, i) == gen(V, i), 1:dim(V))

    @test dim(V) == length(basis(V))
    @test all(i -> basis(V, i) == basis(V)[i], 1:dim(V))

    @test iszero(zero(V))

    @test sum(v[i] * basis(V, i) for i in 1:dim(V)) == v

    @test v == v
    @test deepcopy(v) == v
    @test hash(deepcopy(v)) == hash(v)
  end

  @testset "parent object call overload" begin
    @test V() == zero(V) == V(zeros(base_ring(V), dim(V)))

    for _ in 1:num_random_tests
      coeffs = rand(-10:10, dim(V))
      v1 = V(coeffs)
      v2 = V(base_ring(V).(coeffs))
      v3 = V(matrix(base_ring(V), 1, dim(V), coeffs))
      v4 = V(sparse_row(matrix(base_ring(V), 1, dim(V), coeffs)))
      v5 = V(v1)
      @test v1 == v2
      @test v1 == v3
      @test v1 == v4
      @test v1 == v5
    end
  end

  @testset "vector space axioms" begin
    for _ in 1:num_random_tests
      v = V(rand(-10:10, dim(V)))
      w = V(rand(-10:10, dim(V)))
      w2 = V(rand(-10:10, dim(V)))

      @test v + w == w + v
      @test v + (w + w2) == (v + w) + w2

      @test v + zero(V) == v
      @test zero(V) + v == v

      @test -v + v == zero(V)
      @test v + (-v) == zero(V)

      @test v - w == v + (-w)

      @test v * 0 == zero(V)
      @test 0 * v == zero(V)

      @test 2 * v == v + v
      @test v * 2 == v + v
      @test base_ring(V)(2) * v == v + v
      @test v * base_ring(V)(2) == v + v
    end
  end

  @testset "lie algebra action axioms" begin
    for _ in 1:num_random_tests
      x = L(rand(-10:10, dim(L)))
      y = L(rand(-10:10, dim(L)))
      v = V(rand(-10:10, dim(V)))
      w = V(rand(-10:10, dim(V)))

      @test (x * v) isa elemT
      @test parent(x * v) == parent(v)

      @test (x + y) * v == x * v + y * v
      @test x * (v + w) == x * v + x * w

      @test bracket(x, y) * v == x * (y * v) - y * (x * v)
    end
  end
end

@testset "PBWDeformations.LieAlgebraModule" begin
  sc = Matrix{SRow{elem_type(QQ)}}(undef, 3, 2)
  sc[1, 1] = sparse_row(QQ, [1, 2], [0, 0])
  sc[1, 2] = sparse_row(QQ, [1, 2], [1, 0])
  sc[2, 1] = sparse_row(QQ, [1, 2], [0, 1])
  sc[2, 2] = sparse_row(QQ, [1, 2], [0, 0])
  sc[3, 1] = sparse_row(QQ, [1, 2], [1, 0])
  sc[3, 2] = sparse_row(QQ, [1, 2], [0, -1])

  @testset "conformance tests" begin
    @testset "V of sl_2(QQ) using structure constants" begin
      L = special_linear_lie_algebra(QQ, 2)
      V = abstract_module(L, 2, sc)
      lie_algebra_module_conformance_test(L, V)
    end

    @testset "λ = [1,1,0] of sl_4(QQ)" begin
      L = special_linear_lie_algebra(QQ, 4)
      V = highest_weight_module(L, [1, 1, 0])
      lie_algebra_module_conformance_test(L, V)
    end

    # does not work (some GAP thing)
    # @testset "λ = [1,1,0,0] of so_4(QQ)" begin
    #   L = special_orthogonal_lie_algebra(QQ, 4)
    #   V = highest_weight_module(L, [1, 1, 0, 0])
    #   lie_algebra_module_conformance_test(L, V)
    # end

    @testset "λ = [0,1] of A_2(QQ)" begin
      L = lie_algebra(QQ, ('A', 2))
      V = highest_weight_module(L, [0, 1])
      lie_algebra_module_conformance_test(L, V)
    end

    @testset "λ = [0,1,1] of B_3(QQ)" begin
      L = lie_algebra(QQ, ('B', 3))
      V = highest_weight_module(L, [0, 1, 1])
      lie_algebra_module_conformance_test(L, V)
    end

    @testset "V of sl_4(QQ)" begin
      L = special_linear_lie_algebra(QQ, 4)
      V = standard_module(L)
      lie_algebra_module_conformance_test(L, V)
    end

    @testset "V of so_4(CL(4))" begin
      L = special_orthogonal_lie_algebra(cyclotomic_field(4)[1], 4)
      V = standard_module(L)
      lie_algebra_module_conformance_test(L, V)
    end

    @testset "⋀^2 S^2 V of so_4(QQ)" begin
      L = special_orthogonal_lie_algebra(QQ, 4)
      V = exterior_power(symmetric_power(standard_module(L), 2), 2)
      lie_algebra_module_conformance_test(L, V)
    end

    @testset "⋀^2 V of so_4(CL(4))" begin
      L = special_orthogonal_lie_algebra(cyclotomic_field(4)[1], 4)
      V = exterior_power(standard_module(L), 3)
      lie_algebra_module_conformance_test(L, V)
    end

    @testset "S^2 ⋀^2 V of so_4(QQ)" begin
      L = special_orthogonal_lie_algebra(QQ, 4)
      V = tensor_power(symmetric_power(standard_module(L), 2), 2)
      lie_algebra_module_conformance_test(L, V)
    end

    @testset "S^2 V of so_4(CL(4))" begin
      L = special_orthogonal_lie_algebra(cyclotomic_field(4)[1], 4)
      V = symmetric_power(standard_module(L), 2)
      lie_algebra_module_conformance_test(L, V)
    end

    @testset "T^2 ⋀^2 V of sl_4(QQ)" begin
      L = special_linear_lie_algebra(QQ, 4)
      V = tensor_power(exterior_power(standard_module(L), 2), 2)
      lie_algebra_module_conformance_test(L, V)
    end

    @testset "T^2 V of sl_4(CL(4))" begin
      L = special_linear_lie_algebra(cyclotomic_field(4)[1], 4)
      V = tensor_power(standard_module(L), 2)
      lie_algebra_module_conformance_test(L, V)
    end
  end

  @testset "standard_module" begin
    L = special_orthogonal_lie_algebra(QQ, 4)
    V = standard_module(L)
    @test dim(V) == 4

    x = L(rand(-10:10, dim(L)))
    v = V(rand(-10:10, dim(V)))
    @test x * v == V(transpose(matrix_repr(x) * transpose(Generic._matrix(v))))
  end

  @testset "exterior_power" begin
    L = special_orthogonal_lie_algebra(QQ, 4)
    V = symmetric_power(standard_module(L), 2)
    pow_V = exterior_power(V, 2)

    @test dim(pow_V) == binomial(dim(V), 2)
    a = gen(V, 1)
    b = gen(V, 2)
    @test !iszero(pow_V([a, b]))
    @test iszero(pow_V([a, b]) + pow_V([b, a]))
    @test !iszero(pow_V([a, b]) - pow_V([b, a]))
    @test iszero(pow_V([a, a]))
  end

  @testset "symmetric_power" begin
    L = special_orthogonal_lie_algebra(QQ, 4)
    V = exterior_power(standard_module(L), 2)
    pow_V = symmetric_power(V, 2)

    @test dim(pow_V) == binomial(dim(V) + 2 - 1, 2)
    a = gen(V, 1)
    b = gen(V, 2)
    @test !iszero(pow_V([a, b]))
    @test !iszero(pow_V([a, b]) + pow_V([b, a]))
    @test iszero(pow_V([a, b]) - pow_V([b, a]))
    @test !iszero(pow_V([a, a]))
  end

  @testset "tensor_power" begin
    L = special_orthogonal_lie_algebra(QQ, 4)
    V = standard_module(L)
    pow_V = tensor_power(V, 2)

    @test dim(pow_V) == dim(V)^2
    a = gen(V, 1)
    b = gen(V, 2)
    @test !iszero(pow_V([a, b]))
    @test !iszero(pow_V([a, b]) + pow_V([b, a]))
    @test !iszero(pow_V([a, b]) - pow_V([b, a]))
    @test !iszero(pow_V([a, a]))
  end
end
