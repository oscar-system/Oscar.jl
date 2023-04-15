
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

    @test base_lie_algebra(V) === L

    # this block stays only as long as `ngens` and `gens` are not specialized for Lie algebra modules
    @test dim(V) == ngens(V)
    @test basis(V) == gens(V)
    @test all(i -> basis(V, i) == gen(V, i), 1:dim(V))

    @test dim(V) == length(basis(V))
    @test all(i -> basis(V, i) == basis(V)[i], 1:dim(V))

    @test iszero(zero(V))

    @test coefficients(v) == [coeff(v, i) for i in 1:dim(V)]
    @test all(i -> coeff(v, i) == v[i], 1:dim(V))
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

      @test (x * y) * v == x * (y * v) - y * (x * v)
    end
  end
end

@testset "LieAlgebras.LieAlgebraModule" begin
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

    @testset "λ = [1,1,0,0] of so_4(CF(4))" begin
      L = special_orthogonal_lie_algebra(cyclotomic_field(4)[1], 4)
      V = highest_weight_module(L, [1, 1, 0, 0])
      lie_algebra_module_conformance_test(L, V)
    end

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

    @testset "V of so_4(CF(4))" begin
      L = special_orthogonal_lie_algebra(cyclotomic_field(4)[1], 4)
      V = standard_module(L)
      lie_algebra_module_conformance_test(L, V)
    end

    @testset "⋀^2 S^2 V of so_4(QQ)" begin
      L = special_orthogonal_lie_algebra(QQ, 4)
      V = exterior_power(symmetric_power(standard_module(L), 2), 2)
      lie_algebra_module_conformance_test(L, V)
    end

    @testset "⋀^2 V of so_4(CF(4))" begin
      L = special_orthogonal_lie_algebra(cyclotomic_field(4)[1], 4)
      V = exterior_power(standard_module(L), 3)
      lie_algebra_module_conformance_test(L, V)
    end

    @testset "S^2 ⋀^2 V of so_4(QQ)" begin
      L = special_orthogonal_lie_algebra(QQ, 4)
      V = tensor_power(symmetric_power(standard_module(L), 2), 2)
      lie_algebra_module_conformance_test(L, V)
    end

    @testset "S^2 V of so_4(CF(4))" begin
      L = special_orthogonal_lie_algebra(cyclotomic_field(4)[1], 4)
      V = symmetric_power(standard_module(L), 2)
      lie_algebra_module_conformance_test(L, V)
    end

    @testset "T^2 ⋀^2 V of sl_4(QQ)" begin
      L = special_linear_lie_algebra(QQ, 4)
      V = tensor_power(exterior_power(standard_module(L), 2), 2)
      lie_algebra_module_conformance_test(L, V)
    end

    @testset "T^2 V of sl_4(CF(4))" begin
      L = special_linear_lie_algebra(cyclotomic_field(4)[1], 4)
      V = tensor_power(standard_module(L), 2)
      lie_algebra_module_conformance_test(L, V)
    end
  end

  module_type_bools(V) =
    (is_standard_module(V), is_exterior_power(V), is_symmetric_power(V), is_tensor_power(V))

  @testset "standard_module" begin
    L = special_orthogonal_lie_algebra(QQ, 4)
    V = standard_module(L)
    @test dim(V) == 4
    @test length(repr(V)) < 10^4 # outputs tend to be excessively long due to recursion

    @test is_standard_module(V)
    @test !is_dual(V)
    @test !is_exterior_power(V)
    @test !is_symmetric_power(V)
    @test !is_tensor_power(V)

    x = L(rand(-10:10, dim(L)))
    v = V(rand(-10:10, dim(V)))
    @test x * v == V(matrix_repr(x) * coefficients(v))
  end

  @testset "dual" begin
    L = special_orthogonal_lie_algebra(QQ, 4)
    V = symmetric_power(standard_module(L), 2)
    type_V = module_type_bools(V)

    dual_V = dual(V)
    @test type_V == module_type_bools(V) # construction of dual_V should not change type of V
    @test base_module(dual_V) === V
    @test dim(dual_V) == dim(V)
    @test length(repr(dual_V)) < 10^4 # outputs tend to be excessively long due to recursion

    @test !is_standard_module(dual_V)
    @test is_dual(dual_V)
    @test !is_exterior_power(dual_V)
    @test !is_symmetric_power(dual_V)
    @test !is_tensor_power(dual_V)

    dual_dual_V = dual(dual_V)
    @test dim(dual_dual_V) == dim(V)
    @test all(
      i ->
        Oscar.LieAlgebras.transformation_matrix(dual_dual_V, i) ==
        Oscar.LieAlgebras.transformation_matrix(V, i),
      1:dim(L),
    )
  end

  @testset "exterior_power" begin
    L = special_orthogonal_lie_algebra(QQ, 4)
    V = symmetric_power(standard_module(L), 2)
    type_V = module_type_bools(V)

    for k in 1:3
      pow_V = exterior_power(V, k)
      @test type_V == module_type_bools(V) # construction of pow_V should not change type of V
      @test base_module(pow_V) === V
      @test dim(pow_V) == binomial(dim(V), k)
      @test length(repr(pow_V)) < 10^4 # outputs tend to be excessively long due to recursion

      @test !is_standard_module(pow_V)
      @test !is_dual(V)
      @test is_exterior_power(pow_V)
      @test !is_symmetric_power(pow_V)
      @test !is_tensor_power(pow_V)

      if k == 1
        x = L(rand(-10:10, dim(L)))
        a = V(rand(-10:10, dim(V)))
        @test pow_V([x * a]) == x * pow_V([a])
      elseif k == 2
        a = V(rand(-10:10, dim(V)))
        b = V(rand(-10:10, dim(V)))
        @test !iszero(pow_V([a, b]))
        @test iszero(pow_V([a, b]) + pow_V([b, a]))
        @test !iszero(pow_V([a, b]) - pow_V([b, a]))
        @test iszero(pow_V([a, a]))
      end
    end
  end

  @testset "symmetric_power" begin
    L = special_orthogonal_lie_algebra(QQ, 4)
    V = exterior_power(standard_module(L), 2)
    type_V = module_type_bools(V)

    for k in 1:3
      pow_V = symmetric_power(V, k)
      @test type_V == module_type_bools(V) # construction of pow_V should not change type of V
      @test base_module(pow_V) === V
      @test dim(pow_V) == binomial(dim(V) + k - 1, k)
      @test length(repr(pow_V)) < 10^4 # outputs tend to be excessively long due to recursion

      @test !is_standard_module(pow_V)
      @test !is_dual(V)
      @test !is_exterior_power(pow_V)
      @test is_symmetric_power(pow_V)
      @test !is_tensor_power(pow_V)

      if k == 1
        x = L(rand(-10:10, dim(L)))
        a = V(rand(-10:10, dim(V)))
        @test pow_V([x * a]) == x * pow_V([a])
      elseif k == 2
        a = V(rand(-10:10, dim(V)))
        b = V(rand(-10:10, dim(V)))
        @test !iszero(pow_V([a, b]))
        @test !iszero(pow_V([a, b]) + pow_V([b, a]))
        @test iszero(pow_V([a, b]) - pow_V([b, a]))
        @test !iszero(pow_V([a, a]))
      end
    end
  end

  @testset "tensor_power" begin
    L = special_orthogonal_lie_algebra(QQ, 4)
    V = standard_module(L)
    type_V = module_type_bools(V)

    for k in 1:3
      pow_V = tensor_power(V, k)
      @test type_V == module_type_bools(V) # construction of pow_V should not change type of V
      @test base_module(pow_V) === V
      @test dim(pow_V) == dim(V)^k
      @test length(repr(pow_V)) < 10^4 # outputs tend to be excessively long due to recursion

      @test !is_standard_module(pow_V)
      @test !is_dual(V)
      @test !is_exterior_power(pow_V)
      @test !is_symmetric_power(pow_V)
      @test is_tensor_power(pow_V)

      if k == 1
        x = L(rand(-10:10, dim(L)))
        a = V(rand(-10:10, dim(V)))
        @test pow_V([x * a]) == x * pow_V([a])
      elseif k == 2
        a = V(rand(-10:10, dim(V)))
        b = V(rand(-10:10, dim(V)))
        @test !iszero(pow_V([a, b]))
        @test !iszero(pow_V([a, b]) + pow_V([b, a]))
        @test !iszero(pow_V([a, b]) - pow_V([b, a]))
        @test !iszero(pow_V([a, a]))
      end
    end
  end

  @testset "so_n correctness regression" begin
    function lie_algebra_module_struct_const(
      L::LieAlgebra{C}, V::LieAlgebraModule{C}
    ) where {C<:RingElement}
      R = base_ring(L)
      dimL = dim(L)
      dimV = dim(V)
      struct_const_V = Matrix{Vector{Tuple{elem_type(R),Int}}}(undef, dimL, dimV)
      for (i, xi) in enumerate(basis(L)), (j, vj) in enumerate(basis(V))
        struct_const_V[i, j] = [
          (c, k) for (k, c) in enumerate(coefficients(xi * vj)) if !iszero(c)
        ]
      end
      return struct_const_V
    end

    L = special_orthogonal_lie_algebra(QQ, 3)

    struct_const_V = Matrix{Vector{Tuple{QQFieldElem,Int64}}}(undef, 3, 3)
    struct_const_V[1, :] = Vector{Tuple{QQFieldElem,Int64}}[[(QQ(-1), 2)], [(QQ(1), 1)], []]
    struct_const_V[2, :] = Vector{Tuple{QQFieldElem,Int64}}[[(QQ(-1), 3)], [], [(QQ(1), 1)]]
    struct_const_V[3, :] = Vector{Tuple{QQFieldElem,Int64}}[[], [(QQ(-1), 3)], [(QQ(1), 2)]]
    @test lie_algebra_module_struct_const(L, standard_module(L)) == struct_const_V

    struct_const_V = Matrix{Vector{Tuple{QQFieldElem,Int64}}}(undef, 3, 3)
    struct_const_V[1, :] = Vector{Tuple{QQFieldElem,Int64}}[[(QQ(-1), 2)], [(QQ(1), 1)], []]
    struct_const_V[2, :] = Vector{Tuple{QQFieldElem,Int64}}[[(QQ(-1), 3)], [], [(QQ(1), 1)]]
    struct_const_V[3, :] = Vector{Tuple{QQFieldElem,Int64}}[[], [(QQ(-1), 3)], [(QQ(1), 2)]]
    @test lie_algebra_module_struct_const(L, symmetric_power(standard_module(L), 1)) ==
      struct_const_V

    struct_const_V = Matrix{Vector{Tuple{QQFieldElem,Int64}}}(undef, 3, 6)
    struct_const_V[1, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(-2), 2)],
      [(QQ(1), 1), (QQ(-1), 4)],
      [(QQ(-1), 5)],
      [(QQ(2), 2)],
      [(QQ(1), 3)],
      [],
    ]
    struct_const_V[2, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(-2), 3)],
      [(QQ(-1), 5)],
      [(QQ(1), 1), (QQ(-1), 6)],
      [],
      [(QQ(1), 2)],
      [(QQ(2), 3)],
    ]
    struct_const_V[3, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [],
      [(QQ(-1), 3)],
      [(QQ(1), 2)],
      [(QQ(-2), 5)],
      [(QQ(1), 4), (QQ(-1), 6)],
      [(QQ(2), 5)],
    ]
    @test lie_algebra_module_struct_const(L, symmetric_power(standard_module(L), 2)) ==
      struct_const_V

    struct_const_V = Matrix{Vector{Tuple{QQFieldElem,Int64}}}(undef, 3, 10)
    struct_const_V[1, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(-3), 2)],
      [(QQ(1), 1), (QQ(-2), 4)],
      [(QQ(-2), 5)],
      [(QQ(2), 2), (QQ(-1), 7)],
      [(QQ(1), 3), (QQ(-1), 8)],
      [(QQ(-1), 9)],
      [(QQ(3), 4)],
      [(QQ(2), 5)],
      [(QQ(1), 6)],
      [],
    ]
    struct_const_V[2, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(-3), 3)],
      [(QQ(-2), 5)],
      [(QQ(1), 1), (QQ(-2), 6)],
      [(QQ(-1), 8)],
      [(QQ(1), 2), (QQ(-1), 9)],
      [(QQ(2), 3), (QQ(-1), 10)],
      [],
      [(QQ(1), 4)],
      [(QQ(2), 5)],
      [(QQ(3), 6)],
    ]
    struct_const_V[3, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [],
      [(QQ(-1), 3)],
      [(QQ(1), 2)],
      [(QQ(-2), 5)],
      [(QQ(1), 4), (QQ(-1), 6)],
      [(QQ(2), 5)],
      [(QQ(-3), 8)],
      [(QQ(1), 7), (QQ(-2), 9)],
      [(QQ(2), 8), (QQ(-1), 10)],
      [(QQ(3), 9)],
    ]
    @test lie_algebra_module_struct_const(L, symmetric_power(standard_module(L), 4)) ==
      struct_const_V

    struct_const_V = Matrix{Vector{Tuple{QQFieldElem,Int64}}}(undef, 3, 3)
    struct_const_V[1, :] = Vector{Tuple{QQFieldElem,Int64}}[[(QQ(-1), 2)], [(QQ(1), 1)], []]
    struct_const_V[2, :] = Vector{Tuple{QQFieldElem,Int64}}[[(QQ(-1), 3)], [], [(QQ(1), 1)]]
    struct_const_V[3, :] = Vector{Tuple{QQFieldElem,Int64}}[[], [(QQ(-1), 3)], [(QQ(1), 2)]]
    @test lie_algebra_module_struct_const(L, exterior_power(standard_module(L), 2)) ==
      struct_const_V

    struct_const_V = Matrix{Vector{Tuple{QQFieldElem,Int64}}}(undef, 3, 3)
    struct_const_V[1, :] = Vector{Tuple{QQFieldElem,Int64}}[[], [(QQ(-1), 3)], [(QQ(1), 2)]]
    struct_const_V[2, :] = Vector{Tuple{QQFieldElem,Int64}}[[(QQ(1), 3)], [], [(QQ(-1), 1)]]
    struct_const_V[3, :] = Vector{Tuple{QQFieldElem,Int64}}[[(QQ(-1), 2)], [(QQ(1), 1)], []]
    @test lie_algebra_module_struct_const(L, exterior_power(standard_module(L), 3)) ==
      struct_const_V

    L = special_orthogonal_lie_algebra(QQ, 4)

    struct_const_V = Matrix{Vector{Tuple{QQFieldElem,Int64}}}(undef, 6, 4)
    struct_const_V[1, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(-1), 2)], [(QQ(1), 1)], [], []
    ]
    struct_const_V[2, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(-1), 3)], [], [(QQ(1), 1)], []
    ]
    struct_const_V[3, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(-1), 4)], [], [], [(QQ(1), 1)]
    ]
    struct_const_V[4, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [], [(QQ(-1), 3)], [(QQ(1), 2)], []
    ]
    struct_const_V[5, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [], [(QQ(-1), 4)], [], [(QQ(1), 2)]
    ]
    struct_const_V[6, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [], [], [(QQ(-1), 4)], [(QQ(1), 3)]
    ]
    @test lie_algebra_module_struct_const(L, standard_module(L)) == struct_const_V

    struct_const_V = Matrix{Vector{Tuple{QQFieldElem,Int64}}}(undef, 6, 4)
    struct_const_V[1, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(-1), 2)], [(QQ(1), 1)], [], []
    ]
    struct_const_V[2, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(-1), 3)], [], [(QQ(1), 1)], []
    ]
    struct_const_V[3, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(-1), 4)], [], [], [(QQ(1), 1)]
    ]
    struct_const_V[4, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [], [(QQ(-1), 3)], [(QQ(1), 2)], []
    ]
    struct_const_V[5, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [], [(QQ(-1), 4)], [], [(QQ(1), 2)]
    ]
    struct_const_V[6, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [], [], [(QQ(-1), 4)], [(QQ(1), 3)]
    ]
    @test lie_algebra_module_struct_const(L, symmetric_power(standard_module(L), 1)) ==
      struct_const_V

    struct_const_V = Matrix{Vector{Tuple{QQFieldElem,Int64}}}(undef, 6, 4)
    struct_const_V[1, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(-1), 2)], [(QQ(1), 1)], [], []
    ]
    struct_const_V[2, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(-1), 3)], [], [(QQ(1), 1)], []
    ]
    struct_const_V[3, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(-1), 4)], [], [], [(QQ(1), 1)]
    ]
    struct_const_V[4, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [], [(QQ(-1), 3)], [(QQ(1), 2)], []
    ]
    struct_const_V[5, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [], [(QQ(-1), 4)], [], [(QQ(1), 2)]
    ]
    struct_const_V[6, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [], [], [(QQ(-1), 4)], [(QQ(1), 3)]
    ]
    @test lie_algebra_module_struct_const(L, exterior_power(standard_module(L), 1)) ==
      struct_const_V

    struct_const_V = Matrix{Vector{Tuple{QQFieldElem,Int64}}}(undef, 6, 6)
    struct_const_V[1, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [], [(QQ(-1), 4)], [(QQ(-1), 5)], [(QQ(1), 2)], [(QQ(1), 3)], []
    ]
    struct_const_V[2, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(1), 4)], [], [(QQ(-1), 6)], [(QQ(-1), 1)], [], [(QQ(1), 3)]
    ]
    struct_const_V[3, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(1), 5)], [(QQ(1), 6)], [], [], [(QQ(-1), 1)], [(QQ(-1), 2)]
    ]
    struct_const_V[4, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(-1), 2)], [(QQ(1), 1)], [], [], [(QQ(-1), 6)], [(QQ(1), 5)]
    ]
    struct_const_V[5, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(-1), 3)], [], [(QQ(1), 1)], [(QQ(1), 6)], [], [(QQ(-1), 4)]
    ]
    struct_const_V[6, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [], [(QQ(-1), 3)], [(QQ(1), 2)], [(QQ(-1), 5)], [(QQ(1), 4)], []
    ]
    @test lie_algebra_module_struct_const(L, exterior_power(standard_module(L), 2)) ==
      struct_const_V

    struct_const_V = Matrix{Vector{Tuple{QQFieldElem,Int64}}}(undef, 6, 4)
    struct_const_V[1, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [], [], [(QQ(-1), 4)], [(QQ(1), 3)]
    ]
    struct_const_V[2, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [], [(QQ(1), 4)], [], [(QQ(-1), 2)]
    ]
    struct_const_V[3, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(-1), 4)], [], [], [(QQ(1), 1)]
    ]
    struct_const_V[4, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [], [(QQ(-1), 3)], [(QQ(1), 2)], []
    ]
    struct_const_V[5, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(1), 3)], [], [(QQ(-1), 1)], []
    ]
    struct_const_V[6, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(-1), 2)], [(QQ(1), 1)], [], []
    ]
    @test lie_algebra_module_struct_const(L, exterior_power(standard_module(L), 3)) ==
      struct_const_V
  end
end
