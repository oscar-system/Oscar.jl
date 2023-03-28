@testset "PBWDeformations.LieAlgebraAbstractModule" begin
  R = QQ
  sc = Matrix{SRow{elem_type(R)}}(undef, 3, 2)
  sc[1, 1] = sparse_row(R, [1, 2], [0, 0])
  sc[1, 2] = sparse_row(R, [1, 2], [1, 0])
  sc[2, 1] = sparse_row(R, [1, 2], [0, 1])
  sc[2, 2] = sparse_row(R, [1, 2], [0, 0])
  sc[3, 1] = sparse_row(R, [1, 2], [1, 0])
  sc[3, 2] = sparse_row(R, [1, 2], [0, -1])

  @testset "conformance tests" begin
    @testset "V of sl_2(QQ) using structure constants" begin
      L = special_linear_liealgebra(QQ, 2)
      V = abstract_module(L, 2, sc)
      liealgebra_module_conformance_test(
        L,
        V,
        LieAlgebraAbstractModule{QQFieldElem},
        LieAlgebraAbstractModuleElem{QQFieldElem},
      )
    end

    @testset "λ = [1,1,0] of sl_4(QQ)" begin
      L = special_linear_liealgebra(QQ, 4)
      V = highest_weight_module(L, [1, 1, 0])
      liealgebra_module_conformance_test(
        L,
        V,
        LieAlgebraAbstractModule{QQFieldElem},
        LieAlgebraAbstractModuleElem{QQFieldElem},
      )
    end

    # does not work (some GAP thing)
    # @testset "λ = [1,1,0,0] of so_4(QQ)" begin
    #   L = special_orthogonal_liealgebra(QQ, 4)
    #   V = highest_weight_module(L, [1, 1, 0, 0])
    #   liealgebra_module_conformance_test(L, V, LieAlgebraAbstractModule{QQFieldElem}, LieAlgebraAbstractModuleElem{QQFieldElem})
    # end

    @testset "λ = [0,1] of A_2(QQ)" begin
      L = liealgebra(QQ, ('A', 2))
      V = highest_weight_module(L, [0, 1])
      liealgebra_module_conformance_test(
        L,
        V,
        LieAlgebraAbstractModule{QQFieldElem},
        LieAlgebraAbstractModuleElem{QQFieldElem},
      )
    end

    @testset "λ = [0,1,1] of B_3(QQ)" begin
      L = liealgebra(QQ, ('B', 3))
      V = highest_weight_module(L, [0, 1, 1])
      liealgebra_module_conformance_test(
        L,
        V,
        LieAlgebraAbstractModule{QQFieldElem},
        LieAlgebraAbstractModuleElem{QQFieldElem},
      )
    end

    @testset "V of sl_4(QQ)" begin
      L = special_linear_liealgebra(QQ, 4)
      V = standard_module(L)
      liealgebra_module_conformance_test(
        L,
        V,
        LieAlgebraAbstractModule{QQFieldElem},
        LieAlgebraAbstractModuleElem{QQFieldElem},
      )
    end

    @testset "V of so_4(CL(4))" begin
      L = special_orthogonal_liealgebra(cyclotomic_field(4)[1], 4)
      V = standard_module(L)
      liealgebra_module_conformance_test(
        L, V, LieAlgebraAbstractModule{nf_elem}, LieAlgebraAbstractModuleElem{nf_elem}
      )
    end

    @testset "⋀^2 S^2 V of so_4(QQ)" begin
      L = special_orthogonal_liealgebra(QQ, 4)
      V = exterior_power(symmetric_power(standard_module(L), 2), 2)
      liealgebra_module_conformance_test(
        L,
        V,
        LieAlgebraAbstractModule{QQFieldElem},
        LieAlgebraAbstractModuleElem{QQFieldElem},
      )
    end

    @testset "⋀^2 V of so_4(CL(4))" begin
      L = special_orthogonal_liealgebra(cyclotomic_field(4)[1], 4)
      V = exterior_power(standard_module(L), 3)
      liealgebra_module_conformance_test(
        L, V, LieAlgebraAbstractModule{nf_elem}, LieAlgebraAbstractModuleElem{nf_elem}
      )
    end

    @testset "S^2 ⋀^2 V of so_4(QQ)" begin
      L = special_orthogonal_liealgebra(QQ, 4)
      V = tensor_power(symmetric_power(standard_module(L), 2), 2)
      liealgebra_module_conformance_test(
        L,
        V,
        LieAlgebraAbstractModule{QQFieldElem},
        LieAlgebraAbstractModuleElem{QQFieldElem},
      )
    end

    @testset "S^2 V of so_4(CL(4))" begin
      L = special_orthogonal_liealgebra(cyclotomic_field(4)[1], 4)
      V = symmetric_power(standard_module(L), 2)
      liealgebra_module_conformance_test(
        L, V, LieAlgebraAbstractModule{nf_elem}, LieAlgebraAbstractModuleElem{nf_elem}
      )
    end

    @testset "T^2 ⋀^2 V of sl_4(QQ)" begin
      L = special_linear_liealgebra(QQ, 4)
      V = tensor_power(exterior_power(standard_module(L), 2), 2)
      liealgebra_module_conformance_test(
        L,
        V,
        LieAlgebraAbstractModule{QQFieldElem},
        LieAlgebraAbstractModuleElem{QQFieldElem},
      )
    end

    @testset "T^2 V of sl_4(CL(4))" begin
      L = special_linear_liealgebra(cyclotomic_field(4)[1], 4)
      V = tensor_power(standard_module(L), 2)
      liealgebra_module_conformance_test(
        L, V, LieAlgebraAbstractModule{nf_elem}, LieAlgebraAbstractModuleElem{nf_elem}
      )
    end
  end

  @testset "standard_module" begin
    L = special_orthogonal_liealgebra(QQ, 4)
    V = standard_module(L)
    @test dim(V) == 4

    x = L(rand(-10:10, dim(L)))
    v = V(rand(-10:10, dim(V)))
    @test x * v == V(transpose(matrix_repr(x) * transpose(Generic._matrix(v))))
  end

  @testset "exterior_power" begin
    L = special_orthogonal_liealgebra(QQ, 4)
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
    L = special_orthogonal_liealgebra(QQ, 4)
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
    L = special_orthogonal_liealgebra(QQ, 4)
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
