@testset "LieAlgebras.LinearLieAlgebra" begin
  @testset "constructors for R=$R, n=$n" for R in [QQ, cyclotomic_field(4)[1]], n in 1:5
    L = general_linear_lie_algebra(R, n)
    @test dim(L) == n^2

    L = special_linear_lie_algebra(R, n)
    @test dim(L) == n^2 - 1

    L = special_orthogonal_lie_algebra(R, n)
    @test dim(L) == div(n^2 - n, 2)
  end

  @testset "conformance tests" begin
    @testset "gl_4(QQ)" begin
      L = general_linear_lie_algebra(QQ, 4)
      lie_algebra_conformance_test(
        L, LinearLieAlgebra{QQFieldElem}, LinearLieAlgebraElem{QQFieldElem}
      )
    end

    @testset "sl_4(QQ)" begin
      L = special_linear_lie_algebra(QQ, 4)
      lie_algebra_conformance_test(
        L, LinearLieAlgebra{QQFieldElem}, LinearLieAlgebraElem{QQFieldElem}
      )
    end

    @testset "so_4(QQ)" begin
      L = special_orthogonal_lie_algebra(QQ, 4)
      lie_algebra_conformance_test(
        L, LinearLieAlgebra{QQFieldElem}, LinearLieAlgebraElem{QQFieldElem}
      )
    end

    @testset "gl_4(CF(4))" begin
      L = general_linear_lie_algebra(cyclotomic_field(4)[1], 4)
      lie_algebra_conformance_test(
        L, LinearLieAlgebra{nf_elem}, LinearLieAlgebraElem{nf_elem}
      )
    end

    @testset "sl_4(CF(4))" begin
      L = special_linear_lie_algebra(cyclotomic_field(4)[1], 4)
      lie_algebra_conformance_test(
        L, LinearLieAlgebra{nf_elem}, LinearLieAlgebraElem{nf_elem}
      )
    end

    @testset "so_4(CF(4))" begin
      L = special_orthogonal_lie_algebra(cyclotomic_field(4)[1], 4)
      lie_algebra_conformance_test(
        L, LinearLieAlgebra{nf_elem}, LinearLieAlgebraElem{nf_elem}
      )
    end
  end

  @testset "so_n correctness regression" begin
    function lie_algebra_struct_const(L::LieAlgebra{C}) where {C<:RingElement}
      R = coefficient_ring(L)
      dimL = dim(L)
      struct_const_L = Matrix{Vector{Tuple{elem_type(R),Int}}}(undef, dimL, dimL)
      for (i, xi) in enumerate(basis(L)), (j, xj) in enumerate(basis(L))
        struct_const_L[i, j] = [
          (c, k) for (k, c) in enumerate(coefficients(xi * xj)) if !iszero(c)
        ]
      end
      return struct_const_L
    end

    L = special_orthogonal_lie_algebra(QQ, 3)
    struct_const_L = Matrix{Vector{Tuple{QQFieldElem,Int64}}}(undef, 3, 3)
    struct_const_L[1, :] = Vector{Tuple{QQFieldElem,Int64}}[[], [(QQ(-1), 3)], [(QQ(1), 2)]]
    struct_const_L[2, :] = Vector{Tuple{QQFieldElem,Int64}}[[(QQ(1), 3)], [], [(QQ(-1), 1)]]
    struct_const_L[3, :] = Vector{Tuple{QQFieldElem,Int64}}[[(QQ(-1), 2)], [(QQ(1), 1)], []]
    @test lie_algebra_struct_const(L) == struct_const_L

    L = special_orthogonal_lie_algebra(QQ, 4)
    struct_const_L = Matrix{Vector{Tuple{QQFieldElem,Int64}}}(undef, 6, 6)
    struct_const_L[1, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [], [(QQ(-1), 4)], [(QQ(-1), 5)], [(QQ(1), 2)], [(QQ(1), 3)], []
    ]
    struct_const_L[2, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(1), 4)], [], [(QQ(-1), 6)], [(QQ(-1), 1)], [], [(QQ(1), 3)]
    ]
    struct_const_L[3, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(1), 5)], [(QQ(1), 6)], [], [], [(QQ(-1), 1)], [(QQ(-1), 2)]
    ]
    struct_const_L[4, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(-1), 2)], [(QQ(1), 1)], [], [], [(QQ(-1), 6)], [(QQ(1), 5)]
    ]
    struct_const_L[5, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(-1), 3)], [], [(QQ(1), 1)], [(QQ(1), 6)], [], [(QQ(-1), 4)]
    ]
    struct_const_L[6, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [], [(QQ(-1), 3)], [(QQ(1), 2)], [(QQ(-1), 5)], [(QQ(1), 4)], []
    ]
    @test lie_algebra_struct_const(L) == struct_const_L

    L = special_orthogonal_lie_algebra(QQ, 5)
    struct_const_L = Matrix{Vector{Tuple{QQFieldElem,Int64}}}(undef, 10, 5)
    struct_const_L[1, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(-1), 2)], [(QQ(1), 1)], [], [], []
    ]
    struct_const_L[2, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(-1), 3)], [], [(QQ(1), 1)], [], []
    ]
    struct_const_L[3, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(-1), 4)], [], [], [(QQ(1), 1)], []
    ]
    struct_const_L[4, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(-1), 5)], [], [], [], [(QQ(1), 1)]
    ]
    struct_const_L[5, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [], [(QQ(-1), 3)], [(QQ(1), 2)], [], []
    ]
    struct_const_L[6, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [], [(QQ(-1), 4)], [], [(QQ(1), 2)], []
    ]
    struct_const_L[7, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [], [(QQ(-1), 5)], [], [], [(QQ(1), 2)]
    ]
    struct_const_L[8, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [], [], [(QQ(-1), 4)], [(QQ(1), 3)], []
    ]
    struct_const_L[9, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [], [], [(QQ(-1), 5)], [], [(QQ(1), 3)]
    ]
    struct_const_L[10, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [], [], [], [(QQ(-1), 5)], [(QQ(1), 4)]
    ]
    @test lie_algebra_struct_const(L) == struct_const_L
  end
end
