@testset ExtendedTestSet "All LinearLieAlgebra.jl tests" begin

    @testset "constructors for R=$R, n=$n" for R in [QQ, cyclotomic_field(4)[1]], n in 1:5
        L = general_linear_liealgebra(R, n)
        @test dim(L) == n^2

        L = special_linear_liealgebra(R, n)
        @test dim(L) == n^2 - 1

        L = special_orthogonal_liealgebra(R, n)
        @test dim(L) == div(n^2 - n, 2)
    end

    @testset "conformance tests for $desc" for (desc, L, parentT, elemT) in [
        ("sl_4(QQ)", special_linear_liealgebra(QQ, 4), LinearLieAlgebra{fmpq}, LinearLieAlgebraElem{fmpq}),
        ("so_4(QQ)", special_orthogonal_liealgebra(QQ, 4), LinearLieAlgebra{fmpq}, LinearLieAlgebraElem{fmpq}),
        (
            "sl_4(CF(4))",
            special_linear_liealgebra(cyclotomic_field(4)[1], 4),
            LinearLieAlgebra{nf_elem},
            LinearLieAlgebraElem{nf_elem},
        ),
        (
            "so_4(CF(4))",
            special_orthogonal_liealgebra(cyclotomic_field(4)[1], 4),
            LinearLieAlgebra{nf_elem},
            LinearLieAlgebraElem{nf_elem},
        ),
    ]
        liealgebra_conformance_test(L, parentT, elemT)
    end

    @testset "so_n correctness regression" begin
        function liealgebra_struct_const(L::LieAlgebra{C}) where {C <: RingElement}
            R = base_ring(L)
            dimL = dim(L)
            struct_const_L = Matrix{Vector{Tuple{elem_type(R), Int}}}(undef, dimL, dimL)
            for (i, xi) in enumerate(gens(L)), (j, xj) in enumerate(gens(L))
                struct_const_L[i, j] = [(c, k) for (k, c) in enumerate(PD._matrix(bracket(xi, xj))) if !iszero(c)]
            end
            return struct_const_L
        end

        function liealgebra_module_struct_const(L::LieAlgebra{C}, V::LieAlgebraModule{C}) where {C <: RingElement}
            R = base_ring(L)
            dimL = dim(L)
            dimV = dim(V)
            struct_const_V = Matrix{Vector{Tuple{elem_type(R), Int}}}(undef, dimL, dimV)
            for (i, xi) in enumerate(gens(L)), (j, vj) in enumerate(gens(V))
                struct_const_V[i, j] = [(c, k) for (k, c) in enumerate(PD._matrix(xi * vj)) if !iszero(c)]
            end
            return struct_const_V
        end

        L = special_orthogonal_liealgebra(QQ, 3)
        @test repr("text/plain", liealgebra_struct_const(L)) ==
              "3×3 Matrix{Vector{Tuple{fmpq, Int64}}}:\n []         [(-1, 3)]  [(1, 2)]\n [(1, 3)]   []         [(-1, 1)]\n [(-1, 2)]  [(1, 1)]   []"
        @test repr("text/plain", liealgebra_module_struct_const(L, symmetric_power(standard_module(L), 1))) ==
              "3×3 Matrix{Vector{Tuple{fmpq, Int64}}}:\n [(-1, 2)]  [(1, 1)]   []\n [(-1, 3)]  []         [(1, 1)]\n []         [(-1, 3)]  [(1, 2)]"
        @test repr("text/plain", liealgebra_module_struct_const(L, symmetric_power(standard_module(L), 2))) ==
              "3×6 Matrix{Vector{Tuple{fmpq, Int64}}}:\n [(-2, 2)]  [(1, 1), (-1, 4)]  [(-1, 5)]          [(2, 2)]   [(1, 3)]           []\n [(-2, 3)]  [(-1, 5)]          [(1, 1), (-1, 6)]  []         [(1, 2)]           [(2, 3)]\n []         [(-1, 3)]          [(1, 2)]           [(-2, 5)]  [(1, 4), (-1, 6)]  [(2, 5)]"
        @test repr("text/plain", liealgebra_module_struct_const(L, symmetric_power(standard_module(L), 3))) ==
              "3×10 Matrix{Vector{Tuple{fmpq, Int64}}}:\n [(-3, 2)]  [(1, 1), (-2, 4)]  [(-2, 5)]          [(2, 2), (-1, 7)]  [(1, 3), (-1, 8)]  [(-1, 9)]           [(3, 4)]   [(2, 5)]           [(1, 6)]            []\n [(-3, 3)]  [(-2, 5)]          [(1, 1), (-2, 6)]  [(-1, 8)]          [(1, 2), (-1, 9)]  [(2, 3), (-1, 10)]  []         [(1, 4)]           [(2, 5)]            [(3, 6)]\n []         [(-1, 3)]          [(1, 2)]           [(-2, 5)]          [(1, 4), (-1, 6)]  [(2, 5)]            [(-3, 8)]  [(1, 7), (-2, 9)]  [(2, 8), (-1, 10)]  [(3, 9)]"
        @test repr("text/plain", liealgebra_module_struct_const(L, exterior_power(standard_module(L), 1))) ==
              "3×3 Matrix{Vector{Tuple{fmpq, Int64}}}:\n [(-1, 2)]  [(1, 1)]   []\n [(-1, 3)]  []         [(1, 1)]\n []         [(-1, 3)]  [(1, 2)]"
        @test repr("text/plain", liealgebra_module_struct_const(L, exterior_power(standard_module(L), 2))) ==
              "3×3 Matrix{Vector{Tuple{fmpq, Int64}}}:\n []         [(-1, 3)]  [(1, 2)]\n [(1, 3)]   []         [(-1, 1)]\n [(-1, 2)]  [(1, 1)]   []"

        L = special_orthogonal_liealgebra(QQ, 4)
        @test repr("text/plain", liealgebra_struct_const(L)) ==
              "6×6 Matrix{Vector{Tuple{fmpq, Int64}}}:\n []         [(-1, 4)]  [(-1, 5)]  [(1, 2)]   [(1, 3)]   []\n [(1, 4)]   []         [(-1, 6)]  [(-1, 1)]  []         [(1, 3)]\n [(1, 5)]   [(1, 6)]   []         []         [(-1, 1)]  [(-1, 2)]\n [(-1, 2)]  [(1, 1)]   []         []         [(-1, 6)]  [(1, 5)]\n [(-1, 3)]  []         [(1, 1)]   [(1, 6)]   []         [(-1, 4)]\n []         [(-1, 3)]  [(1, 2)]   [(-1, 5)]  [(1, 4)]   []"
        @test repr("text/plain", liealgebra_module_struct_const(L, symmetric_power(standard_module(L), 1))) ==
              "6×4 Matrix{Vector{Tuple{fmpq, Int64}}}:\n [(-1, 2)]  [(1, 1)]   []         []\n [(-1, 3)]  []         [(1, 1)]   []\n [(-1, 4)]  []         []         [(1, 1)]\n []         [(-1, 3)]  [(1, 2)]   []\n []         [(-1, 4)]  []         [(1, 2)]\n []         []         [(-1, 4)]  [(1, 3)]"
        @test repr("text/plain", liealgebra_module_struct_const(L, exterior_power(standard_module(L), 1))) ==
              "6×4 Matrix{Vector{Tuple{fmpq, Int64}}}:\n [(-1, 2)]  [(1, 1)]   []         []\n [(-1, 3)]  []         [(1, 1)]   []\n [(-1, 4)]  []         []         [(1, 1)]\n []         [(-1, 3)]  [(1, 2)]   []\n []         [(-1, 4)]  []         [(1, 2)]\n []         []         [(-1, 4)]  [(1, 3)]"
        @test repr("text/plain", liealgebra_module_struct_const(L, exterior_power(standard_module(L), 2))) ==
              "6×6 Matrix{Vector{Tuple{fmpq, Int64}}}:\n []         [(-1, 4)]  [(-1, 5)]  [(1, 2)]   [(1, 3)]   []\n [(1, 4)]   []         [(-1, 6)]  [(-1, 1)]  []         [(1, 3)]\n [(1, 5)]   [(1, 6)]   []         []         [(-1, 1)]  [(-1, 2)]\n [(-1, 2)]  [(1, 1)]   []         []         [(-1, 6)]  [(1, 5)]\n [(-1, 3)]  []         [(1, 1)]   [(1, 6)]   []         [(-1, 4)]\n []         [(-1, 3)]  [(1, 2)]   [(-1, 5)]  [(1, 4)]   []"
        @test repr("text/plain", liealgebra_module_struct_const(L, exterior_power(standard_module(L), 3))) ==
              "6×4 Matrix{Vector{Tuple{fmpq, Int64}}}:\n []         []         [(-1, 4)]  [(1, 3)]\n []         [(1, 4)]   []         [(-1, 2)]\n [(-1, 4)]  []         []         [(1, 1)]\n []         [(-1, 3)]  [(1, 2)]   []\n [(1, 3)]   []         [(-1, 1)]  []\n [(-1, 2)]  [(1, 1)]   []         []"

        L = special_orthogonal_liealgebra(QQ, 5)
        @test repr("text/plain", liealgebra_module_struct_const(L, exterior_power(standard_module(L), 1))) ==
              "10×5 Matrix{Vector{Tuple{fmpq, Int64}}}:\n [(-1, 2)]  [(1, 1)]   []         []         []\n [(-1, 3)]  []         [(1, 1)]   []         []\n [(-1, 4)]  []         []         [(1, 1)]   []\n [(-1, 5)]  []         []         []         [(1, 1)]\n []         [(-1, 3)]  [(1, 2)]   []         []\n []         [(-1, 4)]  []         [(1, 2)]   []\n []         [(-1, 5)]  []         []         [(1, 2)]\n []         []         [(-1, 4)]  [(1, 3)]   []\n []         []         [(-1, 5)]  []         [(1, 3)]\n []         []         []         [(-1, 5)]  [(1, 4)]"
    end
end
