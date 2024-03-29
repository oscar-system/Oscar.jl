using ..Oscar
using ..Oscar: GAPWrap, IntegerUnion, _is_weighted
using Test

include("../src/VectorDemazure.jl")

@testset "Test demazure_vw" begin
    L = lie_algebra(:A, 2)
    highest_weight = [ZZ(1), ZZ(0)]
    chevalley_basis = chevalley_basis_gap(L)
    operators = chevalley_basis[1]
    matrices_of_operators = matrices_of_operators_gap(L, highest_weight, operators)
    reduced_expression = [1, 2, 1]

    result_vw, result_extremal_weight = demazure_vw(L, reduced_expression, highest_weight, matrices_of_operators)

    expected_vw = sparse_row(ZZ, Int64[], ZZRingElem[])
    expected_extremal_weight = ZZRingElem[0, -1]

    @test isequal(result_vw, expected_vw)
    @test isequal(result_extremal_weight, expected_extremal_weight)
end