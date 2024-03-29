using ..Oscar
using ..Oscar: GAPWrap, IntegerUnion, _is_weighted
using Test

include("../src/LieAlgebras.jl")

@testset "Test matrices_of_operators_gap" begin
    highest_weight = [ZZ(1), ZZ(0)]
    L = lie_algebra(:A, 2)
    chevalley_basis = chevalley_basis_gap(L)
    operators = chevalley_basis[1]

    result = matrices_of_operators_gap(L, highest_weight, operators)

    mat1 = sparse_matrix(ZZ, 0, 3)
    push!(mat1, sparse_row(ZZ, Int64[2], ZZRingElem[ZZ(1)]))
    push!(mat1, sparse_row(ZZ, Int64[], ZZRingElem[]))
    push!(mat1, sparse_row(ZZ, Int64[], ZZRingElem[]))

    mat2 = sparse_matrix(ZZ, 0, 3)
    push!(mat2, sparse_row(ZZ, Int64[], ZZRingElem[]))
    push!(mat2, sparse_row(ZZ, Int64[3], ZZRingElem[ZZ(-1)]))
    push!(mat2, sparse_row(ZZ, Int64[], ZZRingElem[]))

    mat3 = sparse_matrix(ZZ, 0, 3)
    push!(mat3, sparse_row(ZZ, Int64[3], ZZRingElem[ZZ(1)]))
    push!(mat3, sparse_row(ZZ, Int64[], ZZRingElem[]))
    push!(mat3, sparse_row(ZZ, Int64[], ZZRingElem[]))

    expected = [mat1, mat2, mat3]

    @test isequal(result, expected)
end

@testset "Test weight" begin
    L = lie_algebra(:A, 2)
    chevalley_basis = chevalley_basis_gap(L)
    operators = chevalley_basis[1]
    operator = operators[1]

    result = weight(L, operator)

    expected = ZZRingElem[2, -1]
    @test isequal(result, expected)
end