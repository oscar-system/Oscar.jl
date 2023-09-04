using Oscar
using Test

include("../src/NewMonomial.jl")

@testset "Test NewMonomial" begin
    ZZx, _ = PolynomialRing(ZZ, 2)
    x = gens(ZZx)
    mon1 = ZZx(1)
    mon2 = x[1]^2 * x[2]
    weights = [[ZZ(1), ZZ(1)], [ZZ(2), ZZ(1)]]
    A = sparse_matrix(ZZ, 2, 2) # [0, 1; 2, 1]
    setindex!(A, sparse_row(ZZ, [2], [ZZ(1)]), 1)
    setindex!(A, sparse_row(ZZ, [1, 2], [ZZ(2), ZZ(1)]), 2)
    B = sparse_matrix(ZZ, 2, 2) # [1, 0; 2, 0]
    setindex!(B, sparse_row(ZZ, [1], [ZZ(1)]), 1)
    setindex!(B, sparse_row(ZZ, [2], [ZZ(2)]), 2)
    matrices_of_operators = [A, B]
    v0 =  sparse_row(ZZ, [1], [1])::SRow{ZZRingElem} # [1, 0]
    calc_monomials = Dict{ZZMPolyRingElem, Tuple{SRow{ZZRingElem}, Vector{Int}}}(ZZx(1) => (v0, [0, 0])) 

    mon2_vec = sparse_row(ZZ, [1, 2], [2, 2])::SRow{ZZRingElem}

    @testset "calc_weight" begin
        @test isequal(calc_weight(mon1, weights), [ZZ(0), ZZ(0)])
        @test isequal(calc_weight(mon2, weights), [ZZ(4), ZZ(3)])
    end

    @testset "calc_vec" begin
        @test isequal(calc_vec(v0, mon1, matrices_of_operators), v0)
        @test isequal(calc_vec(v0, mon2, matrices_of_operators), mon2_vec)
    end

    @testset "highest_calc_sub_monomial" begin
        @test isequal(highest_calc_sub_monomial(x, mon2, calc_monomials), ZZx(1))
    end
end
