using ..Oscar
using ..Oscar: GAPWrap, IntegerUnion, _is_weighted
using Test

include("../src/LieAlgebras.jl")
include("../src/WeylPolytope.jl")

@testset "Test orbit_weylgroup" begin
    L = lie_algebra(:A, 2)
    weight = [ZZ(1), ZZ(1)]
    result = orbit_weylgroup(L, weight)

    expected = Vector{ZZRingElem}[[1, 1], [-1, -1], [-1, 2], [-2, 1], [2, -1], [1, -2]]
    @test isequal(result, expected)
end

@testset "Test get_dim_weightspace" begin
    L = lie_algebra(:A, 2)
    weight = [ZZ(1), ZZ(1)]
    result = get_dim_weightspace(L, weight)

    expected = Dict{Vector{ZZRingElem}, Int64}([2, -1] => 1, [1, 1] => 2, 
                                               [0, 0] => 1, [0, 3] => 1, 
                                               [3, 0] => 1, [2, 2] => 1, 
                                               [-1, 2] => 1)
    @test isequal(result, expected)
end

@testset "Test convert_lattice_points_to_monomials" begin
    ZZx, _ = polynomial_ring(ZZ, 2)
    lattice_points_weightspace = [[ZZ(1), ZZ(2)], [ZZ(5), ZZ(0)]]
    result = convert_lattice_points_to_monomials(ZZx, lattice_points_weightspace)

    expected = [ZZx[1]*ZZx[2]^2, ZZx[1]^5]
    @test isequal(result, expected)
end

@testset "Test get_lattice_points_of_weightspace" begin
    root_weights = Vector{QQFieldElem}[[1, 0], [0, 1], [1, 0]]
    weight = QQFieldElem[1, 0]
    zero_coordinates = Int64[]
    result = get_lattice_points_of_weightspace(root_weights, weight, zero_coordinates)

    expected = Vector{ZZRingElem}[[0, 0, 1], [1, 0, 0]]
    @test isequal(result, expected)
end

@testset "Test compute_zero_coordinates" begin
    highest_weight = [ZZ(1), ZZ(0)]

    # bir_sequence
    L = lie_algebra(:A, 2)
    chevalley_basis = chevalley_basis_gap(L)
    operators = chevalley_basis[1]
    weights_w = [weight(L, op) for op in operators] # weights of the operators
    weights_alpha = [w_to_alpha(L, weight_w) for weight_w in weights_w] # other root system
    asVec(v) = GAP.gap_to_julia(GAPWrap.ExtRepOfObj(v)) # TODO
    bir_sequence = BirationalSequence(
      operators, [asVec(v) for v in operators], weights_w, weights_alpha
    )

    result = compute_zero_coordinates(bir_sequence, highest_weight)

    expected = [2]
    @test isequal(result, expected)
end