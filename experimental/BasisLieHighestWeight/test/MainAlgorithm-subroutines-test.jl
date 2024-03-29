using ..Oscar
using ..Oscar: GAPWrap, IntegerUnion, _is_weighted
using Test

include("../src/MainAlgorithm.jl")

@testset "Test operators_by_index" begin
  L = lie_algebra(:A, 2)
  chevalley_basis = chevalley_basis_gap(L)
  birational_sequence = [1, 2, 1]

  result = operators_by_index(L, chevalley_basis, birational_sequence)

  expected = [chevalley_basis[1][1], chevalley_basis[1][2], chevalley_basis[1][1]]
  @test isequal(result, expected)
end

@testset "Test operators_by_simple_roots" begin
  L = lie_algebra(:A, 2)
  chevalley_basis = chevalley_basis_gap(L)
  birational_sequence = [[1, 0], [0, 1], [1, 0]]

  result = operators_by_simple_roots(L, chevalley_basis, birational_sequence)

  expected = [chevalley_basis[1][1], chevalley_basis[1][2], chevalley_basis[1][1]]
  @test isequal(result, expected)
end

@testset "Test operators_lusztig" begin
  L = lie_algebra(:A, 2)
  chevalley_basis = chevalley_basis_gap(L)
  birational_sequence = [1, 2, 1]

  result = operators_lusztig(L, chevalley_basis, birational_sequence)

  expected = [chevalley_basis[1][1], chevalley_basis[1][3], chevalley_basis[1][2]]
  @test isequal(result, expected)
end

@testset "Test operators_lusztig_indices" begin
  L = lie_algebra(:A, 2)
  birational_sequence = [1, 2, 1]

  result = operators_lusztig_indices(L, birational_sequence)

  expected = [1, 3, 2]
  @test isequal(result, expected)
end

@testset "Test compute_sub_weights" begin
  highest_weight = [ZZ(2), ZZ(2)]

  result = compute_sub_weights(highest_weight)

  expected = Any[ZZRingElem[1, 0], ZZRingElem[0, 1], 
                ZZRingElem[1, 1], ZZRingElem[2, 0], 
                ZZRingElem[0, 2], ZZRingElem[2, 1], 
                ZZRingElem[1, 2]]
  @test isequal(result, expected)
end