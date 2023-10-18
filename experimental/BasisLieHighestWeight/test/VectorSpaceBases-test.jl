using Oscar
using Test

include("../src/VectorSpaceBases.jl")

@testset "Test VectorSpaceBases" begin
  a = sparse_row(ZZ, [1, 3, 6], [3, 2, 4])::SRow{ZZRingElem} # [3, 0, 2, 0, 0, 4]
  b = sparse_row(ZZ, [3], [2])::SRow{ZZRingElem} # [0, 0, 2, 0, 0, 0]
  c = sparse_row(ZZ, [1, 6], [4, 3])::SRow{ZZRingElem} # [4, 0, 0, 0, 0, 3]
  d = sparse_row(ZZ, [1, 3], [4, 3])::SRow{ZZRingElem} # [6, 0, 4, 0, 0, 0]
  sparse_vector_space_basis = SparseVectorSpaceBasis([a, b], [1, 3])

  @testset "reduce_col" begin
    a_reduced_b_3 = sparse_row(ZZ, [1, 6], [6, 8])::SRow{ZZRingElem} # [6, 0, 0, 0, 0, 8]
    a_reduced_c_1 = sparse_row(ZZ, [3, 6], [8, 7])::SRow{ZZRingElem} # [0, 0, 8, 0, 0, 7]
    @test isequal(reduce_col(a, b, 3), a_reduced_b_3)
    @test isequal(reduce_col(a, c, 1), a_reduced_c_1)
  end

  @testset "normalize" begin
    a_normalized = sparse_row(ZZ, [1, 3, 6], [3, 2, 4])::SRow{ZZRingElem} # [3, 0, 2, 0, 0, 4]
    b_normalized = sparse_row(ZZ, [3], [1])::SRow{ZZRingElem} # [0, 0, 1, 0, 0, 0]
    @test isequal(normalize(a), (a_normalized, 1))
    @test isequal(normalize(b), (b_normalized, 3))
  end

  @testset "add_and_reduce!" begin
    add_and_reduce!(sparse_vector_space_basis, c)
    c_reduced = sparse_row(ZZ, [6], [1])::SRow{ZZRingElem} # [0, 0, 0, 0, 0, 1] 
    @test isequal(sparse_vector_space_basis.basis_vectors, [a, b, c_reduced])
    @test isequal(sparse_vector_space_basis.pivot, [1, 3, 6])
    add_and_reduce!(sparse_vector_space_basis, d) # d = 2*a, therefore basis should not change
    @test isequal(sparse_vector_space_basis.basis_vectors, [a, b, c_reduced])
    @test isequal(sparse_vector_space_basis.pivot, [1, 3, 6])
  end
end
