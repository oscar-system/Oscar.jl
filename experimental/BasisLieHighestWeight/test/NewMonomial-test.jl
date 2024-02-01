@testset "Test NewMonomial" begin
  weight = BasisLieHighestWeight.weight
  calc_vec = BasisLieHighestWeight.calc_vec

  ZZx, _ = polynomial_ring(ZZ, 2)
  x = gens(ZZx)
  mon1 = ZZx(1)
  mon2 = x[1]^2 * x[2]
  weights = [[ZZ(1), ZZ(1)], [ZZ(2), ZZ(1)]]
  A = sparse_matrix(ZZ, 2, 2) # [0, 2; 1, 1]
  setindex!(A, sparse_row(ZZ, [2], [ZZ(2)]), 1)
  setindex!(A, sparse_row(ZZ, [1, 2], [ZZ(1), ZZ(1)]), 2)
  B = sparse_matrix(ZZ, 2, 2) # [1, 0; 0, 2]
  setindex!(B, sparse_row(ZZ, [1], [ZZ(1)]), 1)
  setindex!(B, sparse_row(ZZ, [2], [ZZ(2)]), 2)
  matrices_of_operators = [A, B]
  v0 = sparse_row(ZZ, [1], [1]) # [1, 0]

  mon2_vec = sparse_row(ZZ, [1, 2], [2, 2])

  @testset "weight" begin
    @test isequal(weight(mon1, weights), [ZZ(0), ZZ(0)])
    @test isequal(weight(mon2, weights), [ZZ(4), ZZ(3)])
  end

  @testset "calc_vec" begin
    @test isequal(calc_vec(v0, mon1, matrices_of_operators), v0)
    @test isequal(calc_vec(v0, mon2, matrices_of_operators), mon2_vec)
  end
end
