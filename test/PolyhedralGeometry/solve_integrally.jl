@testset "solve_integrally" begin
  @testset "solve_mixed" begin
    A = ZZMatrix([1 1])
    b = zero_matrix(ZZ, 1, 1)
    b[1, 1] = 7
    C = ZZMatrix([1 0; 0 1])
    d = zero_matrix(ZZ, 2, 1)
    d[1, 1] = 2
    d[2, 1] = 3
    S0 = solve_mixed(A, b, C, d)
    S1 = solve_mixed(ZZMatrix, A, b, C, d)
    S2 = solve_mixed(SubObjectIterator{PointVector{ZZRingElem}}, A, b, C, d)
    @test S0 isa ZZMatrix
    @test S1 isa ZZMatrix
    @test S2 isa SubObjectIterator{PointVector{ZZRingElem}}
    @test nrows(S0) == 3
    @test nrows(S1) == 3
    @test length(S2) == 3

    S3 = solve_mixed(A, b, C)
    S4 = solve_mixed(ZZMatrix, A, b, C)
    S5 = solve_mixed(SubObjectIterator{PointVector{ZZRingElem}}, A, b, C)
    @test S3 isa ZZMatrix
    @test S4 isa ZZMatrix
    @test S5 isa SubObjectIterator{PointVector{ZZRingElem}}
    @test nrows(S3) == 8
    @test nrows(S4) == 8
    @test length(S5) == 8
  end

  @testset "solve_ineq" begin
    A = ZZMatrix([1 0; 0 1; -1 0; 0 -1])
    b = zero_matrix(ZZ, 4, 1)
    b[1, 1] = 1
    b[2, 1] = 1
    b[3, 1] = 0
    b[4, 1] = 0
    S0 = solve_ineq(A, b)
    S1 = solve_ineq(ZZMatrix, A, b)
    S2 = solve_ineq(SubObjectIterator{PointVector{ZZRingElem}}, A, b)
    @test S0 isa ZZMatrix
    @test S1 isa ZZMatrix
    @test S2 isa SubObjectIterator{PointVector{ZZRingElem}}
    @test nrows(S0) == 4
    @test nrows(S1) == 4
    @test length(S2) == 4
  end

  @testset "solve_non_negative" begin
    A = ZZMatrix([1 1])
    b = zero_matrix(ZZ, 1, 1)
    b[1, 1] = 3
    S0 = solve_non_negative(A, b)
    S1 = solve_non_negative(ZZMatrix, A, b)
    S2 = solve_non_negative(SubObjectIterator{PointVector{ZZRingElem}}, A, b)
    @test S0 isa ZZMatrix
    @test S1 isa ZZMatrix
    @test S2 isa SubObjectIterator{PointVector{ZZRingElem}}
    @test nrows(S0) == 4
    @test nrows(S1) == 4
    @test length(S2) == 4
  end
end
