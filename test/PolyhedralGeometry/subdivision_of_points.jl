@testset "SubdivisionOfPoints" begin
  C = cube(2)
  square_weights = [0, 0, 1, 2]
  square_max_cells = [[1, 2, 3], [2, 3, 4]]
  square_incidence = incidence_matrix(square_max_cells)

  square_by_weights = subdivision_of_points(C, square_weights)
  square_by_cells = subdivision_of_points(C, square_max_cells)
  square_by_incidence = subdivision_of_points(C, square_incidence)

  # input from vectors of vectors
  @test subdivision_of_points(vertices(C), square_weights) isa SubdivisionOfPoints
  @test subdivision_of_points([[0, 0], [1, 0], [0, 1], [1, 1]], square_weights) isa
    SubdivisionOfPoints

  @testset "alternative inputs" begin
    @test issetequal(maximal_cells(square_by_incidence),
                     maximal_cells(square_by_weights))
    @test issetequal(square_max_cells, maximal_cells(subdivision_of_points(C,min_weights(square_by_cells))))
  end

  moaepts = [4 0 0; 0 4 0; 0 0 4; 2 1 1; 1 2 1; 1 1 2]
  fulldim_moaepts = moaepts[:, 2:3]
  moaeimreg0 = incidence_matrix(1, 6)
  moaeimreg0[1, :] = [1, 1, 1, 0, 0, 0]
  moaeimnonreg0 = incidence_matrix([
    [4, 5, 6], [1, 4, 2], [2, 4, 5], [2, 3, 5], [3, 5, 6], [1, 3, 6], [1, 4, 6]
  ])

  MOAE = subdivision_of_points(moaepts, moaeimnonreg0)
  fulldim_MOAE = subdivision_of_points(fulldim_moaepts, moaeimnonreg0)
  SOP0 = subdivision_of_points(moaepts, moaeimreg0)

  CMOAE = secondary_cone(MOAE)

  SOP1 = subdivision_of_points(moaepts, [1, 2, 3, 4, 5, 6])
  C1 = secondary_cone(SOP1)

  @testset "core functionality" begin
    @test !is_regular(MOAE)
    @test is_regular(SOP0)
    @test number_of_maximal_cells(MOAE) == 7
    @test number_of_maximal_cells(SOP0) == 1
    @test number_of_maximal_cells(SOP1) == 1
    @test min_weights(SOP1) == [0, 0, 0, 1, 1, 1]
    @test dim(C1) == 6
    @test dim(CMOAE) == 4
    @test _check_im_perm_rows(moaeimnonreg0, maximal_cells(IncidenceMatrix, MOAE))
    @test number_of_points(MOAE) == 6
    @test length(points(MOAE)) == 6
    @test [0, 0, 4] in points(MOAE)
    @test gkz_vector(fulldim_MOAE) == [9, 9, 9, 7, 7, 7]
    @test_throws ArgumentError subdivision_of_points(matrix(QQ,[[0, 0], [1, 0], [1, 0], [1, 1]]), [1, 2, 2, 4])
    @test_throws ArgumentError subdivision_of_points([[0, 0], [1, 0], [1, 0], [1, 1]], square_weights)
  end
end
