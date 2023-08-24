@testset "SubdivisionOfPoints" begin

   C = cube(2)
   square_weights = [0,0,1,2]
   square_max_cells = [[1,2,3],[2,3,4]]
   square_incidence = IncidenceMatrix(square_max_cells)

   square_by_weights = subdivision_of_points(C,square_weights)
   square_by_cells = subdivision_of_points(C,square_max_cells)
   square_by_incidence = subdivision_of_points(C,square_incidence)

   @testset "alternative inputs" begin
      @test collect(maximal_cells(square_by_incidence)) == collect(maximal_cells(square_by_weights))
      @test min_weights(square_by_cells) == min_weights(square_by_weights)
   end


   moaepts = [4 0 0; 0 4 0; 0 0 4; 2 1 1; 1 2 1; 1 1 2]
   fulldim_moaepts = moaepts[:,2:3]
   moaeimreg0 = IncidenceMatrix(1,6)
   moaeimreg0[1,:] = [1,1,1,0,0,0]
   moaeimnonreg0 = IncidenceMatrix([[4,5,6],[1,4,2],[2,4,5],[2,3,5],[3,5,6],[1,3,6],[1,4,6]])

   MOAE = subdivision_of_points(moaepts, moaeimnonreg0)
   fulldim_MOAE = subdivision_of_points(fulldim_moaepts,moaeimnonreg0)
   SOP0 = subdivision_of_points(moaepts, moaeimreg0)

   CMOAE = secondary_cone(MOAE)

   SOP1 = subdivision_of_points(moaepts, [1,2,3,4,5,6])
   C1 = secondary_cone(SOP1)

    @testset "core functionality" begin
        @test !is_regular(MOAE)
        @test is_regular(SOP0)
        @test n_maximal_cells(MOAE) == 7
        @test n_maximal_cells(SOP0) == 1
        @test n_maximal_cells(SOP1) == 1
        @test min_weights(SOP1) == [0,0,0,1,1,1]
        @test dim(C1) == 6
        @test dim(CMOAE) == 4
        @test moaeimnonreg0 == maximal_cells(IncidenceMatrix, MOAE)
        @test npoints(MOAE) == 6
        @test length(points(MOAE)) == 6
        @test collect(points(MOAE))[3] == [0,0,4]
        @test gkz_vector(fulldim_MOAE) == [9,9,9,7,7,7]
    end

end
