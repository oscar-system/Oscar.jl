@testset "SubdivisionOfPoints" begin

   moaepts = [4 0 0; 0 4 0; 0 0 4; 2 1 1; 1 2 1; 1 1 2]
   moaeimreg0 = IncidenceMatrix([[1,2,3]])
   moaeimnonreg0 = IncidenceMatrix([[4,5,6],[1,4,2],[2,4,5],[2,3,5],[3,5,6],[1,3,6],[1,4,6]])

   MOAE = SubdivisionOfPoints(moaepts, moaeimnonreg0)
   SOP0 = SubdivisionOfPoints(moaepts, moaeimreg0)

   CMOAE = secondary_cone(MOAE)

   SOP1 = SubdivisionOfPoints(moaepts, [1,2,3,4,5,6])
   C1 = secondary_cone(SOP1)

    @testset "core functionality" begin
        @test !isregular(MOAE)
        @test isregular(SOP0)
        @test nmaximal_cells(MOAE) == 7
        @test nmaximal_cells(SOP0) == 1
        @test nmaximal_cells(SOP1) == 1
        @test min_weights(SOP1) == [0,0,0,1,1,1]
        @test dim(C1) == 6
        @test dim(CMOAE) == 4
        @test moaeimnonreg0 == maximal_cells_as_incidence_matrix(MOAE)
    end

end
