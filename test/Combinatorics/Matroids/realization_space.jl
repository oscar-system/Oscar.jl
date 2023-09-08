@testset "Matroid Strata in the Grassmannian" begin

    M1 = uniform_matroid(2,5)
    M2 = fano_matroid()
    M3 = cycle_matroid(complete_graph(4))
    M4 = pappus_matroid(
    mb6 = matroid_from_bases([['a','b'], ['a','c'], ['a','d'], ['b','c'], ['b','d']], ['a','b','c','d']))
    
    mat_M1 = Oscar.realization_space_matrix(bases(mb1), [1,2]);
    mat_M2 = Oscar.bases_matrix_coordinates(bases(mb2), [1,2,4]);
    
    @testset "realization_space_matrix" begin
	
        @test bmc1 == [[1,1],[2,1],[1,2],[2,2],[1,3],[2,3]]
        @test bmc2 == [[1,1],[2,1],[1,2],[3,2],[2,3],[3,3],[1,4],[2,4],[3,4]]
            
    end


    
end

