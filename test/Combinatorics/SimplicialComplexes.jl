@testset "SimplicialComplex" begin
    
    @testset "properties" begin
        
        sphere = SimplicialComplex([[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4]])
        
        sphere2 = SimplicialComplex(IncidenceMatrix([[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4]]))
        
        @test sphere isa SimplicialComplex
        @test sphere2 isa SimplicialComplex
        @test vertexindices(sphere) == [1, 2, 3, 4]
        @test nvertices(sphere) == 4
        @test facets(sphere) == Set{Int}.([[1, 3, 2], [1, 4, 2], [1, 3, 4], [3, 4, 2]])
        @test facets(sphere2) == facets(sphere)
        @test dim(sphere) == 2
        @test f_vector(sphere) == [4, 6, 4]
        @test h_vector(sphere) == [1, 1, 1, 1]
        @test betti_numbers(sphere) == [0, 0, 1]
        @test euler_characteristic(sphere) == 1
        @test minimal_nonfaces(sphere) == [Set{Int}([1, 2, 3, 4])]
        P, _ = PolynomialRing(ZZ, 4)
        R, _ = PolynomialRing(ZZ, ["a", "x", "i_7", "n"])
        @test stanley_reisner_ideal(sphere) == ideal([P([1], [[1, 1, 1, 1]])])
        @test stanley_reisner_ideal(R, sphere) == ideal([R([1], [[1, 1, 1, 1]])])
        @test isisomorphic(fundamental_group(sphere)[1], free_group())[1]
        @test string.(gens(fundamental_group(sphere)[1])) == string.(gens(free_group(3)))
        
    end
    
    @testset "standard examples" begin
        
        for (SC, fv, bn) in ((torus(), [7, 21, 14], [0, 2, 1]),
                            (klein_bottle(), [9, 27, 18], [0, 1, 0]),
                            (real_projective_plane(), [6, 15, 10], [0, 0, 0]),
                            (complex_projective_plane(), [9, 36, 84, 90, 36], [0, 0, 1, 0, 1]))
            
            @test SC isa SimplicialComplex
            @test f_vector(SC) == fv
            @test betti_numbers(SC) == bn
            
        end
        
    end
    
end
