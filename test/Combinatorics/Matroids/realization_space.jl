@testset "matroid realization spaces" begin


    M1 = fano_matroid()
    M2 = pappus_matroid()
    M3 = vamos_matroid()
    M4 = matroid_from_bases([['a','b'], ['a','c'], ['a','d'], ['b','c'], ['b','d']], ['a','b','c','d']))
    
    R1, mat_M1 = Oscar.realization_space_matrix(M1, [1,2,4], GF(2));
    R2, mat_M2 = Oscar.realization_space_matrix(M2, [2,3,8], QQ);
    R3, mat_M3 = Oscar.realization_space_matrix(M3, [2,4,5,7], ZZ);
    
    x1 = gens(R1)
    x2 = gens(R2)
    x3 = gens(R3)
      
    @testset "realization_space_matrix" begin
        
        @test coefficient_ring(R1) == GF(2)
        @test coefficient_ring(R2) == QQ
        @test coefficient_ring(R3) == ZZ
        
        @test length(x1) == 3
        @test length(x2) == 8
        @test length(x3) == 9
	
        @test mat_M1 ==  matrix(R1, [1 0 1 0 1 0 1; 0 1 1 0 0 1 x1[2]; 0 0 0 1 1 x1[1] x1[3]])
	@test mat_M2 ==  matrix(R2, [1 1 0 0 1 1 1 0 1; 1 0 1 1 x2[1] x2[3] x2[5] 0 x2[7]; 0 0 0 1 x2[2] x2[4] x2[6] 1 x2[8]])
	@test mat_M3 ==  matrix(R3, [1 1 1 0 0 1 0 1; 1 0 x3[1] 1 0 x3[4] 0 x3[7]; 1 0 x3[2] 0 1 x3[5] 0 x3[8]; 1 0 x3[3] 0 0 x3[6] 1 x3[9]])
	
    end
    
    B1 = find_good_basis_heuristically(M1)
    B2 = find_good_basis_heuristically(M2)
    B3 = find_good_basis_heuristically(M3)
    
    @testset "find_good_basis_heuristically" begin
        @test B1 == [2,6,7]
        @test B2 == [1,2,4]
        @test B3 == [1,2,3,5]
    end
    
    R, (x,y,z) = QQ["x", "y", "z"]
    f = 2*(x^2+y)^2*(x+y*z)^4
    
    @testset "poly_2_factors" begin
        @test poly_2_factors(f) == [x^2+y, x+y*z]
    end

    Sgens = [2*(x^2+y)^2*(x+y*z)^4, 3*x^2*(x+y*z)^5] 

    @testset "gens_2_factors" begin
        @test gens_2_factors(Sgens) == [x^2+y, x+y*z, x]
    end
    
    I = ideal(R, [x^2*(y+z), y^3*(y+z), z*(y+z)])
    Sgens = [x,y,z]
    
    @testset "stepwise_saturation" begin
        @test stepwise_saturation(I,Sgens) == ideal(R, [y+z])
    end
    
    
     
    
end

