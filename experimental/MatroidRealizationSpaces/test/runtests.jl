@testset "matroid realization spaces" begin

    M1 = fano_matroid()
    M2 = pappus_matroid()
    M3 = vamos_matroid()
    M4 = matroid_from_bases([['a','b'], ['a','c'], ['a','d'], ['b','c'], ['b','d']], ['a','b','c','d'])
    
    R1, mat_M1 = Oscar.realization_space_matrix(M1, [1,2,4], GF(2));
    R2, mat_M2 = Oscar.realization_space_matrix(M2, [2,3,8], QQ);
    R3, mat_M3 = Oscar.realization_space_matrix(M3, [2,4,5,7], ZZ);
    
    x1 = gens(R1)
    x2 = gens(R2)
    x3 = gens(R3)

    @testset "is_realizable" begin
        @test is_realizable(M1) == true
        @test is_realizable(M2) == true
        @test is_realizable(M3) == false
        @test is_realizable(M4) == true

        @test is_realizable(M1, char = 2) == true
        @test is_realizable(M1, char = 3) == false
        @test is_realizable(M1, q = 2) == true
        @test is_realizable(M1, q = 3) == false

        @test is_realizable(M2, char = 2) == true
        @test is_realizable(M2, char = 3) == true
        @test is_realizable(M2, q = 2) == false
        @test is_realizable(M2, q = 7) == true
    end
      
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
    
    B1 = Oscar.find_good_basis_heuristically(M1)
    B2 = Oscar.find_good_basis_heuristically(M2)
    B3 = Oscar.find_good_basis_heuristically(M3)
    
    @testset "find_good_basis_heuristically" begin
        @test B1 == [2,6,7]
        @test B2 == [1,2,4]
        @test B3 == [1,2,3,5]
    end
    
    R, (x,y,z) = QQ[:x, :y, :z]
    f = 2*(x^2+y)^2*(x+y*z)^4
    
    @testset "poly_2_prime_divisors" begin
      @test all(in([x^2+y, x+y*z]), Oscar.poly_2_prime_divisors(f))
    end

    Sgens = [2*(x^2+y)^2*(x+y*z)^4, 3*x^2*(x+y*z)^5] 

    @testset "gens_2_prime_divisors" begin
      @test all(in(Oscar.gens_2_prime_divisors(Sgens)), [x^2+y, x+y*z, x])
    end
    
    I = ideal(R, [x^2*(y+z), y^3*(y+z), z*(y+z)])
    Sgens = [x,y,z]
    
    @testset "stepwise_saturation" begin
        @test Oscar.stepwise_saturation(I,Sgens) == ideal(R, [y+z])
    end

    Igens = [y^3-x^2+1, x^2*y-z^3+6]
    Sgens = [x,y,z]

    @testset "find_solution_v" begin
      @test Oscar.find_solution_v(y, Igens, Sgens, R) == (true, (z^3-6)//x^2)
    end

    phi = Oscar.sub_map(y, (z^3-6)//x^2, R, [x,y,z])
    
    @testset "sub_map" begin
        @test phi.([x,y,z]) == [x,(z^3-6)//x^2,z]
    end

    t = (z^3-6)//x^2
    f = y^3-x^2+1

    @testset "sub_v" begin
        @test Oscar.sub_v(y, t, f, R, [x,y,z]) == (z^3-6)^3 + x^6*(-x^2 + 1)
    end
    
    f = x*y^2*z^3*(x^2+y-x*z^3)
    Sgens = [x,y,z]
    
    @testset "clean" begin
        @test Oscar.clean(f, R, Sgens) == x^2+y-x*z^3
    end
    
    Igens = [x+y, x^2]
    
    @testset "ideal_vars" begin
        @test Oscar.ideal_vars(Igens) == [x,y]
    end

    newSgens = Oscar.n_new_Sgens(y,t,[y^3-x^2+1, y+x^3 ],R,[x,y,z])

    @testset "n_new_Sgens" begin
        @test newSgens == [x^8 - x^6 - z^9 + 18*z^6 - 108*z^3 + 216, x^5 + z^3 - 6]
    end

    @testset "n_new_Igens" begin
        @test Oscar.n_new_Igens(y,t,[y^3-x-2, y+z^3], newSgens ,R,[x,y,z]) == [x^7 + 2*x^6 - z^9 + 18*z^6 - 108*z^3 + 216, x^2*z^3 + z^3 - 6]
    end

    X = matrix(fraction_field(R), [x//(z-1) y//(x+1) z; -y+1 x//z 2*x*y ]) 

    @testset "matrix_clear_den_in_col" begin
        @test Oscar.matrix_clear_den_in_col(X, 2) == matrix(fraction_field(R), [x//(z-1) y*z z; -y+1 x^2+x 2*x*y]) 
    end

    @testset "matrix_clear_den" begin
        @test Oscar.matrix_clear_den(X) == matrix(fraction_field(R), [x y*z z; (-y+1)*(z-1) x^2+x 2*x*y]) 
    end


end

@testset "realization spaces as schemes" begin
  X = realization_space(pappus_matroid())
  @test X isa AbsAffineScheme
  @test !isdefined(X, :underlying_scheme)
  R = OO(X)
  @test R isa Oscar.MPolyQuoLocRing
  f = sum(gens(R))
  U = PrincipalOpenSubset(X, f)
  @test U isa AbsAffineScheme
end

@testset "saturation of realization space" begin
    s = "0******0******0**********0********0*******0*********************0*****************0*************0***********0***************************0*******************0*****************************0************0**0************0****"
    MD = matroid_from_revlex_basis_encoding(s,3,12)
    RS = realization_space(MD, char = 0, saturate = true)
    @test is_reduced(RS) == false
end

@testset "selfprojecting realization spaces" begin
    M1 = fano_matroid() #not realizable over char 0
    M2 = uniform_matroid(3,6)
    M3 = matroid_from_bases([[1,2,7,8],[3,4,6,8],[2,4,7,8],[2,5,7,8],[3,4,6,7],[3,4,7,8],[1,5,6,7],[1,4,5,8],[1,4,5,7],[1,4,5,6],[4,5,6,8],[2,5,6,7],[4,5,7,8],[1,3,6,8],[3,4,5,7],[3,4,5,6],[1,5,6,8],[1,3,6,7],[4,5,6,7],[3,5,7,8],[2,5,6,8],[3,4,5,8],[1,3,4,7],[1,3,4,6],[1,3,4,8],[2,4,5,7],[2,4,5,6],[1,3,5,8],[1,3,5,7],[2,3,5,7],[2,3,4,8],[1,3,5,6],[2,3,5,6],[2,3,5,8],[2,4,5,8],[1,2,3,7],[1,2,3,6],[1,2,3,8],[2,3,4,7],[2,3,4,6],[1,3,7,8],[1,2,6,8],[1,2,6,7],[1,5,7,8],[2,3,7,8],[2,4,6,7],[1,2,5,8],[1,2,5,7],[1,2,5,6],[1,4,6,7],[1,2,4,7],[1,2,4,6],[3,5,6,7],[1,4,7,8],[2,3,6,8],[3,5,6,8],[2,4,6,8],[2,3,6,7],[1,2,4,8],[1,4,6,8]],8) #isomorphic to number 9 in the database (4,8), R not equal to S
    M4 = vamos_matroid() # not selfprojecting
    M5 = matroid_from_bases([[1,5,8],[4,6,7],[3,7,8],[5,7,8],[1,2,8],[3,4,8],[1,2,7],[1,3,8],[1,2,6],[4,7,8],[4,6,8],[1,3,7],[1,3,6],[3,4,7],[2,7,8],[3,4,6],[4,5,7],[4,5,6],[1,4,8],[1,4,7],[1,4,6],[2,5,7],[2,5,6],[4,5,8],[2,5,8],[3,5,7],[5,6,7],[3,5,6],[2,4,7],[2,4,6],[3,6,8],[3,5,8],[3,6,7],[2,6,8],[2,3,8],[5,6,8],[2,4,8],[2,6,7],[1,5,7],[1,5,6],[2,3,7],[2,3,6]],8) #isomorphic to number 52 in database (3,8); R = S
    M6 = matroid_from_bases([[3,4,8,9],[2,3,5,7],[1,4,6,7],[2,3,5,6],[2,3,5,9],[1,4,6,9],[1,4,6,8],[1,5,8,9],[2,5,7,9],[2,5,7,8],[2,6,7,9],[2,6,7,8],[4,5,8,9],[4,6,8,9],[3,5,8,9],[1,5,6,7],[2,5,8,9],[2,4,6,9],[2,4,6,8],[5,7,8,9],[1,5,6,9],[1,5,6,8],[1,6,8,9],[5,6,8,9],[2,5,6,7],[1,3,6,9],[1,6,7,9],[1,3,5,9],[2,5,6,9],[1,3,6,8],[1,6,7,8],[1,3,5,8],[2,5,6,8],[1,2,7,9],[1,2,7,8],[1,2,6,9],[1,3,5,6],[1,2,6,8],[1,3,6,7],[6,7,8,9],[1,3,4,7],[1,2,6,7],[1,3,4,6],[1,3,4,5],[4,5,7,9],[1,3,4,9],[1,3,4,8],[2,3,7,9],[2,3,7,8],[2,4,8,9],[1,2,3,5],[1,2,3,7],[1,2,3,6],[1,2,3,9],[1,2,3,8],[2,4,7,9],[2,4,7,8],[4,6,7,9],[1,4,7,9],[1,5,7,9],[4,6,7,8],[1,4,7,8],[1,5,7,8],[2,6,8,9],[3,5,6,7],[3,6,7,9],[3,6,7,8],[3,5,6,9],[3,5,6,8],[4,7,8,9],[3,5,7,9],[3,5,7,8],[3,4,5,7],[3,4,5,6],[3,4,7,9],[2,3,6,9],[3,4,7,8],[2,3,6,8],[3,4,5,8],[2,3,6,7],[1,2,4,7],[1,2,4,6],[1,2,4,5],[2,4,5,7],[2,4,5,6],[1,2,4,9],[1,2,4,8],[2,3,4,9],[2,3,4,8],[1,3,8,9],[5,6,7,9],[2,4,5,9],[5,6,7,8],[2,4,5,8],[3,7,8,9],[2,3,4,5],[1,2,8,9],[2,3,4,7],[2,3,4,6],[3,4,6,9],[3,4,6,8],[1,7,8,9],[1,3,7,9],[1,3,7,8],[4,5,6,9],[4,5,6,8],[1,4,5,9],[1,4,5,8],[1,2,5,9],[1,2,5,8],[3,4,6,7],[1,2,5,7],[4,5,6,7],[1,4,5,7],[1,4,5,6],[2,3,8,9]],9)# isomorphic to numer 2375 in the database (4,9), R not empty, S empty

    R1, mat_M1 = Oscar.selfprojecting_realization_matrix(M1, [1,2,4]); #this should be (nothing, nothing)
    R2, mat_M2 = Oscar.selfprojecting_realization_matrix(M2, [1,2,3]);
    R3, mat_M3 = Oscar.selfprojecting_realization_matrix(M3, [4,5,6,7]);
    @test_throws ArgumentError Oscar.selfprojecting_realization_matrix(M4, [1,2,3,5]); #this should give an error because the matroid is not selfprojecting
    R5, mat_M5 = Oscar.selfprojecting_realization_matrix(M5, [1,5,6]);
    R6, mat_M6 = Oscar.selfprojecting_realization_matrix(M6, [1,2,5,7]);#this should give (Ring, nothing)
    
    x2 = gens(R2)
    x3 = gens(R3)
    x5 = gens(R5)
    x6 = gens(R6)

    @testset "is_selfprojecting" begin
        @test is_selfprojecting(M1) == true
        @test is_selfprojecting(M2) == true
        @test is_selfprojecting(M3) == true
        @test is_selfprojecting(M4) == false
        @test is_selfprojecting(M5) == true
        @test is_selfprojecting(M6) == true
    end
    @testset "selfprojecting_realization_matrix" begin
        
        @test coefficient_ring(R2) == QQ
        @test coefficient_ring(R3) == QQ
        @test coefficient_ring(R5) == QQ
        
        @test length(x2) == 4
        @test length(x3) == 5
        @test length(x5) == 3
        @test length(x6) == 3

        @test isnothing(mat_M1)
        @test mat_M2 ==  matrix(R2, [1 0 0 1 1 1;0 1 0 1 x2[1] x2[3];0 0 1 1 x2[2] x2[4]])
        @test mat_M3 ==  matrix(R3,  [1 1 1 1 0 0 0 0; 1 x3[1] x3[3] 0 1 0 0 0; 1 x3[2] x3[4] 0 0 1 0 1; 1 x3[2] x3[4] 0 0 0 1 x3[5]])
        @test mat_M5 ==  matrix(R5, [1 1 1 1 0 0 1 1; 0 1 x5[1] x5[2] 1 0 0 0; 0 0 0 0 0 1 1 x5[3]])
        @test isnothing(mat_M6)
    end
    @testset "selfprojecting_realization_space_ideal" begin
        #R2 and R3 are quotient rings. For the selfprojecting realization ideal we need the underlying polynomial ring
        y2 = gens(base_ring(R2))
        y3 = gens(base_ring(R3))
        RR1 = base_ring(defining_ideal(realization_space(M1,char=0)))

        #add a test for a completely not realizable matroid? Haven't found an example yet... 
        @test selfprojecting_realization_ideal(M1) == ideal(RR1,[1]) 
        @test selfprojecting_realization_ideal(M2) == ideal(base_ring(R2),[y2[1]*y2[2]*y2[3] - y2[1]*y2[2]*y2[4] - y2[1]*y2[3]*y2[4] + y2[1]*y2[4] + y2[2]*y2[3]*y2[4] - y2[2]*y2[3]])
        @test selfprojecting_realization_ideal(M3) ==ideal(base_ring(R3),[y3[1]*y3[2]*y3[3] - y3[1]*y3[2]*y3[4] - y3[1]*y3[3]*y3[4] + y3[1]*y3[4] + y3[2]*y3[3]*y3[4] - y3[2]*y3[3]])
        @test_throws ArgumentError selfprojecting_realization_ideal(M4) #this should give an error because M4 is not selfprojecting
        @test selfprojecting_realization_ideal(M5) == ideal(R5,[0])
        @test selfprojecting_realization_ideal(M6) == ideal(R6,[1]) 
    end
    @testset "dimension" begin
        @test dimension(selfprojecting_realization_space(M2)) == 3
        @test dimension(selfprojecting_realization_space(M3)) == 4
        @test dimension(selfprojecting_realization_space(M5)) == 3
        @test dimension(selfprojecting_realization_space(M6)) == -inf
    end
end