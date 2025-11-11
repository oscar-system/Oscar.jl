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

@testset "realization spaces with loops" begin
  r = matroid_from_bases([[3,5],[3,4],[2,5],[2,4]],5)
  X = realization_space(r,char=0)
end

