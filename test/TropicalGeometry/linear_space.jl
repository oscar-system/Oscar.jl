@testset "src/TropicalGeometry/linear_space.jl" begin

    @testset "det(A::Generic.MatSpaceElem{<:TropicalSemiringElem})" begin
        M = matrix(QQ,[0 -3 0 -3 0 3 9 9 4 -9 -10 -6;
                       -8 -7 -10 -10 0 -4 0 -2 -8 0 1 0;
                       -2 7 0 1 -5 0 -5 0 0 -10 -3 6;
                       9 2 -6 5 0 -6 -5 -6 -2 -2 4 2;
                       -1 0 5 7 9 1 -5 2 -6 0 0 0;
                       0 -9 8 3 -10 -9 4 -2 8 -9 -6 2])

        TropL = tropical_linear_space(M)
        @test dim(TropL) == rank(M)
        @test ambient_dim(TropL) == ncols(M)
        @test lineality_dim(TropL) == 1
    end

    @testset "tropical linear space attributes" begin
        R,(w,x,y,z) = GF(2)[:w, :x, :y, :z]
        I = ideal([w+x+y+z])
        TropL = tropical_linear_space(I)
        @test has_attribute(TropL, :tropical_pluecker_vector)

        plueckerVector = algebraic_pluecker_vector(TropL)
        delete!(TropL.__attrs, :algebraic_pluecker_vector)
        @test issetequal(plueckerVector,algebraic_pluecker_vector(TropL))

        TT = tropical_semiring()
        M = matrix(TT,[1 2 3; 2 3 4])
        TropL = tropical_linear_space(M)
        plueckerVector = tropical_pluecker_vector(TropL)
        delete!(TropL.__attrs, :tropical_pluecker_vector)
        @test issetequal(plueckerVector,tropical_pluecker_vector(TropL))
    end

    @testset "tropical linear space from graphs" begin
        G = complete_graph(3)
        TropG = tropical_linear_space(G)
        @test f_vector(TropG) == [0,1,3]
    end

end
