@testset "Markov basis" begin
    M =  [-1  -1  1   1   1   0  0   0   -1; -1  0   0   0   0   2  0   2   -3;
          -1  0   0   0   2   0  0   0   -1; -1  0   0   1   0   1  0   0   -1;
          -1  0   0   2   0   0  -1  0   0; -1  0   1   0   0   0  1   0   -1;
          -1  1   -1  0   1   1  0   0   -1; -1  1   0   0   0   0  0   1   -1;
          -1  1   0   0   0   0  1   -1  0;  -1  1   0   1   -1  0  0   0   0;
          -1  2   -1  0   0   0  0   0   0; 0   -1  0   0   0   2  0   1   -2;
          0   -1  0   0   1   1  0   0   -1; 0   -1  0   1   0   1  0   -1  0;
          0   -1  0   1   1   0  -1  0   0; 0   -1  1   0   0   0  0   1   -1;
          0   -1  1   1   -1  0  0   0   0; 0   0   -1  0   0   2  0   0   -1;
          0   0   0   -1  0   1  0   2   -2; 0   0   0   -1  0   1  1   0   -1;
          0   0   0   -1  1   0  0   1   -1; 0   0   0   -1  1   0  1   -1  0;
          0   0   0   0   -1  1  0   1   -1; 0   0   0   0   0   0  -1  2   -1;
          0   1   -1  0   -1  1  0   0   0; 1   -1  0   0   -1  1  0   0   0;
          1   0   -1  -1  0   1  0   0   0]

    S = [homogenize(v, 1) for v in lattice_points(cube(2))]
    T = transpose(matrix(ZZ, Vector{Vector{fmpz}}(S)))
    L = kernel(T)[2]

    @test M == markov_basis(cube(2))
    @test M == markov_basis(S)
    @test M == markov_basis(T)
    @test M == markov_basis(L; use_kernel=false)
end
