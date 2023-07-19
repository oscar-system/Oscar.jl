@testset "LieAlgebras.RootSystem" begin
<<<<<<< HEAD
    R = RootSystem(:F,4)
  @test cartan_matrix(R) ==
    matrix(QQ, 4, 4, [2, -1, 0, 0, -1, 2, -2, 0, 0, -1, 2, -1, 0, 0, -1, 2])
    @test number_of_roots(R) == 48
    @test number_of_roots("F4") == number_of_roots(R)
end

=======
    R = RootSystem("F4")
    @test cartan_matrix(R) == MatrixSpace(QQ, 4, 4)([2, -1, 0, 0, -1, 2, -2, 0, 0, -1, 2, -1, 0, 0, -1, 2])
<<<<<<< HEAD
    
end
>>>>>>> changes to LieAlgebras export, added tests
=======
    @test size(R,1) == 48
    @test size("F4") = size(R)
end
>>>>>>> updated LieAlgebras.jl, root_systems.jl, simple_lie_algebra.jlm root_systems-test.jl, and runtests.jl
