@testset "modules over nc rings" begin
    E,x = exterior_algebra(QQ, 3)
    M = free_module(E, 3)
    @test ngens(M) == 3
    N = Oscar.SubModuleOfFreeModule(M, gens(M)[1:2])
    @test gens(N) == gens(M)[1:2]
    H = hom(M,M,gens(M))
    #fails here
    #W = free_module(E,2)
    #G = hom(W,M, gens(M)[1:2])
    #this is becuase the line below fails
    #@test parent(gens(M[1])) == M
end