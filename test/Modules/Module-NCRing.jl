@testset "modules over nc rings" begin
    E,x = exterior_algebra(QQ, 3)
    M = free_module(E, 3)
    @test ngens(M) == 3
    N = Oscar.SubModuleOfFreeModule(M, gens(M)[1:2])
    @test gens(N) == gens(M)[1:2]
end