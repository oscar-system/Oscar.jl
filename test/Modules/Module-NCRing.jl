@testset "modules over nc rings" begin
    E,x = exterior_algebra(QQ, 3)
    M = FreeMod(E, 2)
    @test ngens(M) == 2
end