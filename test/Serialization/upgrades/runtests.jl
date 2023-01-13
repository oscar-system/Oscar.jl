@testset "Serialization.Upgrades" begin
    L = ones(fmpq, 15)
    R, x = QQ["x"]
    p = R(L)
    loaded_p = load(joinpath(@__DIR__, "file_version<=0.11.2.json"); parent=R);

    @test p == loaded_p
end
