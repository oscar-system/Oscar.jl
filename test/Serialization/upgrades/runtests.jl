@testset "Serialization.Upgrades" begin
    R, x = QQ["x"]
    p = x^2 - x + fmpq(1, 2)
    loaded_p = load(joinpath(@__DIR__, "file_version<=0.11.1.json"); parent=R);

    @test p == loaded_p
end
