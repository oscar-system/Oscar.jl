@testset "Serialization.Upgrades" begin
  @testset "< 0.11.3 Upgrade" begin
    L = ones(QQFieldElem, 15)
    R, x = QQ["x"]
    p = R(L)
    loaded_p = load(joinpath(@__DIR__, "file_version<=0.11.2.json"); params=R);
    @test p == loaded_p
  end

  @testset "< 0.12.0 Upgrade" begin
    Zt, t = polynomial_ring(residue_ring(ZZ, 2), "t")
    Fin, d = finite_field(t^2 + t + 1)
    Rx, x = Fin["x"]
    p = x^2 + d * x + 1
    loaded_p =  load(joinpath(@__DIR__, "file_version<=0.12.0.json"); params=Rx);
    @test p == loaded_p
  end
end
