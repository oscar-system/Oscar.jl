@testset "Serialization.Upgrades" begin
  @testset "< 0.11.3 Upgrade" begin
    # test loading
    load(joinpath(@__DIR__, "file_version_less_than_0.11.2.json"))

    L = ones(QQFieldElem, 15)
    R, x = QQ[:x]
    p = R(L)
    loaded_p = load(joinpath(@__DIR__, "file_version_less_than_0.11.2.json"); params=R);
    @test p == loaded_p
  end

  @testset "< 0.12.0 Upgrade" begin
    # test loading
    load(joinpath(@__DIR__, "file_version_less_than_0.12.0.json"));
    
    Zt, t = polynomial_ring(residue_ring(ZZ, 2)[1], :t)
    Fin, d = Nemo.Native.finite_field(t^2 + t + 1)
    Rx, x = Fin[:x]
    p = x^2 + d * x + 1
    loaded_p =  load(joinpath(@__DIR__, "file_version_less_than_0.12.0.json"); params=Rx);
    @test p == loaded_p
  end

  @testset "< 1.2.0 Upgrade" begin
    # test loading
    loaded = load(joinpath(@__DIR__, "file_version_less_than_1.2.0.json"));
    @test loaded isa Dict
  end

  @testset "< 1.3.0 Upgrade" begin
    load(joinpath(@__DIR__, "GF_2_2.json"));
    load(joinpath(@__DIR__, "GF_2.json"));
    Oscar.reset_global_serializer_state()
    load(joinpath(@__DIR__, "poly1.0.5.json"));
  end

  @testset "< 1.4.0 Upgrade" begin
    test_1_4_0_upgrade(; exclude=[
      # upgrading the following is tested in experimental/LieAlgebras/test/Serialization-upgrade-test.jl
      "AbstractLieAlgebra", "AbstractLieAlgebraElem",
      "DirectSumLieAlgebra", "DirectSumLieAlgebraElem",
      "LieAlgebraModule", "LieAlgebraModuleElem",
      "LinearLieAlgebra", "LinearLieAlgebraElem",
      "Vector" => 49:97,
      # upgrading the following is tested in experimental/QuadFormAndIsom/test/runtests.jl
      "ZZLatWithIsom",
    ])
  end
end
