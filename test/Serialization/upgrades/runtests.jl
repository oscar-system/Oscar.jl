@testset "Serialization.Upgrades" begin
  @testset "< 0.11.3 Upgrade" begin
    # test loading
    path = joinpath(Main.serialization_upgrade_test_path, "version_0_11_0", "QQPolyRingElem.mrdi")
    @test load(path) isa QQPolyRingElem

    L = ones(QQFieldElem, 15)
    R, x = QQ[:x]
    p = R(L)
    loaded_p = load(path; params=R);
    @test p == loaded_p

    loaded_container1 = load(joinpath(Main.serialization_upgrade_test_path, "version_0_11_0", "Containers.mrdi"))
    loaded_container2 = load(joinpath(Main.serialization_upgrade_test_path, "version_0_11_3", "Containers.mrdi"))
    @test loaded_container1 == loaded_container2
    @test loaded_container1 == (r = QQFieldElem(1, 2), m = QQFieldElem[1//2 1; 0 1], t = (1, 2, 3))             
  end

  @testset "< 0.12.0 Upgrade" begin
    # test loading
    path = joinpath(Main.serialization_upgrade_test_path, "version_0_11_3", "fqPolyRepPolyRingElem.mrdi")
    @test load(path) isa fqPolyRepPolyRingElem
    
    Zt, t = polynomial_ring(residue_ring(ZZ, 2)[1], :t)
    Fin, d = Nemo.Native.finite_field(t^2 + t + 1)
    Rx, x = Fin[:x]
    p = x^2 + d * x + 1
    loaded_p = load(path; params=Rx);
    @test p == loaded_p
  end

  @testset "< 1.2.0 Upgrade" begin
    # test loading
    loaded = load(joinpath(Main.serialization_upgrade_test_path, "version_1_1_0", "Dict.mrdi"))
    @test loaded isa Dict
    @test length(loaded) == 7
  end

  @testset "< 1.3.0 Upgrade" begin
    loaded = load(joinpath(Main.serialization_upgrade_test_path, "version_1_2_0", "FqField-1.mrdi"))
    @test loaded isa FqField
    @test order(loaded) == 2
    loaded = load(joinpath(Main.serialization_upgrade_test_path, "version_1_2_0", "FqField-2.mrdi"))
    @test loaded isa FqField
    @test order(loaded) == 4
    Oscar.reset_global_serializer_state()
    # This is the polynomial in the "A FAIR File Format for Mathematical Software" paper
    loaded = load(joinpath(Main.serialization_upgrade_test_path, "version_1_0_5", "FqMPolyRingElem.mrdi"))
    @test loaded isa FqMPolyRingElem
  end

  @testset "< 1.4.0 Upgrade" begin
    test_1_4_0_upgrade(; exclude=[
      # upgrading the following is tested in experimental/LieAlgebras/test/Serialization-upgrade-test.jl
      "AbstractLieAlgebra", "AbstractLieAlgebraElem",
      "DirectSumLieAlgebra", "DirectSumLieAlgebraElem",
      "LieAlgebraModule", "LieAlgebraModuleElem",
      "LinearLieAlgebra", "LinearLieAlgebraElem",
      "Vector" => 49:97,
    ])
  end
end
