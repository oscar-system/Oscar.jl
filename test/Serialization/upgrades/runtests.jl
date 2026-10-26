@testset "Serialization.Upgrades" begin
  @testset "load containers serialized with 0.11.3" begin
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

  @testset "load polynomial serialized with 0.11.3" begin
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

  @testset "load files serialized with 1.1.0" begin
    # test loading
    loaded = load(joinpath(Main.serialization_upgrade_test_path, "version_1_1_0", "Dict.mrdi"))
    @test loaded isa Dict
    @test length(loaded) == 7
  end

  @testset "load files serialized with 1.2.0" begin
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

  @testset "load files serialized with 1.3.0" begin
    test_upgrade_folder("version_1_3_0"; exclude=[
      # upgrading the following is tested in experimental/LieAlgebras/test/Serialization-upgrade-test.jl
      "AbstractLieAlgebra", "AbstractLieAlgebraElem",
      "DirectSumLieAlgebra", "DirectSumLieAlgebraElem",
      "LieAlgebraModule", "LieAlgebraModuleElem",
      "LinearLieAlgebra", "LinearLieAlgebraElem",
      "Vector" => 49:97,
    ])
  end

  @testset "load files serialized with 1.5.0" begin
    test_upgrade_folder("version_1_5_0")
  end

  @testset "load files serialized with 1.6.0" begin
    test_upgrade_folder("version_1_6_0")
  end

  @testset "load files serialized with 1.7.0" begin
    test_upgrade_folder("version_1_7_0"; exclude=[
      # upgrading the following is tested in experimental/AlgebraicStatistics/test/serialization-upgrades.jl
      "GroupBasedPhylogeneticModel",
    ])

    # The main point is to test the relation between loaded objects.
    loaded = load(joinpath(Main.serialization_upgrade_test_path, "version_1_7_0", "Tuple-3", "0.mrdi"))
    @test length(loaded) == 6
    (F, U, G, H, P, S) = map(parent, loaded)

    @test full_group(F)[1] === F
    @test full_group(U)[1] === F
    @test free_group(G) === F
    @test full_group(G)[1] === G
    @test full_group(H)[1] === G

    @test full_group(S)[1] === P

    @test GAP.Globals.FamilyObj(GapObj(F)) === GAP.Globals.FamilyObj(GapObj(U))
    @test GAP.Globals.FamilyObj(GapObj(G)) === GAP.Globals.FamilyObj(GapObj(H))
    @test GAP.Globals.FamilyObj(GapObj(F)) !== GAP.Globals.FamilyObj(GapObj(G))
    @test GAP.Globals.FamilyObj(GapObj(P)) === GAP.Globals.FamilyObj(GapObj(S))
  end
  @testset "load files serialized with 1.8.2" begin
    test_upgrade_folder("version_1_8_2")
  end

  @testset "load groups serialized via GapObj" begin
    # full groups viewed as `SubFPGroup` or `SubPcGroup`, which are created by
    # `as_sub_fp_group` and `as_sub_pc_group`, and free groups of infinite rank
    loadold(name) = load(joinpath(Main.serialization_upgrade_test_path, "version_1_9_0_dev", name))

    Oscar.reset_global_serializer_state()
    x = loadold("SubFPGroupElem-1/0.mrdi")
    @test x isa SubFPGroupElem
    G = full_group(parent(x))[1]
    @test G isa FPGroup
    @test order(G) == 60
    @test order(parent(x)) == 60
    @test order(x) == 5

    Oscar.reset_global_serializer_state()
    G, H, x = loadold("Tuple-3/0.mrdi")
    @test G isa FPGroup
    @test full_group(H)[1] === G
    @test full_group(parent(x))[1] === G
    @test order(parent(x)) == 60

    Oscar.reset_global_serializer_state()
    y = loadold("SubFPGroupElem-1/1.mrdi")
    @test y isa SubFPGroupElem
    F = full_group(parent(y))[1]
    @test F isa FPGroup
    @test free_group(F) === F
    @test ngens(parent(y)) == 2

    Oscar.reset_global_serializer_state()
    F, U, y = loadold("Tuple-3/1.mrdi")
    @test full_group(U)[1] === F
    @test full_group(parent(y))[1] === F

    Oscar.reset_global_serializer_state()
    z = loadold("SubPcGroupElem-2/0.mrdi")
    @test z isa SubPcGroupElem
    @test full_group(parent(z))[1] isa PcGroup
    @test order(parent(z)) == 24

    Oscar.reset_global_serializer_state()
    P, P1, z = loadold("Tuple-3/2.mrdi")
    @test P isa PcGroup
    @test full_group(P1)[1] === P
    @test full_group(parent(z))[1] === P

    Oscar.reset_global_serializer_state()
    @test_throws "infinite rank" loadold("FPGroup-4/0.mrdi")
    @test_throws "infinite rank" loadold("FPGroupElem-5/0.mrdi")
    Oscar.reset_global_serializer_state()
  end

end
