@testset "Serialization.Upgrades" begin
  @testset "load Lie GroupBasedPhylogeneticModel serialized wtih 1.7.0" begin
    test_upgrade_folder("version_1_7_0";
      only=[
        "GroupbasedPhylogeneticModel",
      ],
    )
  end
end
