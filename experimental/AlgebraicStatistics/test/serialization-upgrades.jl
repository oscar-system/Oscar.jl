@testset "Serialization.Upgrades" begin
  @testset "load GroupBasedPhylogeneticModel serialized wtih 1.7.0" begin
    test_upgrade_folder("version_1_7_0";
      only=[
        "GroupBasedPhylogeneticModel",
      ],
    )
  end
  # TODO: add PhylogeneticModel and GroupBasedPhylogeneticModel example files to
  # oscar-system/serialization-upgrade-tests repo under version_1_8_0+1/ before enabling this test
  # @testset "load PhylogeneticModel and GroupBasedPhylogeneticModel serialized with 1.8.0+1" begin
  #   test_upgrade_folder("version_1_8_0+1";
  #     only=[
  #       "GroupBasedPhylogeneticModel",
  #       "PhylogeneticModel",
  #     ],
  #   )
  # end
end
