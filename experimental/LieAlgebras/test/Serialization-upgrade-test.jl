@testset "Serialization.Upgrades" begin
  @testset "< 1.4.0 Upgrade" begin
    test_upgrade_folder("version_1_3_0";
      only=[
        "AbstractLieAlgebra",
        "AbstractLieAlgebraElem",
        "DirectSumLieAlgebra",
        "DirectSumLieAlgebraElem",
        "LieAlgebraModule",
        "LieAlgebraModuleElem",
        "LinearLieAlgebra",
        "LinearLieAlgebraElem",
        "Vector" => 49:97,
      ],
    )
  end
end
