@testset "Serialization.Upgrades" begin
  @testset "load Lie Algebra files serialized wtih 1.3.0" begin
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
