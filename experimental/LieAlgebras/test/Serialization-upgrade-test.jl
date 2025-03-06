@testset "Serialization.Upgrades" begin
  @testset "< 1.4.0 Upgrade" begin
    test_1_4_0_upgrade(;
      only=[
        "AbstractLieAlgebra",
        "AbstractLieAlgebra",
        "DirectSumLieAlgebra",
        "DirectSumLieAlgebraElem",
        "LieAlgebraModule",
        "LieAlgebraModuleElem",
        "LinearLieAlgebra",
        "LinearLieAlgebraElem",
      ],
    )
  end
end
