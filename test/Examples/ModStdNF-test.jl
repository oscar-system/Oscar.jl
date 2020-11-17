@testset "Examples.ModStdNF" begin
  I = Oscar.DerejeGB.example_1()
  @test length(groebner_basis(I)) == 3
end


