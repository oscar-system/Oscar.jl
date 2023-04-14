@testset "Experimental template tests" begin
  S = oscar.ExampleStruct(5)
  @test 5 == oscar.content(S)
end
