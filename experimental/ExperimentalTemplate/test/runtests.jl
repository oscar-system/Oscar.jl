@testset "Experimental template tests" begin
  S = Oscar.ExampleStruct(5)
  @test 5 == Oscar.my_access_func(S)
end
