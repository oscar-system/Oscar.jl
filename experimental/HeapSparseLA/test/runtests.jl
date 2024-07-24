@testset "heap based sparse linear algebra" begin
  S = Oscar.ExampleStruct(5)
  @test 5 == Oscar.my_access_func(S)
end
