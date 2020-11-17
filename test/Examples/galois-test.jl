@testset "Examples.galois_group" begin

  Zx, x = ZZ["x"]
  k, a = number_field(x^5-2)
  G, _ = galois_group(k)
  @test transitive_group_identification(G) == 3
end


