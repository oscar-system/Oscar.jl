@testset "Experimental.gmodule" begin

  G = small_group(7*3, 1)
  z = Oscar.RepPc.reps(abelian_closure(QQ)[1], G)
  @test length(z) == 5

  z = irreducible_modules(G)
  @test length(z) == 5

  z = irreducible_modules(ZZ, G)
  @test length(z) == 5

  l = irreducible_modules(AnticNumberField, small_group(48, 17), minimal_degree = true)
end

