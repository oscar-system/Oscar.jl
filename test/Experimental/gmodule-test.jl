@testset "Experimental.gmodule" begin

  G = small_group(7*3, 1)
  z = Oscar.RepPc.reps(abelian_closure(QQ)[1], G)
  @test length(z) == 5

  z = irreducible_modules(G)
  @test length(z) == 5

  z = irreducible_modules(ZZ, G)
  @test length(z) == 5

end

