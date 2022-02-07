@testset "Experimental.gmodule" begin

  G = small_group(5*11, 1)
  z = Oscar.RepPc.reps(abelian_closure(QQ)[1], G)
  @test length(z) == 7

  z = irreducible_modules(G)
  @test length(z) == 7

  z = irreducible_modules(ZZ, G)
  @test length(z) == 7

end

