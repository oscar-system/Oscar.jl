@testset "Experimental.gmodule" begin

  G = small_group(7*3, 1)
  z = Oscar.RepPc.reps(abelian_closure(QQ)[1], G)
  @test length(z) == 5

  z = irreducible_modules(G)
  @test length(z) == 5

  z = irreducible_modules(ZZ, G)
  @test length(z) == 5

  l = irreducible_modules(AnticNumberField, small_group(48, 17), minimal_degree = true)
  ds = degree.(base_ring.(l))
  @test length(l) == 12
  @test count(isequal(1), ds) == 8
  @test count(isequal(2), ds) == 4

  l = irreducible_modules(AnticNumberField, small_group(48, 29), minimal_degree = true)
  ds = degree.(base_ring.(l))
  @test length(l) == 8
  @test count(isequal(1), ds) == 6
  @test count(isequal(2), ds) == 2

  G = SL(2, 3)
  @test length(Oscar.RepPc.reps(abelian_closure(QQ)[1], G)) == 7
end

